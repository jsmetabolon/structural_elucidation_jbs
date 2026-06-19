"""Elucidate MLIMS unknown (X) compounds from stored MS2 spectra.

Standalone CLI alternative to the full discovery-to-elucidation pipeline in
elucidator.py. Used in spec2mol v21 inference to batch-annotate library
unknowns and export summary spreadsheets with structure predictions.
"""

import argparse

import elucidation.chemFacts as cf
import numpy as np
from elucidation.binning import bin_for_mass
from elucidation.feature import Feature
from elucidation.library_services.mlimsLibraryService import MlimsLibraryService as LibraryService
from elucidation.ms1Analyzer import MS1Analyzer
from elucidation.ms2Consensus import MS2Consensus
from elucidation.ms2Ion import MS2Ion
from elucidation.postgresSmilesDb import PostgresSmilesDb
from elucidation.rdkitUtils import to_excel_with_molecules
from elucidation.conveniences import *
from elucidation.tableLoaders import get_config_var
from elucidation.userFeatureReporter import UserFeatureReporter
from mdsutils.easyQuery import EasyQuery
from elucidation.config import BASE_DIR, DELIVERY_DIR


def get_arguments(argv):
    """
    Use argparse get arguments from command line.
    :param argv: strings from command line.
    :return: a Namespace object with option values.
    """
    parser = argparse.ArgumentParser(description="""Library Unknown Elucidation.""")
    parser.add_argument("-v", "--verbosity", action="count", help="increase output verbosity")
    parser.add_argument(
        "-l",
        "--library",
        choices=[209, 400, 402, 305],
        help="library, only 209, 400, 402, 305 are allowed",
        type=int,
    )
    parser.add_argument(
        "-c",
        "--comp_id",
        action="append",
        help="elucidate this comp_id. May be repeated. (Might not be an X!)",
        default=[],
    )

    parser.add_argument(
        "-C",
        "--comp_id_file",
        type=str,
        help="elucidate the comp_ids named in this file",
        default="",
    )
    return parser.parse_args(argv[1:])


def get_unknowns_from_mlims(lib_id):
    df = LibraryService.get_instance().get_unnamed_ms2(lib_id)
    df = df[(df.mass_ratio >= 1) & (df.comp_id < 2000000)].copy()
    d = dict(df.groupby("comp_id").apply(lambda tab: list(zip(tab.mass, tab.mass_ratio))))
    df = df[
        qw("comp_id lib_id parent_mass unknown_name ri created_on comments putative_details")
    ].drop_duplicates()
    df["raw_spectrum"] = df.comp_id.apply(d.get)
    return df


def build_fake_spectrum(spectrum, ri, precursor):
    masses, intensities = lzip(*spectrum)
    return dict(
        masses=np.array(masses),
        intensities=np.array(intensities),
        ri=ri,
        precursor=precursor,
        scan_num=0,
        rt=0,
        max_int=max(intensities),  # should be 1
        sample_name="dummy",
    )


def analyze_feature(d):
    d["mz"] = d["parent_mass"]
    d["report_name"] = d["unknown_name"]
    d["ms1_spectra"] = []
    d["ms2_spectra"] = []
    d["peak_df"] = pd.DataFrame(
        [
            {
                "mass": d["mz"],
                "ri": d["ri"],
                "sample_name": "dummy",
                "area": 1,
                "apex_rt": 0,
                "task_id": -1,
            }
        ]
    )
    d["ms2_counts"] = {"dummy": 1}
    d["ms2_spectra"] = [build_fake_spectrum(d["raw_spectrum"], d["ri"], d["mz"])]
    d["ri_mz_tuples"] = []
    d["study_id"] = "X_compounds"
    d["serial"] = f"{d['report_name']}: comp_id {d['comp_id']}"
    d["discovery_time"] = d["created_on"]
    d["discoverer_version"] = "IonTracker"

    comp_id = d["comp_id"]
    f = Feature(d, opts)
    f.serial = d["report_name"]
    f.comments = [d["comments"], d["putative_details"]]
    f.study_id = "Library Unknowns"
    f.get_candidate_fmlas()
    ms1 = MS1Analyzer(f, [])
    f.ms1 = ms1
    ms1.library_unknown_extract_buddies(d["comp_id"], d["lib_id"], ppm=3)
    ms1.analyze_isotopes()
    ms1.build_buddy_library_df()
    f.iso_count = ms1.iso_count
    f.complete_candidate_fmlas()

    f.create_fmla_consistency_reports()
    f.create_library_report()
    pgsmiles = PostgresSmilesDb.get_instance()
    df = pgsmiles.get_records_for_neutral_mass(cf.neutralize_precursor_mz(f.mz, f.lib_id))
    f.smilesdb_df = df[qw("smiles name source formula exact_mass")].copy()

    ms2 = MS2Consensus(f)
    f.consensus = ms2
    spectrum = d["raw_spectrum"]
    mzs, ints = lzip(*spectrum)

    ms2.ms2_ion_dict = dict()
    for mz in mzs:
        b = bin_for_mass(mz)
        ion = MS2Ion(b, f.consensus)
        ion.mz = mz
        ion.lib_id = f.lib_id
        ms2.ms2_ion_dict[b] = ion
        ion.build_structure_assignment_reports()

    ms2.ms2_ions = lmap(ms2.ms2_ion_dict.get, sorted(ms2.ms2_ion_dict.keys()))
    ms2.consensus_df = pd.DataFrame(
        {"mz": mzs, "norm_int": ints, "bucket": bin_for_mass(np.array(mzs))}
    )
    ms2.build_ms2_library_report()
    ms2.build_ms2_smiles_prediction_report()
    ms2.build_chemont_prediction_report()
    # ms2.build_general_ms2_smiles_match_report()
    # ms2.build_zapped_ms2_smiles_match_report()
    ms2.build_ms2_library_report()
    # ms2.build_zapped_ms2_library_report()
    ms2.build_structure_prediction_report()
    ms2.build_predicted_spectrum_match_report()
    d.update(f.as_dict())
    fr = UserFeatureReporter(d)
    d = fr.write_report(opts)
    d["comp_id"] = comp_id
    return d


def add_links(df):
    """
    Adds a column of report hyperlinks to a dataframe.
    :param df: dataframe receiving links.
    :return: the updated dataframe.
    """
    df.serial = [
        f'=HYPERLINK("file:{file}","{feature}")' for file, feature in zip(df.report, df.serial)
    ]
    return df.drop("report", axis=1)


def format_summary_spreadsheet(summaries, opts):
    """
    create and reorganize spreadsheet summarizing all features.
    :param summaries: list of dicts, one for each feature.
    :return: summary spreadsheet.
    """
    summaries = pd.DataFrame(summaries)
    summaries = add_links(summaries)
    k = len(summaries)
    summaries["date"] = [opts.tag] * k
    summaries["lib_id"] = [opts.library] * k
    summaries["study"] = [opts.study] * k
    area_df = pd.DataFrame(list(summaries["areas"]))

    first_columns = qw(
        "date study lib_id serial mz ri cand_fmlas "
        + "max_harmony num_library_matches isotopes "
        + "max_smiles_db_score max_smiles_db largest_smiles_prediction fraction_predicted"
    )
    claimed = set(first_columns)
    leftover_columns = [c for c in summaries.columns if c not in claimed and c != "areas"]
    cols = first_columns + leftover_columns
    cols = [c for c in summaries if c in summaries.columns]
    return summaries[cols].drop(qw("internal_standard num_peaks num_spectra"), axis=1)


def get_ms2_df(opts):
    df = get_unknowns_from_mlims(
        opts.library
    )  # TODO It would be useful to be able to elucidate named compounds for debugging.
    # If required, get individual comp_ids according to command line options.
    comp_ids = []
    if opts.comp_id:
        comp_ids.extend(opts.comp_id)
    if opts.comp_id_file:
        with open(opts.comp_id_file) as f:
            comp_ids.extend(map(str.strip, f.readlines()))
    if comp_ids:
        comp_ids = set(map(int, comp_ids))
        df = df[df.comp_id.isin(comp_ids)]
        missing = comp_ids - set(df.comp_id)
        if missing:
            print("Couldn't find MS2s for these comp_ids:", missing)
        if len(df) == 0:
            print("Couldn't find MS2s for any compounds")
            exit(-1)
    return df


if __name__ == "__main__":
    local_base_dir = BASE_DIR
    delivery_base_dir = DELIVERY_DIR

    opts = get_arguments(sys.argv)
    opts.gitbranch = get_git_branch()
    lib_id = opts.library
    eq = EasyQuery(db="mlims")
    df = get_ms2_df(opts)

    opts.study = f"X-{lib_id}"
    opts.now = now_string()
    opts.out_dir = f"{local_base_dir}/X_{opts.library}/{opts.now}"
    opts.feature_dir = f"{opts.out_dir}/features"
    os.makedirs(opts.feature_dir, exist_ok=True)
    opts.tag = f"X_{opts.library}_{opts.now}"
    opts.annot_df = None

    summaries = []
    for j, row in enumerate(df.iterrows()):
        print(f"unknown {j+1} / {len(df)}", now_string())
        summaries.append(analyze_feature(dict(row[1])))
        # pprint(summaries[-1])
        # print()

    summary_df = format_summary_spreadsheet(summaries, opts)
    columns = qw(
        "date study lib_id comp_id serial mz ri cand_fmlas max_harmony "
        + "isotopes ion_names inferred_masses lib_buddies "  # MS1 stuff
        + "largest_smiles_prediction fraction_predicted num_smiles_preds "  # MS2SMILES
        + "best_spectrum_prediction num_best_predictions max_spectrum_score "  # SMILES2MS2
        + "chemont_preds"  # MS2 => Chemical Ontology
    )
    summary_df = summary_df[columns]
    summary_df.to_csv(f"{opts.out_dir}/feature_summaries_{opts.tag}.csv", index=False)
    to_excel_with_molecules(
        summary_df,
        f"{opts.out_dir}/debug_feature_summaries_{opts.tag}.xlsx",
        qw("best_spectrum_prediction largest_smiles_prediction"),
        size=(150, 100),
    )
    summary_df.serial = [s.replace(local_base_dir, delivery_base_dir) for s in summary_df.serial]
    to_excel_with_molecules(
        summary_df,
        f"{opts.out_dir}/feature_summaries_{opts.tag}.xlsx",
        qw("best_spectrum_prediction largest_smiles_prediction"),
        size=(150, 100),
    )
    print(f"finished; output in {opts.out_dir}")
