import argparse
import gc
import logging
import os
import subprocess
import sys
import time
import json
from collections import defaultdict
from functools import reduce
from glob import glob
from operator import add
from traceback import print_exception

exec("import elucidation.initPackageOptions")  # pycharm won't optimize this away

import pandas as pd

import regex as re
from elucidation.binning import within_20ppm
from elucidation.conveniences import (
    qw,
    resource_string,
    pickle_load,
    get_git_branch,
    pickle_dump,
    now_string,
    ram_gb_in_use,
)
from elucidation.feature import Feature
from elucidation.htmlUtils import heading, boxit, NullReporter
from elucidation.ms1Clusterer import Ms1ClusterAnalyzer, Ms1Cluster, Ms1ClusterReporter
from elucidation.userFeatureReporter import UserFeatureReporter as FeatureReporter
from elucidation.reformatter import reformat
from elucidation.tableLoaders import system_base_dir
from elucidation.libraryService import LibraryService
from elucidation.rdkitUtils import to_excel_with_molecules
from elucidation.config import BASE_DIR, DELIVERY_DIR
from mdsutils.easyQuery import EasyQuery


def elapsed_time(comment, run_report):
    """
    prints and logs the current elapsed time.
    :param comment: annotation next to the time, e.g., what just finished.
    :return:
    """
    s = resource_string()
    logging.info(f"{s}: {comment}")
    run_report.comment(f"<h4>{s}: {comment}</h4>")


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


def make_argv_from_opts(opts):
    """check an argparse namespace against a module's get_args method.
    Ideally, there would be something built in to argparse, but no such luck.
    This tries to reconstruct the arg list that argparse.parse_args would expect
    """
    arg_list = [[k, v] for k, v in sorted(vars(opts).items())]
    argparse_formatted_list = []
    for l in arg_list:
        ####  deal with flag arguments (store true/false)
        if l[1] == True:
            argparse_formatted_list.append("--{}".format(l[0]))
        elif l[1] == False or l[1] is None:
            pass  # dont add this arg
        elif l[1] == []:
            pass
        elif l[1] == "":
            pass
        elif l[0] == "workers":
            pass
        # add the named arguments
        else:
            argparse_formatted_list.append("--{}".format(l[0]))
            argparse_formatted_list.append(str(l[1]))
    return argparse_formatted_list


def get_arguments():
    """
    Use argparse get arguments from command line.
    :return: a Namespace object with option values.
    """
    parser = argparse.ArgumentParser(description="""Feature Elucidation.""")
    parser.add_argument("-v", "--verbosity", action="count", help="increase output verbosity")
    parser.add_argument("-s", "--study", "--ssid", help="sample set id", type=str)
    parser.add_argument(
        "-l",
        "--library",
        choices=[209, 400, 402, 305],
        help="library, only 209, 305, 400, 402 are allowed",
        type=int,
    )
    parser.add_argument("--ri_lo", help="elucidate features with RI >= this", type=int, default=0)
    parser.add_argument(
        "--ri_hi", help="elucidate features with RI <= this", type=int, default=1000000
    )
    parser.add_argument(
        "--mz",
        "-m",
        help="only features within 20ppm of this mass are considered. May be repeated",
        action="append",
        type=float,
        default=[],
    )

    parser.add_argument(
        "--min_ms2_intensity",
        help="minimum relative intensity for retaining an MS2 ion in consensus",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--out_dir",
        help="directory to which output will be delivered",
        type=str,
        default="",
    )
    parser.add_argument(
        "-f",
        "--feature",
        action="append",
        help="elucidate this feature. May be repeated.  May include wild cards: 'F0123', 'F1*1', 'F1?1?'",
        default=[],
    )

    parser.add_argument(
        "-F",
        "--feature_file",
        type=str,
        help="elucidate the features named in this file",
        default="",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="if set, any exception will abort the program.  If not, program continues with next feature.",
    )

    parser.add_argument(
        "--discovery_dir",
        help="folder containing discovered features to elucidate. "
        + "Default: most recent of form discoveries_????-??-??-????",
        type=str,
        default="",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=7,  # seems to be only slightly slower than 7.
        help="number of workers. default=7. 0=>do it yourself",
    )

    parser.add_argument(
        "--share",
        type=str,
        help="This is used by the main program to start the workers.",
        default="",
    )

    return parser.parse_args()


def update_trace_spreadsheet(summaries, trace_csv):
    """
    Update the intermediate summary spreadsheet.  Called very few (10) features.
    :param summaries: A list of summary rows for features.
    :return:
    """
    try:
#        pd.DataFrame(summaries).to_csv(f"{out_dir}/trace_features{worker_number}.csv", index=False)  # current directory
        pd.DataFrame(summaries).to_csv(trace_csv, index=False)  # current directory
    except Exception as e:  # Maybe the spreadsheet is open.
        logging.warning(e)


def print_brief_feature_summary(summary, run_report):  # TODO push this down into Feature.
    """
    print a brief summary of feature on terminal and log it in run_report
    :param summary:
    :return:log
    """
    brief = {
        k: str(v)
        for k, v in summary.items()
        if not re.match(r"[A-Z]", k) and k not in ["report", "areas"]
    }
    for k, v in brief.items():
        logging.info(f" {k} => {v}")
    df = pd.DataFrame([brief])
    run_report.htmls(reformat(df))


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
    summaries["ms1_cluster"] = [""] * k
    area_df = pd.DataFrame(list(summaries["areas"]))
    return summaries, area_df


def perform_ms1_clustering(summary_df, run_report, opts):
    """
    Identify clusters of features based on RI and MS1 inferred mass and creates a report for each cluster.
    @param summary_df: The proto-summary spreadsheet
    @type summary_df: dataframe
    @return: summary_df with a new column added with links to the MS1 cluster reports.
    @rtype: dataframe
    """
    df = summary_df
    cpa = Ms1ClusterAnalyzer(df)
    cl_objs = []
    f_to_clink_dict = defaultdict(str)
    for i, cl in enumerate(cpa.clusters):
        id = f"C{i + 1:03d}"
        s = f"{i + 1}/{len(cpa.clusters)}: {id:} {cl}"
        logging.info(s)
        run_report.htmls(heading(s))
        features = [pickle_load(glob(f"{opts.elucidations_dir}/e{f}.pkl")[0]) for f in cl]

        cl_obj = Ms1Cluster(id, features).get_me()
        if cl_obj is None:
            continue
        cr = Ms1ClusterReporter(cl_obj)
        cr.write_report(opts)
        cl_objs.append(cl_obj)
        link = f'=HYPERLINK("file:{cr.report_file}", "{id}")'
        for f in cl:
            f_to_clink_dict[f] = link

    df.ms1_cluster = [f_to_clink_dict[re.search(r"[FGHJ]\d\d\d\d", s).group(0)] for s in df.serial]
    return df


def get_files_for_feature_opt(f):
    """
    Finds the discovery pickle file matching the given feature specification, which may include wildcards.
    @param f: feature specification, such as "F0013"  or "F1*3"
    @type f: str
    @return: list of all matching files
    @rtype: list of str
    """
    return glob(f"{opts.discovery_dir}/{f}.pkl")


def get_feature_files(opts):
    all_files = sorted(glob(f"{opts.discovery_dir}/*.pkl"))
    if (
        opts.feature == []
        and opts.feature_file == ""
        and opts.ri_lo == 0
        and opts.ri_hi == 1000000
        and opts.mz == []
    ):
        return all_files
    if opts.feature_file != "":
        with open(opts.feature_file) as f:
            features = set(map(str.strip, f.readlines()))

        def is_a_feature(fn):
            m = re.search(r"([FGHJ]\d{4})", fn)
            return m is not None and m.group(1) in features

        our_files = [fn for fn in all_files if is_a_feature(fn)]
        return our_files

    elif opts.ri_lo == 0 and opts.ri_hi == 1000000 and opts.mz == []:
        our_files = reduce(add, map(get_files_for_feature_opt, opts.feature))
        return sorted(set(our_files))
    else:  # If there are mass and/or RI restrictions, we will have to scan all the files.
        # This isn't going to work for Lund because the RI didn't get put into the object.
        our_files = reduce(add, map(get_files_for_feature_opt, opts.feature))
        for fn in all_files:
            f = pickle_load(fn)
            if opts.mz == [] or any([within_20ppm(f["mz"], mz) for mz in opts.mz]):
                our_files.append(fn)
    return sorted(set(our_files))


def process_a_pickle(i, fn, features, run_report):
    f = Feature(pickle_load(fn), opts)
    s = f"{i + 1}/{len(features)}: {f.serial} mz {f.mz:.4f} ri {int(f.ri)}"
    logging.info(s)
    print(s)
    run_report.htmls(heading(s))
    f.get_candidate_fmlas()
    f.create_library_report()
    f.perform_ms1_analysis()
    f.perform_ms2_analysis(opts)
    f_pkl = f"{opts.elucidations_dir}/e{f.serial}.pkl"
    pickle_dump(f_pkl, f.as_dict())
    reporter = FeatureReporter(f_pkl)
    summary = reporter.write_report(opts)
    if i == 0 and opts.is_worker and opts.offset == 0:  # Yeah, a kludge.
        # document software versions for this run.  Going to plunk all this in a database as part of deliver.py. #TODO
        versions = reporter.get_software_versions()
        versions.to_csv(f"{opts.out_dir}/software_status_{opts.tag}.csv", index=False)
    return summary


def process_features(features, run_report):

    # We will build a spreadsheet with a summary row for each feature.  This is a list of rows. See if a trace file exists, and if it does, load it.
    summaries = list()
    if os.path.isfile(opts.trace_csv):
        summaries = pd.read_csv(opts.trace_csv, converters={'areas': lambda x: json.loads(x.replace("'", "\""))}).to_dict(orient="records")
        existing_serials = [d['serial'] for d in summaries]
        features_to_process = [f for f in features if (os.path.splitext(os.path.basename(f))[0]) not in existing_serials]
    else:
        features_to_process = features
    failures = []

    # get dataframe of sample information.  From study module for now, later from spreadsheet.
    elapsed_time(f"About to elaborate {len(features_to_process)} features.", run_report)

    # Now go through the features and generate reports.
    for i, fn in enumerate(features_to_process):
        if ram_gb_in_use() > 2.5:
            elapsed_time(f"GC stats before:  {gc.get_stats()}", run_report)
            gc.collect()
            elapsed_time(f"GC stats after:  {gc.get_stats()}", run_report)
        s = ""
        try:
            # Print headers on terminal and in run report HTML file.
            logging.info("=" * 100)
            summary = process_a_pickle(i, fn, features_to_process, run_report)
            if summary:  # If something goes awry, nothing goes to spreadsheet.
                summaries.append(summary)
            print_brief_feature_summary(summary, run_report)

            # And periodically update the summary spreadsheet as the program runs.
            if i % 10 == 0:
                update_trace_spreadsheet(summaries, opts.trace_csv)
        except Exception as e:
            logging.error(f"{fn} failed {e}", exc_info=e, stack_info=True)
            print_exception(*sys.exc_info())
            run_report.htmls(f"<h1>{fn} failed {e}</h1>")
            failures.append(fn)
            if opts.debug:
                raise e
        elapsed_time(s, run_report)
    summaries_cucumber = format_summary_spreadsheet(summaries, opts)
    pickle_dump(f"{opts.out_dir}/feature_summaries_{opts.tag}_worker{opts.offset}.pkl", summaries_cucumber)

    if failures:
        run_report.htmls("<h1>Failures</h1>", str(failures))
        logging.error(f"Failures: {failures}")
    elapsed_time("Finished All Feature Reports", run_report)


def start_workers(opts_original, opts):
    opts_original.out_dir = '"' + opts.out_dir + '"'
    queue = []
    for i in range(opts.workers):
        opts_original.share = f"{opts.workers},{i},{opts.tag}"
        arg_list = make_argv_from_opts(opts_original)
        queue.append([sys.executable] + [__file__] + arg_list)
    # We have to pass the tag in to be sure it matches the master; otherwise the worker could be 1 min later :(
    workers = []
    start_times = {}
    while queue:
        # Wait for an opening.
        while len(workers) >= opts.workers:
            time.sleep(15)
            workers = [p for p in workers if p.poll() is None]
        # Go!
        argv, *queue = queue
        logging.info(f"Starting {' '.join(argv)}")
        worker = subprocess.Popen(argv)
        start_times[worker] = time.monotonic()
        workers.append(worker)

    while workers:
        logging.info(f"Waiting for {len(workers)} workers")
        time.sleep(15)
        # Check if any workers have been running for more than 4 hours
        for worker in workers:
            elapsed_time = time.monotonic() - start_times[worker]
            if elapsed_time >= 4 * 60 * 60:  # 4 hours in seconds
                logging.info(f"Restarting worker with pid {worker.pid}")
                worker.terminate()
                worker.wait()
                argv = worker.args
                new_worker = subprocess.Popen(argv)
                start_times[new_worker] = time.monotonic()
                workers.append(new_worker)

        # Remove workers that have completed
        for worker in workers[:]:
            if worker.poll() is not None:
                workers.remove(worker)
                del start_times[worker]


    logging.info(f"workers finished.")


def init_logging(opts):
    opts.log_level = logging.DEBUG if opts.debug else logging.INFO
    format = "%(asctime)s - %(levelname)s - %(message)s"
    if opts.is_worker:
        filename = f"{opts.out_dir}/log_{opts.tag}_worker{opts.offset}.log"
    else:
        filename = f"{opts.out_dir}/log_{opts.tag}_master.log"

    logging.basicConfig(
        level=opts.log_level,
        format=format,
        filename=filename,
        filemode="w",
    )
    console = logging.StreamHandler()
    console.setLevel(opts.log_level)
    console.setFormatter(logging.Formatter(format))
    # if opts.echo:
    #     print("echoing")
    #     logging.getLogger("").addHandler(console)
    logging.info(f"Starting elucidator for study {opts.study}, lib {opts.library} features.")
    logging.info(LibraryService.get_instance().get_study_title(opts.study))
    for k, v in vars(opts).items():
        logging.info(f"option {k}: {v}")
    mpl = logging.getLogger("matplotlib")
    mpl.setLevel(logging.WARNING)


def run_nme(opts):
    local_base_dir = BASE_DIR
    delivery_base_dir = DELIVERY_DIR
    # Save original opts for later
    opts_original = argparse.Namespace(**vars(opts))
    # Master or worker?
    opts.is_worker = opts.share != ""
    if opts.is_worker:
        opts.multiplier, opts.offset, opts.tag = opts.share.split(",")
        opts.multiplier = int(opts.multiplier)
        opts.offset = int(opts.offset)
    else:
        opts.multiplier, opts.offset, opts.tag = opts.workers, 0, ""

    # get lib_id, study_id, and set names of output directories.
    lib_id = opts.library
    study_id = opts.study
    opts.out_dir = opts.out_dir.replace('"', "")
    if opts.out_dir == "":
        opts.now = now_string()  # Build directory from the date.
        opts.out_dir = f"{local_base_dir}\\{opts.study}\\{opts.library}\\{opts.now}"
    elif study_id == "REGRESSION_TESTING":
        opts.now = now_string()
    else:
        opts.now = opts.out_dir[-15:]  # Extract the date from the directory.

    opts.trace_csv = f"{opts.out_dir}/trace_features{opts.offset}.csv"
    opts.feature_dir = f"{opts.out_dir}/features"
    opts.elucidations_dir = f"{opts.out_dir}/epkl"
    os.makedirs(opts.feature_dir, exist_ok=True)
    os.makedirs(opts.elucidations_dir, exist_ok=True)

    # Find the discovery directory that contains the inputs.
    if opts.tag == "":  # This would happen for the master only
        opts.tag = f"{opts.study}_{opts.library}_{opts.now}"
    if opts.discovery_dir == "":
        glb = glob(f"{system_base_dir}/{study_id}/{lib_id}/discoveries_*????-??-??-????/")
        if glb:
            opts.discovery_dir = max(glb)
        else:
            opts.discovery_dir = f"{system_base_dir}/{study_id}/{lib_id}/clusters/"

    #
    init_logging(opts)
    logging.info(f"logging for {opts.offset}")

    # Exit early, so we don't accidentally lock our user by messing up the login
    # 1,000 times :)
    try:
        mlims_credentials = dict(db=os.environ["DB_NAME"], user=os.environ["DB_USER"], )
        eq = EasyQuery(**mlims_credentials)
        eq._sql_alchemy_engine.connect()
    except Exception as err:
        logging.error("Was not able to setup DB connection")
        raise err

    # This is a file of annotations (e.g., statistics against metadata) that can be merged in to final spreadsheet.
    annot_file = f"{opts.discovery_dir}/stats_{study_id}_{lib_id}.csv"
    try:
        opts.annot_df = pd.read_csv(annot_file)
    except Exception as e:
        logging.info(f"Not found (or not loaded): {annot_file}")
        opts.annot_df = None

    features = get_feature_files(opts)
    opts.workers = min([opts.workers, len(features)])

    init_logging(opts)

    s = f"Starting elucidator for study {study_id}, lib {opts.library}, {len(features)} features."
    logging.info(s)
    print(s)
    print(f"output in {opts.out_dir}")
    study_title = LibraryService.get_instance().get_study_title(study_id)

    # Open the HTML run report file.
    # run_report = FileReporter(f"{opts.out_dir}/run_report_{opts.tag}_{opts.offset}.html")
    run_report = NullReporter()
    run_report.htmls(boxit(f"{study_title} -- Library {lib_id}", tag="h1"))

    lib_id = opts.library
    opts.gitbranch = get_git_branch()
    # Begin run report.
    logging.info(study_title)
    run_report.htmls(heading(study_title))
    run_report.comment(time.asctime())
    run_report.comment(opts.gitbranch)
    run_report.flush()

    if opts.is_worker:  # workers will process features and stop.
        features = features[opts.offset :: opts.multiplier]
        process_features(features, run_report)
        elapsed_time(f"worker {opts.offset} DONE", run_report)
        run_report.close_report()
        sys.exit()

    #  Master:   Manage workers, then collect results.
    start_workers(opts_original, opts)

    # Pretty up and save summary dataframe as a spreadsheet.

    summary_dfs, area_dfs = zip(
        *[
            pickle_load(f"{opts.out_dir}/feature_summaries_{opts.tag}_worker{i}.pkl")
            for i in range(opts.workers)
        ]
    )
    summary_df = pd.concat(summary_dfs, ignore_index=True).sort_values(qw("ri mz"))
    columns = qw(
        "date study lib_id serial mz ri internal_standard num_peaks num_spectra cand_fmlas max_harmony "
        + "isotopes ion_names inferred_masses lib_buddies ms1_cluster wins_ms1 "  # MS1 stuff
        + "largest_smiles_prediction fraction_predicted num_smiles_preds "  # MS2SMILES
        + "best_spectrum_prediction num_best_predictions max_spectrum_score "  # SMILES2MS2
        + "chemont_preds rediscovery_of"  # MS2 => Chemical Ontology
    )
    summary_df = summary_df[columns]

    area_df = pd.concat(area_dfs, ignore_index=True).sort_values("serial")
    area_df = area_df[["serial"] + list(area_df.columns)[:-1]]

    summary_df = perform_ms1_clustering(summary_df, run_report, opts)
    elapsed_time("Finished MS1 Clustering", run_report)

    # Write summary CSV for debugging, then write local excel, also for debugging, I guess you'd say.
    summary_df.to_csv(f"{opts.out_dir}/feature_summaries_{opts.tag}.csv", index=False)
    to_excel_with_molecules(
        summary_df,
        f"{opts.out_dir}/local_feature_summaries_{opts.tag}.xlsx",
        qw("best_spectrum_prediction largest_smiles_prediction"),
        size=(150, 100),
    )

    # Rewrite hyperlinks for excel file for delivery to the L: drive.
    summary_df.serial = [s.replace(local_base_dir, delivery_base_dir) for s in summary_df.serial]
    summary_df.ms1_cluster = [
        s.replace(local_base_dir, delivery_base_dir) for s in summary_df.ms1_cluster
    ]
    to_excel_with_molecules(
        summary_df,
        f"{opts.out_dir}/feature_summaries_{opts.tag}.xlsx",
        qw("best_spectrum_prediction largest_smiles_prediction"),
        size=(150, 100),
    )

    # Rewrite a spreadsheet of areas. This is also being done in discovery now, but we might need it if
    # rerunning an older study.
    area_df.to_csv(f"{opts.out_dir}/feature_areas_{opts.tag}.csv", index=False)
    elapsed_time("Finished Writing Spreadsheets", run_report)
    run_report.close_report()
    elapsed_time(f"Results in {opts.out_dir}", run_report)


if __name__ == "__main__":
    opts = get_arguments()
    run_nme(opts)
