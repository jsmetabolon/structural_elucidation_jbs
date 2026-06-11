from dotenv import load_dotenv
load_dotenv()
import discovery.discoverer as nmd
import elucidation.elucidator as nme
import NCDPKLLoader.PKLLoader as loader
from argparse import Namespace, ArgumentParser
import sys


def make_opts_for_nmd(ssid, libid):
    opts = Namespace(
        cache=False,
        debug=False,
        echo=False,
        ignore_missing=False,
        library=libid,
        max_samples=1000000000.0,
        no_rediscovery=False,
        out_dir="L:/New Compound Discovery/",
        path="L:/New Compound Discovery/{}/thermo".format(ssid),
        sample_type="BMSRC",
        study=ssid.split('+'),
        verbosity=None,
    )
    return opts


def make_opts_for_nme(ssid, libid):
    opts = Namespace(
        debug=False,
        discovery_dir="",
        feature=[],
        feature_file="",
        library=libid,
        min_ms2_intensity=0.01,
        mz=[],
        out_dir="",
        ri_hi=1000000,
        ri_lo=0,
        share="",
        workers=7,
        study='+'.join(sorted(ssid.split('+'))),
        verbosity=None,
    )
    return opts


def make_opts_for_loader(elucidations_dir):
    opts = Namespace(db="prod", filepath=elucidations_dir)
    return opts


def run_pipeline(ssid):
    libs = [209, 305, 400, 402]
    for libid in libs:
        opts_nmd = make_opts_for_nmd(ssid, libid)
        nmd.main(opts_nmd)
        opts_nme = make_opts_for_nme(ssid, libid)
        nme.run_nme(opts_nme)
        count_query = (
            f"select count(*) from chro_lib_entry where lib_id = {libid} and origin = 'NCD{ssid}'"
        )
        lib_and_origin_count = loader.run_query(count_query, loader.dbs["prod"])
        if lib_and_origin_count[0][0] == 0:
            opts_load = make_opts_for_loader(opts_nme.elucidations_dir)
            loader.load_pkls(opts_load)


def get_arguments(argv):
    """
    Use argparse get arguments from command line.
    :param argv: strings from command line.
    :return: a Namespace object with option values.
    """
    parser = ArgumentParser(
        description="""Pipeline runner for discovery, elucidation, and loading."""
    )
    parser.add_argument("-s", "--study", "--ssid", help="sample set id, if you are combining multiple sets, separate by a plus sign, e.g. 11111+22222+33333", type=str)
    return parser.parse_args(argv[1:])


if __name__ == "__main__":
    opts = get_arguments(sys.argv)
    run_pipeline(opts.study)
