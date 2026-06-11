"""
This program is to fetch the most recent reference files from the L: drive.
"""
import argparse
import os
import shutil
import sys
from glob import glob

# default; override on command line
_local_path = "c:/cdtemp/references"
_archive_path = "l:/New Compound Discovery/references"


def pickle_updater(file, local_path, archive_path, force):
    """
    Compares dates of file on two paths, and makes sure the most recent is in both.
    Note that no files are overwritten or deleted; the new one is just put in both directories.
    @param file: prefix of file to synchronized, e.g. ms2smiles_models_209
    @type file: str
    @param local_path: a path; not necessarily anywhere in particular.
    @type local_path: str
    @param archive_path: a second path; not necessarily anywhere in particular.
    @type archive_path: str
    @param force: boolean: if False, ask before syncing.
    @type force: bool
    @return: None
    @rtype: NoneType
    """
    local_filepath = max(glob(f"{os.path.join(local_path, file)}*.pkl") + ["<none>"])
    archive_filepath = max(glob(f"{os.path.join(archive_path, file)}*.pkl") + ["<none>"])
    local_file = local_filepath.split("\\")[-1]
    archive_file = archive_filepath.split("\\")[-1]
    print()
    print(file, local_file, archive_file)
    if local_file == archive_file:
        print(f"{local_path} == {archive_path}")
    elif local_file < archive_file or local_file == "<none>":
        print(f"Copying {archive_filepath} to supersede {local_filepath}")
        if force or input("Proceed? [y/n] ").lower()[0] == "y":
            shutil.copy(archive_filepath, local_path)
            print(f"Copied {archive_filepath} to {local_path}")
    elif local_file > archive_file or archive_file == "<none>":
        print(f"Copying {local_filepath} to supersede {archive_filepath}")
        if force or input("Proceed? [y/n] ").lower()[0] == "y":
            shutil.copy(local_filepath, archive_path)
            print(f"Copied {local_filepath} to {archive_path}")


def get_arguments(argv):
    """
    Use argparse get arguments from command line.
    :param argv: strings from command line.
    :return: a Namespace object with option values.
    """
    parser = argparse.ArgumentParser(description="""Synchronize model and table files.""")
    parser.add_argument(
        "-f", "--force", action="store_true", help="Copy files without asking.",
    )
    parser.add_argument(
        "-l", "--local", type=str, help="local (or any) path", default=_local_path,
    )
    parser.add_argument(
        "-a", "--archive", type=str, help="archive (or any) path", default=_archive_path,
    )

    return parser.parse_args(argv[1:])


if __name__ == "__main__":
    opts = get_arguments(sys.argv)

    for pat in [
        "fmla_finder_cache",
        "harmony",
        "ff_dbe",
        "library_fragments",
        "fragment_masses",
    ]:
        pickle_updater(pat, opts.local, opts.archive, opts.force)

    for pat in ["ms2smiles_models_%d", "smiles2ms2_models_%d", "chemont_models_%d"]:
        for lib_id in (209, 400, 402):
            pickle_updater(pat % lib_id, opts.local, opts.archive, opts.force)

#        chemont_definition_loader = csv_loader("chemont_defs")
