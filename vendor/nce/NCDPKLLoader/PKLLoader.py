# ---------------------------------------------------------
#
# Name: New Compound Discovery PKL Loader
# Author: Duncan McPherson
# Description: This script will load PKL files output by the elucidator.py script into the Oracle database
#
# ---------------------------------------------------------
import sys

user = "DSRWUSER"
env_pw_key = user.upper() + "_PASSWORD"
dbs = {
    "prod": "10.0.0.38/mlims.world",
    "preprod": "10.0.0.60/mlims.world",
    "beta": "rtpbtaodb001.metabolon.com/beta.metabolon.com",
}

libs = {209: "-", 305: "-", 400: "+", 402: "+"}

import argparse
from argparse import Namespace
import numpy as np
import os
import cx_Oracle
from glob import glob
from operator import itemgetter
import pickle

password = os.environ[env_pw_key]


def run_query(statement, connect_string):
    """
    :param statement: an Oracle Query to be run
    :return: The result of the query
    """
    with cx_Oracle.connect(user, password, connect_string) as con:
        cur = con.cursor()
        cur.execute(statement)
        return cur.fetchall()


def is_part_of_load(file):
    with open(file, "rb") as f:
        data = pickle.load(f)
    is_load = not bool(data.get("rediscovery_of"))
    # ri_windows as requested in #18326
    if is_load:
        lib_id = data.get("lib_id")
        ri = data.get("ri")
        if lib_id not in (209, 305, 400, 402):
            is_load = False
        if lib_id == 209 and ri < 750:
            is_load = False
        if lib_id == 305 and ri < 700:
            is_load = False
        if lib_id == 400 and ri < 500:
            is_load = False
        if lib_id == 402 and ri < 700:
            is_load = False
    print("Checking " + file + ": " + str(is_load))
    return is_load


class BulkCompoundUpload(object):
    def __init__(self, files, connect_string):
        """
        :param files: A list object containing the names of all the PKL files to load
        """
        self.files = files
        self.connect_string = connect_string
        comp_id = run_query("select max(comp_id)+1 as comp_id from limsuser.chro_comp_entry", self.connect_string)
        comp_id = int(comp_id[0][0])
        if comp_id < 2000000:
            comp_id = 2000000

        unknown_id = run_query(
            "select name from (select * from limsuser.chro_comp_entry "
            "where name like 'NMD-%' order by comp_id desc) "
            "where rownum = 1", self.connect_string
        )
        if unknown_id == []:
            unknown_id = "2000000"
        else:
            unknown_id = unknown_id[0][0]
            unknown_id = unknown_id.replace("NMD-", "")
        if int(unknown_id) <= 2000000:
            unknown_id = 2000000
        else:
            unknown_id = int(unknown_id) + 1
        list1 = list(range(comp_id, comp_id + len(files)))
        list2 = list(range(unknown_id, unknown_id + len(files)))
        list3 = ["NMD-" + str(s) for s in list2]
        list4 = ["NMD-" + str(s) for s in list2]
        self.comp_list = list(map(lambda x, y, z: (x, y, z), list1, list3, list4))
        self.lib_list = []
        # Get last chro_lib_entry_id
        c = run_query("select max(chro_lib_entry_id)+1 as pkey from limsuser.chro_lib_entry", self.connect_string)
        self.lib_entry_id = int(c[0][0])
        if self.lib_entry_id < 1000000:
            self.lib_entry_id = 1000000
        self.masses_list = []
        # Get last masses_id, set to seq
        d = run_query("select max(chro_lib_masses_id)+1 as pkey from limsuser.chro_lib_masses", self.connect_string)
        self.masses_seq = int(d[0][0])
        if self.masses_seq < 100000000:
            self.masses_seq = 100000000

    def get_lib_list(self, num, data):
        """
        :param num: An iterator
        :return: None
        """

        # Calculate RT
        rt = np.median(list(filter(lambda x: x == x, data.peak_df["apex_rt"])))
        # Calculate Target Mass
        if libs[data.lib_id] == "-":
            target_mass = data.mz + 1.007276466
        else:
            target_mass = data.mz - 1.007276466
        # Create RI Window as requested in #18326
        if data.lib_id == 209:
            ri_window = 125 if data.ri < 5800 else 450
        if data.lib_id == 400:
            ri_window = 150 if data.ri < 3400 else 250
            if data.ri > 2400 and data.ri < 2600:
                ri_window = 165
        if data.lib_id == 402:
            ri_window = round(data.ri / 10)
        if data.lib_id == 305:
            a = [
                (0, 150),
                (800, 175),
                (900, 200),
                (1000, 250),
                (2000, 300),
                (3000, 400),
                (5000, 750),
            ]
            for x in a:
                if data.ri > x[0]:
                    ri_window = x[1]
        # Create Comments
        comments = "Study %s Feature %s Discovered at %s by %s" % (
            data.study_id,
            data.serial,
            data.discovery_time,
            data.discoverer_version,
        )
        comments = comments.replace("'", "''")
        # Create Putative Details
        pds = ""
        if data.fmla_df is not None:
            if "neutral" in list(data.fmla_df.keys()):
                if len(data.fmla_df["neutral"]) == 1:
                    pds = data.fmla_df["neutral"].iloc[0] + "?"
        if data.chemont_predictions is not None:
            if "chem_ont_name" in list(data.chemont_predictions.keys()):
                if len(data.chemont_predictions["chem_ont_name"]) > 0:
                    ismax = np.argmax(data.chemont_predictions["prob"])
                    if len(pds) > 0:
                        pds = pds + " "
                    pds = pds + data.chemont_predictions["chem_ont_name"][ismax] + "?"
        mis = []
        for i in data.buddy_df["ratio"].keys():
            if data.buddy_df["ratio"][i] > 1.01:
                mis.append(str(data.buddy_df["mz_buddy"][i]))
        if len(mis) > 0:
            if len(pds) > 0:
                pds = pds + " "
            pds = pds + "Molecular Ion may be " + ", ".join(mis)
        pds = pds.replace("'", "''")
        row = (
            self.lib_entry_id + num,
            data.lib_id,
            self.comp_list[num][0],
            self.comp_list[num][0],
            self.comp_list[num][0],
            rt,
            data.ri,
            ri_window,
            "NCD" + data.study_id,
            "Study_" + data.study_id,
            target_mass,
            libs[data.lib_id],
            comments,
            pds,
        )
        self.lib_list.append(row)

    def get_ms1(self, num, data):
        """
        :param num: An iterator
        :return: None
        """
        # Generate row
        row = (
            self.lib_list[num][0],
            self.comp_list[num][0],
            data.lib_id,
            "m" + libs[data.lib_id] + "H",
            round(data.mz, 5),
            0.001 if data.mz < 100 else round(data.mz / 100000, 5),
            100,
            "Y",
            "Y",
            100,
            libs[data.lib_id] + "i",
            self.masses_seq,
            None,
        )
        self.masses_list.append(row)
        self.ms1 = len(self.masses_list) - 1
        self.parent_id = self.masses_seq
        self.masses_seq = self.masses_seq + 1

    def get_buddies(self, num, data):
        """
        :param num: An iterator
        :return: None
        """
        # Construct Rows
        rows = []
        for i in data.buddy_df["buddy"].keys():
            if data.buddy_df["feature"][i] in ("x", "m" + libs[data.lib_id] + "H"):
                row = (
                    self.lib_list[num][0],
                    self.comp_list[num][0],
                    data.lib_id,
                    data.buddy_df["buddy"][i].replace("x", "m" + libs[data.lib_id] + "H"),
                    round(data.buddy_df["mz_buddy"][i], 5),
                    0.005,
                    round(float(data.buddy_df["ratio"][i]) * 100, 1),
                    None,
                    None,
                    round(float(data.buddy_df["ratio"][i]) * 100, 1),
                    libs[data.lib_id] + "i",
                    self.masses_seq,
                    None,
                )
                rows.append(row)
                self.masses_seq = self.masses_seq + 1
        # Remove buddies with the same mass
        old_rows = rows
        rows = sorted(rows, key=lambda x: len(x[3]), reverse=True)
        d = {t[4]: t for t in rows}
        rows = list(d.values())
        dropped = [t for t in old_rows if t not in rows]
        print("Dropped: ", dropped)
        self.masses_list.extend(rows)

    def get_ms2(self, num, data):
        """
        :param num: An iterator
        :return: None
        """
        # Construct Rows
        rows = []
        if data.consensus_df is not None:
            normalizer = max(data.consensus_df["norm_int"])
            for i in data.consensus_df["mz"].keys():
                row = [
                    self.lib_list[num][0],
                    self.comp_list[num][0],
                    data.lib_id,
                    None,
                    round(data.consensus_df["mz"][i], 5),
                    0.005,
                    round(data.consensus_df["norm_int"][i] / normalizer * 100, 1),
                    None,
                    None,
                    round(data.consensus_df["norm_int"][i] / normalizer * 100, 1),
                    libs[data.lib_id] + "i",
                    self.masses_seq,
                    self.parent_id,
                ]
                if row[6] == 100:
                    row[7] = "Y"
                    row[8] = "Y"
                row = tuple(row)
                if row[5] > 0:
                    rows.append(row)
                    self.masses_seq = self.masses_seq + 1
            # Purge Dupes
            old_rows = rows
            rows = sorted(rows, key=itemgetter(5))
            d = {(str(t[1]) + str(t[4])): t for t in rows}
            rows = list(d.values())
            dropped = [t for t in old_rows if t not in rows]
            print("Dropped: ", dropped)
            self.masses_list.extend(rows)

    def purge_mass_dupes(self):
        i = self.ms1
        while i < len(self.masses_list) - 1:
            j = i + 1
            while j < len(self.masses_list):
                if self.masses_list[i][4] == self.masses_list[j][4]:
                    print("Dropped "),
                    print(self.masses_list.pop(j))
                    j -= 1
                j += 1
            i += 1

    def comp_execute(self, cur):
        # Execute Comp Query
        query = (
            "insert into limsuser.chro_comp_entry (comp_id, name, report_name, type) "
            "values (:1, :2, :3, 'UNKNOWN')"
        )
        cur.bindarraysize = len(self.comp_list)
        cur.setinputsizes(int, 20, 20)
        cur.executemany(query, self.comp_list)

    def lib_execute(self, cur):
        # Execute Lib Query
        query = (
            "insert into limsuser.chro_lib_entry ("
            "chro_lib_entry_id,"
            "lib_id,"
            "comp_id,"
            "comp_id_set,"
            "comp_id_alt,"
            "version,"
            "version_current,"
            "rt,"
            "rt_window,"
            "ri,"
            "ri_window,"
            "exp_function,"
            "status,"
            "confidence,"
            "origin,"
            "group_name,"
            "display_color,"
            "target_mass,"
            "pos_neg_factor,"
            "quant_stream,"
            "comments,"
            "putative_details,"
            "created_on,"
            "modified_on)"
            "values("
            ":1,"
            ":2,"
            ":3,"
            ":4,"
            ":5,"
            "1,"
            "'T',"
            ":6,"
            "0.05,"
            ":7,"
            ":8,"
            "'Sum',"
            "'silver',"
            "100,"
            ":9,"
            ":10,"
            "'DarkViolet',"
            ":11,"
            "1,"
            ":12,"
            ":13,"
            ":14,"
            "sysdate,"
            "sysdate)"
        )
        # row = (self.lib_entry_id+num, lib_id, comp_id, comp_id, comp_id, rt, data.ri, data.study_id, target_mass, comments, pds)
        cur.bindarraysize = len(self.lib_list)
        cur.setinputsizes(
            int, int, int, int, int, float, float, float, 14, 17, float, 1, 4000, 4000
        )
        cur.executemany(query, self.lib_list)

    def mass_execute(self, cur):
        # Execute Masses Query
        query = (
            "insert into limsuser.chro_lib_masses ("
            "chro_lib_entry_id,"
            "comp_id,"
            "lib_id,"
            "name,"
            "mass,"
            "mass_window,"
            "mass_ratio,"
            "mass_normalizer,"
            "quant_mass,"
            "weighting,"
            "data_stream,"
            "timestamp,"
            "chro_lib_masses_id,"
            "parent_id)"
            "values ("
            ":1,"
            ":2,"
            ":3,"
            ":4,"
            ":5,"
            ":6,"
            ":7,"
            ":8,"
            ":9,"
            ":10,"
            ":11,"
            "sysdate,"
            ":12,"
            ":13)"
        )
        # row = (entry_id, comp_id, lib_id, name, mass, mass_ratio, mass_normalizer, quant_mass, weighting, mass_id, parent_id)
        cur.bindarraysize = len(self.masses_list)
        cur.setinputsizes(int, int, int, 40, float, float, float, 4, 4, float, 2, int, int)
        cur.executemany(query, self.masses_list)


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Loads PKL files output by elucidator.py into the mLIMS Oracle DB"
    )
    parser.add_argument("filepath", help="PKL file or Directory of PKL files to load to database")
    parser.add_argument(
        "--db", help="Database to write features to", choices=["prod", "preprod", "beta"]
    )
    return parser.parse_args()


def load_pkls(opts):
    file_path = opts.filepath
    if opts.db:
        connect_string = dbs[opts.db]
    else:
        connect_string = "rtpbtaodb001.metabolon.com/beta.metabolon.com"

    if os.path.isdir(file_path):
        print("We are using the directory!")
        files = glob(f"{file_path}/*.pkl")
    elif os.path.isfile(file_path) and ".pkl" in file_path:
        print("We are using the PKL file!")
        files = [file_path]
    else:
        print("Input filepath is not valid!")
        return

    # Adding a filter for rediscoveries here as in the creation of the BulkCompoundUpload, a comp_id is created for every tuple in the files list
    load_files = [x for x in files if is_part_of_load(x)]
    if not load_files:
        print("No files present that are to be loaded!")
        return

    thisrun = BulkCompoundUpload(load_files, connect_string)

    for x, file in enumerate(load_files):
        with open(file, "rb") as f:
            data = pickle.load(f)
            data = Namespace(**data)

        print("Starting " + file)

        thisrun.get_lib_list(x, data)
        thisrun.get_ms1(x, data)
        thisrun.get_buddies(x, data)
        thisrun.get_ms2(x, data)
        thisrun.purge_mass_dupes()

    # Open connection and execute
    con = cx_Oracle.connect(user, password, connect_string)
    cur = con.cursor()
    thisrun.comp_execute(cur)
    thisrun.lib_execute(cur)
    thisrun.mass_execute(cur)

    con.commit()
    cur.close()
    con.close()


# ---------------------------------------------------------
#
# Main Loop
#
# ---------------------------------------------------------

if __name__ == "__main__":

    opts = get_arguments()
    load_pkls(opts)
