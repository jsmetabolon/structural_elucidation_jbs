# -*- coding: utf-8 -*-
"""Shared helpers for NCE elucidation runs.

Pickle I/O, timing/memory reporting, git labels, plotting resets, and
database credential loaders used across elucidation modules in spec2mol
v21 inference.
"""

import datetime
import logging
import os
import pickle
import subprocess
import sys
import time
from functools import wraps

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import psutil
from dotenv import load_dotenv


def get_postgres_password():
    load_dotenv()
    return os.environ["POSTGRES_PASSWORD"]


def get_oracle_configuration():
    load_dotenv()
    return {k: v for k, v in os.environ.items() if k in qw("DB_USER DB_NAME DB_PASSWORD")}


start_time = time.time()


def reset_start_time():
    global start_time
    start_time = time.time()


def ram_gb_in_use():
    return psutil.Process(os.getpid()).memory_info().rss / 2 ** 30


def resource_string():
    return "%.2f min %.3f Gb" % ((time.time() - start_time) / 60, ram_gb_in_use())


def get_git_branch():
    try:
        label = (
            subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]).decode().strip()
        )
        label1 = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()
        return f"{label}@{label1}"
    except:
        return "Unknown git branch"


def leng(x):
    return len(x) if x is not None else 0


def describe(lis, percentiles=None):
    """
    Builds a dataframe with a statistical summary of the items in the list (iterable)
    :param lis: The list or other iterable containing the items to summarize
    :param percentiles:  A list of quantiles to include in the summary
    :return: dataframe of statistics just like in pandas.
    """
    return pd.Series(list(lis)).describe(percentiles=percentiles)


def qw(s):
    """
    Takes a list and splits it into words, like Perl qw//.  Saves a lot of keystrokes when selecting columns.
    :param s: The string to split.
    :return: A list of the words.
    """
    return s.split()


def lmap(*args):
    """
    Like map, but returns a list of results instead of a map object.
    :param args: the arguments normally passed to map, i.e., a function and an iterable.
    :return: List of results.
    """
    return list(map(*args))


def lzip(*args):
    """
    Like zip, but returns a list of results instead of a zip object.
    :param args: the arguments normally passed to zip, i.e., a function and an iterable.
    :return: List of results.
    """
    return list(zip(*args))


def lfilter(*args):
    """
    Like filter, but returns a list of results instead of a filter object.
    :param args: the arguments normally passed to filter, i.e., a function and an iterable.
    :return: List of results.
    """
    return list(filter(*args))


def now_string():
    return time.strftime("%Y-%m-%d-%H%M")


def today_string():
    """
    A string for today's date.
    :return: A string for today's date.
    """
    today = datetime.date.today()
    today = "%4d-%02d-%02d" % (today.year, today.month, today.day)
    return today


def pickle_dump(file, tuple_):
    """
    dumps a tuple to a pickle file.  This must exist somewhere else, though.
    :param file: file name
    :param tuple_: Stuff to dump
    :return:
    """
    with open(file, "wb") as pkl:
        pickle.dump(tuple_, pkl)


def pickle_load(file):
    """
    reads a pickle file and returns contents.
    :param file: file name
    :return: python object from pickle file.
    """
    with open(file, "rb") as pkl:
        return pickle.load(pkl)


def identity(x):
    """
    The identity function.  Useful with lfilter and on some other occasions.
    :param x: anything at all
    :return: the same anything
    """
    return x


# Save matplotlib default parameters.
my_params = mpl.rcParams


def plt_reset():
    """
    restore matplotlib to defaults.  Sometimes seaborn will change stuff.
    :return: None
    """
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rc("font", size=15)


def add_to_path(path):
    """
    Adds paths to sys.path.  I've learned how to not have to do this.
    :param path: path to add
    :return: None
    """
    if path not in sys.path:
        sys.path.append(path)


def is_empty_df(df):
    """
    True if a dataframe is not interesting (None or no rows)
    :param df: the dataframe (or None)
    :return: boolean
    """
    return df is None or len(df) == 0


def reorganize_columns(df, cols, drop_leftovers=False):
    """
    Reorganizes columns of a dataframe, moving some columns to the left.
    :param df: the dataframe
    :param cols:  the columns to move to the left.  If some are absent, it's OK.
    :param drop_leftovers: whether to drop or include the other columnns.
    :return: reorganized dataframe.
    """
    new_cols = [col for col in cols if col in set(df.columns)]
    if not drop_leftovers:
        new_cols += [col for col in df.columns if col not in cols]
    return df[new_cols]


def sort_df_by_vector(df, vec):
    """
    Sorts a data frame with keys from a vector.  Puts the keys in the dataframe, sorts, and removes them.
    :param df: dataframe to sort.
    :param vec: vector of keys
    :return: sorted dataframe
    """
    col = "__sort_key__"
    df[col] = vec
    df = df.sort_values(col).drop(col, axis=1)
    return df


def set_pandas_display_options(
    max_columns=1000, max_rows=1000, max_colwidth=199, width=None
) -> None:
    """
    Sets pandas options to display more than default amount of dataframes.
    :return:
    """
    # Ref: https://stackoverflow.com/a/52432757/
    display = pd.options.display

    display.max_columns = max_columns
    display.max_rows = max_rows
    display.max_colwidth = max_colwidth
    display.width = width
    # display.precision = 2  # set as needed


def start_logging(level=logging.INFO):
    print("Running start_logging")
    # logger = logging.getLogger()

    # logging.basicConfig(
    #     format="%(asctime)s %(levelname)s %(module)s %(funcName)s: %(message)s",
    #     level=level,
    #     datefmt="%Y-%m-%d %H:%M:%S %p",
    #     )
    # configDict = {
    #     'version': 1,
    #     'formatters': {
    #         'std': {
    #             'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    #             'datafmt': "%Y-%m-%d %H:%M %p",
    #            # 'class': logging.Formatter
    #         }
    #
    #     },
    #     'handlers': {
    #         'console': {
    #             'class': 'logging.StreamHandler',
    #             'formatter': 'std',
    #             'level': 'INFO',
    #         }
    #     },
    #     'loggers': {
    #         '': {
    #             'handlers': ['console'],
    #             'formatter': 'std',
    #             'propagate': True
    #         },
    #          'compound_discovery': {
    #              'level': 'INFO',
    #              'handlers': ['console'],
    #              'propogate': False
    #          },
    #     },
    #     'disable_existing_loggers': True
    # }
    # dictConfig(configDict)
    # Modifying root logger
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    logger = root.getChild("compound_discovery")
    logger.setLevel(logging.INFO)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)
    # logger = logging.getLogger('compound_discovery')
    logging.info("Logging started")
    print("logger", logger)


def funcTime(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        m0 = ram_gb_in_use()
        t0 = time.time()
        p0 = time.process_time()
        result = function(*args, **kwargs)
        t1 = time.time()
        p1 = time.process_time()
        m1 = ram_gb_in_use()
        logging.debug(
            f"{function.__name__} memory growth: {m1-m0:.3f} wall time: {t1 - t0}  process time: {p1 - p0} seconds"
        )
        # print(f"{function.__name__} memory growth: {m1 - m0:.3f} wall time: {t1 - t0}  process time: {p1 - p0} seconds")
        return result

    return function_timer
