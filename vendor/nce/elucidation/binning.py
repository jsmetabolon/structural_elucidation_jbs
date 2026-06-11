# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 12:17:24 2019

@author: rdwyer
"""

from collections import defaultdict

import numpy as np

ten_ppm_factor = 1 - 10 / 1e6
log_ten_ppm_factor = np.log(ten_ppm_factor)


def mass_interval(m):  # if m is negative, this is a loss.
    if abs(m) > 100:
        return abs(m * 5e-6)
    else:
        return 100 * 5e-6


def plus_5ppm(m):
    """
    :param m: a mass
    :return: the mass + 5ppm
    """
    return (1 + 5e-6) * m


def minus_5ppm(m):
    """
    :param m: a mass
    :return: the mass - 5ppm
    """
    return (1 - 5e-6) * m


def within_5ppm(m1, m2):
    """
    :param m1: a mass
    :param m2: another mass!
    :return: True iff m1 and m2 are within 5ppm.
    """
    # ppm with respect to midpoint; this makes the function symmetric.
    return np.fabs(2 * (m1 - m2) / (m1 + m2)) <= 5e-6


def within_10ppm(m1, m2):
    """
    :param m1: a mass
    :param m2: another mass!
    :return: True iff m1 and m2 are within 5ppm.
    """
    # ppm with respect to midpoint; this makes the function symmetric.
    return np.fabs(2 * (m1 - m2) / (m1 + m2)) <= 10e-6


def within_20ppm(m1, m2):  # Used for debugging feature detection.
    """
    :param m1: a mass
    :param m2: another mass!
    :return: True iff m1 and m2 are within 5ppm.
    """
    # ppm with respect to midpoint; this makes the function symmetric.
    return np.fabs(2 * (m1 - m2) / (m1 + m2)) <= 20e-6


# A bin is 10 ppm wide
def bin_for_mass(mass):
    """
    Finds the bin into which a mass falls.
    :param mass: a mass!
    :return: bin number, 0 <= bin < 300000
    """
    # if np.isnan(mass):
    # return np.nan
    return 299999 + np.int32(  # This puts 0 little below 50
        np.log(1000.0 / mass) / log_ten_ppm_factor
    )


def bin_base(b):
    """
    Returns the mass at the low end of a bin.  Satisfies bin_for_mass(bin_base(bin)) == bin
    :param b: a bin number
    :return: a mass
    """
    bb = 300000 - b
    return 1000 * ten_ppm_factor ** bb


# In: list of bins.  Goal is to combine adjacent bins into buckets.  Out:map bin# to bucket#.
def create_bucket_map(bin_list, ignore_anomalies=False):
    """
    Creates a dictionary mapping bin numbers to bucket numbers by merging adjacent bins. (So to speak...)
    :param bin_list: A list of bin numbers of non-empty bins to be mapped into buckets
    :return: A dictionary mapping adjacent bins to the same bucket number, and a list of "anomalies".
    """

    # Map bin to interesting bucket.
    # this will map adjacent bin numbers to the smaller value, the index of a bucket that contains it.
    # This won't be correct if three buckets in a row are occupied, but this should be rare, I think.
    bin_list = sorted(list(bin_list))  # in case we are passed keys, array, etc.

    # Form runs of consecutive bin numbers.
    groups = []
    group = [bin_list[0]]
    for b in bin_list[1:]:
        if b != group[-1] + 1:  # gap > 1; start a new group.
            groups.append(group)
            group = [b]
        else:
            group.append(b)
    groups.append(group)
    # We really hope all groups are 1 or 2 bins.  If groups are bigger, then masses >10ppm apart are being grouped
    # together.  All we do at the moment is to go ahead and group, and then report the long run as an "anomaly" to advise
    # caution in interpretation.
    # TODO fix this.  Or fix the (two) places where the results are used.
    bin_anomalies = [g for g in groups if len(g) > 2]
    groups_to_bin = groups if ignore_anomalies else [g for g in groups if len(g) <= 2]
    bin_map = defaultdict(int)
    for g in groups:
        bucket = max(g)
        bin_map.update({b: bucket for b in g})
    return bin_map, bin_anomalies
