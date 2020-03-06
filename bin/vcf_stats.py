#!/usr/bin/env python

import argparse

import re
import sys
import gzip

from os.path import join as pjoin

import pandas as pd

import matplotlib
matplotlib.use('agg')

from matplotlib import pyplot as plt
import seaborn as sns

def plotter(vals, field, path, logx=False):
    fig, ax = plt.subplots()
    sns.kdeplot(vals, ax=ax)
    ax.set_xlabel(field)
    ax.set_ylabel("density")

    if logx:
        ax.set_xscale('log')
    fig.savefig(path)
    plt.close()
    print("Saved plot to", path)
    return

def main(infile, prefix):

    df = pd.read_table(infile, na_values="NA")

    for field in df.columns:
        vals = df[field].dropna()

        if field == "FS":
            logx = True
        else:
            logx = False

        outpath = prefix + field + ".pdf"
        plotter(vals, field, outpath, logx)

    targets = [
        ("qual_by_depth", "QD"),
        ("fisher_strand_bias", "FS"),
        ("strand_odds_ratio", "SOR"),
        ("rms_mapping_quality", "MQ"),
        ("mapping_quality_rank_sum", "MQRankSum"),
        ("read_position_rank_sum", "ReadPosRankSum"),
        ("base_quality_rank_sum", "BaseQRankSum"),
        ]

    return

if __name__== '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "-i", "--infile",
        required=True,
        help="files"
        )
    arg_parser.add_argument(
        "-o", "--prefix",
        dest="prefix",
        default="out",
        help="prefix"
        )

    args = arg_parser.parse_args()

    main(**args.__dict__)

