#!/usr/bin/env python
# coding: utf-8

import argparse
import git
import pandas as pd
import os
from urllib.request import urlretrieve

TRACKED_FILE = "data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv"
CASES_URL = "https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/commitid/data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv"
SEQUENCE_URL = "https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/commitid/data/omicron-countries-split/omicron-countries-split_location-variant-sequence-counts.tsv"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Getting versioned data files from rt-from-frequency-dynamics"
    )
    parser.add_argument(
        "--repo_location", default="../rt-from-frequency-dynamics"
    )
    parser.add_argument("--export_location", default="../data")
    args = parser.parse_args()

    # Obtaining commits for case counts
    repo_location = args.repo_location
    g = git.Git(repo_location)
    commits = list(
        g.log(
            "--follow",
            "--pretty=%H",
            "--",
            TRACKED_FILE,
        ).split("\n")
    )
    dates = list(
        g.log(
            "--follow",
            "--pretty=%as",
            "--",
            TRACKED_FILE,
        ).split("\n")
    )

    # Loop over dates and commits
    for date, commit in zip(dates, commits):
        # Make folder for this date
        path = args.export_location + date
        if not os.path.exists(path):
            os.mkdir(path)

        # get case, seq url
        commit_url_case = CASES_URL.replace("commitid", commit)
        commit_url_seq = SEQUENCE_URL.replace("commitid", commit)

        # save both in folder
        # cases_df = pd.read_csv(commit_url_case, sep="\t")
        # cases_df.to_csv(
        #     f"{path}/case-counts_{date}.tsv", index=False, sep="\t"
        # )

        # seq_df = pd.read_csv(commit_url_seq, sep="\t")
        # seq_df.to_csv(f"{path}/seq-counts_{date}.tsv", index=False, sep="\t")
        urlretrieve(commit_url_case, f"{path}/case-counts_{date}.tsv")
        urlretrieve(commit_url_seq, f"{path}/seq-counts_{date}.tsv")
