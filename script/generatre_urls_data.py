#!/usr/bin/env python
# coding: utf-8

# In[10]:


import fileinput
import requests
import git
from git import repo
import pandas as pd
import os
import argparse


# In[11]:


BASE_URL_CASE = "https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/commitid/data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv"
BASE_URL_SEQ = "https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/commitid/data/omicron-countries-split/omicron-countries-split_location-variant-sequence-counts.tsv"

parser = argparse.ArgumentParser()
parser.add_argument("--location", type=str, help="Repo-location")
parser.add_argument("--output", type =str, help="output-location")
args = parser.parse_args()

####obtaining commits for case counts and seqs
g = git.Git(args.location)
commits = list(g.log("--follow","--pretty=%H", "--", "data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv").split('\n'))
dates = list(g.log("--follow","--pretty=%as", "--", "data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv").split('\n'))

#obtaining datasets and date_stamping output files
for date, commit in zip(dates, commits):
    #making folder
    path = args.output + date
    if not os.path.exists(path):
        os.mkdir(path)

    # get case url
    commit_url_case = BASE_URL_CASE.replace('commitid',commit)
    # get seq url
    commit_url_seq = BASE_URL_SEQ.replace('commitid',commit)
    
    # save both in folder
    case_df = pd.read_csv(commit_url_case, sep="\t")
    case_df.to_csv(f"{path}/case-counts_{date}.tsv", index=False, sep="\t")
    seq_df = pd.read_csv(commit_url_seq, sep="\t")
    seq_df.to_csv(f"{path}/seq-counts_{date}.tsv", index=False, sep="\t")
    

