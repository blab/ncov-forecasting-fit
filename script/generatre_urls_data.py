#!/usr/bin/env python
# coding: utf-8

# In[34]:


import fileinput
import requests
import git
from git import repo
import pandas as pd
import os


# In[ ]:


####obtaining commits for case counts

g = git.Git("../rt-from-frequency-dynamics")
full_commit = list(g.log("--follow","--pretty=short", "--", "data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv").split('\n'))
commits = list(g.log("--follow","--pretty=%H", "--", "data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv").split('\n'))
dates = list(g.log("--follow","--pretty=%as", "--", "data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv").split('\n'))


# In[ ]:


BASE_URL_case = "https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/commitid/data/omicron-countries-split/omicron-countries-split_location-case-counts.tsv"
BASE_URL_seq = "https://raw.githubusercontent.com/blab/rt-from-frequency-dynamics/commitid/data/omicron-countries-split/omicron-countries-split_location-variant-sequence-counts.tsv"


# In[50]:


for date, commit in zip(dates, commits):
    # Make folder for this date
    path = "data/" + date
    if not os.path.exists(path):
        os.mkdir(path)
    # get case url
    commit_url_case = BASE_URL_case.replace('commitid',commit)
    # get seq url
    commit_url_seq = BASE_URL_seq.replace('commitid',commit)
    
    # save both in folder
    this_df = pd.read_csv(commit_url_case, sep="\t")
    this_df.to_csv(f"./data/{date}/case-counts_{date}.tsv", index=False, sep="\t")
    other_df = pd.read_csv(commit_url_seq, sep="\t")
    other_df.to_csv(f"./data/{date}/seq-counts_{date}.tsv", index=False, sep="\t")
    
    
    

