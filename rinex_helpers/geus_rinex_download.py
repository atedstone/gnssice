# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Download GEUS Russell GNSS base station data 
#
# Data are obtained from draft GEUS Dataverse dataset. We have access to this via a Dataverse key provided especially.
#
#
# Andrew Tedstone, September 2025

# %%
import json
import requests
import os

# %%
SAVE_TO = '/work/atedstone/gnss_2024_2025/lrhp/'
headers = {
    'X-Dataverse-key':'INSERT_KEY_HERE'
}

response = requests.get('https://dataverse.geus.dk/api/datasets/:persistentId/versions/:draft?persistentId=doi:10.22008/FK2/QIGJCR', 
                       headers=headers)
listing = json.loads(response.content)


# %%
files = listing['data']['files']

for f in files:
    # Get the Dataverse file ID (which is not the same thing as the filename)
    fid = f['dataFile']['id']
    # Get the filename, which we will use to save the file locally
    fname = f['label']
    # Skip the NAV files, we only want the compressed OBS files
    if fname[-1] != 'D':
        continue
        
    print(fname)
    
    with open(os.path.join(SAVE_TO, fname), 'wb') as fh:
        response = requests.get(f'https://dataverse.geus.dk/api/access/datafile/{fid}', headers=headers)
        fh.write(response.content)
