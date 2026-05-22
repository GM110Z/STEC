#!/usr/bin/env python3

"""
Retrieve metadata for GCA accessions from NCBI Assembly/BioSample

"""

import sys
import time
import requests
import pandas as pd
from xml.etree import ElementTree as ET



EMAIL = "your@email.com"

API_KEY = None

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"




def ncbi_get(url, params=None, retries=5):

    if params is None:
        params = {}

    params["email"] = EMAIL

    if API_KEY:
        params["api_key"] = API_KEY

    for attempt in range(retries):

        try:

            r = requests.get(url, params=params, timeout=30)

            # success
            if r.status_code == 200:
                return r.text

            # rate limit
            elif r.status_code == 429:

                wait = 2 ** attempt

                print(f"429 received. Waiting {wait}s...")
                time.sleep(wait)

            else:
                print(f"HTTP error {r.status_code}")
                time.sleep(2)

        except requests.exceptions.RequestException as e:

            print(f"Request failed: {e}")
            time.sleep(2)

    raise RuntimeError("Failed after retries")




def assembly_to_biosample(gca):

    search_url = BASE + "esearch.fcgi"

    txt = ncbi_get(
        search_url,
        {
            "db": "assembly",
            "term": gca,
            "retmode": "json"
        }
    )

    uid_list = requests.models.complexjson.loads(txt)["esearchresult"]["idlist"]

    if not uid_list:
        return None

    uid = uid_list[0]

    summary_url = BASE + "esummary.fcgi"

    txt = ncbi_get(
        summary_url,
        {
            "db": "assembly",
            "id": uid,
            "retmode": "json"
        }
    )

    data = requests.models.complexjson.loads(txt)

    biosample = data["result"][uid].get("biosampleaccn")

    return biosample




def fetch_biosample_metadata(biosample):

    fetch_url = BASE + "efetch.fcgi"

    xml_text = ncbi_get(
        fetch_url,
        {
            "db": "biosample",
            "id": biosample,
            "retmode": "xml"
        }
    )

    root = ET.fromstring(xml_text)

    attrs = {}

    for attr in root.findall(".//Attribute"):

        name = attr.attrib.get("attribute_name", "").lower()
        value = attr.text

        if value:
            attrs[name] = value

    fields = {

        "geo_location": (
            attrs.get("geo_loc_name")
            or attrs.get("geographic location")
            or attrs.get("geographic location (country and/or sea)")
        ),

        "isolation_source": (
            attrs.get("isolation_source")
        ),

        "host": (
            attrs.get("host")
        ),

        "disease": (
            attrs.get("disease")
            or attrs.get("host_disease")
        ),

        "collection_date": (
            attrs.get("collection_date")
        ),

        "strain": (
            attrs.get("strain")
        )
    }

    return fields


if len(sys.argv) != 3:

    print("Usage:")
    print("python fetch_metadata.py accessions.txt output.tsv")
    sys.exit()

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as f:
    gcas = [x.strip() for x in f if x.strip()]

results = []

for i, gca in enumerate(gcas, start=1):

    print(f"[{i}/{len(gcas)}] Processing {gca}")

    try:

        biosample = assembly_to_biosample(gca)

        if biosample is None:

            results.append({
                "GCA": gca,
                "BioSample": None
            })

            continue

        metadata = fetch_biosample_metadata(biosample)

        row = {
            "GCA": gca,
            "BioSample": biosample,
            **metadata
        }

        results.append(row)

        time.sleep(0.5)

    except Exception as e:

        print(f"Failed for {gca}: {e}")

        results.append({
            "GCA": gca,
            "error": str(e)
        })

df = pd.DataFrame(results)

df.to_csv(output_file, sep="\t", index=False)

print(f"\nSaved results to: {output_file}")
