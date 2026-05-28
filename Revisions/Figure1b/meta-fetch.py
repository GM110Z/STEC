#!/usr/bin/env python3

import argparse
import json
import time
import requests
import pandas as pd
from xml.etree import ElementTree as ET


BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

UNKNOWN_VALUES = {
    "", " ", "missing", "unknown", "not applicable", "not_applicable",
    "not applicable.", "n/a", "na", "none", "not provided",
    "not provided.", "not collected", "not collected.",
    "not determined", "not determined.", "nd", "null", "nan", "-", "--"
}


def clean_value(x):
    if x is None:
        return "Unknown"
    x = str(x).strip()
    if x.lower() in UNKNOWN_VALUES:
        return "Unknown"
    return x


def ncbi_get(endpoint, params, email, api_key=None, retries=8, sleep=0.8):
    params = dict(params)
    params["email"] = email

    if api_key:
        params["api_key"] = api_key

    url = BASE + endpoint

    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=90)

            if r.status_code == 200:
                time.sleep(sleep)
                return r.text

            if r.status_code == 429:
                wait = max(1.0, sleep) * (2 ** attempt)
                print(f"429 rate limit. Waiting {wait:.1f}s...")
                time.sleep(wait)
                continue

            print(f"HTTP {r.status_code}: {r.text[:300]}")
            time.sleep(max(1.0, sleep) * (attempt + 1))

        except requests.exceptions.RequestException as e:
            print(f"Request failed: {e}")
            time.sleep(max(1.0, sleep) * (attempt + 1))

    raise RuntimeError(f"NCBI request failed after retries: {endpoint}")


def esearch_assembly_uid(gca, email, api_key=None, sleep=0.8):
    """
    Try several Assembly searches because NCBI field matching can be annoying.
    """

    terms = [
        f"{gca}[Assembly Accession]",
        f'"{gca}"[Assembly Accession]',
        f"{gca}",
        f'"{gca}"'
    ]

    for term in terms:
        txt = ncbi_get(
            "esearch.fcgi",
            {
                "db": "assembly",
                "term": term,
                "retmode": "json",
                "retmax": 5
            },
            email=email,
            api_key=api_key,
            sleep=sleep
        )

        data = json.loads(txt)
        ids = data.get("esearchresult", {}).get("idlist", [])

        if ids:
            return ids[0], term

    return "Unknown", "No assembly search term matched"


def assembly_summary(assembly_uid, email, api_key=None, sleep=0.8):
    txt = ncbi_get(
        "esummary.fcgi",
        {
            "db": "assembly",
            "id": assembly_uid,
            "retmode": "json"
        },
        email=email,
        api_key=api_key,
        sleep=sleep
    )

    data = json.loads(txt)
    rec = data.get("result", {}).get(str(assembly_uid), {})

    assembly_accession = clean_value(
        rec.get("assemblyaccession")
        or rec.get("assemblyAccession")
        or rec.get("synonym", {}).get("genbank")
    )

    biosample = clean_value(
        rec.get("biosampleaccn")
        or rec.get("biosample")
        or rec.get("biosampleid")
    )

    return assembly_accession, biosample


def assembly_elink_to_biosample_uid(assembly_uid, email, api_key=None, sleep=0.8):
    """
    Fallback: link Assembly UID directly to BioSample UID.
    """

    txt = ncbi_get(
        "elink.fcgi",
        {
            "dbfrom": "assembly",
            "db": "biosample",
            "id": assembly_uid,
            "retmode": "json"
        },
        email=email,
        api_key=api_key,
        sleep=sleep
    )

    data = json.loads(txt)

    linksets = data.get("linksets", [])
    for linkset in linksets:
        for db in linkset.get("linksetdbs", []):
            links = db.get("links", [])
            if links:
                return str(links[0])

    return "Unknown"


def biosample_acc_to_uid(biosample_acc, email, api_key=None, sleep=0.8):
    biosample_acc = clean_value(biosample_acc)

    if biosample_acc == "Unknown":
        return "Unknown"

    terms = [
        f"{biosample_acc}[Accession]",
        f'"{biosample_acc}"[Accession]',
        biosample_acc,
        f'"{biosample_acc}"'
    ]

    for term in terms:
        txt = ncbi_get(
            "esearch.fcgi",
            {
                "db": "biosample",
                "term": term,
                "retmode": "json",
                "retmax": 5
            },
            email=email,
            api_key=api_key,
            sleep=sleep
        )

        data = json.loads(txt)
        ids = data.get("esearchresult", {}).get("idlist", [])

        if ids:
            return ids[0]

    return "Unknown"


def normalise_key(x):
    return (
        str(x)
        .lower()
        .strip()
        .replace("_", "")
        .replace("-", "")
        .replace(" ", "")
        .replace("/", "")
        .replace("(", "")
        .replace(")", "")
        .replace("[", "")
        .replace("]", "")
    )


def get_attr(attrs, names):
    # exact lowercase
    for name in names:
        key = name.lower().strip()
        if key in attrs:
            return clean_value(attrs[key])

    # normalised
    norm_attrs = {normalise_key(k): v for k, v in attrs.items()}

    for name in names:
        key = normalise_key(name)
        if key in norm_attrs:
            return clean_value(norm_attrs[key])

    return "Unknown"


def parse_biosample_xml(xml_text):
    root = ET.fromstring(xml_text)

    bs = root.find(".//BioSample")
    if bs is None:
        raise ValueError("No BioSample element found in XML")

    biosample_acc = clean_value(bs.attrib.get("accession"))
    biosample_uid = clean_value(bs.attrib.get("id"))

    attrs = {}

    for attr in bs.findall(".//Attribute"):
        value = clean_value(attr.text)

        # Important: SAMEA/ENA records may use different attribute labels.
        for key in ["attribute_name", "harmonized_name", "display_name", "name"]:
            k = attr.attrib.get(key)
            if k:
                attrs[k.lower().strip()] = value

    ids = {}
    for id_elem in bs.findall(".//Id"):
        value = clean_value(id_elem.text)

        for key in ["db", "db_label", "namespace"]:
            k = id_elem.attrib.get(key)
            if k:
                ids[k.lower().strip()] = value

    title = "Unknown"
    title_elem = bs.find(".//Description/Title")
    if title_elem is not None:
        title = clean_value(title_elem.text)

    organism = "Unknown"
    org_elem = bs.find(".//Organism")
    if org_elem is not None:
        organism = clean_value(
            org_elem.attrib.get("taxonomy_name")
            or org_elem.attrib.get("organism_name")
        )

    meta = {
        "BioSample_refetched": biosample_acc,
        "BioSample_UID": biosample_uid,

        "geo_location": get_attr(attrs, [
            "geo_loc_name",
            "geographic location",
            "geographic location (country and/or sea)",
            "geographic location (country)",
            "country"
        ]),

        "isolation_source": get_attr(attrs, [
            "isolation_source",
            "isolation source",
            "sample source",
            "source"
        ]),

        "host": get_attr(attrs, [
            "host",
            "host scientific name"
        ]),

        "disease": get_attr(attrs, [
            "disease",
            "host_disease",
            "host disease",
            "disease state"
        ]),

        "host_health_state": get_attr(attrs, [
            "host health state",
            "host_health_state",
            "health state"
        ]),

        "collection_date": get_attr(attrs, [
            "collection_date",
            "collection date",
            "collection-date"
        ]),

        "strain": get_attr(attrs, [
            "strain",
            "isolate",
            "sample name",
            "sample_name"
        ]),

        "sample_name": get_attr(attrs, [
            "sample name",
            "sample_name"
        ]),

        "submitter_id": clean_value(
            ids.get("submitter id")
            or ids.get("external id")
            or ids.get("sample name")
            or get_attr(attrs, ["submitter id", "submitter_id"])
        ),

        "organism": organism,
        "biosample_title": title,

        # This is useful for debugging if a row still looks wrong.
        "all_biosample_attributes": "; ".join(
            f"{k}={v}" for k, v in sorted(attrs.items())
        )
    }

    return meta


def fetch_biosample_metadata(biosample_uid, email, api_key=None, sleep=0.8):
    xml_text = ncbi_get(
        "efetch.fcgi",
        {
            "db": "biosample",
            "id": biosample_uid,
            "retmode": "xml"
        },
        email=email,
        api_key=api_key,
        sleep=sleep
    )

    return parse_biosample_xml(xml_text)


def make_empty_row(gca):
    return {
        "GCA": gca,
        "Assembly_UID": "Unknown",
        "Assembly_search_term": "Unknown",
        "Assembly_accession": "Unknown",
        "BioSample": "Unknown",
        "BioSample_refetched": "Unknown",
        "BioSample_UID": "Unknown",
        "geo_location": "Unknown",
        "isolation_source": "Unknown",
        "host": "Unknown",
        "disease": "Unknown",
        "host_health_state": "Unknown",
        "collection_date": "Unknown",
        "strain": "Unknown",
        "sample_name": "Unknown",
        "submitter_id": "Unknown",
        "organism": "Unknown",
        "biosample_title": "Unknown",
        "all_biosample_attributes": "",
        "error": ""
    }


def main():
    parser = argparse.ArgumentParser(
        description="Fetch metadata from GCA accessions via NCBI Assembly and BioSample."
    )

    parser.add_argument("accessions", help="Text file with one GCA accession per line")
    parser.add_argument("output", help="Output TSV")
    parser.add_argument("--email", required=True, help="Email for NCBI E-utilities")
    parser.add_argument("--api_key", default=None, help="Optional NCBI API key")
    parser.add_argument("--sleep", type=float, default=0.8, help="Sleep between requests")
    parser.add_argument("--resume", action="store_true", help="Resume from existing output file")

    args = parser.parse_args()

    with open(args.accessions) as f:
        gcas = [x.strip() for x in f if x.strip()]

    # preserve order but remove duplicates
    seen = set()
    gcas_unique = []
    for g in gcas:
        if g not in seen:
            seen.add(g)
            gcas_unique.append(g)

    rows = []
    done = set()

    if args.resume:
        try:
            old = pd.read_csv(args.output, sep="\t", dtype=str).fillna("")
            rows = old.to_dict(orient="records")
            done = set(old["GCA"].astype(str).str.strip())
            print(f"Resuming. Already done: {len(done)}")
        except FileNotFoundError:
            pass

    for i, gca in enumerate(gcas_unique, start=1):
        if gca in done:
            continue

        print(f"[{i}/{len(gcas_unique)}] Processing {gca}")

        row = make_empty_row(gca)

        try:
            assembly_uid, search_term = esearch_assembly_uid(
                gca,
                email=args.email,
                api_key=args.api_key,
                sleep=args.sleep
            )

            row["Assembly_UID"] = clean_value(assembly_uid)
            row["Assembly_search_term"] = clean_value(search_term)

            if assembly_uid == "Unknown":
                row["error"] = "Assembly accession not resolved"
                rows.append(row)
                pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
                continue

            assembly_acc, biosample_acc = assembly_summary(
                assembly_uid,
                email=args.email,
                api_key=args.api_key,
                sleep=args.sleep
            )

            row["Assembly_accession"] = assembly_acc
            row["BioSample"] = biosample_acc

            biosample_uid = "Unknown"

            if biosample_acc != "Unknown":
                biosample_uid = biosample_acc_to_uid(
                    biosample_acc,
                    email=args.email,
                    api_key=args.api_key,
                    sleep=args.sleep
                )

            # Fallback if BioSample accession is missing or cannot be resolved
            if biosample_uid == "Unknown":
                biosample_uid = assembly_elink_to_biosample_uid(
                    assembly_uid,
                    email=args.email,
                    api_key=args.api_key,
                    sleep=args.sleep
                )

            row["BioSample_UID"] = clean_value(biosample_uid)

            if biosample_uid == "Unknown":
                row["error"] = "No BioSample found from Assembly summary or elink"
                rows.append(row)
                pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
                continue

            meta = fetch_biosample_metadata(
                biosample_uid,
                email=args.email,
                api_key=args.api_key,
                sleep=args.sleep
            )

            row.update(meta)
            row["error"] = ""

        except Exception as e:
            print(f"FAILED {gca}: {e}")
            row["error"] = str(e)

        rows.append(row)

        # save after every accession
        pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)

    final = pd.DataFrame(rows)
    final.to_csv(args.output, sep="\t", index=False)

    print()
    print("Done.")
    print(f"Wrote: {args.output}")
    print(f"Rows: {len(final)}")
    print(f"Assembly resolved: {(final['Assembly_UID'] != 'Unknown').sum()}")
    print(f"BioSample resolved: {(final['BioSample_UID'] != 'Unknown').sum()}")
    print(f"With geo_location: {(final['geo_location'] != 'Unknown').sum()}")
    print(f"With isolation_source: {(final['isolation_source'] != 'Unknown').sum()}")
    print(f"With host: {(final['host'] != 'Unknown').sum()}")
    print(f"With collection_date: {(final['collection_date'] != 'Unknown').sum()}")


if __name__ == "__main__":
    main()
