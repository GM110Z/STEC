#!/usr/bin/env python3

import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt


CLASSES = [
    "Human/clinical",
    "Animal/livestock/wildlife",
    "Food",
    "Environmental/water",
    "Lab/reference/cell culture",
    "Unknown"
]


COLOURS = {
    "Human/clinical": "#d73027",
    "Animal/livestock/wildlife": "#fc8d59",
    "Food": "#1a9850",
    "Environmental/water": "#4575b4",
    "Lab/reference/cell culture": "#984ea3",
    "Unknown": "#f0f0f0"
}


UNKNOWN_VALUES = {
    "",
    " ",
    "missing",
    "unknown",
    "not applicable",
    "not_applicable",
    "n/a",
    "na",
    "none",
    "not provided",
    "not collected",
    "not determined",
    "nd",
    "null",
    "nan",
    "-",
    "--"
}


def clean(x):
    if pd.isna(x):
        return "Unknown"
    x = str(x).strip()
    if x.lower() in UNKNOWN_VALUES:
        return "Unknown"
    return x


def clean_gca(x):
    return str(x).replace("\ufeff", "").strip().replace(" ", "")


def clean_phylogroup(x):
    x = clean(x)
    if x == "Unknown":
        return "Unknown"
    return x.strip()


def has_useful_metadata(text):
    low = text.lower()

    useful_keys = [
        "host=",
        "host scientific name=",
        "host disease=",
        "host_disease=",
        "host health state=",
        "isolation source=",
        "isolation_source=",
        "source type=",
        "source_type=",
        "sample type=",
        "sample_type=",
        "env_medium=",
        "environmental medium=",
        "attribute_package=",
        "ifsac",
        "ontological term=",
    ]

    if not any(k in low for k in useful_keys):
        return False

    informative_words = [
        "homo sapiens", "human", "stool", "faeces", "feces", "clinical",
        "cattle", "bovine", "cow", "chicken", "pig", "swine", "deer",
        "animal", "livestock", "food", "flour", "lettuce", "sprout",
        "water", "pond", "soil", "environmental", "cell culture",
        "atcc", "biomaterial", "reference material"
    ]

    return any(w in low for w in informative_words)


def classify(row):
    parts = []

    for col in [
        "all_biosample_attributes",
        "isolation_source",
        "host",
        "disease",
        "host_health_state",
        "biosample_title",
        "organism",
        "strain",
        "sample_name"
    ]:
        if col in row.index:
            parts.append(clean(row[col]))

    text = " ; ".join(parts)
    low = text.lower()

    if not has_useful_metadata(low):
        return "Unknown", "no useful source metadata"

    if re.search(r"source type=human|source_type=human", low):
        return "Human/clinical", "source_type=human"

    if re.search(r"source type=animal|source_type=animal", low):
        return "Animal/livestock/wildlife", "source_type=animal"

    if re.search(r"source type=food|source_type=food", low):
        return "Food", "source_type=food"

    if re.search(r"source type=environmental|source_type=environmental", low):
        return "Environmental/water", "source_type=environmental"

    if re.search(r"attribute_package=clinical/host-associated", low):
        if re.search(r"homo sapiens|\bhuman\b|clinical|patient|stool|faeces|feces|urine|blood|diarr|hus", low):
            return "Human/clinical", "attribute_package clinical + human/clinical terms"
        return "Animal/livestock/wildlife", "attribute_package clinical/host-associated"

    if re.search(r"attribute_package=environmental/food/other", low):
        if re.search(r"food|flour|lettuce|sprout|meat|milk|retail|supermarket", low):
            return "Food", "attribute_package environmental/food/other + food terms"
        return "Environmental/water", "attribute_package environmental/food/other"

    human_patterns = [
        r"homo sapiens",
        r"\bhuman\b",
        r"clinical",
        r"patient",
        r"stool",
        r"faeces",
        r"feces",
        r"fecal",
        r"faecal",
        r"human faeces",
        r"human feces",
        r"human intestinal",
        r"intestinal microflora",
        r"urine",
        r"blood",
        r"diarr",
        r"\bhus\b",
        r"hemorrhagic colitis",
        r"haemorrhagic colitis",
        r"hemolytic uremic",
        r"haemolytic uremic",
        r"gastroenteritis",
        r"bloody diarr",
    ]

    animal_patterns = [
        r"host=bos taurus",
        r"host=bovine",
        r"\bbovine\b",
        r"\bcattle\b",
        r"\bcow\b",
        r"sus scrofa",
        r"\bswine\b",
        r"\bpig\b",
        r"\bporcine\b",
        r"chicken",
        r"gallus",
        r"deer",
        r"cervus",
        r"animal manure",
        r"livestock",
        r"canine",
        r"feline",
        r"rabbit",
        r"duck",
        r"horse",
        r"equus",
        r"wild deer",
        r"wildlife",
    ]

    food_patterns = [
        r"\bfood\b",
        r"flour",
        r"lettuce",
        r"alfalfa sprout",
        r"sprout",
        r"meat",
        r"retail",
        r"supermarket",
        r"milk",
        r"seeded vegetables",
        r"foodon",
    ]

    environmental_patterns = [
        r"environment",
        r"environmental",
        r"water",
        r"wastewater",
        r"waste water",
        r"pond",
        r"soil",
        r"envo",
        r"env_medium",
        r"environmental medium",
        r"forest biome",
        r"terrestrial biome",
    ]

    lab_patterns = [
        r"sample type=cell culture",
        r"sample_type=cell culture",
        r"cell culture",
        r"culture collection",
        r"culture_collection",
        r"biomaterial provider",
        r"reference material",
        r"type strain",
        r"\batcc\b",
        r"\bnctc\b",
        r"\bdsmz\b",
        r"\bjcm\b",
        r"\bccug\b",
        r"\bfdaargos\b",
    ]

    if any(re.search(p, low) for p in human_patterns):
        return "Human/clinical", "matched human/clinical terms"

    if any(re.search(p, low) for p in animal_patterns):
        return "Animal/livestock/wildlife", "matched animal/host terms"

    if any(re.search(p, low) for p in food_patterns):
        return "Food", "matched food terms"

    if any(re.search(p, low) for p in environmental_patterns):
        return "Environmental/water", "matched environmental/water terms"

    if any(re.search(p, low) for p in lab_patterns):
        return "Lab/reference/cell culture", "matched lab/reference/cell culture terms"

    return "Unknown", "metadata present but no recognised class"


def main():
    parser = argparse.ArgumentParser(
        description="Classify BioSample source and plot source class by phylogroup."
    )

    parser.add_argument(
        "input",
        help="Merged metadata/tox TSV containing Phylogroup and all_biosample_attributes"
    )

    parser.add_argument(
        "--prefix",
        default="source_by_phylogroup",
        help="Output prefix"
    )

    parser.add_argument(
        "--count_mode",
        choices=["isolates", "rows"],
        default="isolates",
        help="Count unique isolates or rows. Default: isolates."
    )

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")
    df.columns = (
        df.columns.astype(str)
        .str.replace("\ufeff", "", regex=False)
        .str.strip()
    )

    if "Phylogroup" not in df.columns:
        raise ValueError("Input needs column: Phylogroup")

    if "all_biosample_attributes" not in df.columns:
        raise ValueError("Input needs column: all_biosample_attributes")

    df["Phylogroup"] = df["Phylogroup"].apply(clean_phylogroup)
    df = df[df["Phylogroup"] != "Unknown"].copy()

    classified = df.apply(classify, axis=1, result_type="expand")
    df["source_class"] = classified[0]
    df["source_reason"] = classified[1]

    accession_col = None
    for c in ["GCA_clean", "GCA", "Bacterial accession", "Assembly_accession"]:
        if c in df.columns:
            accession_col = c
            break

    if args.count_mode == "isolates" and accession_col is not None:
        df["_isolate_id"] = df[accession_col].apply(clean_gca)
        counted = df.drop_duplicates(
            subset=["_isolate_id", "Phylogroup", "source_class"]
        ).copy()
        count_label = "Number of isolates"
    else:
        counted = df.copy()
        count_label = "Number of entries"

    counts = pd.crosstab(counted["Phylogroup"], counted["source_class"])

    for c in CLASSES:
        if c not in counts.columns:
            counts[c] = 0

    counts = counts[CLASSES]

    # Order phylogroups by natural-ish order if present, then extras
    preferred_order = [
        "A", "B1", "B2", "B2-2", "C", "D", "D2",
        "E", "E1", "F", "G", "Shigella"
    ]

    present = [x for x in preferred_order if x in counts.index]
    extra = sorted([x for x in counts.index if x not in present])
    counts = counts.loc[present + extra]

    # Save outputs
    df.to_csv(f"{args.prefix}.classified_all_rows.tsv", sep="\t", index=False)
    counted.to_csv(f"{args.prefix}.counted_rows.tsv", sep="\t", index=False)
    counts.to_csv(f"{args.prefix}.counts.tsv", sep="\t")

    # Horizontal stacked barplot
    ax = counts.plot(
        kind="barh",
        stacked=True,
        figsize=(8, max(5, 0.5 * len(counts))),
        color=[COLOURS[c] for c in CLASSES],
        edgecolor="black",
        linewidth=0.8
    )

    ax.set_xlabel(count_label)
    ax.set_ylabel("Phylogroup")
    ax.legend(
        title="Source class",
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    ax.invert_yaxis()
    plt.tight_layout()

    plt.savefig(f"{args.prefix}.stacked_bar.png", dpi=300)
    plt.savefig(f"{args.prefix}.stacked_bar.pdf")
    plt.close()

    print("Done.")
    print(f"Rows classified: {len(df)}")
    print(f"Counting: {count_label}")
    if accession_col:
        print(f"Accession column used: {accession_col}")

    print("\nSource class counts in all rows:")
    print(df["source_class"].value_counts().to_string())

    print("\nWrote:")
    print(f"  {args.prefix}.classified_all_rows.tsv")
    print(f"  {args.prefix}.counted_rows.tsv")
    print(f"  {args.prefix}.counts.tsv")
    print(f"  {args.prefix}.stacked_bar.png")
    print(f"  {args.prefix}.stacked_bar.pdf")


if __name__ == "__main__":
    main()
