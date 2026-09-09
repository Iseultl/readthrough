import pandas as pd

# Read annotations.tsv (has headers)
annotations = pd.read_csv(
    "/Users/iseult/gitlab/ExtractGenes/annotations.tsv",
    sep="\t"
)

# Read model_annotations.tsv (appears to have no header)
model_annotations = pd.read_csv(
    "model_annotations.tsv",
    sep="\t",
    header=None
)

# Assign column names based on the structure shown
model_annotations.columns = [
    "annotation_id",
    "assembly",
    "assembly_name",
    "Taxon name",
    "taxid",
    "source",
    "annotation_source",
    "model_gff",
    "gff_path",
    "gff_index"
]

# Merge using annotation_id
merged = model_annotations.merge(
    annotations,
    on="annotation_id",
    how="left"
)

# Select the desired columns
output = merged[
    [
        "annotation_id",
        "taxid",
        "Taxon name",
        "annotation_url",
        "assembly_url"
    ]
].copy()

# Remove repeated TaxIDs, keeping the first occurrence
output = output.drop_duplicates(
    subset="taxid",
    keep="first"
)

# Save
output.to_csv(
    "merged_annotations.tsv",
    sep="\t",
    index=False
)

# Summary statistics
print(f"Total model annotations: {len(model_annotations)}")
print(f"Unique TaxIDs retained: {len(output)}")
print(f"Duplicate TaxID rows removed: {len(merged) - len(output)}")
print(f"Matched annotations: {merged['annotation_url'].notna().sum()}")
print(f"Missing annotations: {merged['annotation_url'].isna().sum()}")

# Show unmatched annotation IDs
missing = merged[merged["annotation_url"].isna()]
print("\nFirst 20 missing annotation IDs:")
print(missing["annotation_id"].head(20))