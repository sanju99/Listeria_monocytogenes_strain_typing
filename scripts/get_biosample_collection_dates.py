
import argparse
import subprocess
import time
import xml.etree.ElementTree as ET
import pandas as pd


def fetch_biosample_dates(accessions, batch_size=100, delay=0.4):
    """
    Fetch collection dates from NCBI BioSample in batches.

    Parameters
    ----------
    accessions : iterable
        BioSample accessions, e.g. SAMN23491780.
    batch_size : int
        Number of accessions to request per batch.
    delay : float
        Seconds to wait between requests.

    Returns
    -------
    dict
        Mapping of BioSample accession to collection date.
        Missing dates are represented by pd.NA.
    """
    accessions = list(dict.fromkeys(
        str(x).strip()
        for x in accessions
        if pd.notna(x) and str(x).strip()
    ))

    dates = {}

    for start in range(0, len(accessions), batch_size):
        batch = accessions[start:start + batch_size]
        print(
            f"Fetching batch {start // batch_size + 1}: "
            f"{len(batch)} BioSamples"
        )

        result = subprocess.run(
            [
                "efetch",
                "-db", "biosample",
                "-id", ",".join(batch),
                "-format", "xml",
            ],
            capture_output=True,
            text=True,
            check=True,
        )

        root = ET.fromstring(result.stdout)

        # Handle XML elements regardless of namespace.
        for sample in root.iter():
            if sample.tag.split("}")[-1] != "BioSample":
                continue

            accession = sample.attrib.get("accession")
            if not accession:
                continue

            collection_date = pd.NA

            for attribute in sample.iter():
                if attribute.tag.split("}")[-1] != "Attribute":
                    continue

                name = attribute.attrib.get(
                    "attribute_name", ""
                ).strip().lower().replace("_", " ")

                harmonized_name = attribute.attrib.get(
                    "harmonized_name", ""
                ).strip().lower().replace("_", " ")

                if (
                    name == "collection date"
                    or harmonized_name == "collection date"
                ):
                    value = (attribute.text or "").strip()

                    if value and value.lower() != "missing":
                        collection_date = value
                    break

            dates[accession] = collection_date

        time.sleep(delay)

    # Include accessions that were not returned by NCBI.
    for accession in accessions:
        dates.setdefault(accession, pd.NA)

    return dates


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve collection dates from NCBI BioSample."
    )
    parser.add_argument("input_csv", help="Input CSV file")
    parser.add_argument("output_csv", help="Output CSV file")
    parser.add_argument(
        "--id-column",
        default="BioSample",
        help="Column containing BioSample accessions "
             "(default: BioSample)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Number of accessions per request (default: 100)",
    )

    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)

    if args.id_column not in df.columns:
        raise ValueError(
            f"Column '{args.id_column}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    dates = fetch_biosample_dates(
        df[args.id_column],
        batch_size=args.batch_size,
    )

    # Map dates back to the original dataframe.
    df["collection_date"] = (
        df[args.id_column]
        .astype("string")
        .str.strip()
        .map(dates)
    )

    df.to_csv(args.output_csv, index=False)

    print(f"Saved results to {args.output_csv}")
    print(f"Total rows: {len(df)}")
    print(f"Dates found: {df['collection_date'].notna().sum()}")
    print(f"Dates missing: {df['collection_date'].isna().sum()}")


if __name__ == "__main__":
    main()