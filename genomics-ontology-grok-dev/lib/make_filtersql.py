#!/usr/bin/env python3
"""
make_filtersql.py
Emit an OpenCRAVAT `--filtersql` WHERE-clause that selects exactly the
actionable variant UIDs, so native OC reporters (excel/vcf/tsv) can export the
same variant set with full annotation columns.
Prints the SQL to stdout (empty string if no records).
"""
import argparse
import json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-json", required=True)
    ap.add_argument("--schema", required=True)
    args = ap.parse_args()
    data = json.load(open(args.in_json))
    uids = [r.get("uid") for r in data.get("records", []) if r.get("uid") not in (None, "")]
    if not uids:
        print("")
        return
    # base__uid is an integer primary key in OpenCRAVAT variant tables.
    ints = ",".join(str(int(u)) for u in uids)
    print(f"base__uid in ({ints})")


if __name__ == "__main__":
    main()
