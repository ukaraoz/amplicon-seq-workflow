#!/usr/bin/env python3
"""Download the fungal trait databases used by fungal-guilds.R and convert them to tab-separated tables.

usage: get-guild-dbs.py <dest_dir>

  FungalTraits 1.2 (Polme et al. 2020, Fungal Diversity 105:1-16; genus-level lifestyles and traits): the
    article's supplementary table (sheet "data") -> FungalTraits_1.2.tsv
  FUNGuild (Nguyen et al. 2016, Fungal Ecology 20:241-248; guilds by taxon with a confidence ranking): the current
    database from the FUNGuild server -> FUNGuild.tsv (the confidence typo "Higly Probable" is corrected)

The originals are kept next to the tables, and README.txt records the sources, download date, checksums and record
counts. Files already in dest_dir are not downloaded again (delete them to update). Standard library only.
"""
import datetime
import hashlib
import json
import os
import re
import sys
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

FT_URL = ("https://static-content.springer.com/esm/art%3A10.1007%2Fs13225-020-00466-2/MediaObjects/"
          "13225_2020_466_MOESM4_ESM.xlsx")
FT_XLSX = "FungalTraits_1.2_ver_16Dec_2020.xlsx"
FT_TSV = "FungalTraits_1.2.tsv"
FG_URL = "http://www.stbates.org/funguild_db_2.php"
FG_RAW = "funguild_db_2.json"
FG_TSV = "FUNGuild.tsv"
FG_FIELDS = ["taxon", "taxonomicLevel", "trophicMode", "guild", "confidenceRanking", "growthForm", "trait",
             "notes", "citationSource"]
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "rel": "http://schemas.openxmlformats.org/package/2006/relationships"}


def fetch(url, path):
    if os.path.exists(path):
        print("present:", path)
        return
    print("downloading", url)
    req = urllib.request.Request(url, headers={"User-Agent": "amplicon-seq-workflow"})
    with urllib.request.urlopen(req, timeout=300) as r, open(path + ".part", "wb") as fh:
        fh.write(r.read())
    os.replace(path + ".part", path)


def clean(v):
    return re.sub(r"[\t\r\n]+", " ", "" if v is None else str(v)).strip()


def xlsx_sheet(path, name):
    """Rows (lists of strings) of one sheet of an .xlsx file."""
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rid = next(s.get("{%s}id" % NS["r"]) for s in wb.find("m:sheets", NS) if s.get("name") == name)
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = next(r.get("Target") for r in rels.findall("rel:Relationship", NS) if r.get("Id") == rid)
    sheet = ET.fromstring(z.read("xl/" + target.lstrip("/").replace("xl/", "", 1)))
    rows = []
    for row in sheet.iter("{%s}row" % NS["m"]):
        cells = {}
        for c in row.findall("m:c", NS):
            col = 0
            for ch in re.match(r"[A-Z]+", c.get("r")).group(0):
                col = col * 26 + ord(ch) - 64
            t, v = c.get("t"), c.find("m:v", NS)
            if t == "s" and v is not None:
                val = shared[int(v.text)]
            elif t == "inlineStr":
                val = "".join(x.text or "" for x in c.iter("{%s}t" % NS["m"]))
            else:
                val = v.text if v is not None else ""
            cells[col - 1] = val
        if cells:
            rows.append([cells.get(i, "") for i in range(max(cells) + 1)])
    return rows


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    dest = sys.argv[1]
    os.makedirs(dest, exist_ok=True)
    p = lambda f: os.path.join(dest, f)

    # FungalTraits: sheet "data", one row per genus; '_template' dropped from the column names
    fetch(FT_URL, p(FT_XLSX))
    rows = xlsx_sheet(p(FT_XLSX), "data")
    hdr = [re.sub(r"_template$", "", clean(h)) for h in rows[0]]
    with open(p(FT_TSV), "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows[1:]:
            r = (r + [""] * len(hdr))[:len(hdr)]
            if clean(r[hdr.index("GENUS")]):
                fh.write("\t".join(clean(x) for x in r) + "\n")
    n_ft = sum(1 for _ in open(p(FT_TSV))) - 1

    # FUNGuild: the server returns an HTML page whose 7th line is the JSON array (as Guilds_v1.1.py reads it)
    if not os.path.exists(p(FG_RAW)):
        tmp = p(FG_RAW + ".html")
        fetch(FG_URL, tmp)
        body = open(tmp, encoding="utf-8").read().split("\n")[6].strip()
        body = body[:-len("</body>")] if body.endswith("</body>") else body
        db = json.loads(body)
        with open(p(FG_RAW), "w") as fh:
            json.dump(db, fh)
        os.remove(tmp)
    db = json.load(open(p(FG_RAW)))
    with open(p(FG_TSV), "w") as fh:
        fh.write("\t".join(FG_FIELDS) + "\n")
        for rec in db:
            rec = dict(rec)
            rec["confidenceRanking"] = clean(rec.get("confidenceRanking")).replace("Higly", "Highly")
            fh.write("\t".join(clean(rec.get(k)) for k in FG_FIELDS) + "\n")

    today = datetime.date.today().isoformat()
    readme = p("README.txt")
    if not os.path.exists(readme):
        with open(readme, "w") as fh:
            fh.write(f"Fungal trait databases for fungal-guilds.R, downloaded {today} by get-guild-dbs.py\n\n")
            fh.write(f"FungalTraits 1.2 (Polme et al. 2020, Fungal Diversity 105:1-16, doi:10.1007/s13225-020-00466-2)\n"
                     f"  source: {FT_URL}\n  {FT_XLSX}  md5 {md5(p(FT_XLSX))}\n"
                     f"  {FT_TSV}: {n_ft} genera (sheet 'data'; '_template' dropped from column names)\n\n")
            fh.write(f"FUNGuild (Nguyen et al. 2016, Fungal Ecology 20:241-248, doi:10.1016/j.funeco.2015.06.006)\n"
                     f"  source: {FG_URL}\n  {FG_RAW}  md5 {md5(p(FG_RAW))}\n"
                     f"  {FG_TSV}: {len(db)} records; taxonomicLevel 3 phylum, 7 order, 9 family, 13 genus, "
                     f"20 species, 21-24 below species; 'Higly Probable' corrected to 'Highly Probable'\n")
    print(f"{FT_TSV}: {n_ft} genera; {FG_TSV}: {len(db)} records -> {dest}")


if __name__ == "__main__":
    main()
