"""Read .xlsx with the Python standard library only (no openpyxl) -> one CSV per sheet."""
import zipfile, re, csv, os, sys
import xml.etree.ElementTree as ET
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def col_idx(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1


def read(path):
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            ss.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    rmap = {r.get("Id"): r.get("Target") for r in rels}
    out = {}
    for sh in wb.find("m:sheets", NS):
        tgt = rmap[sh.get("{%s}id" % NS["r"])].lstrip("/")
        tgt = tgt if tgt.startswith("xl/") else "xl/" + tgt
        rows = []
        for row in ET.fromstring(z.read(tgt)).iter("{%s}row" % NS["m"]):
            vals = {}
            for c in row.findall("m:c", NS):
                v = c.find("m:v", NS); t = c.get("t")
                if t == "inlineStr":
                    x = "".join(e.text or "" for e in c.iter("{%s}t" % NS["m"]))
                elif v is None:
                    continue
                else:
                    x = ss[int(v.text)] if t == "s" else v.text
                vals[col_idx(c.get("r"))] = x
            if vals:
                rows.append([vals.get(i, "") for i in range(max(vals) + 1)])
        out[sh.get("name")] = rows
    return out


if __name__ == "__main__":
    src = sys.argv[1]; od = os.path.join(os.path.dirname(src), "csv"); os.makedirs(od, exist_ok=True)
    for name, rows in read(src).items():
        fn = os.path.join(od, re.sub(r"[^\w.-]+", "_", name) + ".csv")
        with open(fn, "w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerows(rows)
        print(f"{name}: {len(rows)} rows | first: {rows[0][:6] if rows else ''}")
