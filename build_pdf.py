"""manuscript_v2.md → manuscript_v2.html → manuscript_v2.pdf (Edge headless; ไม่ต้องติดตั้งแพ็กเกจ) · รูปวางท้ายเอกสารพร้อมคำบรรยาย"""
import os, re, html, subprocess, sys

D = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(D, sys.argv[1] if len(sys.argv) > 1 else "manuscript_v2.md")
BASE = os.path.splitext(SRC)[0]
FIGS = [("fig1_benchmark", "Figure 1"), ("fig2_knockout_screen", "Figure 2"), ("fig3_why_knockouts_fail", "Figure 3"),
        ("fig4_usnea", "Figure 4"), ("fig5_dnge031_male", "Figure 5")]


def inline(t):
    t = html.escape(t, quote=False)
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<i>\1</i>", t)
    t = re.sub(r"`(.+?)`", r"<code>\1</code>", t)
    t = re.sub(r"(https?://[^\s)<]+)", r'<a href="\1">\1</a>', t)
    return t


def convert(md):
    out, i, L = [], 0, md.splitlines()
    legends = {}
    while i < len(L):
        ln = L[i]
        if ln.startswith("**Figure ") and "Figure legends" in "".join(out[-3:] + [""]) or re.match(r"^\*\*Figure \d\.\*\*", ln):
            m = re.match(r"^\*\*(Figure \d)\.\*\*\s*(.*)", ln)
            if m:
                legends[m.group(1)] = m.group(2); i += 1; continue
        if ln.startswith("|"):
            rows = []
            while i < len(L) and L[i].startswith("|"):
                cells = [c.strip() for c in L[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells):
                    rows.append(cells)
                i += 1
            t = "<table><thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in rows[0]) + "</tr></thead><tbody>"
            t += "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in rows[1:]) + "</tbody></table>"
            out.append(t); continue
        if re.match(r"^\s*[-*] ", ln) or re.match(r"^\s*\d+\. ", ln):
            ordered = bool(re.match(r"^\s*\d+\. ", ln)); tag = "ol" if ordered else "ul"; items = []
            while i < len(L) and (re.match(r"^\s*[-*] ", L[i]) or re.match(r"^\s*\d+\. ", L[i])):
                items.append(re.sub(r"^\s*([-*]|\d+\.)\s+", "", L[i])); i += 1
            out.append(f"<{tag}>" + "".join(f"<li>{inline(x)}</li>" for x in items) + f"</{tag}>"); continue
        if ln.startswith("### "): out.append(f"<h3>{inline(ln[4:])}</h3>")
        elif ln.startswith("## "): out.append(f"<h2>{inline(ln[3:])}</h2>")
        elif ln.startswith("# "): out.append(f"<h1>{inline(ln[2:])}</h1>")
        elif ln.strip() == "---": out.append("<hr>")
        elif ln.strip() == "": pass
        else:
            para = [ln]; i += 1
            while i < len(L) and L[i].strip() and not re.match(r"^(#|\||\s*[-*] |\s*\d+\. |---|\*\*Figure \d)", L[i]):
                para.append(L[i]); i += 1
            out.append("<p>" + inline(" ".join(para)).replace("\\*", "*") + "</p>"); continue
        i += 1
    return "\n".join(out), legends


md = open(SRC, encoding="utf-8").read()
body, legends = convert(md)
figs = "".join(f'<div class="fig"><img src="figures/{f}.png"><p><b>{n}.</b> {inline(legends.get(n, ""))}</p></div>' for f, n in FIGS)
page = f"""<!doctype html><html><head><meta charset="utf-8"><title>FlyBrain preprint</title><style>
@page {{ size: A4; margin: 18mm 17mm; }}
body {{ font-family: "Times New Roman", serif; font-size: 10.5pt; line-height: 1.45; color: #000; }}
h1 {{ font-size: 16pt; line-height: 1.25; margin: 0 0 8px; }} h2 {{ font-size: 12.5pt; margin: 16px 0 6px; }} h3 {{ font-size: 11pt; margin: 12px 0 4px; }}
p {{ margin: 0 0 7px; text-align: justify; }} table {{ border-collapse: collapse; width: 100%; font-size: 8.3pt; margin: 6px 0 10px; }}
td, th {{ border: 0.5pt solid #999; padding: 2px 4px; vertical-align: top; }} th {{ background: #eee; }}
hr {{ border: 0; border-top: 0.5pt solid #bbb; margin: 10px 0; }} code {{ font-size: 9pt; }}
.fig {{ page-break-inside: avoid; margin: 14px 0; }} .fig img {{ width: 100%; }} .fig p {{ font-size: 9pt; }}
a {{ color: #000; }} li {{ margin-bottom: 3px; }}
</style></head><body>{body}<h2>Figures</h2>{figs}</body></html>"""
open(BASE + ".html", "w", encoding="utf-8").write(page)
edge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
subprocess.run([edge, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--print-to-pdf={BASE}.pdf", "file:///" + (BASE + ".html").replace("\\", "/")],
               check=True, timeout=180, capture_output=True)
print("wrote", BASE + ".pdf", os.path.getsize(BASE + ".pdf"), "bytes")
