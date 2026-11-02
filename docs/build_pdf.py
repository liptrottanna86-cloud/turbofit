"""
Build docs/turbofit_user_guide.pdf from the Markdown chapters in docs/.

    python docs/build_pdf.py

Needs: markdown, pygments, playwright (with Chromium: `playwright install chromium`).
The Markdown files stay the single source; this only restyles them
(pink/purple theme) and prints them with Chromium.
"""

import os
import re
import glob
import datetime
import markdown

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_HTML = os.path.join(HERE, '_guide.html')
OUT_PDF = os.path.join(HERE, 'turbofit_user_guide.pdf')

PURPLE, DEEP, PINK, LPINK, LAV = '#7B2D9E', '#3A1450', '#D6458F', '#FCE4F1', '#F6E8FB'

CALLOUTS = {'NOTE': ('Note', PURPLE, LAV), 'TIP': ('Tip', PINK, LPINK),
            'IMPORTANT': ('Important', DEEP, LAV), 'WARNING': ('Warning', '#B0246B', LPINK),
            'CAUTION': ('Caution', '#B0246B', LPINK)}


def chapters():
    files = sorted(glob.glob(os.path.join(HERE, '[0-9][0-9]-*.md')))
    return [os.path.join(HERE, 'README.md')] + files


def anchor(fname):
    return 'ch-' + os.path.splitext(os.path.basename(fname))[0].lower()


def preprocess(text, fname):
    # GitHub callouts -> styled divs
    def callout(m):
        kind = m.group(1)
        body = re.sub(r'^> ?', '', m.group(2), flags=re.M).strip()
        title, col, bg = CALLOUTS[kind]
        inner = markdown.markdown(body, extensions=['tables', 'fenced_code'])
        return (f'\n<div class="callout" style="border-color:{col};background:{bg}">'
                f'<div class="ctitle" style="color:{col}">{title}</div>{inner}</div>\n')
    text = re.sub(r'^> \[!(\w+)\]\s*\n((?:>.*\n?)*)', callout, text, flags=re.M)
    # links to other chapters -> in-document anchors
    text = re.sub(r'\]\((\d\d-[\w-]+)\.md\)', lambda m: f'](#{anchor(m.group(1))})', text)
    text = re.sub(r'\]\(README\.md\)', '](#ch-readme)', text)
    # images -> absolute paths
    text = re.sub(r'\]\(img/', '](file://' + HERE + '/img/', text)
    text = text.replace('img/workflow.svg', 'img/workflow.png')   # crisper in print
    return text


def build_html():
    exts = ['tables', 'fenced_code', 'codehilite', 'toc', 'attr_list', 'sane_lists']
    cfg = {'codehilite': {'guess_lang': False, 'css_class': 'hl'}}
    parts, toc = [], []
    for i, f in enumerate(chapters()):
        txt = preprocess(open(f, encoding='utf-8').read(), f)
        if i == 0:
            # the docs index: drop its H1 (the cover replaces it) and its contents table
            txt = re.sub(r'^# .*\n', '', txt, count=1)
            txt = re.sub(r'^\*Anna Liptrott.*\*\n', '', txt, count=1, flags=re.M)
            txt = re.sub(r'## Contents\n.*?(?=\n## )', '', txt, flags=re.S)
            txt = '# Overview\n' + txt
        html = markdown.markdown(txt, extensions=exts, extension_configs=cfg)
        title = re.search(r'<h1[^>]*>(.*?)</h1>', html).group(1)
        h2 = re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', html)
        a = anchor(f)
        # make heading ids unique per chapter
        html = re.sub(r'id="([^"]+)"', lambda m: f'id="{a}-{m.group(1)}"', html)
        html = html.replace('<h1', f'<h1 class="chapter" data-a="{a}"', 1)
        parts.append(f'<section id="{a}">{html}</section>')
        keep_subs = i > 0 and not f.endswith('api-reference.md')
        toc.append((a, title, [(f'{a}-{j}', t) for j, t in h2] if keep_subs else []))

    toc_html = '<div class="toc"><h1 class="chapter">Contents</h1><ul class="toc1">'
    for a, t, subs in toc:
        toc_html += f'<li><a href="#{a}">{t}</a>'
        if subs:
            toc_html += '<ul>' + ''.join(f'<li><a href="#{i}">{s}</a></li>' for i, s in subs) + '</ul>'
        toc_html += '</li>'
    toc_html += '</ul></div>'

    date = datetime.date.today().strftime('%B %Y')
    cover = f"""
<div class="cover">
  <div class="band"></div>
  <div class="ctext">
    <div class="kicker">User guide</div>
    <div class="title">turbofit</div>
    <div class="sub">Turbospectrum synthesis, FERRE grids, NLTE<br>and line-by-line abundances
      with CNO model atmospheres</div>
    <div class="author">Anna Liptrott</div>
    <div class="inst">Astrophysics Research Institute, Liverpool John Moores University</div>
    <div class="meta">Version 0.1 &middot; {date}</div>
  </div>
  <div class="foot">Synthesis &middot; NLTE &middot; LSFs &middot; FERRE &middot; radial velocities
    &middot; automatic line masks &middot; line-by-line abundances</div>
</div>"""

    css = f"""
@page {{ size: A4; margin: 18mm; }}
@page :first {{ margin: 0; }}
body {{ margin: 0; font-family: 'DejaVu Serif', Georgia, serif; font-size: 10.3pt; line-height: 1.45;
        color: #241430; }}
h1, h2, h3 {{ font-family: 'DejaVu Sans', Helvetica, sans-serif; }}
h1.chapter {{ page-break-before: always; color: {DEEP}; font-size: 22pt; margin: 0 0 4mm 0;
             padding-bottom: 3mm; border-bottom: 3px solid {PINK}; }}
h2 {{ color: {PURPLE}; font-size: 14pt; margin-top: 7mm; border-left: 4px solid {PINK};
      padding-left: 3mm; page-break-after: avoid; }}
h3 {{ color: {PINK}; font-size: 11.5pt; margin-top: 5mm; page-break-after: avoid; }}
a {{ color: {PURPLE}; text-decoration: none; border-bottom: 1px dotted {PINK}; }}
code {{ font-family: 'DejaVu Sans Mono', monospace; font-size: 8.8pt; background: {LAV};
        color: {DEEP}; padding: 0 2px; border-radius: 2px; }}
.hl {{ background: #FBF6FD; border-left: 4px solid {PURPLE}; border-radius: 3px;
       padding: 2mm 3mm; margin: 3mm 0; page-break-inside: avoid; }}
.hl pre {{ margin: 0; white-space: pre-wrap; }}
.hl code {{ background: none; padding: 0; font-size: 8.3pt; color: #2b1838; }}
.hl .k, .hl .kn, .hl .kc, .hl .ow {{ color: {PURPLE}; font-weight: bold; }}
.hl .s, .hl .s1, .hl .s2, .hl .sd {{ color: {PINK}; }}
.hl .mi, .hl .mf, .hl .m {{ color: #B0246B; }}
.hl .c, .hl .c1, .hl .cm {{ color: #8a7a96; font-style: italic; }}
.hl .nb, .hl .bp {{ color: #5B2380; }}
.hl .nf, .hl .nc {{ color: {DEEP}; font-weight: bold; }}
table {{ border-collapse: collapse; width: 100%; margin: 3mm 0; font-size: 9pt;
         page-break-inside: avoid; }}
th {{ background: {PURPLE}; color: white; text-align: left; padding: 1.5mm 2mm;
      font-family: 'DejaVu Sans', sans-serif; font-weight: 600; }}
td {{ padding: 1.3mm 2mm; border-bottom: 1px solid #ead7f2; vertical-align: top; }}
tr:nth-child(even) td {{ background: #FDF5FA; }}
img {{ max-width: 100%; display: block; margin: 3mm auto; }}
.callout {{ border-left: 5px solid; border-radius: 3px; padding: 2mm 4mm; margin: 4mm 0;
            page-break-inside: avoid; }}
.callout .ctitle {{ font-family: 'DejaVu Sans', sans-serif; font-weight: bold; font-size: 9pt;
                    text-transform: uppercase; letter-spacing: 0.06em; margin-bottom: 1mm; }}
.callout p {{ margin: 1mm 0; }}
.cover {{ height: 296mm; position: relative; margin: 0; }}
.cover .band {{ position: absolute; top: 0; left: 0; right: 0; bottom: 0;
               background: linear-gradient(160deg, {DEEP} 0%, {PURPLE} 50%, {PINK} 100%); }}
.cover .ctext {{ position: absolute; top: 88mm; left: 24mm; right: 24mm; color: white; }}
.cover .foot {{ position: absolute; bottom: 28mm; left: 24mm; right: 24mm; color: white;
               font-family: 'DejaVu Sans', sans-serif; font-size: 9.5pt; opacity: 0.85;
               border-top: 1px solid rgba(255,255,255,0.5); padding-top: 3mm; }}
.cover .kicker {{ font-family: 'DejaVu Sans', sans-serif; letter-spacing: 0.3em;
                 text-transform: uppercase; font-size: 11pt; opacity: 0.85; }}
.cover .title {{ font-family: 'DejaVu Sans', sans-serif; font-size: 60pt; font-weight: bold;
                margin: 4mm 0 2mm 0; }}
.cover .sub {{ font-size: 14pt; line-height: 1.35; }}
.cover .author {{ margin-top: 14mm; font-family: 'DejaVu Sans', sans-serif; font-size: 15pt;
                 font-weight: bold; }}
.cover .inst {{ font-size: 11pt; opacity: 0.92; margin-top: 1mm; }}
.cover .meta {{ margin-top: 8mm; font-size: 10pt; opacity: 0.9; }}
.cover .cfig {{ position: absolute; top: 125mm; left: 0; width: 100%; }}
.toc ul.toc1 {{ padding-left: 0; list-style: none; }}
.toc ul.toc1 > li {{ margin: 2.5mm 0; font-family: 'DejaVu Sans', sans-serif;
                     font-weight: bold; color: {PURPLE}; }}
.toc ul {{ font-family: 'DejaVu Serif', serif; font-weight: normal; font-size: 9.5pt;
           list-style: none; padding-left: 5mm; margin: 1mm 0; columns: 2; }}
.toc ul li {{ margin: 0.5mm 0; }}
.toc a {{ border: none; }}
"""
    html = (f'<!doctype html><html><head><meta charset="utf-8"><style>{css}</style></head>'
            f'<body>{cover}{toc_html}{"".join(parts)}</body></html>')
    with open(OUT_HTML, 'w', encoding='utf-8') as f:
        f.write(html)
    return OUT_HTML


def build_pdf():
    from playwright.sync_api import sync_playwright
    html = build_html()
    footer = (f'<div style="width:100%;font-size:8px;font-family:DejaVu Sans,sans-serif;'
              f'color:{PURPLE};padding:0 18mm;display:flex;justify-content:space-between">'
              f'<span>turbofit user guide</span><span class="pageNumber"></span></div>')
    with sync_playwright() as p:
        b = p.chromium.launch()
        pg = b.new_page()
        pg.goto('file://' + html)
        pg.wait_for_timeout(500)
        pg.pdf(path=OUT_PDF, format='A4', print_background=True, display_header_footer=True,
               header_template='<div></div>', footer_template=footer,
               prefer_css_page_size=True)
        b.close()
    os.remove(html)
    print('wrote', OUT_PDF)


if __name__ == '__main__':
    build_pdf()
