"""
Step 8 — Convert manuscript_comparative_v2.md → PDF
Uses pandoc (HTML) + Chrome headless (same method as genesis project).

Usage:
  python analysis/08_convert_to_pdf.py [--input reports/manuscript_comparative_v2.md]
"""
import argparse
import base64
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path('e:/miniconda3/envs/llama-env/comparative_amr_project')
PANDOC  = 'e:/miniconda3/bin/pandoc.exe'
CHROME  = 'C:/Program Files/Google/Chrome/Application/chrome.exe'

CSS = '''
@page { size: A4; margin: 2.5cm 2.5cm 2.5cm 3.0cm; }
* { box-sizing: border-box; }
body {
  font-family: "Times New Roman", Times, serif;
  font-size: 11pt; line-height: 1.5; color: #000; margin: 0; padding: 0;
}
h1 { font-size: 13pt; font-weight: bold; text-align: center;
     margin: 0 0 0.5em 0; line-height: 1.35; }
.author-block { text-align: center; margin-bottom: 1em; }
.author-block .name  { font-size: 11pt; font-weight: bold; }
.author-block .affil { font-size: 10pt; font-style: italic; }
.author-block .corresp { font-size: 10pt; }
h2 { font-size: 12pt; font-weight: bold; margin: 1.2em 0 0.4em 0;
     border-bottom: 1px solid #444; padding-bottom: 2pt; }
h3 { font-size: 11pt; font-weight: bold; font-style: italic;
     margin: 0.9em 0 0.3em 0; }
h4 { font-size: 11pt; font-weight: bold; margin: 0.8em 0 0.2em 0; }
p { margin: 0 0 0.5em 0; text-align: justify; }
.abstract-section {
  background: #dce6ff; padding: 0.35cm 0.5cm; margin: 0.4cm 0;
  font-size: 10pt; line-height: 1.5;
}
table { border-collapse: collapse; width: 100%; font-size: 9pt;
        margin: 0.7em 0 1em 0; page-break-inside: avoid; }
th { background: #dce6ff; border: 1px solid #666; padding: 4pt 6pt;
     font-weight: bold; text-align: center; }
td { border: 1px solid #888; padding: 3pt 5pt;
     text-align: center; vertical-align: middle; }
ul, ol { font-size: 10pt; margin: 0.3em 0 0.3em 1.4em; }
li { margin-bottom: 2pt; }
hr { border: none; border-top: 1px solid #aaa; margin: 0.7em 0; }
.references p { font-size: 9pt; text-indent: -1.4em; padding-left: 1.4em;
                margin-bottom: 3pt; text-align: left; }
'''


def convert_md_to_html(md_text: str) -> str:
    with tempfile.NamedTemporaryFile(suffix='.md', delete=False,
                                     mode='w', encoding='utf-8') as f:
        f.write(md_text)
        tmp = f.name
    try:
        r = subprocess.run(
            [PANDOC, tmp, '--from=gfm', '--to=html5', '--no-highlight'],
            capture_output=True, text=True, encoding='utf-8', check=True)
        return r.stdout
    finally:
        os.unlink(tmp)


def post_process(html: str) -> str:
    # Unicode superscripts → <sup>
    for ch, d in {'⁰':'0','¹':'1','²':'2','³':'3','⁴':'4',
                  '⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁻':'-'}.items():
        html = html.replace(ch, f'<sup>{d}</sup>')
    # Abstract block styling
    html = re.sub(
        r'(<h3[^>]*>(?:Background|Methods|Results|Conclusions)</h3>\s*)(.*?)(<h3|<h2|<hr)',
        lambda m: m.group(1) + '<div class="abstract-section">' +
                  m.group(2) + '</div>' + m.group(3),
        html, flags=re.DOTALL | re.IGNORECASE)
    # References styled div
    html = re.sub(
        r'(<h2[^>]*>References</h2>\s*)(.*?)($)',
        lambda m: m.group(1) + '<div class="references">' + m.group(2) + '</div>',
        html, flags=re.DOTALL | re.IGNORECASE)
    return html


def assemble_html(body: str, title: str) -> str:
    author_html = '''
<div class="author-block">
  <div class="name">ZJY</div>
  <div class="affil">Independent Researcher</div>
  <div class="corresp">Correspondence: jiayu6954@gmail.com</div>
</div>'''
    title_m = re.search(r'<h1[^>]*>(.*?)</h1>', body, re.DOTALL)
    if title_m:
        body = body.replace(title_m.group(0),
                            title_m.group(0) + author_html, 1)
    return f'''<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<title>{title}</title>
<style>{CSS}</style>
</head><body><div class="main-content">{body}</div></body></html>'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default=str(
        PROJECT_ROOT / 'reports' / 'manuscript_comparative_v2.md'))
    args = ap.parse_args()

    md_path  = Path(args.input).resolve()
    pdf_path = md_path.with_suffix('.pdf')

    print(f'Input : {md_path}')
    print(f'Output: {pdf_path}')

    text     = md_path.read_text(encoding='utf-8')
    html_frag = convert_md_to_html(text)
    html_frag = post_process(html_frag)
    full_html = assemble_html(html_frag,
                              'IS Element Ecology Across CR-GNB Pathogens')

    with tempfile.NamedTemporaryFile(suffix='.html', delete=False,
                                     mode='w', encoding='utf-8') as f:
        f.write(full_html)
        tmp_html = f.name

    print('Chrome headless rendering…')
    r = subprocess.run([
        CHROME,
        '--headless', '--disable-gpu', '--disable-extensions',
        '--no-sandbox', '--run-all-compositor-stages-before-draw',
        '--virtual-time-budget=8000',
        f'--print-to-pdf={pdf_path}',
        '--no-pdf-header-footer',
        f'file:///{Path(tmp_html).as_posix()}',
    ], capture_output=True)   # bytes mode — avoids GBK/UTF-8 decoding errors

    os.unlink(tmp_html)

    if r.returncode != 0:
        err = r.stderr.decode('utf-8', errors='replace')[:400]
        print('Chrome error:', err)
        sys.exit(1)

    kb = round(pdf_path.stat().st_size / 1024)
    print(f'\nPDF saved : {pdf_path}')
    print(f'Size      : {kb} KB')


if __name__ == '__main__':
    main()
