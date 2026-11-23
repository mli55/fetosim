"""Put the block-by-block section into docs/index.html: blocks.css before </head>, section.html before the circuit
section, traces.json and blocks.js before </body>, each between marker comments (replaced on every build).
Run after docs/blocks/run.py: python docs/blocks/build.py [page.html]"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAGE = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / 'index.html'


def put(html, name, body, anchor):
    """Replace the block between <!-- blocks:name --> markers, or insert it before `anchor` (first build)."""
    a, b = f'<!-- blocks:{name} -->', f'<!-- /blocks:{name} -->'
    block = f'{a}\n{body}\n{b}'
    if a in html:
        return re.sub(re.escape(a) + '.*?' + re.escape(b), lambda _: block, html, flags=re.S)
    assert html.count(anchor) == 1, anchor
    return html.replace(anchor, block + '\n' + anchor)


html = PAGE.read_text()
html = put(html, 'css', '<style>\n' + (HERE / 'blocks.css').read_text() + '</style>', '</head>')
html = put(html, 'section', (HERE / 'section.html').read_text().strip(), '<section id="circuit-sec"')
data = (HERE / 'traces.json').read_text().replace('</', '<\\/')
html = put(html, 'script', f'<script type="application/json" id="bk-traces">{data}</script>\n<script>\n'
           + (HERE / 'blocks.js').read_text() + '</script>', '</body>')
PAGE.write_text(html)
print(PAGE, len(html) // 1024, 'KB')
