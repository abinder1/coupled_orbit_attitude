"""Build an interactive HTML report from MATLAB test and coverage results.

Reads the JUnit test results and Cobertura coverage produced by
``buildtool test`` and writes one self-contained HTML file (no network or
third-party packages needed).

Usage (from the project root):
    python test_results/test_report.py                 # write report.html
    python test_results/test_report.py --open          # ...and open it
    python test_results/test_report.py --run --open    # run buildtool test first
    python test_results/test_report.py --watch --open  # live: the page updates
                                                       # after every test run
"""

import argparse
import json
import subprocess
import sys
import time
import webbrowser
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def parse_tests(path):
    """Return suites with per-test status, time, and failure details."""
    suites = []
    for suite in ET.parse(path).getroot().iter("testsuite"):
        cases = []
        for case in suite.findall("testcase"):
            status, detail = "passed", ""
            for kind, label in (("failure", "failed"), ("error", "error"),
                                ("skipped", "skipped")):
                node = case.find(kind)
                if node is not None:
                    status = label
                    detail = (node.get("message") or "") + "\n" + \
                        (node.text or "")
                    break
            cases.append({"name": case.get("name", ""), "status": status,
                          "time": float(case.get("time") or 0),
                          "detail": detail.strip()})
        suites.append({"name": suite.get("name", ""),
                       "time": float(suite.get("time") or 0),
                       "cases": cases})
    return suites


def read_source(path):
    try:
        return path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return None


def parse_coverage(path):
    """Return per-file line coverage and the source text for each file."""
    root = ET.parse(path).getroot()
    source_dirs = [Path(s.text.strip()) for s in root.iter("source")
                   if s.text and s.text.strip()]
    files = []
    for package in root.iter("package"):
        for cls in package.iter("class"):
            relative = cls.get("filename", "")
            hits = {int(line.get("number")): int(line.get("hits"))
                    for line in cls.iter("line")}
            source = None
            for directory in source_dirs:
                source = read_source(directory / relative)
                if source is not None:
                    break
            files.append({"name": cls.get("name", relative),
                          "path": relative, "package": package.get("name", ""),
                          "hits": hits, "source": source})
    return files


def build_data(results_dir):
    tests_path = results_dir / "test-results.xml"
    coverage_path = results_dir / "coverage.xml"
    data = {"generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "tests": None, "coverage": None}
    if tests_path.exists():
        data["tests"] = parse_tests(tests_path)
        data["tests_timestamp"] = datetime.fromtimestamp(
            tests_path.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
    if coverage_path.exists():
        data["coverage"] = parse_coverage(coverage_path)
    return data


PAGE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Test and coverage report</title>
<style>
:root { --pass:#2e7d32; --fail:#c62828; --skip:#ef6c00; --bg:#fafafa;
        --line:#e0e0e0; --muted:#666; }
* { box-sizing:border-box; }
body { font:14px/1.45 system-ui,Segoe UI,sans-serif; margin:0;
       background:var(--bg); color:#222; }
header { background:#1f2933; color:#fff; padding:14px 24px; }
header h1 { margin:0; font-size:18px; }
header small { color:#aab; }
nav { display:flex; gap:4px; padding:0 24px; background:#fff;
      border-bottom:1px solid var(--line); }
nav button { border:0; background:none; padding:12px 16px; cursor:pointer;
             font-size:14px; border-bottom:3px solid transparent; }
nav button.active { border-color:#1976d2; font-weight:600; }
main { padding:20px 24px; max-width:1200px; }
.cards { display:flex; gap:12px; flex-wrap:wrap; margin-bottom:18px; }
.card { background:#fff; border:1px solid var(--line); border-radius:6px;
        padding:10px 16px; min-width:130px; }
.card b { display:block; font-size:22px; }
.card span { color:var(--muted); font-size:12px; }
.toolbar { display:flex; gap:10px; margin-bottom:12px; flex-wrap:wrap; }
input[type=search], select { padding:6px 8px; border:1px solid #bbb;
                             border-radius:4px; font-size:13px; }
table { border-collapse:collapse; width:100%; background:#fff;
        border:1px solid var(--line); }
th, td { text-align:left; padding:6px 10px;
         border-bottom:1px solid var(--line); }
th { background:#f0f2f4; cursor:pointer; user-select:none; }
tr.clickable { cursor:pointer; }
tr.clickable:hover { background:#f5f9ff; }
.pill { display:inline-block; padding:0 8px; border-radius:10px; color:#fff;
        font-size:12px; }
.passed { background:var(--pass); } .failed, .error { background:var(--fail); }
.skipped { background:var(--skip); }
.bar { background:#eee; border-radius:3px; height:10px; width:120px;
       display:inline-block; vertical-align:middle; overflow:hidden; }
.bar i { display:block; height:100%; }
.suite { background:#fff; border:1px solid var(--line); border-radius:6px;
         margin-bottom:8px; }
.suite > summary { padding:8px 12px; cursor:pointer; font-weight:600; }
.suite table { border:0; }
pre.detail { margin:0; padding:8px 12px; background:#fff3f3; color:#7a1010;
             white-space:pre-wrap; font-size:12px; }
.code { background:#fff; border:1px solid var(--line); font:12px/1.5 Consolas,
        monospace; overflow:auto; }
.code div { display:flex; white-space:pre; }
.code .n { width:50px; text-align:right; padding-right:8px; color:#999;
           flex:none; }
.code .h { width:54px; text-align:right; padding-right:8px; color:#555;
           flex:none; }
.code .hit { background:#e8f5e9; } .code .miss { background:#ffebee; }
.muted { color:var(--muted); } a.back { cursor:pointer; color:#1976d2; }
</style></head><body>
<header><h1>Test and coverage report</h1>
<small id="meta"></small></header>
<nav><button id="tab-tests">Unit tests</button>
<button id="tab-cov">Code coverage</button></nav>
<main id="view"></main>
<script>
let DATA = __DATA__;
const WATCH = __WATCH__;
let currentFile;
const view = document.getElementById('view');
const sortState = {};

function el(tag, attrs, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs || {})) {
    if (key === 'class') node.className = value;
    else if (key.startsWith('on')) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value);
  }
  for (const child of children.flat(Infinity)) {
    if (child === null || child === undefined) continue;
    node.append(child.nodeType ? child : document.createTextNode(child));
  }
  return node;
}
const pct = value => (100 * value).toFixed(1) + '%';
const color = rate => rate >= 0.9 ? '#2e7d32' : rate >= 0.75 ? '#ef6c00'
                                                             : '#c62828';
function card(value, label) {
  return el('div', {class: 'card'}, el('b', {}, value), el('span', {}, label));
}
function fileStats(file) {
  const lines = Object.values(file.hits);
  const covered = lines.filter(hit => hit > 0).length;
  return {covered, valid: lines.length,
          rate: lines.length ? covered / lines.length : 1};
}
function bar(rate) {
  return el('span', {class: 'bar'}, el('i', {style:
    `width:${100 * rate}%;background:${color(rate)}`}));
}

function renderTests() {
  const suites = DATA.tests;
  if (!suites) return el('p', {}, 'No test-results.xml found. Run buildtool test.');
  const all = suites.flatMap(s => s.cases.map(c => ({...c, suite: s.name})));
  const count = status => all.filter(c => c.status === status).length;
  const bad = count('failed') + count('error');
  const root = el('div');
  root.append(el('div', {class: 'cards'},
    card(all.length, 'tests'), card(count('passed'), 'passed'),
    card(bad, 'failed / errors'), card(count('skipped'), 'skipped'),
    card(suites.length, 'suites'),
    card(all.reduce((a, c) => a + c.time, 0).toFixed(1) + ' s', 'total time')));

  const search = el('input', {type: 'search', placeholder: 'Filter tests...'});
  const status = el('select', {}, ['all', 'passed', 'failed', 'skipped']
    .map(v => el('option', {value: v}, v)));
  const list = el('div');
  const slow = el('div');
  root.append(el('div', {class: 'toolbar'}, search, status), list,
              el('h3', {}, 'Slowest tests'), slow);

  function draw() {
    const text = search.value.toLowerCase();
    const wanted = c => (status.value === 'all' ||
      (status.value === 'failed' ? c.status === 'failed' || c.status === 'error'
                                 : c.status === status.value)) &&
      (c.name + ' ' + c.suite).toLowerCase().includes(text);
    list.replaceChildren();
    for (const suite of suites) {
      const cases = suite.cases.filter(c => wanted({...c, suite: suite.name}));
      if (!cases.length) continue;
      const failing = suite.cases.some(c => c.status === 'failed' ||
                                            c.status === 'error');
      const details = el('details', {class: 'suite'},
        el('summary', {}, suite.name + ' ',
          el('span', {class: 'pill ' + (failing ? 'failed' : 'passed')},
             `${cases.length}/${suite.cases.length}`),
          el('span', {class: 'muted'}, `  ${suite.time.toFixed(2)} s`)),
        el('table', {}, cases.map(c => [
          el('tr', {}, el('td', {}, el('span', {class: 'pill ' + c.status},
            c.status)), el('td', {}, c.name),
            el('td', {class: 'muted'}, c.time.toFixed(3) + ' s')),
          c.detail ? el('tr', {}, el('td', {colspan: 3},
            el('pre', {class: 'detail'}, c.detail))) : null])));
      if (failing || text || status.value !== 'all') details.open = true;
      list.append(details);
    }
    if (!list.children.length) list.append(el('p', {class: 'muted'}, 'No matches.'));
  }
  const slowest = [...all].sort((a, b) => b.time - a.time).slice(0, 10);
  slow.append(el('table', {}, slowest.map(c => el('tr', {},
    el('td', {}, c.suite), el('td', {}, c.name),
    el('td', {}, c.time.toFixed(3) + ' s')))));
  search.addEventListener('input', draw);
  status.addEventListener('change', draw);
  draw();
  return root;
}

function renderCoverage(selected) {
  const files = DATA.coverage;
  if (!files) return el('p', {}, 'No coverage.xml found. Run buildtool test.');
  if (selected !== undefined) return renderSource(files[selected]);
  const stats = files.map(fileStats);
  const covered = stats.reduce((a, s) => a + s.covered, 0);
  const valid = stats.reduce((a, s) => a + s.valid, 0);
  const root = el('div');
  root.append(el('div', {class: 'cards'},
    card(pct(valid ? covered / valid : 1), 'line coverage'),
    card(`${covered}/${valid}`, 'lines covered'),
    card(valid - covered, 'lines missed'), card(files.length, 'files')));

  const rows = files.map((file, index) => ({index, file, ...stats[index]}));
  const columns = {name: r => r.file.name, rate: r => r.rate,
                   covered: r => r.covered, missed: r => r.valid - r.covered};
  const body = el('tbody');
  function draw() {
    const {key, dir} = sortState.cov || {key: 'rate', dir: 1};
    rows.sort((a, b) => {
      const x = columns[key](a), y = columns[key](b);
      return (x < y ? -1 : x > y ? 1 : 0) * dir;
    });
    body.replaceChildren(...rows.map(r => el('tr', {class: 'clickable',
        onclick: () => show(r.index)},
      el('td', {}, r.file.name), el('td', {}, bar(r.rate), ' ', pct(r.rate)),
      el('td', {}, `${r.covered}/${r.valid}`), el('td', {}, r.valid - r.covered))));
  }
  const header = el('tr', {}, [['name', 'File'], ['rate', 'Coverage'],
    ['covered', 'Covered'], ['missed', 'Missed']].map(([key, label]) =>
    el('th', {onclick: () => {
      const old = sortState.cov;
      sortState.cov = {key, dir: old && old.key === key ? -old.dir : 1};
      draw();
    }}, label)));
  root.append(el('table', {}, el('thead', {}, header), body));
  draw();
  return root;
}

function renderSource(file) {
  const stats = fileStats(file);
  const root = el('div');
  root.append(el('p', {}, el('a', {class: 'back', onclick: () => show()},
    '< All files')),
    el('h3', {}, file.name + '  ', el('span', {class: 'muted'},
      `${pct(stats.rate)} (${stats.covered}/${stats.valid} lines)`)));
  if (!file.source) {
    root.append(el('p', {class: 'muted'}, 'Source file not found on disk.'));
    return root;
  }
  const code = el('div', {class: 'code'});
  file.source.forEach((text, i) => {
    const hits = file.hits[i + 1];
    code.append(el('div', {class: hits === undefined ? '' : hits > 0 ? 'hit' : 'miss'},
      el('span', {class: 'n'}, i + 1),
      el('span', {class: 'h'}, hits === undefined ? '' : hits + 'x'),
      el('span', {}, text || ' ')));
  });
  root.append(code);
  return root;
}

let current = 'tests';
function show(file) {
  currentFile = file;
  document.getElementById('tab-tests').className = current === 'tests' ? 'active' : '';
  document.getElementById('tab-cov').className = current === 'cov' ? 'active' : '';
  view.replaceChildren(current === 'tests' ? renderTests() : renderCoverage(file));
  window.scrollTo(0, 0);
}
document.getElementById('tab-tests').onclick = () => { current = 'tests'; show(); };
document.getElementById('tab-cov').onclick = () => { current = 'cov'; show(); };
function updateMeta() {
  document.getElementById('meta').textContent = 'Generated ' + DATA.generated +
    (DATA.tests_timestamp ? ' from results written ' + DATA.tests_timestamp : '') +
    (WATCH ? ' - live: updates automatically' : '');
}
updateMeta();
show();

// In watch mode, reload report_data.js (a script tag works over file://)
// and redraw in place when the results change.
function poll() {
  const script = document.createElement('script');
  const next = () => { script.remove(); setTimeout(poll, 2000); };
  script.src = 'report_data.js?' + Date.now();
  script.onload = () => {
    const fresh = window.REPORT_DATA;
    if (fresh && fresh.signature !== DATA.signature) {
      DATA = fresh;
      if (currentFile !== undefined && !(DATA.coverage || [])[currentFile])
        currentFile = undefined;
      updateMeta();
      show(currentFile);
    }
    next();
  };
  script.onerror = next;
  document.head.append(script);
}
if (WATCH) poll();
</script></body></html>
"""


RESULT_FILES = ("test-results.xml", "coverage.xml")


def signature(results_dir):
    """Identify the current results by file modification times."""
    parts = []
    for name in RESULT_FILES:
        path = results_dir / name
        parts.append(f"{path.stat().st_mtime_ns}" if path.exists() else "-")
    return "|".join(parts)


def write_report(results_dir, output, watch=False):
    data = build_data(results_dir)
    data["signature"] = signature(results_dir)
    payload = json.dumps(data).replace("</", "<\\/")
    page = PAGE.replace("__WATCH__", "true" if watch else "false")
    output.write_text(page.replace("__DATA__", payload), encoding="utf-8")
    # The open page reloads this file to pick up new results in watch mode.
    output.with_name("report_data.js").write_text(
        "window.REPORT_DATA = " + payload + ";", encoding="utf-8")
    return data


def summarize(data, output):
    tests = [c for s in data["tests"] or [] for c in s["cases"]]
    failed = sum(c["status"] in ("failed", "error") for c in tests)
    print(f"{datetime.now():%H:%M:%S} Wrote {output} "
          f"({len(tests)} tests, {failed} failing)", flush=True)


def watch_results(results_dir, output, interval=1.0):
    """Regenerate the report whenever the result files change."""
    last = None
    print(f"Watching {results_dir} (Ctrl+C to stop)...", flush=True)
    try:
        while True:
            current = signature(results_dir)
            if current != last:
                try:
                    summarize(write_report(results_dir, output, True), output)
                    last = current
                except (ET.ParseError, OSError):
                    pass  # files are mid-write; retry on the next poll
            time.sleep(interval)
    except KeyboardInterrupt:
        print("Stopped.")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--results", type=Path, default=HERE,
                        help="folder with test-results.xml and coverage.xml")
    parser.add_argument("--output", type=Path,
                        help="HTML file to write (default: <results>/report.html)")
    parser.add_argument("--run", action="store_true",
                        help="run 'buildtool test' in MATLAB first")
    parser.add_argument("--open", action="store_true",
                        help="open the report in the default browser")
    parser.add_argument("--watch", action="store_true",
                        help="keep running; the open page updates itself "
                             "whenever new results are written")
    args = parser.parse_args()

    if args.run:
        completed = subprocess.run(["matlab", "-batch", "buildtool test"],
                                   cwd=ROOT)
        if completed.returncode != 0:
            print("buildtool reported failures; building report anyway.")

    output = args.output or args.results / "report.html"
    if args.watch:
        # The page can show a placeholder until the first results arrive.
        if args.open:
            write_report(args.results, output, True)
            webbrowser.open(output.resolve().as_uri())
        watch_results(args.results, output)
        return

    if not any((args.results / name).exists() for name in RESULT_FILES):
        sys.exit(f"No results found in {args.results}. Run buildtool test "
                 "or use --run.")
    summarize(write_report(args.results, output), output)
    if args.open:
        webbrowser.open(output.resolve().as_uri())


if __name__ == "__main__":
    main()
