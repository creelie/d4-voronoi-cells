#!/usr/bin/env python3
"""verify_dois.py -- check every DOI in the bibliography of D4.tex against its registry.

For each \\bibitem with a doi.org link, the DOI is resolved through Crossref
(api.crossref.org) or, for DataCite prefixes (arXiv 10.48550, Zenodo 10.5281,
4TU 10.4121, TU Delft 10.4233, Preprints.org 10.20944 and others not in
Crossref), through DataCite (api.datacite.org). The registered title is compared
with the title in the bibitem after both are reduced to lower-case letters and
digits. Standard library only.

Usage:  python3 verify_dois.py [path/to/D4.tex]
Exit status 0 when every DOI resolves and every title matches.
"""
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from difflib import SequenceMatcher

TEX = sys.argv[1] if len(sys.argv) > 1 else 'D4.tex'
UA = 'd4-voronoi-cells-doi-check/1.0 (mailto:unknown@example.org)'


def norm(s):
    s = re.sub(r'\\[a-zA-Z]+\s*', ' ', s)          # drop TeX commands
    s = re.sub(r'[{}$\\"\'`^~]', '', s)
    return re.sub(r'[^a-z0-9]+', ' ', s.lower()).strip()


def bibitems(tex):
    body = tex[tex.index(r'\begin{thebibliography}'):tex.index(r'\end{thebibliography}')]
    for chunk in re.split(r'\\bibitem', body)[1:]:
        key = re.match(r'(?:\[[^\]]*\])?\{([^}]*)\}', chunk).group(1)
        doi = re.search(r'https://doi\.org/([^}\s]+)\}', chunk)
        title = re.search(r'\\textit\{((?:[^{}]|\{[^{}]*\})*)\}', chunk)
        yield key, doi.group(1) if doi else None, title.group(1) if title else ''


def fetch(url):
    req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def registered_title(doi):
    q = urllib.parse.quote(doi, safe='')
    try:
        m = fetch('https://api.crossref.org/works/' + q)['message']
        return 'Crossref', (m.get('title') or [''])[0]
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
    a = fetch('https://api.datacite.org/dois/' + q)['data']['attributes']
    return 'DataCite', (a.get('titles') or [{}])[0].get('title', '')


def main():
    tex = open(TEX, encoding='utf-8').read()
    bad = 0
    rows = list(bibitems(tex))
    for key, doi, title in rows:
        if not doi:
            print('%-10s  no DOI' % key)
            continue
        try:
            reg, rt = registered_title(doi)
        except Exception as e:                      # network or registry failure
            bad += 1
            print('%-10s  %-45s  UNRESOLVED (%s)' % (key, doi, e))
            continue
        ratio = SequenceMatcher(None, norm(title), norm(rt)).ratio()
        ok = ratio > 0.8 or norm(rt) in norm(title) or norm(title) in norm(rt)
        bad += not ok
        print('%-10s  %-45s  %-8s  %s  %.2f  %s' % (key, doi, reg, 'OK      ' if ok else 'MISMATCH', ratio, rt[:70]))
        time.sleep(0.3)
    print('\n%d entries, %d with a DOI, %d problems' % (len(rows), sum(1 for r in rows if r[1]), bad))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
