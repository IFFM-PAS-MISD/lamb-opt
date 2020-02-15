@echo off
SETLOCAL
set windowtext="Bibtex"

start %windowtext% /B bibtex.exe monograph_a.aux
start %windowtext% /B bibtex.exe monograph_b.aux
start %windowtext% /B bibtex.exe journal_c.aux
start %windowtext% /B bibtex.exe journal_d.aux
start %windowtext% /B bibtex.exe conference_e.aux
start %windowtext% /B bibtex.exe conference_f.aux
start %windowtext% /B bibtex.exe report_g.aux
start %windowtext% /B bibtex.exe submitted_h.aux
ENDLOCAL