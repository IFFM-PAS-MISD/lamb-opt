powershell -Command "(gc SPIE2020_LR_wdiff.tex) -replace '{ ', '{' | Out-File -encoding UTF8 SPIE2020_L_wdiff.tex"
powershell -Command "(gc SPIE2020_L_wdiff.tex) -replace ' }', '}' | Out-File -encoding UTF8 SPIE2020_wdiff.tex"
del "SPIE2020_L_wdiff.tex" /f /q
del "SPIE2020_LR_wdiff.tex" /f /q