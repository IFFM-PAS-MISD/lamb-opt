powershell -Command "(gc SPIE2020.tex) -replace '{', '{ ' | Out-File -encoding UTF8 SPIE2020_L.tex"
powershell -Command "(gc SPIE2020_L.tex) -replace '}', ' }' | Out-File -encoding UTF8 SPIE2020_LR.tex"
del "SPIE2020_L.tex" /f /q
powershell -Command "(gc SPIE2020old.tex) -replace '{', '{ ' | Out-File -encoding UTF8 SPIE2020old_L.tex"
powershell -Command "(gc SPIE2020old_L.tex) -replace '}', ' }' | Out-File -encoding UTF8 SPIE2020old_LR.tex"
del "SPIE2020old_L.tex" /f /q