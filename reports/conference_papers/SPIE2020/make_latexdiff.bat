SETLOCAL
set latexdiffpath="C:\Strawberry\perl\bin\latexdiff\latexdiff.pl"
set original=SPIE2020.tex
set revised=SPIE2020_R1.tex
set output=diff.tex
set preamblefile=%~dp0\..\..\latexdiff_helper\diffpreamble.txt

:: %latexdiffpath% --math-markup=0 %original% %revised%>%output%
:: %latexdiffpath% --math-markup=1 --flatten %original% %revised%>%output%

:: preferable style
:: %latexdiffpath% --math-markup=1 --type=CFONT %original% %revised%>%output%

:: %latexdiffpath% --math-markup=1 --type=CULINECHBAR %original% %revised%>%output%
:: %latexdiffpath% --math-markup=1 --type=CFONTCHBAR %original% %revised%>%output%

:: another preferable style
::%latexdiffpath% --math-markup=1 --preamble=%preamblefile% %original% %revised%>%output%

::%latexdiffpath% --type=CFONT SPIE2020.bbl SPIE2020_R1.bbl>diff.bbl
::%latexdiffpath% --math-markup=1 --graphics-markup=0 --preamble=%preamblefile% --disable-citation-markup --flatten --allow-spaces --config="PICTUREENV=(?:picture|DIFnomarkup|align)[\w\d*@]*" %original% %revised%>%output%
%latexdiffpath% --math-markup=1 --preamble=%preamblefile% --disable-citation-markup --flatten --allow-spaces --config="PICTUREENV=(?:picture|DIFnomarkup|align)[\w\d*@]*" %original% %revised%>%output%
ENDLOCAL
