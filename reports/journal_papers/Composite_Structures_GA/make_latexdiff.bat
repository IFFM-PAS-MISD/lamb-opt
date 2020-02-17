
SETLOCAL
set latexdiffpath="C:\Strawberry\perl\bin\latexdiff\latexdiff.pl"
set original=Identification_GA.tex
set revised=Identification_GA_R1.tex
set output=diff.tex

%latexdiffpath% --math-markup=0 %original% %revised%>%output%

ENDLOCAL
