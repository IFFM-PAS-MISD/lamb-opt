
SETLOCAL
set latexdiffpath="C:\Strawberry\perl\bin\latexdiff\latexdiff.pl"
set original=Identification_GA.tex
set revised=Identification_GA_R1.tex
set output=diff.tex
set preamblefile=%~dp0\..\..\latexdiff_helper\diffpreamble.txt

::%latexdiffpath% %original% %revised%>%output%
::%latexdiffpath% --type=CFONT --math-markup=0 --graphics-markup=0 --allow-spaces %original% %revised%>%output%
::%latexdiffpath% --math-markup=1 --preamble=%preamblefile% %original% %revised%>%output%
::%latexdiffpath% --math-markup=0 --preamble=%preamblefile% --type=CFONT --exclude-safecmd="cite,equation,align,figure,label,enumerate" %original% %revised%>%output%
%latexdiffpath% --math-markup=0 --graphics-markup=0 --preamble=%preamblefile% --type=CFONT --exclude-safecmd="cite,equation,align,figure,label,enumerate,modelname" --allow-spaces --config="PICTUREENV=(?:picture|DIFnomarkup|align)[\w\d*@]*" --config "MATHENV=(?:displaymath|equation*)" %original% %revised%>%output%


ENDLOCAL
