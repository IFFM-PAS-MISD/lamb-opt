SETLOCAL
set original=SPIE2020.tex
::set output=detex-out.tex
set output="\\tsclient\E\grammarly_web\detex-out.txt"
set detex=%~dp0\..\..\..\bin\external\opendetex\detex.exe

%detex% -l %original%>%output%

ENDLOCAL
