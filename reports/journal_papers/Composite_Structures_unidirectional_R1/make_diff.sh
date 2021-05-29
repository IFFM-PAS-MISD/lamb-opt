#!/bin/bash
OLD_FILE="$HOME/work/projects/opus15/lamb-opt/reports/journal_papers/Composite_Structures_unidirectional/main_Identification_unidirectional.tex"
NEW_FILE="$HOME/work/projects/opus15/lamb-opt/reports/journal_papers/Composite_Structures_unidirectional_R1/main_Identification_unidirectional_R1.tex"
latexdiff $OLD_FILE $NEW_FILE > diff.tex  
