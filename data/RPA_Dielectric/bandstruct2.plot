#!/usr/bin/gnuplot -persist
set xtics ( "Gamma" 0,  "Mprime" 60 )
unset key
nRows = real(system("awk '$1==\"kpoint\" {nRows++} END {print nRows}' bandstruct2.kpoints"))
nCols = real(system("wc -c < bandstruct2.eigenvals")) / (8*nRows)
formatString = system(sprintf("echo '' | awk 'END { str=\"\"; for(i=0; i<%d; i++) str = str \"%%\" \"lf\"; print str}'", nCols))
plot for [i=1:nCols] "bandstruct2.eigenvals" binary format=formatString u 0:i w l
