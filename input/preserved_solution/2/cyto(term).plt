set key font "Times Italic,10" tmargin

set ytics nomirror
set xtics nomirror
#set xrange[-0.5:38]
#set yrange[-0.05:0.89]
set border 3

set grid back lt 0
set ylabel "Cytokines (cy)"
plot	'd:/Diploma/source/input/preserved_solution/2/Cy.txt' \
			using 2:3 w lp ps 0.01 lt rgb 'black' title 'solution',\
		'd:/Diploma/source/input/PreservedSol/CY.txt' \
			using 2:3 w lp ps 0.007 lt rgb 'brown' title 'V-spline',\
		'd:/Diploma/source/output/SPL/SPL_CY.txt' \
			using 1:2 w l dt 2 lc rgb 'brown' title 'spline interpolation (Fouda)',\
		'd:/Diploma/source/input/exp/cyto4SPL.txt' \
			using 1:2 with points pointtype 3 lt rgb 'brown' title 'experiment (Fouda)'

