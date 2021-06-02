set key font "Times Italic,10" tmargin

set ytics nomirror
set xtics nomirror
set xrange[-0.5:38]
set yrange[-0.05:0.89]
set border 3

set grid back lt 0
set ylabel "Necrotic cells density (N)"
plot 	'Necr.txt' \
			using 2:3 w lp ps 0.01 lt rgb 'black' title 'solution',\
		'd:/Diploma/source/input/PreservedSol/Necr.txt' \
			using 2:3 w lp ps 0.007 lt rgb 'brown' title 'V-spline',\
		'd:/Diploma/source/input/exp/accurate_necr.txt' \
			using 1:2 with points pointtype 3 lt rgb 'brown' title 'experiment (Garcia)',\
		'd:/Diploma/source/output/SPL/SPL_Necr.txt' \
			using 1:2 w l dt 2  lc rgb 'brown' title 'spline interpolation (Garcia)'