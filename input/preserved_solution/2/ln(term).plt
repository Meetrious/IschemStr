set key font "Times Italic,10" tmargin

set ytics nomirror
set xtics nomirror
#set xrange[-0.5:38]
#set yrange[-0.05:0.89]
set border 3

set grid back lt 0
set ylabel "LeuNeutrophils (ln)"
plot	'd:/Diploma/source/input/preserved_solution/2/Ln.txt' \
			using 2:3 w lp ps 0.01 lt rgb 'black' title 'solution',\

