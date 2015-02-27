set arrow from 1,1.07 to 4,1.07 nohead lt 3 lw 10
set arrow from 5,1.09 to 27,1.09 nohead lt 1 lw 40
set arrow from 28,1.11 to 91,1.11 nohead lt 4 lw 10
set arrow from 92,1.09 to 114,1.09 nohead lt 1 lw 40
set arrow from 115,1.07 to 150,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for gi|298286506|ref|NP_002090.4|"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:150]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_10771/gi_298286506_ref_NP_002090.4_.eps"
plot "./TMHMM_10771/gi_298286506_ref_NP_002090.4_.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
