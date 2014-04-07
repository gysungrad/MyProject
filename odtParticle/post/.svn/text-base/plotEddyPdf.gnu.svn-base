set sty dat lp;

set term png;

#------- planar means of T for given times

set output 'eddyPdf.png';
set xlabel 'Eddy Size (m)';
set ylabel 'PDF';
set log x;
set log y;
plot  \
'eddyPDF.dat' us 1:3 ti 'sampled', 'eddyPDF.dat' us 1:2 ti 'ODT'

#unset output;
#set term x11;
#rep;
#pause -1;
