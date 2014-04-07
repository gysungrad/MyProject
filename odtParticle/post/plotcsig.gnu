set sty dat l;

set term png;

#------- planar sig of T for given times

set output 'T_csig.png';
set yrange [0:1000];
set xrange [0:1];
set xlabel 'Mixture Fraction';
set ylabel 'RMS T (K)';
plot  \
'csig_8.dat' us 1:2 ti "t1", \
'csig_8.dat' us 1:4 ti "t3", \
'csig_8.dat' us 1:6 ti "t5", \
'csig_8.dat' us 1:8 ti "t7", \
'csig_8.dat' us 1:10 ti "t9", \
'csig_8.dat' us 1:12 ti "t11", \
'csig_8.dat' us 1:14 ti "t13"

#unset output;
#set term x11;
#rep;
#pause -1;

#------- planar sig of T for given times

set output 'Chi_csig.png';
set auto y;
set xrange [0:1];
set xlabel 'Mixture Fraction';
set ylabel 'RMS Scalar Dissipation Rate (1/s)';
plot  \
'csig_7.dat' us 1:2 ti "t1", \
'csig_7.dat' us 1:4 ti "t3", \
'csig_7.dat' us 1:6 ti "t5", \
'csig_7.dat' us 1:8 ti "t7", \
'csig_7.dat' us 1:10 ti "t9", \
'csig_7.dat' us 1:12 ti "t11", \
'csig_7.dat' us 1:14 ti "t13"

#unset output;
#set term x11;
#rep;
#pause -1;
