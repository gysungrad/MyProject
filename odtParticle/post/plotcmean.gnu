set sty dat l;

set term png;

#------- planar means of T for given times

set output 'T_cmeans.png';
set yrange [400:2600];
set xrange [0:1];
set xlabel 'Mixture Fraction';
set ylabel 'T (K)';
plot  \
'cmean_8.dat' us 1:2 ti "t1", \
'cmean_8.dat' us 1:4 ti "t3", \
'cmean_8.dat' us 1:6 ti "t5", \
'cmean_8.dat' us 1:8 ti "t7", \
'cmean_8.dat' us 1:10 ti "t9", \
'cmean_8.dat' us 1:12 ti "t11", \
'cmean_8.dat' us 1:14 ti "t13"

#unset output;
#set term x11;
#rep;
#pause -1;

#------- planar means of Chi for given times

set output 'Chi_cmeans.png';
set auto y;
set xrange [0:1];
set xlabel 'Mixture Fraction';
set ylabel 'Scalar Dissipation Rate (1/s)'; 
plot  \
'cmean_7.dat' us 1:2 ti "t1", \
'cmean_7.dat' us 1:4 ti "t3", \
'cmean_7.dat' us 1:6 ti "t5", \
'cmean_7.dat' us 1:8 ti "t7", \
'cmean_7.dat' us 1:10 ti "t9", \
'cmean_7.dat' us 1:12 ti "t11", \
'cmean_7.dat' us 1:14 ti "t13"
                           
#unset output;
#set term x11;
#rep;
#pause -1;
