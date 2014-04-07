set sty dat l;

set term png;

#------- planar means of T for given times

set output 'T_means.png';
set yrange [400:2600];
set xlabel 'Position (m)';
set ylabel 'T (K)';
plot  \
'means_8.dat' us 1:2 ti "t1", \
'means_8.dat' us 1:4 ti "t3", \
'means_8.dat' us 1:6 ti "t5", \
'means_8.dat' us 1:8 ti "t7", \
'means_8.dat' us 1:10 ti "t9", \
'means_8.dat' us 1:12 ti "t11", \
'means_8.dat' us 1:14 ti "t13"

#unset output;
#set term x11;
#rep;
#pause -1;

#------- planar means of mixf for given times

set output 'Mixf_means.png';
set yrange [0:1];
set xlabel 'Position (m)';
set ylabel 'Mixture Fraction';
plot  \
'means_6.dat' us 1:2 ti "t1", \
'means_6.dat' us 1:4 ti "t3", \
'means_6.dat' us 1:6 ti "t5", \
'means_6.dat' us 1:8 ti "t7", \
'means_6.dat' us 1:10 ti "t9", \
'means_6.dat' us 1:12 ti "t11", \
'means_6.dat' us 1:14 ti "t13"

#------- 

#unset output;
#set term x11;
#rep;
#pause -1;

#------- planar means of chi for given times

set output 'Chi_means.png';
set auto y;
set xlabel 'Position (m)';
set ylabel 'Scalar Dissipation Rate (1/s)';
plot  \
'means_7.dat' us 1:2 ti "t1", \
'means_7.dat' us 1:4 ti "t3", \
'means_7.dat' us 1:6 ti "t5", \
'means_7.dat' us 1:8 ti "t7", \
'means_7.dat' us 1:10 ti "t9", \
'means_7.dat' us 1:12 ti "t11", \
'means_7.dat' us 1:14 ti "t13"

#------- 

#unset output;
#set term x11;
#rep;
#pause -1;
