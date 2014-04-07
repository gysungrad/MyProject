set sty dat l;

set term png;

#------- planar sig of T for given times

set output 'T_sig.png';
set yrange [0:1000];
set xlabel 'Position (m)';
set ylabel 'RMS T (K)';
plot  \
'sig_8.dat' us 1:2 ti "t1", \
'sig_8.dat' us 1:4 ti "t3", \
'sig_8.dat' us 1:6 ti "t5", \
'sig_8.dat' us 1:8 ti "t7", \
'sig_8.dat' us 1:10 ti "t9", \
'sig_8.dat' us 1:12 ti "t11", \
'sig_8.dat' us 1:14 ti "t13"

#unset output;
#set term x11;
#rep;
#pause -1;

#------- planar sig of mixf for given times

set output 'Mixf_sig.png';
set yrange [0:1];
set xlabel 'Position (m)';
set ylabel 'RMS Mixture Fraction';
plot  \
'sig_6.dat' us 1:2 ti "t1", \
'sig_6.dat' us 1:4 ti "t3", \
'sig_6.dat' us 1:6 ti "t5", \
'sig_6.dat' us 1:8 ti "t7", \
'sig_6.dat' us 1:10 ti "t9", \
'sig_6.dat' us 1:12 ti "t11", \
'sig_6.dat' us 1:14 ti "t13"

#------- 

#unset output;
#set term x11;
#rep;
#pause -1;

#------- planar sig of chi for given times

set output 'Chi_sig.png';
set auto y;
set xlabel 'Position (m)';
set ylabel 'RMS Scalar Dissipation Rate (1/s)';
plot  \
'sig_7.dat' us 1:2 ti "t1", \
'sig_7.dat' us 1:4 ti "t3", \
'sig_7.dat' us 1:6 ti "t5", \
'sig_7.dat' us 1:8 ti "t7", \
'sig_7.dat' us 1:10 ti "t9", \
'sig_7.dat' us 1:12 ti "t11", \
'sig_7.dat' us 1:14 ti "t13"

#------- 

#unset output;
#set term x11;
#rep;
#pause -1;
