set sty dat lp;

set term png;

#------- planar means of T for given times

set output 'FWHM_CL_mixf.png';
set xlabel 'Time (s)';
set ylabel 'FWHM (m)';
set y2label 'cL Mixture Fraction';
set y2tics;
set ytics nomirror;
plot  \
'fwhm_odt.dat' us 1:2 ti "FWHM ODT", \
'fwhm_odt.dat' us 1:3 axis x1y2 ti "cL ODT"

#unset output;
#set term x11;
#rep;
#pause -1;
