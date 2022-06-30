b = 1000
reset
#set terminal postscript enhanced eps  font 'Helvetica, 20'
#set style data lines
#set terminal postscript enhanced 
set border linewidth 1.5
set grid
#set parametric
set nokey
set key left # bottom 
#set xrange[0.0:1.01]
#set yrange[-0.01:0.05]  
#v= sprintf("%g",a)
#set output 'test_shear_'.v.'.eps'   
set multiplot layout 2, 1 columnsfirst title " " 
set xlabel 'x'
set ylabel "h"
plot './resu/1Dtest.out' u 1:3 every:::a::a  lt rgb 'blue'  title 'Numerical results'
set ylabel "u"
plot './resu/1Dtest.out' u 1:4 every:::a::a  lt rgb 'red'    
#unset multiplot
pause -1
a=a+1
if(a<b) reread


