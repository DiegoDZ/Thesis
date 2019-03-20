dt = 0.004
set nokey
xarrow = 0.15
set arrow from xarrow,-0.001 to xarrow,0.006 nohead lc 1
plot for [i=1:50] 'EigenvaluesCt-WALLS.dat' u ($0*dt):i w l  lt 1  lc i lw 1
pause -1

