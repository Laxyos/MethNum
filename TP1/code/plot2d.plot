set view map
set dgrid3d
set pm3d interpolate 10,10
splot "heat.dat" using 1:2:3 with pm3d t "Température"