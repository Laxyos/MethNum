set terminal png
if (!exists("filename")) filename='default.png'
set output "img/".filename
set view map
set dgrid3d
set pm3d interpolate 100,100
set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)
splot "data/heat.dat" using 1:2:3 with pm3d t "Température"
