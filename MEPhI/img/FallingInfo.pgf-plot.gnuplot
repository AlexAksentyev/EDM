set table "FallingInfo.pgf-plot.table"; set format "%.5f"
set format "%.7e";; set samples 25; set dummy x; plot [x=0:4*pi] exp(-.11*(x))*sin(deg(x));
