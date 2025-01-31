unset border; unset tics; unset key;
set term gif animate delay 10
set output "trace.gif"
set xrange [0:30]
set yrange [0:8]
set size ratio -1
set object 1 circle front at 5.0,4.0 size 1.0 fillcolor rgb "black" lw 2


do for [ts=0:500] {
    plot "particles_".ts.".dat" with points pointtype 7 pointsize 0.3
}
unset output
