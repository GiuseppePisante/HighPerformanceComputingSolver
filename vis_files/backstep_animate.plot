unset border; unset tics; unset key;
set term gif animate delay 10
set output "trace.gif"
set xrange [0:7]
set yrange [0:1.5]
set size ratio -1

set object 1 rect from 0.0,0.0 to 1.0,0.5 lw 5


do for [ts=0:300] {
    plot "particles_".ts.".dat" with points pointtype 7
}
unset output
