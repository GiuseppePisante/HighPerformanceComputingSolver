unset border; unset tics; unset key;
set term gif animate delay 30
set output "trace.gif"
set xrange [0:1]
set yrange [0:1]

do for [ts=0:120] {
    plot "particles_".ts.".dat" with points pointtype 7
}
unset output
