#set terminal pdfcairo enh color font "Helvetica,20" lw 1
set terminal pngcairo enhanced font 'Times-New-Roman,20'

set style line 1 lt 1 lc rgb "black" lw 3 pt 1 ps 1
set style line 2 lt 1 lc rgb "red" lw 3 pt 1 ps 1
set style line 3 lt 1 lc rgb "black" lw 3 pt 1 ps 1
set style arrow 1 head front filled size screen 0.03,15,45 ls 3


start = 1
end = 100 
bulk_gamma = -0.1
bulk_gamma_prime = 0.1

g(x)    = x**2*(1-x)**2

h(x) 	= x**2*(3-2*x)

omega(x,gp,g)  = g(x) + h(x)*gp + h(1-x)*g


set xzeroaxis
set ylabel "Energiedichte"
set ytics 1
set xtics 1

set autoscale noextend

set key bottom center

do for [inum = start:end] {
    set output sprintf("potential%04d.png",inum)
    bulk_gamma_prime = bulk_gamma_prime - 0.002
    bulk_gamma = bulk_gamma + 0.002
    set arrow 1 from 0, 0 to 0, bulk_gamma front as 1
    set arrow 2 from 1.0, 0 to 1.0, bulk_gamma_prime front as 1
    plot [-0.1: 1.1] g(x)  title "g({/Symbol f})" with lines ls 1,\
                     g(x) + omega(x,bulk_gamma_prime,bulk_gamma) title sprintf("g + %.2f*h + %.2f*(1-h)",bulk_gamma_prime, bulk_gamma) with lines ls 2 
		
}

quit
