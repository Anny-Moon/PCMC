reset
set terminal postscript eps enhanced colour "Times-Roman,38";

n=20 #number of intervals
max=7.0 #max value
min=-7.0 #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0

set output "collapsed_phase_tau.eps";
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/10,max
set boxwidth width*0.9
set style fill solid #fillstyle
set tics out nomirror

set xlabel "{/Symbol t}" font ",50";
set ylabel "Probability distribution" 


#count and plot
plot  "results/Configurations/kappaTauN100.dat" u (hist($2,width)):(1.0) smooth freq w boxes lc rgb "steelblue" notitle
