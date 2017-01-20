reset
set terminal postscript eps enhanced colour "Times-Roman,20";

n=200 #number of intervals
max=2.5 #max value
min=-2.5 #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0

set output "collapsed_phase_kappa.eps";
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid #fillstyle
set tics out nomirror
#set ytics 
set xlabel "{/Symbol k}" font ",50";
set ylabel "Probability distribution" 


#count and plot
plot  "results/Conf/kappaTauN100.dat" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb "red" notitle
