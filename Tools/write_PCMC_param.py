import sys

fileIn = sys.argv[1]; 
fileOut = "PCMC_parameters.dat";

with open(fileIn, "r") as fp:
    content = fp.read().splitlines()

lastSite = -1;

with open(fileOut, "w") as fp:
    line = 0;
    print(len(content))
    while(line < len(content)-1):
	print(line);
	print content[line];
	length, cent = content[line].split();
	length = int(length);
	cent = int(cent);
	line+=1;
	
	From = lastSite + 1;
	To = From + cent - 1;
	lastSite = To;
	
	b = content[line];
	b = float(b);
	line+=1;
	c1 = content[line];
	c1 = float(c1);
	line+=1;
	c2 = content[line];
	c2 = float(c2);
	line+=1;
	d = content[line];
	d = float(d);
	line+=1;
	e = content[line];
	e = float(e);
	line+=1;
	q = content[line];
	q = float(q);
	line+=1;
	m1 = content[line];
	m1 = float(m1);
	line+=1;
	m2 = content[line];
	m2 = float(m2);
	line+=1;
	
	fp.write("FROM\t%i\n" % From);
	fp.write("S_HAM_Q\t%f\n" % c1);
	fp.write("S_HAM_M\t%f\n" % m1);
	
	fp.write("S_HAM_A\t%f\n" % (-d));
	fp.write("S_HAM_B\t%f\n" % (q/d));
	
	fp.write("S_HAM_C\t%f\n" % (2.0*e));
	fp.write("S_HAM_D\t%f\n" % (b/e));
	fp.write("TO\t%i\n" %To);
	fp.write("\n");
	
	From = lastSite + 1;
	To = From + length - cent - 1;
	lastSite = To;
	
	fp.write("FROM\t%i\n" % From);
	fp.write("S_HAM_Q\t%f\n" % c2);
	fp.write("S_HAM_M\t%f\n" % m2);
	
	fp.write("S_HAM_A\t%f\n" % (-d));
	fp.write("S_HAM_B\t%f\n" % (q/d));
	
	fp.write("S_HAM_C\t%f\n" % (2.0*e));
	fp.write("S_HAM_D\t%f\n" % (b/e));
	fp.write("TO\t%i\n" %To);
	fp.write("\n");
	