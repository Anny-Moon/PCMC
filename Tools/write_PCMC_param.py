import sys

fileIn = sys.argv[1]; 
fileOut = "PCMC_parameters.dat";

with open(fileIn, "r") as fp:
    content = fp.read().splitlines()

lastSite = -1;

print("Comments in file '%s':" % sys.argv[1]);
print("\t%s\n" %content[0]);

with open(fileOut, "w") as fp:
    
    fp.write("#Generated automaticaly from data from Vladivostok\n")
    fp.write("#2l86\n")
    fp.write("\n\n")
    fp.write("#----------------Polymer---------------------------\n")
    fp.write("NUMBER_OF_MONOMERS  %i\n" % int(content[1]));
    fp.write("MONOMER_LENGTH 3.8  #fixed for current version of the program\n")
    fp.write("\n")
    fp.write("#----------------Hamiltonian-----------------------\n")
    fp.write("HAM_MU	0.0\n")
    fp.write("HAM_ALPHA	1.0\n")
    fp.write("\n")
    fp.write("#---------------Solitons---------------------------\n")

    
    line = 2;
    while(line < len(content)-1):
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
	fp.write("S_HAM_Q\t%e\n" % c1);
	fp.write("S_HAM_M\t%e\n" % m1);
	
	fp.write("S_HAM_A\t%e\n" % (-d));
	fp.write("S_HAM_B\t%e\n" % (q/d));
	
	fp.write("S_HAM_C\t%e\n" % (2.0*e));
	fp.write("S_HAM_D\t%e\n" % (b/e));
	fp.write("TO\t%i\n" %To);
	fp.write("\n");
	
	From = lastSite + 1;
	To = From + length - cent - 1;
	lastSite = To;
	
	fp.write("FROM\t%i\n" % From);
	fp.write("S_HAM_Q\t%e\n" % c2);
	fp.write("S_HAM_M\t%e\n" % m2);
	
	fp.write("S_HAM_A\t%e\n" % (-d));
	fp.write("S_HAM_B\t%e\n" % (q/d));
	
	fp.write("S_HAM_C\t%e\n" % (2.0*e));
	fp.write("S_HAM_D\t%e\n" % (b/e));
	fp.write("TO\t%i\n" %To);
	fp.write("\n");


    fp.write("#----------------Interaction-----------------------\n")
    fp.write("LENNARD_JONES_MIN	1.1\n");
    fp.write("LENNARD_JONES_R_MIN	5.0\n");
    fp.write("\n");
    fp.write("#---------------Monte Carlo------------------------\n");
    fp.write("REGIME	2\n\n");
    fp.write("LOOPS_PER_CORE	1\n");
    fp.write("\n");
    fp.write("MAX_LOG_T 1.0\n");
    fp.write("LOG_T_STEP 0.5\n");
    fp.write("MIN_LOG_T -13.0\n");
    fp.write("\n");
    fp.write("SWEEPS_PER_STEP 505\n");

print("Generated file '%s'." % fileOut);