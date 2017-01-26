#CC := g++
CC:= mpicxx
OPT:= -O3
CFLAGS := -Wall -std=c++11
PROGRAM := pcmc
LIB:= -lm
#Path to GSL if it is not standart (/usr/local/lib )
LIB+= -L/Users/annsi118/gslInstalled/lib -lgsl -lgslcblas

all:
	@echo "Generating executable file..." $(PROGRAM)
	@$(CC) $(OPT) $(CFLAGS) mainMPI_benchmark.cpp -o $(PROGRAM) -lpcmc $(LIB)
#	@$(CC) $(OPT) $(CFLAGS) mainMPI_mac.cpp -o $(PROGRAM) -lpcmc $(LIB)
#	@$(CC) $(OPT) $(CFLAGS) mainAnalyzator.cpp -o $(PROGRAM) -lpcmc $(LIB)
clean:
	rm $(PROGRAM)