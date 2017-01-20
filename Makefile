CC := g++
OPT:= -O3
CFLAGS := -Wall -std=c++11
PROGRAM := pcmc
LIB:= -lm -lgsl -lgslcblas

all:
	@echo "Generating executable file..." $(PROGRAM)
	@$(CC) $(OPT) $(CFLAGS) main.cpp -o $(PROGRAM) -I./PCMC_lib/include  -L./PCMC_lib/build  -lpcmc $(LIB)

clean:
	rm $(PROGRAM)