CC := g++
OPT:= -O3
CFLAGS := -Wall -std=c++11
PROGRAM := converter

all:
	@echo "Generating executable file..." $(PROGRAM)
	@$(CC) $(OPT) $(CFLAGS) main.cpp -o $(PROGRAM)  -lpcmc -lm -lgsl -lgslcblas

clean:
	rm $(PROGRAM)