#compilator name (for example could be mpicc instead)
CC:=g++ -std=c++11

#LIBRARIES (can also include path to these libraries as -L/SOME_PATH)
LIB:=-lm

#Optimization (-O3 or -O2)
OPT:=-O3

#directory for .o files
OBJECT_DIR:=./object

#directory for executables
BIN_DIR:=.

#directory for .cpp files
SOURCE_DIR:=./source

#directory for .h files
INCLUDE_DIR:=./include

polymerMC: main.o Vector.o Utilities.o Polymer.o PolymerMC.o PolymerQuantum.o StepFunctionCalculator.o TrancatedExpCalculator.o YukawaCalculator.o
	$(CC) $(OPT) $(OBJECT_DIR)/main.o $(OBJECT_DIR)/Vector.o $(OBJECT_DIR)/Utilities.o $(OBJECT_DIR)/Polymer.o \
	$(OBJECT_DIR)/PolymerMC.o  $(OBJECT_DIR)/Quantum/PolymerQuantum.o \
	$(OBJECT_DIR)/Quantum/StepFunctionCalculator.o $(OBJECT_DIR)/Quantum/TrancatedExpCalculator.o $(OBJECT_DIR)/Quantum/YukawaCalculator.o \
	-o $(BIN_DIR)/polymerMC $(LIB)
    
main.o: $(SOURCE_DIR)/main.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/main.cpp  -o  $(OBJECT_DIR)/main.o  -I $(INCLUDE_DIR)
    
Vector.o: $(SOURCE_DIR)/Vector.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/Vector.cpp -o $(OBJECT_DIR)/Vector.o -I $(INCLUDE_DIR)
    
Utilities.o: $(SOURCE_DIR)/Utilities.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/Utilities.cpp   -o $(OBJECT_DIR)/Utilities.o -I $(INCLUDE_DIR)

Polymer.o: $(SOURCE_DIR)/Polymer.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/Polymer.cpp -o $(OBJECT_DIR)/Polymer.o -I $(INCLUDE_DIR)

PolymerMC.o: $(SOURCE_DIR)/PolymerMC.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/PolymerMC.cpp -o $(OBJECT_DIR)/PolymerMC.o -I $(INCLUDE_DIR)

PolymerQuantum.o: $(SOURCE_DIR)/Quantum/PolymerQuantum.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/Quantum/PolymerQuantum.cpp -o $(OBJECT_DIR)/Quantum/PolymerQuantum.o -I $(INCLUDE_DIR)/Quantum

StepFunctionCalculator.o: $(SOURCE_DIR)/Quantum/StepFunctionCalculator.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/Quantum/StepFunctionCalculator.cpp -o $(OBJECT_DIR)/Quantum/StepFunctionCalculator.o -I $(INCLUDE_DIR)/Quantum

TrancatedExpCalculator.o: $(SOURCE_DIR)/Quantum/TrancatedExpCalculator.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/Quantum/TrancatedExpCalculator.cpp -o $(OBJECT_DIR)/Quantum/TrancatedExpCalculator.o -I $(INCLUDE_DIR)/Quantum

YukawaCalculator.o: $(SOURCE_DIR)/Quantum/YukawaCalculator.cpp
	$(CC) $(OPT) -c $(SOURCE_DIR)/Quantum/YukawaCalculator.cpp -o $(OBJECT_DIR)/Quantum/YukawaCalculator.o -I $(INCLUDE_DIR)/Quantum


#cleaning
clean:
	rm -rf $(OBJECT_DIR)/*.o $(BIN_DIR)/polymerMC