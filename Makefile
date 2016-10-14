#compilator name (for example could be mpicc instead)
CC:=g++

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

polymerMC: main.o Vector.o Utilities.o Polymer.o PolymerMC.o
	$(CC) $(OPT) $(OBJECT_DIR)/main.o $(OBJECT_DIR)/Vector.o $(OBJECT_DIR)/Utilities.o $(OBJECT_DIR)/Polymer.o $(OBJECT_DIR)/PolymerMC.o -o $(BIN_DIR)/polymerMC $(LIB)
    
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


#cleaning
clean:
	rm -rf $(OBJECT_DIR)/*.o $(BIN_DIR)/polymerMC