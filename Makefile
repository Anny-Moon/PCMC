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


SOURCES = $(wildcard source/*.cpp)\
	  $(wildcard source/Quantum/*.cpp)
	  
OBJECTS = $(SOURCES:.cpp=.o)

#polymerMC : $(OBJECTS)
#	$(CC) $(OPT) -o  $@ $^ $(LIB)

polymerMC : $(OBJECTS)
	$(CC) $(OPT) -o polymerMC $(OBJECTS)

#cleaning
clean:
	rm -rf $(OBJECT_DIR)/*.o $(BIN_DIR)/polymerMC