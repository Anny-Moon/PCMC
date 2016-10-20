#compilator name (for example could be mpicc instead)
CC:= g++

CCFLAGS	:= -std=c++11
#LIBRARIES (can also include path to these libraries as -L/SOME_PATH)
LIB:= -lm

#Optimization (-O3 or -O2)
OPT:= -O3

#directory for executables
BIN_DIR:=.

target = a.out
#directory for .h files
INC1 = ./include/quantum
INC2 = 

INCLUDE_DIRS := -I$(INC1) -I$(INC2)
#List of include files

#directory for .cpp files
SRC_DIR = ./source
SRC_DIR1 = 
#List of source files
SRC = $(wildcard $(SRC_DIR)/*.cpp)
SRC +=$(wildcard $(SRC_DIR1)/*.cpp)
#directory for .o files
OBJECT_DIR:=./object

#List of object files
OBJ = $(patsubst %.c,%.o, $(patsubst %.cpp,%.o, $(SRC))) 
all: $(target)

$(target): $(OBJ)

%.o: %.cpp
    @echo "Compiling: " $(addsuffix .cpp, $(basename $(notdir $@)))
    @$(CC) $(CCFLAGS) -c $< -o $@
    
#cleaning
clean:
	rm -rf $(OBJECT_DIR)/*.o $(BIN_DIR)/polymerMC