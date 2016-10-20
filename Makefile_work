# Compilator name (for example could be mpicc instead)
CC:= g++

# Compilator flags
CFLAGS:= -std=c++11

# Libraries (can also include path to these libraries as -L/SOME_PATH)
LIB:= -lm

# Optimization (-O3 or -O2)
OPT:= -O3

# Directory for executables
BIN_DIR:= .

# Name for exectutable
PROGRAM = polymerMC

# Directories for .h files
INC_DIR1 = ./include
INC_DIR2 = ./include/Quantum

# Includes
LDLIBS:= -I$(INC_DIR1) -I$(INC_DIR2)

# Directories for .cpp files
SRC_DIR1 = ./source
SRC_DIR2 = ./source/Quantum

# Source files
SRC =  $(wildcard $(SRC_DIR1)/*.cpp)
SRC += $(wildcard $(SRC_DIR2)/*.cpp)

# Object files
OBJ = $(SRC:.cpp=.o)


all: $(PROGRAM)

# Linking
$(PROGRAM): $(OBJ)
	@echo "Generating executable file..." $(notdir $(PROGRAM))
	@$(CC) $(OPT) $(CFLAGS) $^ -o $(PROGRAM) $(LIB)

# Compiling
%.o: %.cpp
	@echo "Compiling: " $(addsuffix .cpp, $(basename $@))
	@$(CC) $(OPT) $(CFLAGS) -c $< -o $@ $(LDLIBS)

# Cleaning
clean:
	@echo "Cleaning: "
	rm -rf $(OBJ) $(PROGRAM)


#The original extended version:
#https://github.com/latelee/Makefile_templet/blob/master/Makefile_simple
#