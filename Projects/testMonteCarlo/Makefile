# 
# Run $make list_of_objects to see
# the rule for compilation only onle file.
#
# Anna Sinelnikova
# Uppsala, Sweden 2016

# Compilator name (for example could be mpicc instead)
#CC:= g++
CC:= mpicxx
# Compilator flags
CFLAGS:= -std=c++11 -Wall

# Libraries (can also include path to these libraries as -L/SOME_PATH)
#LIB:= -lm -L/Users/annsi118/gslInstalled/lib -lgsl -lgslcblas -lpcmc
LIB:= -lm -lgsl -lgslcblas -lpcmc
# Optimization (-O3 or -O2)
OPT:= -O3

# Directory for executables
BIN_DIR:= .

# Name for exectutable
PROGRAM = pcmc

# Root directories for .h files
INC_DIR:= ./include  

# Root directory for .cpp files
SRC_DIR:= ./source

# Directory for Objects
OBJ_DIR:= ./objects


#---------------------No need to change the following!-------------------

# Include files are needed fot Includes
#INC_FILES:=  $(wildcard $(INC_DIR)/*.h)
INC_FILES:=  $(foreach inc_dir, $(INC_DIR), $(wildcard $(inc_dir)/*h))
INC_FILES += $(foreach inc_dir, $(INC_DIR), $(wildcard $(inc_dir)/*/*h))
# All include directories are needed for Includes
INC_DIRS:= $(dir $(INC_FILES)) 
#INC_DIRS +=/usr/local/include

# Includes: -I ...
LDLIBS:=$(addprefix -I, $(INC_DIRS))

# Source files
SRC :=  $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*/*.cpp) $(wildcard $(SRC_DIR)/*/*/*.cpp)

# Object files
OBJ := $(addprefix $(OBJ_DIR)/,$(notdir $(SRC:.cpp=.o)))

# Files names of objects without directories name
#OBJ_PURE := $(notdir $(OBJ))

all: $(PROGRAM)

# Linking
$(PROGRAM): $(OBJ)
	@echo "Generating executable file..." $(notdir $(PROGRAM))
	@$(CC) $(CFLAGS) $(OPT)  $(OBJ) -o $(PROGRAM) $(LIB)

# Compiling rule
define app_compile_template
 
 $(1)_OBJ = $$(addprefix $$(OBJ_DIR)/, $$(notdir $$(patsubst %.cpp, %.o, $(1))))
 $(1)_OBJ_PURE = $$(notdir $$($(1)_OBJ))

$$($(1)_OBJ): $(1)
	@echo "Compiling: " $$(addsuffix .cpp, $$(basename $$($(1)_OBJ_PURE)))
	@$$(CC) $$(CFLAGS) $$(OPT) -c $(1) -o $$($(1)_OBJ) $$(LDLIBS)
endef

# Compiling
$(foreach app, $(SRC), $(eval $(call app_compile_template,$(app))))

# Cleaning
clean:
	@echo "Cleaning: "
	rm -rf $(OBJ) $(PROGRAM)
	
list_of_objects:
	@echo $(OBJ)

# For debugging makefile
print:
	@echo $(LDLIBS)
#	@echo $(wildcard $(INC_DIR)/*.h)
#	@echo $(addprefix ./object/, $(notdir $(patsubst %.cpp, %.o, $(SRC))))
