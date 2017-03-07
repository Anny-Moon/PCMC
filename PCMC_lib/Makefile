# 
# Run $make list_of_objects to see
# the rule for compilation only one file.
#
# Anna Sinelnikova
# Uppsala, Sweden 2016

# Compilator name (for example could be mpicc instead)
CC:= g++

# Compilator flags
CFLAGS:= -Wall -std=c++11

# Libraries (can also include path to these libraries as -L/SOME_PATH)
#LIB:= -lm -L/Users/annsi118/gslInstalled -lgsl -lgslcblas

# Optimization (-O3 or -O2)
OPT:= -O3

# Name for exectutable
STATIC_LIB = build/libpcmc.a

# Directories for .h files
INC_DIR = ./include
INC_DIR += ./include/Energy
#INC_DIR += ./include/Quantum
INC_DIR += ./include/Random
INC_DIR += ./include/ReadWriteFiles

# GSL includes (needed only if it is not stadart path: /usr/local/include)
#GSL_LDLIBS = -I/Users/annsi118/gslInstalled/include
GSL_LDLIBS = -I/usr/local/include

# Directory for .cpp files
SRC_DIR = ./source
SRC_DIR += ./source/Energy
#SRC_DIR += ./source/Quantum
SRC_DIR += ./source/Random
SRC_DIR += ./source/ReadWriteFiles
SRC_DIR_CL:= ./source/OpenCL

# Directory for Objects
OBJ_DIR:= ./build


#---------------------No need to change the following!-------------------

# Include files are needed fot Includes
INC_FILES:=  $(foreach inc_dir, $(INC_DIR), $(wildcard $(inc_dir)/*h))

# Includes: -I ...
LDLIBS:=$(addprefix -I, $(INC_DIR))

# Source files
SRC :=  $(foreach src_dir, $(SRC_DIR), $(wildcard $(src_dir)/*cpp))
SRC_CL :=  $(foreach src_dir, $(SRC_DIR_CL), $(wildcard $(src_dir)/*cpp))
#SRC :=  $(wildcard $(SRC_DIR)/*.cpp)
#SCR_CL := $(wildcard $(SRC_DIR_CL)/*.cpp)

# Object files
OBJ := $(addprefix $(OBJ_DIR)/,$(notdir $(SRC:.cpp=.o)))
OBJ_CL := $(addprefix $(OBJ_DIR)/,$(notdir $(SRC_CL:.cpp=.o)))
# Files names of objects without directories name
#OBJ_PURE := $(notdir $(OBJ))

all: $(STATIC_LIB)

# Linking
$(STATIC_LIB): $(OBJ)
	@echo "Generating library..." $(notdir $(STATIC_LIB))
#	@$(CC) $(CFLAGS) $(OPT)  $(OBJ) -o $(PROGRAM) $(LIB)
	@ar rcs $(STATIC_LIB) $(OBJ)


opencl: $(STATIC_LIB) $(OBJ_CL)
	@echo "Generating library..." $(notdir $(STATIC_LIB))
#	@$(CC) $(CFLAGS) $(OPT)  $(OBJ) -o $(PROGRAM) $(LIB)
	@ar rcs $(STATIC_LIB) $(OBJ) $(OBJ_CL)
	
	
# Compiling rule
define app_compile_template
 
 $(1)_OBJ = $$(addprefix $$(OBJ_DIR)/, $$(notdir $$(patsubst %.cpp, %.o, $(1))))
 $(1)_OBJ_PURE = $$(notdir $$($(1)_OBJ))
	

$$($(1)_OBJ): $(1)
	@echo "Compiling: " $$(addsuffix .cpp, $$(basename $$($(1)_OBJ_PURE)))
	@$$(CC) $$(CFLAGS) $$(OPT) -c $(1) -o $$($(1)_OBJ) $$(LDLIBS) $$(GSL_LDLIBS)
#	@$$(CC) $$(CFLAGS) $$(OPT) -c $(1) -o $$(addprefix $$(OBJ_DIR)/, $$(notdir $$(patsubst %.cpp, %.o, $(1)))) $$(LDLIBS) $$(GSL_LDLIBS)
endef

# OpenCL Compiling rule
define app_compile_template_cl

 $(1)_OBJ = $$(addprefix $$(OBJ_DIR)/, $$(notdir $$(patsubst %.cpp, %.o, $(1))))
 $(1)_OBJ_PURE = $$(notdir $$($(1)_OBJ))

$$($(1)_OBJ): $(1)
	@echo "Compiling: " $$(addsuffix .cpp, $$(basename $$($(1)_OBJ_PURE)))
	@$$(CC) $$(CFLAGS) $$(OPT) -c $(1) -o $$($(1)_OBJ) $$(LDLIBS) $$(GSL_LDLIBS)
endef

# Compiling
$(foreach app, $(SRC), $(eval $(call app_compile_template,$(app))))
# OpenCL Compiling
$(foreach app, $(SRC_CL), $(eval $(call app_compile_template_cl,$(app))))

# Cleaning
clean:
	@echo "Cleaning: "
	rm -r build/*

list_of_objects:
	@echo $(OBJ_CL)

# For debugging makefile
print:
#	@echo $(LDLIBS)
	@echo $(SRC_CL)
#	@echo $(addprefix ./object/, $(notdir $(patsubst %.cpp, %.o, $(SRC))))

love:
	@echo "not war?"