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
LIB:= -lm -L/Users/annsi118/gslInstalled -lgsl -lgslcblas

# Optimization (-O3 or -O2)
OPT:= -O3

# Directory for executables
BIN_DIR:= .

# Name for exectutable
STATIC_LIB = build/libpcmc.a

# Root directories for .h files
INC_DIR = ./include

# GSL includes (needed only if it is not stadart path: /usr/local/include)
GSL_LDLIBS = -I/Users/annsi118/gslInstalled/include
#GSL_LDLIBS = -I/usr/local/include

# Root directory for .cpp files
SRC_DIR:= ./source

# Directory for Objects
OBJ_DIR:= ./build


#---------------------No need to change the following!-------------------

# Include files are needed fot Includes
INC_FILES:=  $(foreach inc_dir, $(INC_DIR), $(wildcard $(inc_dir)/*h))
INC_FILES += $(foreach inc_dir, $(INC_DIR), $(wildcard $(inc_dir)/*/*h))

# All include directories are needed for Includes
INC_DIRS = $(dir $(INC_FILES))

# Includes: -I ...
LDLIBS:=$(addprefix -I, $(INC_DIRS))

# Source files
#SRC:=  $(foreach src_dir, $(SRC_DIR), $(wildcard $(src_dir)/*h))
#SRC += $(foreach src_dir, $(SRC_DIR), $(wildcard $(src_dir)/*/*h))
SRC :=  $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*/*.cpp) $(wildcard $(SRC_DIR)/*/*/*.cpp)

# Object files
OBJ := $(addprefix $(OBJ_DIR)/,$(notdir $(SRC:.cpp=.o)))

# Files names of objects without directories name
#OBJ_PURE := $(notdir $(OBJ))

all: $(STATIC_LIB)

# Linking
$(STATIC_LIB): $(OBJ)
	@echo "Generating library..." $(notdir $(STATIC_LIB))
#	@$(CC) $(CFLAGS) $(OPT)  $(OBJ) -o $(PROGRAM) $(LIB)
	@ar rcs $(STATIC_LIB) $(OBJ)

# Compiling rule
define app_compile_template
 
 $(1)_OBJ = $$(addprefix $$(OBJ_DIR)/, $$(notdir $$(patsubst %.cpp, %.o, $(1))))
 $(1)_OBJ_PURE = $$(notdir $$($(1)_OBJ))

$$($(1)_OBJ): $(1)
	@echo "Compiling: " $$(addsuffix .cpp, $$(basename $$($(1)_OBJ_PURE)))
	@$$(CC) $$(CFLAGS) $$(OPT) -c $(1) -o $$($(1)_OBJ) $$(LDLIBS) $$(GSL_LDLIBS)
endef

# Compiling
$(foreach app, $(SRC), $(eval $(call app_compile_template,$(app))))

# Cleaning
clean:
	@echo "Cleaning: "
	rm -r build/*

list_of_objects:
	@echo $(OBJ)

# For debugging makefile
print:
	@echo $(LDLIBS)
#	@echo $(INC_FILES)
#	@echo $(addprefix ./object/, $(notdir $(patsubst %.cpp, %.o, $(SRC))))

love:
	@echo "not war?"