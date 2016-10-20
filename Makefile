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

OBJ_DIR = ./object
# Root directory for .cpp files
SRC_DIR = ./source
#SRC_DIR2 = ./source/Quantum

# Source files
SRC =  $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/*/*.cpp)
#SRC += $(wildcard $(SRC_DIR2)/*.cpp)

# Object files
#OBJ1 = $(SRC:.cpp=.o)
OBJ = $(SRC:.cpp=.o)
OBJ_WD = $(notdir $(OBJ))
#SRC_FILES_NAMES =  $(notdir$(SRC))
OBJ_NAMES = $(addprefix ./object/,$(notdir $(OBJ)))
#OBJ = $(patsubst $(SCR_DIR), $(OBJ_DIR), $(OBJ1))
NAMES =   $(basename $(SRC))

all: $(PROGRAM)

# Linking
$(PROGRAM): $(OBJ_NAMES)
	@echo "Generating executable file..." $(notdir $(PROGRAM))
	@$(CC) $(CFLAGS) $(OPT)  $(OBJ_NAMES) -o $(PROGRAM) $(LIB)

# Compiling
#$(addprefix ./odject/, $(notdir %.o)): $(addprefix ./source/, $(notdir %.cpp))
#%.o: $(addprefix ./source/, $(notdir %.cpp))
#	@echo "Compiling: " $(addsuffix .cpp, $(basename $@))
#	@$(CC) $(OPT) $(CFLAGS) -c $< -o $@ $(LDLIBS)

#define app_compile_template
# $(1)_NAME = $(1)
# $(1)_CPP_NAME = $(1)
# $(1)_OBJ_NAME =$$(addprefix ./object/,($$(notdir($$(patsubst %.o,%.cpp, $(1))))))
#$$(1)_OBJ_NAME: $$($(1)_CPP_NAME)
#	$$(CC) $$(OPT) -c $$($(1)_CPP_NAME)  -o $$($(1)_OBJ_NAME) $$(LDLIBS)
#endef

#$(foreach app, $(SRC), $(eval $(call app_compile_template,$(app))))

define app_compile_template
 $(1)_NAME = $(1)
 $(1)_CPP_NAME = $(1).cpp
 $(1)_OBJ_NAME = $(1).o
 $(1)_OBJ_FULL_NAME = ./object/$$(notdir $$($(1)_OBJ_NAME))

#$$(notdir $$($(1)_OBJ_NAME)): $$(notdir $$($(1)_CPP_NAME))
$$($(1)_OBJ_FULL_NAME): $$($(1)_CPP_NAME)
	$$(CC) $$(CFLAGS) $$(OPT) -c $$($(1)_CPP_NAME) -o $$($(1)_OBJ_FULL_NAME) $$(LDLIBS)
endef

$(foreach app, $(NAMES), $(eval $(call app_compile_template,$(app))))


# Cleaning
clean:
	@echo "Cleaning: "
	rm -rf $(OBJ) $(PROGRAM)
	rm -rf $(OBJ_NAMES) $(PROGRAM)
print:
	@echo $(OBJ_WD)
#	@echo $(addprefix ./object/, $(notdir ($(patsubst %.o,%.cpp, $(SRC)))))
#The original extended version:
#https://github.com/latelee/Makefile_templet/blob/master/Makefile_simple
#