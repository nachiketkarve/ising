CXX ?= g++
CXXFLAGS ?= -O3 -fopenmp -std=c++17
CPPFLAGS ?= -I./include -I./../Libraries/eigen -I./../Libraries/json/include
LDFLAGS ?= -fopenmp

SRC_DIR := src
OBJ_DIR := build/obj
BIN_DIR := build/bin

ifeq ($(OS),Windows_NT)
	EXE_EXT ?= .exe
	MKDIR_BUILD = cmd /C "if not exist build mkdir build"
	MKDIR_OBJ = cmd /C "if not exist build\obj mkdir build\obj"
	MKDIR_BIN = cmd /C "if not exist build\bin mkdir build\bin"
	RM_BUILD = cmd /C "if exist build rmdir /S /Q build"
else
	EXE_EXT ?=
	MKDIR_BUILD = mkdir -p build
	MKDIR_OBJ = mkdir -p $(OBJ_DIR)
	MKDIR_BIN = mkdir -p $(BIN_DIR)
	RM_BUILD = rm -rf build
endif

CORE_SRC := $(SRC_DIR)/ising.cpp
CORE_OBJ := $(OBJ_DIR)/ising.o

PROGRAM_SRCS := $(filter-out $(CORE_SRC),$(wildcard $(SRC_DIR)/*.cpp))
PROGRAMS := $(basename $(notdir $(PROGRAM_SRCS)))

BINS := $(addprefix $(BIN_DIR)/,$(addsuffix $(EXE_EXT),$(PROGRAMS)))
OBJS := $(CORE_OBJ) $(addprefix $(OBJ_DIR)/,$(addsuffix .o,$(PROGRAMS)))

.PHONY: all clean dirs
.SECONDARY: $(OBJS)

all: dirs $(BINS)

dirs:
	$(MKDIR_BUILD)
	$(MKDIR_OBJ)
	$(MKDIR_BIN)

$(CORE_OBJ): $(CORE_SRC) include/ising.hpp include/eigenClasses.hpp | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp include/ising.hpp include/eigenClasses.hpp | dirs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR)/%$(EXE_EXT): $(OBJ_DIR)/%.o $(CORE_OBJ) | dirs
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $^ -o $@

clean:
	$(RM_BUILD)
