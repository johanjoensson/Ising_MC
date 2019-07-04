CXX = g++


# Source files directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

GSLLIBROOT=../GSL-lib
GSLLIBDIR=$(GSLLIBROOT)/lib/GSLpp

WFLAGS = -Wall -Wextra -pedantic -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Weffc++ -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op  -Wuseless-cast
# Flags for the above defined compilers
CXXFLAGS = -std=c++11 $(WFLAGS) -I $(SRC_DIR) -I $(GSLLIBROOT)/include -g #-O3

LDFLAGS = -L$(GSLLIBDIR) -L. -Wl,-rpath=$(GSLLIBDIR) -lGSLpp -lxc -lm -lgsl -g #-O3 -flto

EXE = ising

ISING_OBJ = main.o\


OBJS = $(addprefix $(BUILD_DIR)/, $(ISING_OBJ))
DEPS = $(OBJS:.o=.d)

all: $(EXE)

clean:
	rm -f $(OBJS) $(DEPS)

-include $(DEPS)

$(BUILD_DIR)/%.d: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -MM -MP -MT $(@:.d=.o) $< -MF $@

$(BUILD_DIR)/%.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXE): $(OBJS)
	$(CXX)  $^ -o $@ $(LDFLAGS)

