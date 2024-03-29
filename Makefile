# Source files directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

GSLLIBROOT=../GSL-lib
GSLLIBDIR=$(GSLLIBROOT)/lib/GSLpp

WFLAGS = -Wall -Wextra -pedantic -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference -Wdouble-promotion -Wformat=2 -Weffc++ -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op  -Wuseless-cast
# Flags for the above defined compilers
CXXFLAGS = -std=c++11 $(WFLAGS) -I $(SRC_DIR) -I $(GSLLIBROOT)/include -march=native -Ofast -fopenmp -D_GLIBCXX_PARALLEL

LDFLAGS = -L$(GSLLIBDIR) -L. -Wl,-rpath=$(GSLLIBDIR) -lGSLpp -lm -lgsl -Ofast -flto -fopenmp -D_GLIBCXX_PARALLEL -fuse-ld=gold

EXE = ising

ISING_OBJ = main.o\


OBJS = $(addprefix $(BUILD_DIR)/, $(ISING_OBJ))
DEPS = $(OBJS:.o=.d)

all: $(EXE)

clean:
	@rm -f $(OBJS) $(DEPS)

cleanall : clean
	@rm -f $(EXE)


-include $(DEPS)

$(BUILD_DIR)/%.d: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -MM -MP -MT $(@:.d=.o) $< -MF $@

$(BUILD_DIR)/%.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXE): $(OBJS)
	$(CXX)  $^ -o $@ $(LDFLAGS)

debug : CXXFLAGS = -std=c++11 $(WFLAGS) -I $(SRC_DIR) -I $(GSLLIBROOT)/include -march=native -O0 -g -pg
debug : LDFLAGS = -L$(GSLLIBDIR) -L. -Wl,-rpath=$(GSLLIBDIR) -lGSLpp -lm -lgsl -O0
debug : all

gen-profile : CXXFLAGS += -fprofile-generate
gen-profile : LDFLAGS += -fprofile-generate
gen-profile : all

profile : CXXFLAGS += -fprofile-use
profile : LDFLAGS += -fprofile-use
profile : all
