CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra

# Path to RFL repository root
REPO_ROOT ?= ../..

INCLUDES = -I$(REPO_ROOT)/src/core
LDFLAGS = -L$(REPO_ROOT)/build/src/core
LDLIBS = -lrfl_core -larmadillo

TARGET = main
SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	@if [ ! -f $(REPO_ROOT)/build/src/core/librfl_core.a ]; then \
		echo "Building librfl_core.a first via CMake..."; \
		cmake -B $(REPO_ROOT)/build -S $(REPO_ROOT); \
		cmake --build $(REPO_ROOT)/build --target rfl_core; \
	fi
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) $(LDFLAGS) $(LDLIBS) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) main_gcc $(REPO_ROOT)/build/examples/cpp/main_gcc

.PHONY: all run clean
