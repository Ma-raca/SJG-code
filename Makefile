CXX = g++
CXXFLAGS = -O3 -fopenmp -std=c++17 -Wall
SOURCES = src/main.cpp src/io.cpp src/metric.cpp src/search.cpp src/build.cpp
.PHONY: all clean test
all: sjg
sjg: $(SOURCES) src/sjg.h
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $@
test: sjg
	python3 tests/smoke.py ./sjg
clean:
	rm -f sjg
