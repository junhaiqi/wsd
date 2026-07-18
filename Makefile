CXX = g++
CXXFLAGS = -std=c++11 -static-libstdc++ -Wall -O3 -fopenmp -I./src -w
LDFLAGS = -L./lib -lz 
LIBDIR := -L.

SRCS1 = ./src/wsd.cpp ./src/utils.cpp ./src/main.cpp
OBJS1 = $(SRCS1:.cpp=.o)
TARGET1 = wsd
TEST_TARGET = tests/test_overlap_reconciliation

all: $(TARGET1)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

$(TARGET1): $(OBJS1)
	$(CXX) $(CXXFLAGS) $(OBJS1) -o $(TARGET1) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_TARGET): tests/test_overlap_reconciliation.cpp src/wsd.h src/utils.h src/robin_hood.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean:

	rm -f $(OBJS1) $(OBJS2) $(TARGET1) $(TARGET2) $(TEST_TARGET)
