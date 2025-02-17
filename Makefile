CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3 -fopenmp -I./src -w
LDFLAGS = -L./lib -lz 
LIBDIR := -L.

SRCS1 = ./src/wsd.cpp ./src/utils.cpp ./src/main.cpp
OBJS1 = $(SRCS1:.cpp=.o)
TARGET1 = wsd

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(CXX) $(CXXFLAGS) $(OBJS1) -o $(TARGET1) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS1) $(OBJS2) $(TARGET1) $(TARGET2)