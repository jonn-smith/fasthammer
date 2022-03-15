CXX = g++

.PHONY: all clean

all: fasthammer

fasthammer: 
	$(CXX) -O3 fasthammer.cpp -march=native -lpthread -o fasthammer

debug:
	$(CXX) fasthammer.cpp -march=native -lpthread -o fasthammer

clean:
	rm -f fasthammer

