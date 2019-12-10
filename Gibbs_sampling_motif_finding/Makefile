CXX = g++

CPPFLAGS = -c -std=c++0x -Wall 

SRCS = main.cpp random_motif.cpp SC_random.cpp
OBJS = $(subst .cpp,.o,$(SRCS))

all: motif

motif: $(OBJS)
	$(CXX) $(OBJS) -o motif

main.o: main.cpp
	$(CXX) $(CPPFLAGS) main.cpp

random_motif.o: random_motif.cpp	
	$(CXX) $(CPPFLAGS) random_motif.cpp

SC_random.o: SC_random.cpp
	$(CXX) $(CPPFLAGS) SC_random.cpp
      
.PHONY: clean

clean:
	rm -rf *o motif
