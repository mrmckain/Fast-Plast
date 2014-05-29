OBJS = afin.o

afin: $(OBJS)
	g++ -std=c++11 -o $@ $(OBJS)

afin.o: afin.cpp
	g++ -std=c++11 -c afin.cpp


