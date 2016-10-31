CC := g++
CPPFLAGS := -std=c++11 -g -pthread
TARGET := afin
OBJS := afin.o contig.o process.o read.o revcomp.o mismatch.o log.o gzip.o fusion.o match.o readlist.o extension.o contiglist.o

$(TARGET): $(OBJS)
	$(CC) $(CPPFLAGS) -o $@ $(OBJS) -lz

test: $(OBJS)
	$(CC) $(CPPFLAGS) -o $@ $(OBJS) -lz

%.o : %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@ 

clean:
	/bin/rm -f *o $(TARGET)

