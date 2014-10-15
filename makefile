CC := g++
CPPFLAGS := -std=c++11 -g -pthread
TARGET := afin
OBJS := afin.o contig.o process.o read.o revcomp.o afin_util.o mismatch.o log.o gzip.o

$(TARGET): $(OBJS)
	$(CC) $(CPPFLAGS) -o $@ $(OBJS) -lz

test: $(OBJS)
	$(CC) $(CPPFLAGS) -o $@ $(OBJS) -lz

%.o : %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@ 

clean:
	/bin/rm -f *o $(TARGET)

