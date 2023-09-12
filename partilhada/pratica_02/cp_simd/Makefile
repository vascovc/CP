
CFLAGS   = -Wall
CXXFLAGS = -Wall

.PHONY: all clean

all: sumarray sumarray_copy sumarray_Any

sumarray: sumarray.o
	g++ $(CFLAGS) -o $@ $^

sumarray_copy: sumarray_copy.o
	g++ $(CFLAGS) -o $@ $^

sumarray_Any: sumarray_Any.o
	g++ $(CFLAGS) -o $@ $^

clean:
	rm -f *.o sumarray

