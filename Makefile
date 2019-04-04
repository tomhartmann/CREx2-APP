CC=g++
CFLAGS=-Wall -O3 -g -pg -std=gnu++11 -w

TARGET=crex2

all: $(TARGET)

crex2: crex2.cpp common.o conserved.o costfoo.o counter.o crex.o dstnc.o genom.o helpers.o io.o rearrangements.o rrrmtrand.o sb4type.o tdl.o
	$(CC) $(CFLAGS) -o crex2 common.o conserved.o costfoo.o counter.o crex.o dstnc.o genom.o helpers.o io.o rearrangements.o rrrmtrand.o sb4type.o tdl.o crex2.cpp 

common.o: common.cpp common.hpp
	$(CC) $(CFLAGS) -o common.o -c common.cpp
conserved.o: conserved.cpp conserved.hpp
	$(CC) $(CFLAGS) -o conserved.o -c conserved.cpp
costfoo.o: costfoo.cpp costfoo.hpp
	$(CC) $(CFLAGS) -o costfoo.o -c costfoo.cpp
counter.o: counter.cpp counter.hpp
	$(CC) $(CFLAGS) -o counter.o -c counter.cpp
crex.o: crex.cpp crex.hpp
	$(CC) $(CFLAGS) -o crex.o -c crex.cpp
dstnc.o: dstnc.cpp dstnc.hpp
	$(CC) $(CFLAGS) -o dstnc.o -c dstnc.cpp
genom.o: genom.cpp genom.hpp
	$(CC) $(CFLAGS) -o genom.o -c genom.cpp
helpers.o: helpers.cpp helpers.hpp
	$(CC) $(CFLAGS) -o helpers.o -c helpers.cpp
io.o: io.cpp io.hpp
	$(CC) $(CFLAGS) -o io.o -c io.cpp
rearrangements.o: rearrangements.cpp rearrangements.hpp
	$(CC) $(CFLAGS) -o rearrangements.o -c rearrangements.cpp
rrrmtrand.o: rrrmtrand.cpp rrrmtrand.hpp
	$(CC) $(CFLAGS) -o rrrmtrand.o -c rrrmtrand.cpp
sb4type.o: sb4type.cpp sb4type.hpp
	$(CC) $(CFLAGS) -o sb4type.o -c sb4type.cpp
tdl.o: tdl.cpp tdl.hpp
	$(CC) $(CFLAGS) -o tdl.o -c tdl.cpp

clean:
	rm -f *.o
	rm crex2
