CC = gcc
CFLAGS = -O2 

all: companion

companion: companion.o
	$(CC) $(CFLAGS) -o companion companion.o -lm

companion_pub.o: Makefile

clean:
	-rm -f *.o
