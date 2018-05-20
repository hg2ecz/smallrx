CC=gcc
CFLAGS =-Wall -Ofast -march=native -funroll-loops
#CFLAGS+=-DPRINTFREQZ
LDFLAGS=-lm

OBJS=rx.o
TARGET=rx

all: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
