name = program
src = $(wildcard src/*.cpp)
obj = $(src:.c=.o)

CC = g++

CFLAGS = -std=c++0x -O3

BOOSTDIR = '/usr/include'

all: $(name)

$(name): $(obj)
	$(CC) $(CFLAGS) -o  $@ $^ -I$(BOOSTDIR) -lrt

run:
	./$(name)

clean:
	rm -f $(name)
