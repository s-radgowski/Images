CC=gcc
CFLAGS=-std=c99 -g3 -Wall -O
PROGRAMS=jpeg_crop jpeg_read
all: jpeg_crop jpeg_read

jpeg_crop: jpeg_crop.o jpeg.o huffman.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

jpeg_read: jpeg_read.o jpeg.o huffman.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

test: jpeg_crop
	./jpeg_crop 64 16 400 400 simu.jpeg simu2.jpeg

jpeg.o: jpeg.c jpeg.h
jpeg_crop.o: jpeg_crop.c jpeg.h
huffman.o: huffman.c jpeg.h
jpeg_read.o: jpeg_read.c jpeg.h

clean:
	$(RM) *.o $(PROGRAMS)
	$(RM) simu2.jpeg
