CC=gcc
CFLAGS=-std=c99 -g3 -Wall -lz -O
PROGRAMS=png_crop png_read
all: png_crop png_read

png_crop: png_crop.o png.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

png_read: png_read.o png.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

test: png_crop
	./png_crop 10 10 600 400 simu.png simu2.png

png.o: png.c png.h
png_crop.o: png_crop.c png.h
png_read.o: png_read.c png.h

clean:
	$(RM) *.o $(PROGRAMS)
	$(RM) simu2.png
