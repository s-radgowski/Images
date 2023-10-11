#include <stdio.h>
#include <stdint.h>

// Define PNG Constants
#define IHDR_SIZE (25)
#define MAX_COMPRESSION_RATIO (6)
#define MAX_IDAT_SIZE (16384)
typedef uint8_t BYTE;

/*  PNG File Structure:
    First 8 bytes of a PNG file are always the following values:
    137 80 78 71 13 10 26 10
*/
#define PNGSignatureSize (8)

// Check if a signature is valid
int PNGSignatureValid(FILE *f);

/*  IHDR (first chunk in a PNG) contains:
    Values for width/height in pixels;
    Bit depth: number of bits per sample;
    Color type: 0 for grayscale, 2 for RGB triples,
        3 for palette index with a PLTE chunk,
        4 for grayscale samples with alpha samples,
        6 for RGB triples with alpha samples;
    Compression method: currently can only be 0;
    Filter method; and Interlace method.
*/
typedef struct {
    uint32_t length;
    char type[4];
    uint32_t width;
    uint32_t height;
    BYTE bitDepth;
    BYTE colorType;
    BYTE compressionMethod;
    BYTE filterMethod;
    BYTE interlaceMethod;
    uint32_t crc;
} __attribute__((__packed__))
IHDR;

// The CPU will read the uint32_t in the wrong order
uint32_t reverse_bytes(uint32_t b);

// Check if bit depth is valid for given color
int BDValid(const BYTE bd, const BYTE color);

// Find samples per pixel by color
BYTE samples_by_color(const BYTE color);

// Data structure for holding the results of the crop
struct buffer {
    BYTE *data;
    size_t size;
    size_t capacity;
};
typedef struct buffer *Buffer;

// Add data to the buffer
#define BUFFER_STARTING_SIZE (5)
Buffer create_buffer(void);
void destroy_buffer(Buffer buf);
Buffer write_data(Buffer buf, BYTE *b, uint32_t n);
Buffer write_data_at(Buffer buf, size_t start, BYTE *b, uint32_t n);

// Write first 8 bytes of a PNG to a new buffer
void write_signature(Buffer b);

/* PNG Chunk Layout:
    Length gives the number of bytes in the data field, <= 2^31 -1
    Chunk Type: restricted to uppercase and lowercase ASCII letters;
    Chunk Data: can be of length zero;
    Cyclic Redundancy Check: always present, calculated on all after length.
*/
struct chunk {
    uint32_t length;
    char type[4];
    BYTE *data;
    uint32_t crc;
};
typedef struct chunk *Chunk;

// Create a chunk for storing image data
Chunk create_chunk(void);

// Populate a chunk object from the raw image file
void read_chunk(FILE *f, Chunk c);

// Destroy a chunk and free its memory
void destroy_chunk(Chunk c);

// Determine if a chunk is safe to copy
int safe_to_copy(const Chunk c);

/*  Applies one of five PNG filters to a pixel, based on the 
    preceding pixels in the original image.
    fpx = pointer to the destination filtered pixel
    type = the filter used, as an int (options described below)
    old_data = the original image data as raw bytes
    ind = index of the first byte of the unfiltered pixel in old_data
    status = {pixel index, line index}
    vals = {x, y, width0, bpp}
    h = 1D array of all raw byte values for the image's pixels

    Type   Method
    0      None
    1      Sub
    2      Up
    3      Average
    4      Paeth

    For all filters, bytes to the left of the first pixel or above the
    first row of pixels are all assumed to be 0s.

    Writes the correctly filtered pixel to fpx, returns nothing
*/
void png_filt(BYTE *fpx, const BYTE type, BYTE *old_data, const size_t ind, 
              const int status[2], const size_t vals[4], BYTE *h);
