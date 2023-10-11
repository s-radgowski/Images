#include <stdio.h>
#include <stdint.h>

#define BIT_SET(str, bit) (((str) >> (7 - bit)) & 1)
#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x80 ? '1' : '0'), \
  (byte & 0x40 ? '1' : '0'), \
  (byte & 0x20 ? '1' : '0'), \
  (byte & 0x10 ? '1' : '0'), \
  (byte & 0x08 ? '1' : '0'), \
  (byte & 0x04 ? '1' : '0'), \
  (byte & 0x02 ? '1' : '0'), \
  (byte & 0x01 ? '1' : '0')

// Define JPEG Constants
typedef uint8_t BYTE;
BYTE max_byte2(const BYTE a, const BYTE b);
BYTE max_byte3(const BYTE a, const BYTE b, const BYTE c);

/*  JPEG File Structure:
    First 2 bytes of a JPEG file are always the following values:
    0xff 0xd8
*/
#define JPEGSignatureSize (2)
#define MAX_COMPRESSION_RATIO (6)

// Check if a signature is valid
int JPEGSignatureValid(FILE *f);

// The CPU will read a uint16_t in the wrong order
uint16_t reverse_bytes(uint16_t bytes);

// Set bit i (0 = left, 7 = right) to 1 at address b
void set_bit(BYTE *b, BYTE i);

// Data structure for holding the results of the crop
struct buffer {
    BYTE *data;
    size_t size;
    size_t capacity;
};
typedef struct buffer *Buffer;

// Write first 2 bytes of a JPEG to a new buffer
void write_signature(Buffer b);

// Add data to the buffer, returning 1 if memory error
#define BUFFER_STARTING_SIZE (20)
Buffer create_buffer(void);
BYTE read_bits(BYTE *buf, BYTE buf_offset, BYTE n);
BYTE write_data(Buffer buf, const BYTE *b, uint32_t n);
BYTE write_bits(Buffer buf, BYTE buf_offset, BYTE *b,
                BYTE b_offset, size_t b_size, BYTE b_end);
void change_bits(Buffer buf, const BYTE *b, BYTE n);
void destroy_buffer(Buffer buf);

/*  JPEG APP0 Segment Layout:
    Marker: 2 bytes 0xff 0xe0
    Size: 2 bytes, should be >=16
    Identifier: 5 bytes, 'JFIF'#0
    Major: 1 byte, should be 1
    Minor: 1 byte, usually 0-2
    Units: 1 byte, x/y densities code
    xdens: 2 bytes, nonzero
    ydens: 2 bytes, nonzero
    Thumbnail width: 1 byte
    Thumbnail height: 1 byte
    Thumbnail: height * width * 3 bytes
*/

// Read and copy APP0 segment from file
BYTE read_header(FILE *f, Buffer b);

// Read and skip APP0 segment from file
int skip_header(FILE *f);

/*  JPEG Chunk Layout:
    type: 0xff + 1 additional byte to identify the segment
    size: 2 bytes for the size of the data + size fields
    data: Contents of the segment, max 65533 bytes
    bd: For SOF0, the bits/sample of the image, otherwise N/A
    height0: For SOF0, the height of the image, otherwise N/A
    width0: For SOF0, the width of the image, otherwise N/A
    restart_intr: For DRI, the restart interval
    num_comps: For SOS, the number of components
    table_id: For SOS, 0 for DC scan and 1 for AC scan
    start: For SOS, MCU index of starting pixel for this scan
    end: For SOS, MCU index of final pixel for this scan
    successive_appr: For SOS, the Progressive method of
        approximation used in this scan
*/
struct chunk {
    BYTE type[2];
    uint16_t size;
    BYTE *data;
    uint8_t bd;
    uint16_t height0;
    uint16_t width0;
    uint16_t restart_intr;
    BYTE num_comps;

    // Currently supports max of 3 comps
    BYTE comps[3];
    BYTE table_ids[3];

    BYTE start;
    BYTE end;
    BYTE successive_appr;
};
typedef struct chunk *Chunk;

// Create a chunk for storing image data
Chunk create_chunk(void);

// Determine if a segment is parameterless
BYTE paramless(const Chunk c);

// Populate a chunk object from the raw image file
BYTE read_chunk(FILE *f, Chunk c);

// Destroy a chunk and free its memory
void destroy_chunk(Chunk c);

/*  Huffman Table Layout:
    type: 0xff, 0xc4 to identify DHT marker
    size: The length of the Huffman table
    info: Bit 0-3 = index of table, bit 4 = type of table
        0 = DC table, 1 = AC table; bit 5-7 unused, must be 0
    counts: Number of symbols with codes of length 1-16
    symbols_len: Size of just the symbols array
    symbols: Table containing the symbols, increasing in length
    next: Pointer to attached additional tables
*/
struct huffman {
    BYTE type[2];
    uint16_t size;
    BYTE info;
    BYTE counts[16];
    uint16_t symbols_len;
    BYTE *symbols;
    struct huffman *next;
};
typedef struct huffman *DHT;

// Populate a Huffman Table from the raw image file
DHT read_huffman(FILE *f);

// Insert a new value into a Huffman Table
BYTE add_huffman(BYTE value, DHT d);

// Write a Huffman Table to a results buffer
BYTE write_huffman(Buffer buf, const DHT d);

// Rewrite a corrected Huffman Table to a results buffer
BYTE rewrite_huffman(Buffer buf, size_t start, const DHT d);

// Destroy a Huffman Table and free its memory
void destroy_huffman(DHT d);

// Nodes and methods for Huffman Trees
struct _node {
    BYTE value;
    char *code;
    size_t length;
    BYTE leaf;
    struct _node *left, *right, *parent;
};
typedef struct _node *Node;

// Populate a Huffman Tree from a DHT chunk
Node create_tree(const DHT d);

// If ever necessary, you can print the tree too:
void print_tree(const Node root) __attribute__((unused));

// Recursively destroy the branches of a Huffman Tree
void destroy_tree(Node n);

/*  Full results from decoding a Huffman bitstring
    value: The 1-byte value that was encoded in input
    bytes: The number of full bytes of input covered
    bits: The additional number of bits above those
        bytes covered in finding a match down the tree
*/
typedef struct decoded {
    BYTE value;
    BYTE bytes;
    BYTE bits;
} decoded;

// Decode a series of bits to a size count using a Huffman Tree
decoded decode_size(const BYTE *str, BYTE bit, const Node root);

// Decode value of n bits from b using bitstream table
int decode_value(size_t n, const BYTE *b);

// Feed value bits into decoded_dc byte array (dc) for analysis
void pull_value(BYTE *dc, decoded out, BYTE *place, size_t offset);

/*  Results from recoding a bitstream with a bitstream table,
    or from running a value back through a Huffman tree
    n: The number of bits used in b
    b: Pointer to the first bit of the recoded result
*/
typedef struct recoded {
    size_t n;
    BYTE b[2];
} recoded;

// Recode an integer into a bitstring using a Huffman Tree
void recode_size(size_t n, const Node root, recoded *r);

// Recode value of n bits from b using bitstream table
void recode_value(int value, recoded *r);

// Consolidate recoded values into one "flex" array
void consolidate(recoded re_size, recoded re_value, BYTE *flex);
