#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "jpeg.h"

// Function Definitions for jpeg.h

// Return the maximum value of 2 unsigned bytes
BYTE max_byte2(const BYTE a, const BYTE b)
{
    if (a >= b) return a;
    else return b;
}

// Return the maximum value of 3 unsigned bytes
BYTE max_byte3(const BYTE a, const BYTE b, const BYTE c)
{
    if ((a >= b) && (a >= c)) return a;
    else if ((b >= a) && (b >= c)) return b;
    else return c;
}

// JPEG File Structure: First 2 bytes are always 0xff 0xd8
char JPEGSignature[2] = {0xff, 0xd8};

// Return 1 if a JPEG signature is valid, otherwise 0
int JPEGSignatureValid(FILE *f)
{
    char c;
    for (int i = 0; i < JPEGSignatureSize; i++)
    {
        fread(&c, sizeof(char), 1, f);
        if (c != JPEGSignature[i])
        {
            return 0;
        }
    }

    return 1;
}

// The CPU will read a uint16_t in the wrong order
uint16_t reverse_bytes(uint16_t bytes)
{
    uint16_t aux = 0;
    BYTE b;
    
    for (int i = 0; i < 16; i+=8)
    {
        b = (bytes >> i) & 0xff;
        aux |= b << (8 - i);
    }
    return aux;
}

// Set bit i (0 = left, 7 = right) to 1 at address b
void set_bit(BYTE *b, BYTE i)
{
    *b ^= (-1 ^ *b) & (1UL << (7 - i));
}

// Write first 2 bytes of a JPEG to a new buffer
void write_signature(Buffer b)
{
    for (int i = 0; i < JPEGSignatureSize; i++)
    {
        b->data[i] = JPEGSignature[i];
    }

    b->size = JPEGSignatureSize;
}

// Create a buffer
Buffer create_buffer(void)
{
    Buffer buf = malloc(sizeof(struct buffer));
    assert(buf);
    BYTE *data = malloc(BUFFER_STARTING_SIZE);
    assert(data);

    buf->data = data;
    buf->size = 0;
    buf->capacity = BUFFER_STARTING_SIZE;

    return buf;
}

// Read n bits from a buffer's data, returning value read
BYTE read_bits(BYTE *buf, BYTE buf_offset, BYTE n)
{
    BYTE result = 0;
    if (buf_offset + n < 9)
    {
        result = *buf >> (8 - buf_offset - n);
        if (n < 8) result &= (0xff >> (8 - n));
    }

    else if (buf_offset + n < 17)
    {
        result = (buf[1] >> (16 - buf_offset - n));
        result &= (0xff >> (16 - buf_offset - n));
        result += *buf << (buf_offset + n - 8);
        if (n < 8) result &= (0xff >> (8 - n));
    }

    return result;
}

// Add data to a buffer, returning 1 if memory error
#define SIZE_MULTIPLIER (2)
BYTE write_data(Buffer buf, const BYTE *b, uint32_t n)
{
    // Check if you need to expand the buffer
    while (buf->capacity <= buf->size + n)
    {
        buf->capacity *= SIZE_MULTIPLIER;
        buf->data = realloc(buf->data, buf->capacity);
        if (!buf->data) return 1;
    }

    for (int i = 0; i < n; i++)
    {
        buf->data[buf->size++] = b[i];
    }
    return 0;
}

// Write a certain set of offset bits to a buffer
/*  buf = buffer to which you are writing the bits
    buf_offset = # of bits already used in the last byte of
        buf, with 0 signifying a new byte must be added
    b = pointer to first input byte
    b_offset = # of leading irrelevant bits at the start of b
    b_size = # of bytes in b, rounding up leading/trailing bits
    b_end = # of trailing irrelevant bits at the end of b, 
        with 0 signifying the last byte is complete
*/
BYTE write_bits(Buffer buf, BYTE buf_offset, BYTE *b,
                BYTE b_offset, size_t b_size, BYTE b_end)
{
    int diff = b_offset - buf_offset;
    BYTE *new = malloc(sizeof(BYTE));
    assert(new);
    BYTE *current;
    BYTE err;
    size_t i;

    // If there are more spaces than leading bits
    if (diff > 0)
    {
        // Shift the correct bits in place
        *new = (*b << diff) & (0xff >> buf_offset);

        // Fill remaining spaces with next byte
        if (b_size > 1)
        {
            *new += *(b + 1) >> (8 - diff);
        }

        // Write into buf's latest byte with openings
        if (buf_offset)
        {
            i = buf->size - 1;
            buf->data[i] &= (0xff << (8 - buf_offset));
            buf->data[i] += *new;
        }
        else
        {
            err = write_data(buf, new, 1);
            if (err) return err;
        }
        
        // Write all remaining bytes, shifted over
        current = b + 1;
        while (b_size-- > 2)
        {
            *new = *current << diff;
            *new += *(current + 1) >> (8 - diff);
            err = write_data(buf, new, 1);
            if (err) return err;
            current++;
        }

        // Write final truncated byte up to b_end bits
        if (b_size == 1)
        {
            b_end = b_end ? b_end : 8;
            if (b_end > diff)
            {
                *new = *current << diff;
                *new &= (0xff << (8 - (b_end - diff)));
                err = write_data(buf, new, 1);
                if (err) return err;
            }
            else if (b_end < diff)
            {
                // Cut-off extranneous bits added
                i = buf->size - 1;
                size_t extr = diff - b_end;
                buf->data[i] &= (0xff << extr);
            }

            // Otherwise, correct bits already copied
        }
    }

    // If there are more leading bits than spaces
    else if (diff < 0)
    {
        // Shift the correct bits in place
        *new = (*b >> abs(diff)) & (0xff >> buf_offset);

        // Write into buf's latest byte with openings
        i = buf->size - 1;
        buf->data[i] &= (0xff << (8 - buf_offset));
        buf->data[i] += *new;

        // Write all remaining bytes, shifted over
        current = b;
        while (b_size-- > 1)
        {
            *new = *current << (8 + diff);
            *new += *(current + 1) >> abs(diff);
            err = write_data(buf, new, 1);
            if (err) return err;
            current++;
        }

        // Write final truncated byte up to b_end bits
        if (b_size == 0)
        {
            b_end = b_end ? b_end : 8;
            if (b_end > 8 + diff)
            {
                *new = *current << (8 + diff);
                *new &= (0xff << (16 + diff - b_end));
                err = write_data(buf, new, 1);
                if (err) return err;
            }
            else if (b_end < 8 + diff)
            {
                // Cut-off extranneous bits added
                i = buf->size - 1;
                size_t extr = 8 + diff - b_end;
                buf->data[i] &= (0xff << extr);
            }

            // Otherwise, correct bits already copied
        }
    }

    // If the leading bits fit evenly in the spaces
    else
    {
        // Isolate the correct bits
        *new = *b & (0xff >> b_offset);

        // Write into buf's latest byte with openings
        if (buf_offset)
        {
            i = buf->size - 1;
            buf->data[i] &= (0xff << (8 - buf_offset));
            buf->data[i] += *new;
        }
        else
        {
            err = write_data(buf, new, 1);
            if (err) return err;
        }

        // Original: Write all remaining data
        if (b_size > 2)
        {
            err = write_data(buf, b + 1, b_size - 2);
            if (err) return err;
        }

        // Write the trailing bits
        if (b_size > 1)
        {
            b_end = b_end ? b_end : 8;
            *new = *(b + b_size - 1) & (0xff << (8 - b_end));
            err = write_data(buf, new, 1);
            if (err) return err;
        }
    }

    free(new);
    return 0;
}

// Change the last n bits in the last byte of a buffer
void change_bits(Buffer buf, const BYTE *b, BYTE n)
{
    // Get the relevant bits from byte b
    BYTE bits = *b & (0xff >> (8 - n));

    // Insert these alongside target's remaining bits
    size_t i = buf->size - 1;
    buf->data[i] = (buf->data[i] & (0xff << n)) + bits;
}

// Free up space used by a buffer
void destroy_buffer(Buffer buf)
{
    free(buf->data);
    free(buf);
}

// Read and copy APP0 segment from file
BYTE read_header(FILE *f, Buffer b)
{
    char *marker = calloc(2, sizeof(char));
    fread(marker, sizeof(char), 2, f);

    // Confirm valid APP0
    if (marker[0] != '\xff' || marker[1] != '\xe0')
    {
        fprintf(stderr, "ERROR: Invalid APP0\n");
        free(marker);
        destroy_buffer(b);
        return 1;
    }
    free(marker);

    // Record thumbnail size from segment
    uint16_t *size = malloc(sizeof(uint16_t));
    fread(size, sizeof(uint16_t), 1, f);
    uint16_t true_size = reverse_bytes(*size);
    free(size);

    // Check if you need to expand the buffer
    while (b->capacity <= b->size + true_size + 2)
    {
        b->capacity *= SIZE_MULTIPLIER;
        b->data = realloc(b->data, b->capacity);
        if (!b->data) return 2;
    }

    // Copy everything over
    b->data[b->size++] = '\xff';
    b->data[b->size++] = '\xe0';
    b->data[b->size++] = (BYTE) ((true_size >> 8) & 0xff);
    b->data[b->size++] = (BYTE) true_size & 0xff;
    fread(b->data + b->size, 1, true_size - 2, f);

    // Ensure APP0 fields were all regulation
    char *jfif = "JFIF";
    for (int i = 0; i < 4; i++)
    {
        if (b->data[b->size + i] != jfif[i])
        {
            fprintf(stderr, "ERROR: Invalid APP0\n");
            destroy_buffer(b);
            return 3;
        }
    }
    if (b->data[b->size + 5] != 1)
    {
        fprintf(stderr, "ERROR: Invalid APP0\n");
        destroy_buffer(b);
        return 3;
    }

    b->size += true_size - 2;
    return 0;
}

// Read and skip APP0 segment from file
int skip_header(FILE *f)
{
    char *marker = calloc(2, sizeof(char));
    fread(marker, sizeof(char), 2, f);

    // Confirm valid APP0
    if (marker[0] != '\xff' || marker[1] != '\xe0')
    {
        fprintf(stderr, "ERROR: Invalid APP0\n");
        free(marker);
        return 1;
    }
    free(marker);

    uint16_t *size = malloc(sizeof(uint16_t));
    fread(size, sizeof(uint16_t), 1, f);
    uint16_t true_size = reverse_bytes(*size);
    free(size);
    return 0;
}

// Create a chunk for storing image data
Chunk create_chunk(void)
{
    Chunk c = malloc(sizeof(struct chunk));
    c->data = 0x00;
    c->bd = 0;
    c->height0 = 0;
    c->width0 = 0;
    c->restart_intr = 0;
    return c;
}

// Nine type bytes for segmengs with no length/contents
BYTE tags[10] = {0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 
                     0xd6, 0xd7, 0x00, 0x01};

// Determine if a segment is parameterless
BYTE paramless(const Chunk c)
{
    if (c->type[1] == tags[8] || c->type[1] == tags[9])
    {
        return 1;
    } 
    else if (c->type[1] >= tags[0] && c->type[1] <= tags[7])
    {
        return 1;
    }

    return 0;
}

// Populate a chunk object from the raw image file
BYTE read_chunk(FILE *f, Chunk c)
{
    // Read in type
    if (fread(c->type, sizeof(char), 2, f) != 2)
    {
        return 1;
    }

    // Check if it's a parameterless segment
    if (paramless(c))
    {
        c->size = 0;
        return 0;
    }

    // Read in size
    uint16_t *size = malloc(sizeof(uint16_t));
    assert(size);
    if (!fread(size, sizeof(uint16_t), 1, f))
    {
        free(size);
        return 1;
    }
    c->size = *size;
    uint16_t true_size = reverse_bytes(*size);
    free(size);

    // Check if chunk has already been used, free data
    if (c->data != 0x00) free(c->data);

    // Read in remaining data to allocated array
    size_t z = true_size - 2;
    BYTE *data = malloc(z);
    assert(data);
    if (fread(data, sizeof(BYTE), z, f) != z)
    {
        return 1;
    }
    c->data = data;

    // If it's a SOF, store the image dimensions
    if (c->type[1] == 0xc0)
    {
        // Bit Depth
        c->bd = c->data[0];

        // Height & Width
        c->height0 = (c->data[1] * 256) + c->data[2];
        c->width0 = (c->data[3] * 256) + c->data[4];
    }

    // If it's a DRI, read in the restart interval
    else if (c->type[1] == 0xdd)
    {
        c->restart_intr = (c->data[0] * 256) + c->data[1];
    }

    // If it's a SOS, check scan-specific data
    else if (c->type[1] == 0xda)
    {
        c->num_comps = c->data[0];
        for (int i = 0; i < c->num_comps; i++)
        {
            c->comps[i] = c->data[1 + (2 * i)];
            c->table_ids[i] = c->data[2 + (2 * i)];
        }
        c->start = c->data[1 + (2 * c->num_comps)];
        c->end = c->data[2 + (2 * c->num_comps)];
        c->successive_appr = c->data[3 + (2 * c->num_comps)];
    }

    return 0;
}

// Destroy a chunk and free its memory
void destroy_chunk(Chunk c)
{
    if (c->data != 0x00)
    {
        free(c->data);
    }
    free(c);
}

// Recursively read Huffman Tables while length l > 0
static DHT huffman_internal(FILE *f, DHT parent, uint16_t len)
{
    // Check for impending Segmentation Fault
    if (len < 17)
    {
        fprintf(stderr, "ERROR: Huffman Table wrong size\n");
        return NULL;
    }

    // Set up new table in format of parent
    DHT d = malloc(sizeof(struct huffman));
    assert(d);
    d->type[0] = parent->type[0];
    d->type[1] = parent->type[1];
    d->size = parent->size;

    // Read in info, counts, symbols
    fread(&(d->info), sizeof(char), 1, f);
    fread(d->counts, sizeof(char), 16, f);
    
    // Find actual length of symbols array
    d->symbols_len = 0;
    for (int i = 0; i < 16; i++) 
    {
        d->symbols_len += d->counts[i];
    }

    d->symbols = malloc(d->symbols_len);
    fread(d->symbols, sizeof(BYTE), d->symbols_len, f);

    // Check if there is another Huffman Table next
    len -= 16 + d->symbols_len;
    if (len > 0)
    {
        DHT next = huffman_internal(f, d, len);
        d->next = next;
    }
    else d->next = NULL;

    return d;
}

// Populate a Huffman Table from the raw image file
DHT read_huffman(FILE *f)
{
    DHT d = malloc(sizeof(struct huffman));
    assert(d);

    // Read in type
    fread(d->type, sizeof(char), 2, f);
    
    // Read in size
    uint16_t *size = malloc(sizeof(uint16_t));
    fread(size, sizeof(uint16_t), 1, f);
    d->size = *size;
    uint16_t true_size = reverse_bytes(*size);
    free(size);

    // Read in info and counts
    fread(&(d->info), sizeof(char), 1, f);
    fread(d->counts, sizeof(char), 16, f);

    // Find actual length of symbols array
    d->symbols_len = 0;
    for (int i = 0; i < 16; i++) 
    {
        d->symbols_len += d->counts[i];
    }

    d->symbols = malloc(d->symbols_len);
    fread(d->symbols, sizeof(BYTE), d->symbols_len, f);

    // Check if there is another Huffman Table next
    true_size -= 19 + d->symbols_len;
    if (true_size > 17)
    {
        DHT next = huffman_internal(f, d, true_size);
        d->next = next;
    }
    else if (true_size > 0)
    {
        fprintf(stderr, "ERROR: Huffman Table wrong size\n");
        return NULL;
    }
    else d->next = NULL;

    return d;
}

// Insert a new value into a Huffman Table
// TODO: Handle the next *huffman possibility
BYTE add_huffman(BYTE value, DHT d)
{
    // Check if the table is already full
    size_t occupied = 0;
    for (int i = 0; i < 16; i++)
    {
        occupied += (i + 1) * d->counts[i];
    }
    if (occupied >= 256) return 1;

    // No room for the symbol in the table
    if (d->counts[15]) return 2;

    // Find next available code size
    BYTE code = 15;
    while (!d->counts[code - 1]) code--;

    // Upgrade symbol count and symbols length
    d->counts[code]++;
    d->symbols = realloc(d->symbols, ++(d->symbols_len));
    if (!d->symbols) return 3;
    
    // Add new symbol to the end
    d->symbols[d->symbols_len - 1] = value;
    d->size = reverse_bytes(reverse_bytes(d->size) + 1);
    return 0;
}

// Write a Huffman Table to a results buffer
BYTE write_huffman(Buffer buf, const DHT d)
{
    BYTE err;

    // Write individual components of the Huffman Table
    err = write_data(buf, d->type, 2);
    if (err) return err;

    err = write_data(buf, (BYTE *) &(d->size), 2);
    if (err) return err;

    err = write_data(buf, &(d->info), 1);
    if (err) return err;

    err = write_data(buf, d->counts, 16);
    if (err) return err;

    err = write_data(buf, d->symbols, d->symbols_len);
    if (err) return err;

    DHT c = d;
    while (c->next)
    {
        c = c->next;
        err = write_data(buf, &(c->info), 1);
        if (err) return err;

        err = write_data(buf, c->counts, 16);
        if (err) return err;

        err = write_data(buf, c->symbols, c->symbols_len);
        if (err) return err;
    }

    return 0;
}

// Rewrite a corrected Huffman Table to a results buffer
/* Table d currently begins at index start in buffer buf
   Assumes that the DHT length has only increased by 1 byte
*/
BYTE rewrite_huffman(Buffer buf, size_t start, const DHT d)
{
    // Scroll past type bytes
    start += 2;

    // Update size if space allows
    if (buf->data[start + 1] != 255)
    {
        buf->data[start + 1]++;
    }
    else if (buf->data[start] != 255)
    {
        buf->data[start]++;
        buf->data[start + 1] = 0;
    }
    else return 1;

    // Scroll past size and info
    start += 3;

    // Update counts table
    for (int i = 0; i < 16; i++)
    {
        buf->data[start + i] = d->counts[i];
    }

    // Scroll past counts and old symbols
    // (Here, symbols_len is already +1)
    start += 15 + d->symbols_len;

    // Check if you need to expand the buffer
    while (buf->capacity <= ++(buf->size))
    {
        buf->capacity *= SIZE_MULTIPLIER;
        buf->data = realloc(buf->data, buf->capacity);
        if (!buf->data) return 2;
    }

    for (int i = buf->size - 1; i > start; i--)
    {
        buf->data[i] = buf->data[i - 1];
    }

    // Finally, add the new symbol in
    buf->data[start] = d->symbols[d->symbols_len - 1];
    return 0;
}

static void destroy_huffman_internal(DHT d)
{
    // Destroy bottom table first
    if (d->next) destroy_huffman_internal(d->next);

    free(d->symbols);
    free(d);
}

// Destroy a Huffman Table and free its memory
void destroy_huffman(DHT d)
{
    destroy_huffman_internal(d);
}
