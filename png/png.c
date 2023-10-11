#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "png.h"

// Function Definitions for bmp.h

/*  PNG File Structure:
    First 8 bytes of a PNG file are always the following values:
    137 80 78 71 13 10 26 10
*/
char PNGSignature[8] = {'\x89', 'P', 'N', 'G', '\r', '\n', '\x1a', '\n'};

// Return 1 if a PNG signature is valid, otherwise 0
int PNGSignatureValid(FILE *f)
{
    char c;
    for (int i = 0; i < PNGSignatureSize; i++)
    {
        fread(&c, sizeof(char), 1, f);
        if (c != PNGSignature[i])
        {
            return 0;
        }
    }

    return 1;
}

// The CPU will read the uint32_t in the wrong order
uint32_t reverse_bytes(uint32_t bytes)
{
    uint32_t aux = 0;
    BYTE b;
    
    for (int i = 0; i < 32; i+=8)
    {
        b = (bytes >> i) & 0xff;
        aux |= b << (24 - i);
    }
    return aux;
}

// Return 1 if bit depth is valid for given color, otherwise 0
int BDValid(const BYTE bd, const BYTE color)
{
    if (color == 0)
    {
        int valid[5] = {1, 2, 4, 8, 16};
        for (int i = 0; i < 5; i++)
        {
            if (bd == valid[i])
            {
                return 1;
            }
        }
    }

    else if (color == 2 || color == 4 || color == 6)
    {
        if (bd == 8 || bd == 16)
        {
            return 1;
        }
    }

    else if (color == 3)
    {
        int valid[4] = {1, 2, 4, 6};
        for (int i = 0; i < 4; i++)
        {
            if (bd == valid[i])
            {
                return 1;
            }
        }
    }

    // If none of the above are true, invalid bit depth
    return 0;
}

// Find samples per pixel by input color
BYTE samples_by_color(const BYTE color)
{
    if (color == 0) return 1;
    else if (color == 2 || color == 3) return 3;
    else if (color == 4) return 2;
    else if (color == 6) return 4;

    // If none of the above are true, invalid color
    else return 0;
}

// Create a buffer
Buffer create_buffer(void)
{
    Buffer buf = malloc(sizeof(struct buffer));
    BYTE *data = malloc(BUFFER_STARTING_SIZE);
    buf->data = data;
    buf->size = 0;
    buf->capacity = BUFFER_STARTING_SIZE;

    return buf;
}

// Free up space used by a buffer
void destroy_buffer(Buffer buf)
{
    free(buf->data);
    free(buf);
}

// Add data to a buffer
#define SIZE_MULTIPLIER (2)
Buffer write_data(Buffer buf, BYTE *b, uint32_t n)
{
    // Check if you need to expand the buffer
    while (buf->capacity <= buf->size + n)
    {
        buf->data = realloc(buf->data, buf->capacity * SIZE_MULTIPLIER);
        buf->capacity *= SIZE_MULTIPLIER;
    }

    for (int i = 0; i < n; i++)
    {
        buf->data[buf->size++] = b[i];
    }
    return buf;
}

// Add data to a specific point in a buffer
Buffer write_data_at(Buffer buf, size_t start, BYTE *b, uint32_t n)
{
    buf->size = start;

    // Check if you need to expand the buffer
    while (buf->capacity <= buf->size + n)
    {
        buf->data = realloc(buf->data, buf->capacity * SIZE_MULTIPLIER);
        buf->capacity *= SIZE_MULTIPLIER;
    }

    for (int i = 0; i < n; i++)
    {
        buf->data[buf->size++] = b[i];
    }
    return buf;
}

// Write first 8 bytes of a PNG to a new buffer
void write_signature(Buffer b)
{
    for (int i = 0; i < PNGSignatureSize; i++)
    {
        b->data[i] = PNGSignature[i];
    }

    b->size = PNGSignatureSize;
}

// Create a chunk for storing image data
Chunk create_chunk(void)
{
    Chunk c = malloc(sizeof(struct chunk));
    c->data = 0x00;
    return c;
}

// Populate a chunk object from the raw image file
void read_chunk(FILE *f, Chunk c)
{
    // Read in length
    uint32_t *length = malloc(sizeof(uint32_t));
    fread(length, sizeof(uint32_t), 1, f);
    c->length = *length;
    uint32_t true_length = reverse_bytes(*length);
    free(length);

    // Read in type
    fread(c->type, sizeof(char), 4, f);

    // Read in data to allocated array
    if (c->data != 0x00)
    {
        // Chunk has already been used, free data pointer
        free(c->data);
    }
    BYTE *data = malloc(true_length);
    fread(data, sizeof(BYTE), true_length, f);
    c->data = data;

    // Read in CRC
    uint32_t *crc = malloc(sizeof(uint32_t));
    fread(crc, sizeof(uint32_t), 1, f);
    c->crc = *crc;
    free(crc);
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

// Determine if a chunk is safe to copy
int safe_to_copy(const Chunk c)
{
    char d = c->type[3];
    if (d >= 'a' && d <= 'z')
    {
        return 1;
    }
    else return 0;
}

/* Compute the Paeth Predictor for three bytes.
   Inputs a = left byte, b = upper byte, c = upper left byte.
*/
static BYTE PaethPredictor(BYTE a, BYTE b, BYTE c)
{
    int p = a + b - c;
    size_t pa = abs(p - a);
    size_t pb = abs(p - b);
    size_t pc = abs(p - c);
    if (pa <= pb && pa <= pc)
    {
        return a;
    }
    else if (pb <= pc)
    {
        return b;
    }
    else
    {
        return c;
    }
}

/*  Find the raw integer value of a filtered byte recursively.
    Inputs are the same as for png_filt(), except a solo byte is included 
    instead of a full pixel, and the index within a pixel of this byte
    labeled b (int) is included as a final argument.

    status = {pixel, line}

    While this function is very similar in structure to png_filt,
    it is necessary to split these two functions up for recursion.
    Additional PNG Filtering documentation:
    http://www.libpng.org/pub/png/spec/1.2/PNG-Filters.html

    Returns the unsigned integer value of the raw byte.
*/
static BYTE raw_val(const BYTE type, BYTE byte, const BYTE *old_data, 
                    const size_t ind, const int status[2], 
                    const size_t vals[4], BYTE *h, size_t b)
{
    size_t width0 = vals[2];
    BYTE bpp = (BYTE) vals[3];

    // First, check the hashmap to see if the value is already there
    size_t index = ((status[1] * width0) + status[0]) * bpp + b;
    if (h[index] != 0x00)
    {
        return h[index];
    }

    // Method 0: None
    else if (type == 0)
    {
        h[index] = byte;
        return byte;
    }

    // Method 1: Sub
    /* This method transmits the difference between each byte and the
       corresponding byte of the prior pixel in the same scanline 
    */
    else if (type == 1)
    {
        int raw = 0;
        size_t delta;
        for (int i = 0; i <= status[0]; i++)
        {
            delta = i * bpp;
            raw += old_data[ind + b - delta];
        }
        raw = ((raw % 256) + 256) % 256;

        h[index] = (BYTE) raw;
        return (BYTE) raw;
    }

    // Method 2: Up
    /* This method transmits the difference between each byte and the
       corresponding byte of the pixel immediately above it
    */
    else if (type == 2)
    {
        if (status[1] == 0)
        {
            // Pixel is already raw value
            return byte;
        }

        size_t up_ind = ind - ((width0 * bpp) + 1);
        BYTE up = old_data[up_ind + b];
        size_t type_delta = (status[0] * bpp) + 1;
        BYTE up_type = old_data[up_ind - type_delta];
        int up_status[2] = {status[0], status[1] - 1};

        up = raw_val(up_type, up, old_data, up_ind, up_status, vals, h, b);
        int raw = byte + up;
        raw = ((raw % 256) + 256) % 256;

        h[index] = (BYTE) raw;
        return (BYTE) raw;
    }

    // Method 3: Average
    /* This method uses the average of the two neighboring pixels (left 
       and above) to predict the value of a pixel
    */
    else if (type == 3)
    {
        BYTE left, up;
        if (status[1] == 0)
        {
            up = 0;
        }
        else
        {
            size_t up_ind = ind - ((width0 * bpp) + 1);
            up = old_data[up_ind + b];
            size_t type_delta = (status[0] * bpp) + 1;
            BYTE up_type = old_data[up_ind - type_delta];
            int up_status[2] = {status[0], status[1] - 1};
            up = raw_val(up_type, up, old_data, up_ind, up_status, vals, 
                         h, b);
        }
        if (status[0] == 0)
        {
            left = 0;
        }
        else
        {
            size_t left_ind = ind - bpp;
            left = old_data[left_ind + b];
            int left_status[2] = {status[0] - 1, status[1]};
            left = raw_val(type, left, old_data, left_ind, left_status, 
                           vals, h, b);
        }
        
        int raw = byte + (int) floor((left + up) / 2);
        raw = ((raw % 256) + 256) % 256;

        h[index] = (BYTE) raw;
        return (BYTE) raw;
    }

    // Method 4: Paeth
    /* This method computes a simple linear function of the three 
       neighboring pixels (left, up, upper left), then chooses as 
       predictor the neighboring pixel closest to the computed value
    */
    else
    {
        BYTE left, up, upl;
        if (status[1] == 0)
        {
            up = 0;
            upl = 0;
        }
        else
        {
            size_t up_ind = ind - ((width0 * bpp) + 1);
            up = old_data[up_ind + b];

            size_t type_delta = (status[0] * bpp) + 1;
            BYTE up_type = old_data[up_ind - type_delta];
            int up_status[2] = {status[0], status[1] - 1};
            up = raw_val(up_type, up, old_data, up_ind, up_status, vals, 
                         h, b);

            if (status[0] == 0)
            {
                upl = 0;
            }
            else
            {
                size_t upl_ind = ind - ((width0 + 1) * bpp + 1);
                upl = old_data[upl_ind + b];
                int upl_status[2] = {status[0] - 1, status[1] - 1};
                upl = raw_val(up_type, upl, old_data, upl_ind, upl_status, 
                              vals, h, b);
            }
        }

        if (status[0] == 0)
        {
            left = 0;
        }
        else
        {
            size_t left_ind = ind - bpp;
            left = old_data[left_ind + b];
            int left_status[2] = {status[0] - 1, status[1]};
            left = raw_val(type, left, old_data, left_ind, left_status, 
                           vals, h, b);
        }

        int raw = byte + PaethPredictor(left, up, upl);
        raw = ((raw % 256) + 256) % 256;

        h[index] = (BYTE) raw;
        return (BYTE) raw;
    }
}

/* Applies one of five PNG filters to a pixel, based on the 
   preceding pixels in the original image. 
*/
void png_filt(BYTE *fpx, const BYTE type, BYTE *old_data, const size_t ind, 
              const int status[2], const size_t vals[4], BYTE *h)
{
    size_t x = vals[0];
    size_t y = vals[1];
    size_t width0 = vals[2];
    BYTE bpp = (BYTE) vals[3];

    const BYTE *pix = (const BYTE *) old_data + ind;
    // If we aren't chopping off the top or left, no filter necessary
    if (x == 0 && y == 0)
    {
        for (int i = 0; i < bpp; i++)
        {
            fpx[i] = pix[i];
        }
        return;
    }

    switch(type) {
        // Method 0: None
        case 0:
            for (int i = 0; i < bpp; i++)
            {
                fpx[i] = pix[i];
            }
            break;

        /* Method 1: Sub
           This method transmits the difference between each byte and the value
           of the corresponding byte of the prior pixel in the same scanline
        */
        case 1:
            if (status[0] > x)
            {
                // We've already solved this row by fixing the first value
                for (int i = 0; i < bpp; i++)
                {
                    fpx[i] = pix[i];
                }
            }

            else if (status[0] == x)
            {
                int raw;
                size_t delta;
                for (size_t b = 0; b < bpp; b++)
                {
                    // Determine its raw value
                    raw = 0;
                    for (int i = 0; i <= status[0]; i++)
                    {
                        delta = i * bpp;
                        raw += old_data[ind + b - delta];
                    }
                    raw = ((raw % 256) + 256) % 256;

                    // Because this is the new first pixel, sub(x) = raw(x) - 0
                    fpx[b] = (BYTE) raw;
                }
            }
            
            else
            {
                printf("ERROR: Filtering cropped pixel\n");
            }
            break;

        /* Method 2: Up
           This method transmits the difference between each byte and the value
           of the corresponding byte of the pixel immediately above
        */
        case 2:
            if (y == 0 || status[1] > y)
            {
                // This line is unaffected by the crop
                for (int i = 0; i < bpp; i++)
                {
                    fpx[i] = pix[i];
                }
            }

            else if (status[1] == y)
            {
                BYTE up, up_type;
                size_t up_ind, type_delta;
                int raw;
                for (size_t b = 0; b < bpp; b++)
                {
                    // Determine its raw value
                    up_ind = ind - ((width0 * bpp) + 1);
                    up = old_data[up_ind + b];
                    type_delta = (status[0] * bpp) + 1;
                    up_type = old_data[up_ind - type_delta];
                    int up_status[2] = {status[0], status[1] - 1};

                    up = raw_val(up_type, up, old_data, up_ind, up_status, 
                                 vals, h, b);
                    raw = pix[b] + up;
                    raw = ((raw % 256) + 256) % 256;

                    // Because this is the new first line, up(x) = raw(x) - 0
                    fpx[b] = (BYTE) raw;
                }
            }
            
            else
            {
                printf("ERROR: Filtering cropped pixel\n");
            }
            break;
        
        /* Method 3: Average
           This method uses the average of the two neighboring pixels (left and 
           above) to predict the value of a pixel
        */
        case 3:
            if ((status[0] > x || x == 0) && (status[1] > y || y == 0))
            {
                // This line is unaffected by the crop
                for (int i = 0; i < bpp; i++)
                {
                    fpx[i] = pix[i];
                }
            }

            else if (status[0] == x || status[1] == y)
            {
                // Either top or left side are cropped
                size_t type_delta, left_ind, up_ind;
                BYTE left, up, up_type;
                int avg, raw;
                for (size_t b = 0; b < bpp; b++)
                {
                    // Determine its raw value
                    raw = raw_val(type, pix[b], old_data, ind, status, vals, h, b);

                    // Determine the left pixel's value
                    if (status[0] == x)
                    {
                        left = 0;
                    }
                    else
                    {
                        left_ind = ind - bpp;
                        left = old_data[left_ind + b];
                        int left_status[2] = {status[0] - 1, status[1]};
                        left = raw_val(type, left, old_data, left_ind, left_status, 
                                       vals, h, b);
                    }

                    // Determine the top pixel's value
                    if (status[1] == y)
                    {
                        up = 0;
                    }
                    else
                    {
                        up_ind = ind - ((width0 * bpp) + 1);
                        up = old_data[up_ind + b];
                        type_delta = (status[0] * bpp) + 1;
                        up_type = old_data[up_ind - type_delta];
                        int up_status[2] = {status[0], status[1] - 1};
                        up = raw_val(up_type, up, old_data, up_ind, up_status, 
                                     vals, h, b);
                    }

                    avg = raw - (int) floor((left + up) / 2);
                    avg = ((avg % 256) + 256) % 256;
                    fpx[b] = (BYTE) avg;
                }
            }

            else
            {
                printf("ERROR: Filtering cropped pixel\n");
            }
            break;

        /* Method 4: Paeth
           This method computes a simple linear function of the three neighboring 
           pixels (left, up, upper left), then chooses as a predictor the 
           neighboring pixel closest to the computed value
        */
        case 4:
            if ((status[0] > x || x == 0) && (status[1] > y || y == 0))
            {
                // This line is unaffected by the crop
                for (int i = 0; i < bpp; i++)
                {
                    fpx[i] = pix[i];
                }
            }

            else if (status[0] == x || status[1] == y)
            {
                // Either top or left side are cropped
                size_t up_ind, left_ind, upl_ind, type_delta;
                BYTE up, left, upl, up_type;
                int raw, paeth;
                for (size_t b = 0; b < bpp; b++)
                {
                    // Determine its raw value
                    raw = raw_val(type, pix[b], old_data, ind, status, vals, h, b);

                    // Determine the left pixel's value
                    if (status[0] == x)
                    {
                        left = 0;
                    }

                    else
                    {
                        left_ind = ind - bpp;
                        left = old_data[left_ind + b];
                        int left_status[2] = {status[0] - 1, status[1]};
                        left = raw_val(type, left, old_data, left_ind, left_status, 
                                       vals, h, b);
                    }

                    // Determine the upper and upper left pixel's value
                    if (status[1] == y)
                    {
                        up = 0;
                        upl = 0;
                    }
                    else
                    {
                        up_ind = ind - (width0 * bpp + 1);
                        up = old_data[up_ind + b];
                        type_delta = (status[0] * bpp) + 1;
                        up_type = old_data[up_ind - type_delta];
                        int up_status[2] = {status[0], status[1] - 1};
                        up = raw_val(up_type, up, old_data, up_ind, up_status, 
                                     vals, h, b);

                        if (status[0] == x)
                        {
                            upl = 0;
                        }
                        else
                        {
                            upl_ind = ind - (bpp * (width0 + 1) + 1);
                            upl = old_data[upl_ind + b];
                            int upl_status[2] = {status[0] - 1, status[1] - 1};
                            upl = raw_val(up_type, upl, old_data, upl_ind, 
                                          upl_status, vals, h, b);
                        }
                    }

                    paeth = raw - PaethPredictor(left, up, upl);
                    paeth = ((paeth % 256) + 256) % 256;
                    fpx[b] = (BYTE) paeth;
                }
            }
            
            else
            {
                printf("ERROR: Filtering cropped pixel\n");
            }
            break;
        
        // Other: incorrect input
        default:
            printf("ERROR: Filtering cropped pixel\n");
            break;
    }
}
