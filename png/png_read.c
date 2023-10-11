#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"
#include "png.h"


/*  Read a PNG file out to a pixel map.
    Prints to standard out:
    1. Image Width (in pixels, 32 bit integer)
    2. Image Height (in pixels, 32 bit integer)
    3. Bytes per pixel (8 bit integer)
    4. Color profile (8 bit integer)
    5. Each pixel in order

    Returns status 0 if successful, otherwise returns a positive integer
    signifying the error code.
*/
int
main(int argc, char **argv)
{
    // Check usage
    if (argc != 2)
    {
        fprintf(stderr, "Usage: ./png_read infile\n");
        return 1;
    }

    // Open infile
    FILE *infile = fopen(argv[1], "r");
    if (!infile)
    {
        fprintf(stderr, "ERROR: No infile found\n");
        fclose(infile);
        return 3;
    }

    // Check PNG Signature
    if (!PNGSignatureValid(infile))
    {
        fprintf(stderr, "ERROR: Invalid PNG signature\n");
        fclose(infile);
        return 4;
    }

    // Read in IHDR header
    IHDR *header = malloc(sizeof(IHDR));
    assert(header);
    fread(header, sizeof(IHDR), 1, infile);
    if (strncmp(header->type, "IHDR", 4))
    {
        fclose(infile);
        free(header);
        fprintf(stderr, "ERROR: Missing IHDR chunk\n");
        return 5;
    }

    // Find initial dimensions, converting little -> big endian
    uint32_t width = reverse_bytes(header->width);
    uint32_t height = reverse_bytes(header->height);

    BYTE bd = header->bitDepth;
    if (bd % 8 != 0)
    {
        fclose(infile);
        free(header);
        fprintf(stderr, "ERROR: Bit depth must be multiple of 8\n");
        return 7;
    }

    BYTE bps = (BYTE) bd / 8;
    BYTE color = header->colorType;
    if (!BDValid(bd, color))
    {
        fclose(infile);
        free(header);
        fprintf(stderr, "ERROR: Depth %u forbitten for color %u\n", bd, color);
        return 8;
    }

    BYTE spp = samples_by_color(color);
    BYTE bpp = bps * spp;
    printf("%u%u", width, height);
    printf("%u%u", bpp, color);
    if (header->compressionMethod)
    {
        fclose(infile);
        free(header);
        fprintf(stderr, "ERROR: Unsupported compression method\n");
        return 9;
    }
    if (header->interlaceMethod)
    {
        fclose(infile);
        free(header);
        fprintf(stderr, "ERROR: Unsupported interlace method\n");
        return 10;
    }

    uint32_t crc = 0;
    crc = crc32(crc, (Bytef *) &(header->type), 17);
    crc = reverse_bytes(crc);
    free(header);

    // Set up hash table for raw values
    BYTE *h = calloc(bpp, width * height);
    assert(h);

    // Read through remaining infile contents
    Chunk c = create_chunk();
    read_chunk(infile, c);

    // Before the data
    uint32_t true_l;
    while (strncmp(c->type, "IDAT", 4))
    {
        read_chunk(infile, c);
    }

    // Consolidate data across IDAT chunks for editing
    Buffer raw_buf = create_buffer();
    while (strncmp(c->type, "IEND", 4))
    {
        true_l = reverse_bytes(c->length);
        raw_buf = write_data(raw_buf, (BYTE *) c->data, true_l);
        read_chunk(infile, c);
    }
    fclose(infile);

    // Space for the decompressed data stream
    BYTE *old_data = calloc(MAX_COMPRESSION_RATIO, raw_buf->size);
    assert(old_data);

    // Decompress the full raw data together
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    infstream.avail_in = raw_buf->size;
    infstream.next_in = raw_buf->data;
    infstream.avail_out = raw_buf->size * MAX_COMPRESSION_RATIO;
    infstream.next_out = old_data;

    if (inflateInit(&infstream) != Z_OK)
    {
        fprintf(stderr, "ERROR: Inflation failed\n");
        destroy_chunk(c);
        destroy_buffer(raw_buf);
        free(h);
        free(old_data);
        return 11;
    }
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
    destroy_buffer(raw_buf);

    // Set up important variables for iterating
    size_t vals[4] = {0, 0, (size_t) width, (size_t) bpp};
    size_t front, ind;
    BYTE filter;
    BYTE *fpx = calloc(bpp, sizeof(BYTE));
    assert(fpx);
    int status[2];

    // Iterate through all pixels
    for (int line = 0; line < height; line++)
    {
        front = ((width * bpp) + 1) * line;
        filter = *(old_data + front);
        
        // Filter and copy pixels
        for (int pixel = 0; pixel < width; pixel++)
        {
            // Index of the start of this pixel
            ind = front + (pixel * bpp) + 1;
            status[0] = pixel;
            status[1] = line;
            png_filt(fpx, filter, old_data, ind, status, vals, h);

            for (int b = 0; b < bpp; b++)
            {
                putchar(fpx[b]);
            }
        }
    }

    // Clean up
    free(fpx);
    free(h);
    free(old_data);
    destroy_chunk(c);
    return 0;
}
