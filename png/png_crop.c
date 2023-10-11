#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"
#include "png.h"


/*  Crops a PNG file to the specified dimensions.
    x = first accepted pixel in from the left 
    y = first accepted pixel down from the top
    Then, pixels are accepted up to 'width' right of x and down to 
    'height' below y, non-inclusive. Both rows and columns are 0-indexed.
    infile = path to original uncropped image
    outfile = path to new cropped image

    Returns status 0 if successful, otherwise returns a positive integer
    signifying the error code.
*/
int
main(int argc, char **argv)
{
    // Check usage
    if (argc != 7)
    {
        fprintf(stderr, "Usage: ./png_crop x y width height infile outfile\n");
        return 1;
    }
    size_t x = atoi(argv[1]);
    size_t y = atoi(argv[2]);
    size_t width = atoi(argv[3]);
    size_t height = atoi(argv[4]);

    if (x > width || y > height || !height || !width)
    {
        fprintf(stderr, "ERROR: Invalid crop dimensions\n");
        return 2;
    }

    // Open infile
    FILE *infile = fopen(argv[5], "r");
    if (!infile)
    {
        fprintf(stderr, "ERROR: No infile found\n");
        fclose(infile);
        return 3;
    }
    printf("\nProcessing %s...\n", argv[5]);

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
    uint32_t width0 = reverse_bytes(header->width);
    uint32_t height0 = reverse_bytes(header->height);
    if ((width0 < x + width) | (height0 < y + height))
    {
        fclose(infile);
        free(header);
        fprintf(stderr, "ERROR: Invalid crop dimensions\n");
        return 6;
    }

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
    printf("Dimensions: %ux%u\n", width0, height0);
    printf("Bit Depth: %u\n", bd);
    printf("Color Type: %u\n", color);

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

    // Open new receptacle for results, write the signature
    Buffer results = create_buffer();
    write_signature(results);

    // Write the header
    header->width = reverse_bytes(width);
    header->height = reverse_bytes(height);
    results = write_data(results, (BYTE *) &(header->length), 21);

    uint32_t crc = 0;
    crc = crc32(crc, (Bytef *) &(header->type), 17);
    crc = reverse_bytes(crc);
    results = write_data(results, (BYTE *) &crc, 4);
    free(header);

    // Set up hash table for raw values
    BYTE *h = calloc(bpp, width0 * height0);
    assert(h);

    // Read through remaining infile contents
    Chunk c = create_chunk();
    read_chunk(infile, c);

    // Before the data
    uint32_t true_l;
    while (strncmp(c->type, "IDAT", 4))
    {
        if (safe_to_copy(c))
        {
            // Write individual components of chunk
            results = write_data(results, (BYTE *) &(c->length), 4);
            results = write_data(results, (BYTE *) c->type, 4);
            true_l = reverse_bytes(c->length);
            results = write_data(results, (BYTE *) c->data, true_l);
            results = write_data(results, (BYTE *) &(c->crc), 4);
        }
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
        destroy_buffer(results);
        destroy_buffer(raw_buf);
        free(h);
        free(old_data);
        return 11;
    }
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
    destroy_buffer(raw_buf);

    // Set up important variables for iterating
    Buffer cropped = create_buffer();
    size_t vals[4] = {x, y, (size_t) width0, (size_t) bpp};
    size_t front, ind;
    BYTE filter;
    BYTE *fpx = calloc(bpp, sizeof(BYTE));
    assert(fpx);
    int status[2];

    // Iterate through all relevant pixels
    for (int line = y; line < y + height; line++)
    {
        front = ((width0 * bpp) + 1) * line;
        // Process filter byte at the beginning of scanlines
        filter = *(old_data + front);
        cropped = write_data(cropped, &filter, 1);
        
        // Filter and copy pixels that are staying
        for (int pixel = x; pixel < x + width; pixel++)
        {
            // Index of the start of this pixel
            ind = front + (pixel * bpp) + 1;
            status[0] = pixel;
            status[1] = line;
            png_filt(fpx, filter, old_data, ind, status, vals, h);
            cropped = write_data(cropped, fpx, (uint32_t) bpp);
        }
    }
    free(fpx);
    free(h);
    free(old_data);

    // Recompress
    BYTE *new_data = malloc(cropped->size);
    assert(new_data);

    z_stream defstream;
    defstream.zalloc = Z_NULL;
    defstream.zfree = Z_NULL;
    defstream.opaque = Z_NULL;
    defstream.avail_in = cropped->size;
    defstream.next_in = cropped->data;
    defstream.avail_out = cropped->size;
    defstream.next_out = new_data;

    // Confirm successful deflation
    if (deflateInit(&defstream, Z_DEFAULT_COMPRESSION) != Z_OK)
    {
        fprintf(stderr, "ERROR: Deflation failed\n");
        destroy_chunk(c);
        destroy_buffer(cropped);
        destroy_buffer(results);
        free(new_data);
        return 12;
    }
    deflate(&defstream, Z_FINISH);
    deflateEnd(&defstream);

    size_t new_length = (size_t) defstream.total_out;
    destroy_buffer(cropped);

    // Split into new chunks, calculate CRCs, and add to results
    char *idat = "IDAT";
    Buffer crc_buf = create_buffer();
    crc_buf = write_data(crc_buf, (BYTE *) idat, 4);
    uint32_t left, size, reversed_size, reversed_crc;

    size_t i = 0;
    while (i < new_length)
    {
        left = (uint32_t) (new_length - i);
        size = MAX_IDAT_SIZE < left ? (uint32_t) MAX_IDAT_SIZE : left;

        crc = 0;
        crc_buf = write_data_at(crc_buf, 4, new_data + i, size);
        crc = crc32(crc, (Bytef *) crc_buf->data, size + 4);

        reversed_size = reverse_bytes(size);
        reversed_crc = reverse_bytes(crc);
        results = write_data(results, (BYTE *) &reversed_size, 4);
        results = write_data(results, (BYTE *) idat, 4);
        results = write_data(results, new_data + i, size);
        results = write_data(results, (BYTE *) &reversed_crc, 4);
        i += size;
    }
    destroy_buffer(crc_buf);
    free(new_data);

    // Append the final IEND chunk
    results = write_data(results, (BYTE *) &(c->length), 4);
    results = write_data(results, (BYTE *) c->type, 4);
    true_l = reverse_bytes(c->length);
    results = write_data(results, (BYTE *) c->data, true_l);
    results = write_data(results, (BYTE *) &(c->crc), 4);

    // Write to the output file
    FILE *outfile = fopen(argv[6], "w");
    if (!outfile)
    {
        fprintf(stderr, "ERROR: No outfile found\n");
        destroy_chunk(c);
        destroy_buffer(results);
        fclose(outfile);
        return 13;
    }

    fwrite(results->data, sizeof(BYTE), results->size, outfile);
    printf("\nCrop Complete!\nNew dimensions: %zux%zu\n", width, height);

    // Clean up
    destroy_chunk(c);
    destroy_buffer(results);
    fclose(outfile);
    return 0;
}
