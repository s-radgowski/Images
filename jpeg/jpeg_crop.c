#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jpeg.h"

/*  Crops a JPEG file to the specified dimensions.
    x = first accepted pixel in from the left
    y = first accepted pixel down from the top
    Then, pixels are accepted up to 'width' right of x and down to 
    'height' below y, non-inclusive. Both rows & columns are 0-indexed.
    infile = path to original uncropped image
    outfile = path to new cropped image

    Returns 0 if successful, otherwise returns a positive integer
    signifying the error status code.
*/
int main(int argc, char **argv)
{
    // Check usage
    if (argc != 7)
    {
        fprintf(stderr, "Usage: ./jpeg_crop x y width height infile outfile\n");
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

    // Check JPEG Signature
    if (!JPEGSignatureValid(infile))
    {
        fprintf(stderr, "ERROR: Invalid JPEG signature\n");
        fclose(infile);
        return 4;
    }

    // Open new receptacle for results, write the signature
    Buffer results = create_buffer();
    write_signature(results);

    // Copy APP0 Segment
    BYTE err = read_header(infile, results);
    if (err)
    {
        // Invalid APP0 Segment
        fclose(infile);
        return 5;
    }

    // Quantization tables DQT before SOF
    Chunk c = create_chunk();
    read_chunk(infile, c);
    uint16_t true_l;
    while ((c->type[1] < 0xc0) || (c->type[1] > 0xcf))
    {
        // Write individual components of chunk
        write_data(results, c->type, 2);

        // Check if it's a parameterless segment
        if (!paramless(c))
        {
            write_data(results, (BYTE *) &(c->size), 2);
            true_l = reverse_bytes(c->size) - 2;
            write_data(results, c->data, true_l);
        }
        if (read_chunk(infile, c))
        {
            fprintf(stderr, "ERROR: No SOF found\n");
            destroy_buffer(results);
            destroy_chunk(c);
            fclose(infile);
            return 6;
        }
    }

    // c0 = Baseline, c1 = Extended Sequential
    // c2 = Progressive, c3 = Lossless
    BYTE mode = c->type[1];

    // Currently, only Baseline DCT is supported
    if (mode != 0xc0)
    {
        fprintf(stderr, "ERROR: Only Baseline DCT is supported\n");
        destroy_buffer(results);
        destroy_chunk(c);
        fclose(infile);
        return 7;
    }

    // Interpret data from SOF0 segment
    uint16_t height0 = c->height0;
    uint16_t width0 = c->width0;
    uint16_t restart_intr = c->restart_intr;
    if (!height0 || !width0)
    {
        fprintf(stderr, "ERROR: Invalid starting dimensions\n");
        destroy_buffer(results);
        destroy_chunk(c);
        fclose(infile);
        return 8;
    }
    uint8_t bd = c->bd;
    uint8_t color = (uint8_t) c->data[5];
    if (bd != 8 && bd != 12)
    {
        fprintf(stderr, "ERROR: Unsupported bit depth\n");
        destroy_buffer(results);
        destroy_chunk(c);
        fclose(infile);
        return 9;
    }
    if (color != 1 && color != 3)
    {
        fprintf(stderr, "ERROR: Unsupported color type\n");
        destroy_buffer(results);
        destroy_chunk(c);
        fclose(infile);
        return 10;
    }

    // Screen for YIQ Color mode
    if (c->data[1] == 4 || c->data[1] == 5)
    {
        fprintf(stderr, "ERROR: YIQ Color Mode not supported\n");
        destroy_buffer(results);
        destroy_chunk(c);
        fclose(infile);
        return 11;
    }
    if (color == 3)
    {
        if (c->data[4] == 4 || c->data[4] == 5)
        {
            fprintf(stderr, "ERROR: YIQ Color Mode not supported\n");
            destroy_buffer(results);
            destroy_chunk(c);
            fclose(infile);
            return 12;
        }
        if (c->data[7] == 4 || c->data[7] == 5)
        {
            fprintf(stderr, "ERROR: YIQ Color Mode not supported\n");
            destroy_buffer(results);
            destroy_chunk(c);
            fclose(infile);
            return 12;
        }

        // Confirm unique color component IDs
        if (c->data[1] == c->data[4] || c->data[1] == c->data[7] || 
            c->data[4] == c->data[7])
        {
            fprintf(stderr, "ERROR: Duplicate color component ID\n");
            destroy_buffer(results);
            destroy_chunk(c);
            fclose(infile);
            return 12;
        }
    }

    printf("Dimensions: %ux%u\n", width0, height0);
    printf("Bit Depth: %u\n", bd);
    printf("Color Type: %u\n", color);

    // Determine Chroma Subsampling
    BYTE y_dims;
    BYTE cb_dims = 0;
    BYTE cr_dims = 0;

    y_dims = (BYTE) (c->data[7] >> 4) * (c->data[7] & 0x0f);
    if (color == 3)
    {
        cb_dims = (BYTE) (c->data[10] >> 4) * (c->data[10] & 0x0f);
        cr_dims = (BYTE) (c->data[13] >> 4) * (c->data[13] & 0x0f);
    }

    BYTE x_scales[color];
    BYTE y_scales[color];
    x_scales[0] = c->data[7] >> 4;
    y_scales[0] = c->data[7] & 0x0f;

    BYTE x_scale = x_scales[0];
    BYTE y_scale = y_scales[0];
    if (color == 3)
    {
        x_scales[1] = c->data[10] >> 4;
        x_scales[2] = c->data[13] >> 4;
        y_scales[1] = c->data[10] & 0x0f;
        y_scales[2] = c->data[13] & 0x0f;
        x_scale = max_byte3(x_scale, x_scales[1], x_scales[2]);
        y_scale = max_byte3(y_scale, y_scales[1], y_scales[2]);
    }

    // Adjust crop dimensions for convenience
    BYTE diff = x % (8 * x_scale);
    if (diff)
    {
        // Shift to left-side MCU edge
        if (diff < (5 * x_scale)) x -= diff;
        
        // If not enough room on right, also shift left
        else if (x + width + (8 * x_scale) - diff > width0)
        {
            x -= diff;
        }

        // Shift to right-side MCU edge
        else x += (8 * x_scale) - diff;
    }

    diff = y % (8 * y_scale);
    if (diff)
    {
        // Shift to upper MCU edge
        if (diff < (5 * y_scale)) y -= diff;
        
        // If not enough room on bottom, also shift up
        else if (y + height + (8 * y_scale) - diff > height0)
        {
            y -= diff;
        }

        // Shift to bottom MCU edge
        else y += (8 * y_scale) - diff;
    }

    // Write data from SOF0 to results buffer
    write_data(results, c->type, 2);
    write_data(results, (BYTE *) &(c->size), 2);
    write_data(results, (BYTE *) &(c->bd), 1);

    // Write new crop dimensions to results buffer
    uint16_t height_r = reverse_bytes(height);
    uint16_t width_r = reverse_bytes(width);
    write_data(results, (BYTE *) &height_r, 2);
    write_data(results, (BYTE *) &width_r, 2);
    true_l = reverse_bytes(c->size) - 7;
    write_data(results, c->data + 5, true_l);

    // Huffman Tables DHT before SOS
    DHT dc_huffman_1 = read_huffman(infile);
    Node dc_tree_1 = create_tree(dc_huffman_1);
    DHT ac_huffman_1 = read_huffman(infile);
    Node ac_tree_1 = create_tree(ac_huffman_1);

    DHT dc_huffman_2;
    Node dc_tree_2;
    DHT ac_huffman_2;
    Node ac_tree_2;
    if (color == 3)
    {
        if (dc_huffman_1->next)
        {
            dc_huffman_2 = dc_huffman_1->next;
        }
        else dc_huffman_2 = read_huffman(infile);
        dc_tree_2 = create_tree(dc_huffman_2);

        if (ac_huffman_1->next)
        {
            ac_huffman_2 = ac_huffman_1->next;
        }
        else ac_huffman_2 = read_huffman(infile);
        ac_tree_2 = create_tree(ac_huffman_2);
    }

    // Write Huffman Tables to results
    size_t dc_1_start = results->size;
    size_t dc_2_start;
    write_huffman(results, dc_huffman_1);
    write_huffman(results, ac_huffman_1);
    if (color == 3)
    {
        dc_2_start = results->size;
        write_huffman(results, dc_huffman_2);
        write_huffman(results, ac_huffman_2);
    }

    // Read SOS chunk
    read_chunk(infile, c);
    BYTE num_comps = c->num_comps;
    BYTE *comps = malloc(c->num_comps);
    BYTE *table_ids = malloc(c->num_comps);
    for (int i = 0; i < num_comps; i++)
    {
        comps[i] = c->comps[i];
        table_ids[i] = c->table_ids[i];
    }
    BYTE sos_start = c->start;
    BYTE sos_end = c->end;
    BYTE sa_low = c->successive_appr & 0x0f;

    // Create recepticles for reading in raw image data
    Buffer raw_buf = create_buffer();
    BYTE b1 = 0;
    BYTE b2 = 0;

    // Set up recepticle for saved MCUs
    Buffer cropped = create_buffer();
    size_t cropped_start, cropped_end;
    size_t copied_MCUs = 0;
    BYTE cropped_offset = 0;
    BYTE filler = 0x00;
    BYTE flex[4] = {0x00, 0x00, 0x00, 0x00};

    // Set up important variables for iterating
    BYTE *place = raw_buf->data;
    BYTE *last_cropped = cropped->data;
    BYTE decoded_dc[2] = {0, 0};
    BYTE marker[2] = {0xff, 0xd0};
    BYTE *MCU_start = place;
    BYTE MCU_bit_offset = 0;
    size_t seen_acs, cs_ind, cs_dims, m;
    size_t bit_offset = 0;
    size_t edge_offset = 0;
    size_t mcu_size = 0;

    // Set up new dimensions in terms of MCUs
    size_t MCU_rows = (height0 + 7) / (8 * y_scale);
    size_t MCU_cols = (width0 + 7) / (8 * x_scale);
    size_t min_x_MCU = x / (8 * x_scale);
    size_t min_y_MCU = y / (8 * y_scale);
    size_t max_x_MCU = ((width + 7) / (8 * x_scale)) + min_x_MCU - 1;
    size_t max_y_MCU = ((height + 7) / (8 * y_scale)) + min_y_MCU - 1;

    // New stored difference, and recepticles for Huffman processes
    int new_diff = 0;
    decoded output;
    recoded re_size, re_value;

    // BASELINE: Decompress the full raw data together via Huffman Tables
    if (mode == 0xc0)
    {
        // Don't need these here
        free(comps);
        free(table_ids);

        if (sos_start || sos_end != 63 || sa_low)
        {
            fprintf(stderr, "ERROR: Invalid SOS dimensions\n");
            destroy_buffer(results);
            destroy_chunk(c);
            destroy_huffman(dc_huffman_1);
            destroy_huffman(ac_huffman_1);
            destroy_tree(dc_tree_1);
            destroy_tree(ac_tree_1);
            if (color == 3)
            {
                destroy_huffman(dc_huffman_2);
                destroy_huffman(ac_huffman_2);
                destroy_tree(dc_tree_2);
                destroy_tree(ac_tree_2);
            }
            fclose(infile);
            return 13;
        }

        // Write data from SOS to results buffer
        write_data(results, c->type, 2);
        write_data(results, (BYTE *) &(c->size), 2);
        true_l = reverse_bytes(c->size) - 2;
        write_data(results, c->data, true_l);
        destroy_chunk(c);

        // Read in full image data before Huffman decompression
        BYTE eof_check = 1;
        while ((b1 != 0xff || b2 != 0xd9) && eof_check)
        {
            b1 = b2;
            eof_check = fread(&b2, sizeof(BYTE), 1, infile);

            // Remove any padding by moving to the next byte
            if (b1 != 0xff || (b2 != 0x00 && b2 != 0xff))
            {
                write_data(raw_buf, &b2, 1);
            }
        }

        // If ending mark never came, file is corrupt
        if (!eof_check)
        {
            fprintf(stderr, "ERROR: Ending signature not found\n");
            destroy_buffer(results);
            destroy_buffer(raw_buf);
            destroy_huffman(dc_huffman_1);
            destroy_huffman(ac_huffman_1);
            destroy_tree(dc_tree_1);
            destroy_tree(ac_tree_1);
            if (color == 3)
            {
                destroy_huffman(dc_huffman_2);
                destroy_huffman(ac_huffman_2);
                destroy_tree(dc_tree_2);
                destroy_tree(ac_tree_2);
            }
            fclose(infile);
            return 14;
        }

        // Cut out last two bytes (EOI) before decompression
        if (b1 == 0xff && b2 == 0xd9) raw_buf->size -= 2;
        fclose(infile);

        // Set up important variables for iterating
        BYTE edge_y = 0;
        BYTE edge_cb = 0;
        BYTE edge_cr = 0;

        // Current & last kept actual DC values (not stored differences)
        int current_y_dc = 0;
        int current_cb_dc = 0;
        int current_cr_dc = 0;
        int last_y_dc = 0;
        int last_cb_dc = 0;
        int last_cr_dc = 0;

        for (int m_y = 0; m_y < MCU_rows; m_y++)
        {
            for (int m_x = 0; m_x < MCU_cols; m_x++)
            {
                // Check if you've hit a reset interval in reading
                m = m_x + (m_y * MCU_cols);
                if (restart_intr && m && !(m % restart_intr))
                {
                    if (bit_offset)
                    {
                        bit_offset = 0;
                        place++;
                    }

                    if (place[0] == 0xff && place[1] <= 0xd7)
                    {
                        // Skip marker bytes
                        place += 2;
                    }

                    current_y_dc = 0;
                    current_cb_dc = 0;
                    current_cr_dc = 0;
                }

                // Check if you've hit a reset interval in writing
                if (restart_intr && copied_MCUs && !(copied_MCUs % restart_intr))
                {
                    if ((m_y >= min_y_MCU) && (m_y <= max_y_MCU))
                    {
                        if ((m_x >= min_x_MCU) && (m_x <= max_x_MCU))
                        {
                            marker[1] = 0xd0;
                            marker[1] += ((copied_MCUs / restart_intr) - 1) % 8;

                            // First, slip in last cropped byte if necessary
                            if (cropped_offset)
                            {
                                // Standard says remaining bits should be 1s
                                last_cropped = cropped->data + cropped->size - 1;
                                *last_cropped += (0xff >> cropped_offset);
                                write_data(results, last_cropped, 1);
                                cropped_offset = 0;
                            }

                            // Write directly to results buffer
                            write_data(results, marker, 2);

                            last_y_dc = 0;
                            last_cb_dc = 0;
                            last_cr_dc = 0;
                        }
                    }
                }

                // Update MCU starting pointers
                MCU_start = place;
                MCU_bit_offset = bit_offset;

                // LUMINANCE (Y)
                for (cs_ind = 0; cs_ind < y_dims; cs_ind++)
                {
                    // Find Huffman decoding for Luminance (Y) - DC
                    output = decode_size(place, bit_offset, dc_tree_1);
                    place += output.bytes;
                    bit_offset = output.bits;

                    // Feed value bits into decode_value, update DC
                    pull_value(decoded_dc, output, place, bit_offset);
                    current_y_dc += decode_value(output.value, decoded_dc);

                    // Scroll past value category bits
                    bit_offset += output.value;
                    place += bit_offset / 8;
                    bit_offset %= 8;

                    // If cropping here, write the new value in advance
                    if (m_y >= min_y_MCU && m_y <= max_y_MCU)
                    {
                        if (m_x >= min_x_MCU && m_x <= max_x_MCU)
                        {
                            if (edge_y)
                            {
                                // Find and format the new value as a difference
                                new_diff = current_y_dc - last_y_dc;
                                recode_value(new_diff, &re_value);
                                recode_size(re_value.n, dc_tree_1, &re_size);

                                // Check for red flag
                                if (re_size.n == 100)
                                {
                                    // Insert missing value into Huffman table
                                    add_huffman(re_value.n, dc_huffman_1);
                                    destroy_tree(dc_tree_1);
                                    dc_tree_1 = create_tree(dc_huffman_1);
                                    if (rewrite_huffman(results, dc_1_start, 
                                                        dc_huffman_1))
                                    {
                                        fprintf(stderr, 
                                            "ERROR: Incompatible Huffman Table\n");
                                        destroy_buffer(results);
                                        destroy_buffer(raw_buf);
                                        destroy_buffer(cropped);
                                        destroy_huffman(dc_huffman_1);
                                        destroy_huffman(ac_huffman_1);
                                        destroy_tree(dc_tree_1);
                                        destroy_tree(ac_tree_1);
                                        if (color == 3)
                                        {
                                            destroy_huffman(dc_huffman_2);
                                            destroy_huffman(ac_huffman_2);
                                            destroy_tree(dc_tree_2);
                                            destroy_tree(ac_tree_2);
                                        }
                                        return 15;
                                    }

                                    else
                                    {
                                        dc_2_start++;
                                        recode_size(re_value.n, dc_tree_1, &re_size);
                                    }
                                }

                                // Consolidate re_size and re_value into flex
                                consolidate(re_size, re_value, flex);
                                edge_offset = (re_size.n + re_value.n) % 8;
                                mcu_size = (re_size.n + re_value.n + 7) / 8;

                                cropped_start = cropped->size - 1;
                                if (!cropped_offset) cropped_start++;

                                // Fit relevant bits into cropped buffer
                                write_bits(cropped, cropped_offset, flex, 0, 
                                        mcu_size, edge_offset);
                                cropped_offset += 8 + edge_offset;
                                cropped_offset %= 8;

                                // Last byte isn't finished unless offset = 0
                                cropped_end = cropped->size - 1;
                                if (!cropped_offset) cropped_end++;

                                // Escape 0xff values, copy into results buffer
                                for (int i = cropped_start; i < cropped_end; i++)
                                {
                                    write_data(results, cropped->data + i, 1);
                                    if (cropped->data[i] == 0xff)
                                    {
                                        write_data(results, &filler, 1);
                                    }
                                }

                                // Update MCU starting pointers
                                MCU_start = place;
                                MCU_bit_offset = bit_offset;

                                edge_y = 0;
                            }

                            // Regardless, record all the included Y DC values
                            last_y_dc = current_y_dc;
                        }
                    }

                    // Find Huffman decoding for Luminance (Y) - AC
                    seen_acs = 0;
                    while (seen_acs < 63)
                    {
                        // Decode the run-length pair byte
                        output = decode_size(place, bit_offset, ac_tree_1);
                        place += output.bytes;
                        bit_offset = output.bits;

                        // Number of zeros skipped is the first 4 bits
                        seen_acs += output.value >> 4;
                        
                        // Value category bits from end 4 are skipped
                        bit_offset += output.value & 0x0f;
                        place += bit_offset / 8;
                        bit_offset %= 8;
                        
                        // 0xf0 = another 0 seen, 0x00 = end of MCU
                        seen_acs += output.value & 0xff ? 1 : 63;
                    }
                }

                // If it's inside the crop space, copy outstanding component
                if (m_y >= min_y_MCU && m_y <= max_y_MCU)
                {
                    if (m_x >= min_x_MCU && m_x <= max_x_MCU)
                    {
                        mcu_size = place - MCU_start;
                        if (bit_offset) mcu_size++;

                        cropped_start = cropped->size - 1;
                        if (!cropped_offset) cropped_start++;

                        // Fit relevant bits into cropped buffer
                        write_bits(cropped, cropped_offset, MCU_start, 
                                MCU_bit_offset, mcu_size, bit_offset);
                        cropped_offset += 8 + bit_offset - MCU_bit_offset;
                        cropped_offset %= 8;

                        cropped_end = cropped->size - 1;
                        if (!cropped_offset) cropped_end++;

                        // Escape 0xff values, copy into results buffer
                        for (int i = cropped_start; i < cropped_end; i++)
                        {
                            write_data(results, cropped->data + i, 1);
                            if (cropped->data[i] == 0xff)
                            {
                                write_data(results, &filler, 1);
                            }
                        }

                        edge_y = 0;
                    } else edge_y = 1;
                } else edge_y = 1;

                // CHROMINANCE (Cb)
                if (color == 3)
                {
                    // Update MCU starting pointers
                    MCU_start = place;
                    MCU_bit_offset = bit_offset;

                    for (cs_ind = 0; cs_ind < cb_dims; cs_ind++)
                    {
                        // Find Huffman decoding for Chrominance (Cb) - DC
                        output = decode_size(place, bit_offset, dc_tree_2);
                        place += output.bytes;
                        bit_offset = output.bits;

                        // Feed value bits into decode_value, update DC
                        pull_value(decoded_dc, output, place, bit_offset);
                        current_cb_dc += decode_value(output.value, decoded_dc);

                        // Scroll past value category bits
                        bit_offset += output.value;
                        place += bit_offset / 8;
                        bit_offset %= 8;

                        // If cropping here, write the new value in advance
                        if (m_y >= min_y_MCU && m_y <= max_y_MCU)
                        {
                            if (m_x >= min_x_MCU && m_x <= max_x_MCU)
                            {
                                if (edge_cb)
                                {
                                    // Find and format the new value as a difference
                                    new_diff = current_cb_dc - last_cb_dc;
                                    recode_value(new_diff, &re_value);
                                    recode_size(re_value.n, dc_tree_2, &re_size);

                                    // Check for red flag
                                    if (re_size.n == 100)
                                    {
                                        // Insert missing value into Huffman table
                                        add_huffman(re_value.n, dc_huffman_2);
                                        destroy_tree(dc_tree_2);
                                        dc_tree_2 = create_tree(dc_huffman_2);
                                        if (rewrite_huffman(results, dc_2_start, 
                                                            dc_huffman_2))
                                        {
                                            fprintf(stderr, 
                                                    "ERROR: Incompatible Huffman Table\n");
                                            destroy_buffer(results);
                                            destroy_buffer(raw_buf);
                                            destroy_buffer(cropped);
                                            destroy_huffman(dc_huffman_1);
                                            destroy_huffman(ac_huffman_1);
                                            destroy_tree(dc_tree_1);
                                            destroy_tree(ac_tree_1);
                                            destroy_huffman(dc_huffman_2);
                                            destroy_huffman(ac_huffman_2);
                                            destroy_tree(dc_tree_2);
                                            destroy_tree(ac_tree_2);
                                            return 15;
                                        }

                                        else
                                        {
                                            recode_size(re_value.n, dc_tree_2, &re_size);
                                        }
                                    }

                                    // Consolidate re_size and re_value into flex
                                    consolidate(re_size, re_value, flex);
                                    edge_offset = (re_size.n + re_value.n) % 8;
                                    mcu_size = (re_size.n + re_value.n + 7) / 8;

                                    cropped_start = cropped->size - 1;
                                    if (!cropped_offset) cropped_start++;

                                    // Fit relevant bits into cropped buffer
                                    write_bits(cropped, cropped_offset, flex, 0,
                                            mcu_size, edge_offset);
                                    cropped_offset += 8 + edge_offset;
                                    cropped_offset %= 8;

                                    // Last byte isn't finished unless offset = 0
                                    cropped_end = cropped->size - 1;
                                    if (!cropped_offset) cropped_end++;

                                    // Escape 0xff values, copy into results buffer
                                    for (int i = cropped_start; i < cropped_end; i++)
                                    {
                                        write_data(results, cropped->data + i, 1);
                                        if (cropped->data[i] == 0xff)
                                        {
                                            write_data(results, &filler, 1);
                                        }
                                    }

                                    // Update MCU starting pointers
                                    MCU_start = place;
                                    MCU_bit_offset = bit_offset;

                                    edge_cb = 0;
                                }

                                // Regardless, record all the included Cb DC values
                                last_cb_dc = current_cb_dc;
                            }
                        }

                        // Find Huffman decoding for Chrominance (Cb) - AC
                        seen_acs = 0;
                        while (seen_acs < 63)
                        {
                            // Decode the run-length pair byte
                            output = decode_size(place, bit_offset, ac_tree_2);
                            place += output.bytes;
                            bit_offset = output.bits;

                            // Number of zeros skipped is first 4 bits
                            seen_acs += output.value >> 4;
                            
                            // Value category bits from latter 4 are skipped
                            bit_offset += output.value & 0x0f;
                            place += bit_offset / 8;
                            bit_offset %= 8;
                            
                            // 0xf0 = another 0 seen, 0x00 = end of MCU
                            seen_acs += output.value & 0xff ? 1 : 63;
                        }
                    }

                    // If it's inside the crop space, copy outstanding component
                    if (m_y >= min_y_MCU && m_y <= max_y_MCU)
                    {
                        if (m_x >= min_x_MCU && m_x <= max_x_MCU)
                        {
                            if (bit_offset)
                            {
                                mcu_size = place - MCU_start + 1;
                            }
                            else mcu_size = place - MCU_start;

                            cropped_start = cropped->size - 1;
                            if (!cropped_offset) cropped_start++;

                            // Consolidate relevant bits into cropped buffer
                            write_bits(cropped, cropped_offset, MCU_start, 
                                    MCU_bit_offset, mcu_size, bit_offset);
                            cropped_offset += 8 + bit_offset - MCU_bit_offset;
                            cropped_offset %= 8;

                            cropped_end = cropped->size - 1;
                            if (!cropped_offset) cropped_end++;

                            // Escape 0xff values, copy into results buffer
                            for (int i = cropped_start; i < cropped_end; i++)
                            {
                                write_data(results, cropped->data + i, 1);
                                if (cropped->data[i] == 0xff)
                                {
                                    write_data(results, &filler, 1);
                                }
                            }

                            edge_cb = 0;
                        }
                        else edge_cb = 1;
                    }
                    else edge_cb = 1;
                }

                // CHROMINANCE (Cr)
                if (color == 3)
                {
                    // Update MCU starting pointers
                    MCU_start = place;
                    MCU_bit_offset = bit_offset;

                    for (cs_ind = 0; cs_ind < cr_dims; cs_ind++)
                    {
                        // Find Huffman decoding for Chrominance (Cr) - DC
                        output = decode_size(place, bit_offset, dc_tree_2);
                        place += output.bytes;
                        bit_offset = output.bits;

                        // Feed value bits into decode_value, update dc
                        pull_value(decoded_dc, output, place, bit_offset);
                        current_cr_dc += decode_value(output.value, decoded_dc);

                        // Skip over Value Category bits in reading
                        bit_offset += output.value;
                        place += bit_offset / 8;
                        bit_offset %= 8;

                        // If cropping here, write the new value in advance
                        if (m_y >= min_y_MCU && m_y <= max_y_MCU)
                        {
                            if (m_x >= min_x_MCU && m_x <= max_x_MCU)
                            {
                                if (edge_cr)
                                {
                                    // Find and format the new value as a difference
                                    new_diff = current_cr_dc - last_cr_dc;
                                    recode_value(new_diff, &re_value);
                                    recode_size(re_value.n, dc_tree_2, &re_size);

                                    // Check for red flag
                                    if (re_size.n == 100)
                                    {
                                        // Insert missing value into Huffman table
                                        add_huffman(re_value.n, dc_huffman_2);
                                        destroy_tree(dc_tree_2);
                                        dc_tree_2 = create_tree(dc_huffman_2);
                                        if (rewrite_huffman(results, dc_2_start, 
                                                            dc_huffman_2))
                                        {
                                            fprintf(stderr, 
                                                    "ERROR: Incompatible Huffman Table\n");
                                            destroy_buffer(results);
                                            destroy_buffer(raw_buf);
                                            destroy_buffer(cropped);
                                            destroy_huffman(dc_huffman_1);
                                            destroy_huffman(ac_huffman_1);
                                            destroy_tree(dc_tree_1);
                                            destroy_tree(ac_tree_1);
                                            destroy_huffman(dc_huffman_2);
                                            destroy_huffman(ac_huffman_2);
                                            destroy_tree(dc_tree_2);
                                            destroy_tree(ac_tree_2);
                                            return 15;
                                        }

                                        else
                                        {
                                            recode_size(re_value.n, dc_tree_2, &re_size);
                                        }
                                    }

                                    // Consolidate re_size and re_value into flex
                                    consolidate(re_size, re_value, flex);
                                    edge_offset = (re_size.n + re_value.n) % 8;
                                    mcu_size = (re_size.n + re_value.n + 7) / 8;

                                    cropped_start = cropped->size - 1;
                                    if (!cropped_offset) cropped_start++;

                                    // Fit relevant bits into cropped buffer
                                    write_bits(cropped, cropped_offset, flex, 0,
                                            mcu_size, edge_offset);
                                    cropped_offset += 8 + edge_offset;
                                    cropped_offset %= 8;

                                    // Last byte isn't finished unless offset = 0
                                    cropped_end = cropped->size - 1;
                                    if (!cropped_offset) cropped_end++;

                                    // Escape 0xff values, copy into results buffer
                                    for (int i = cropped_start; i < cropped_end; i++)
                                    {
                                        write_data(results, cropped->data + i, 1);
                                        if (cropped->data[i] == 0xff)
                                        {
                                            write_data(results, &filler, 1);
                                        }
                                    }

                                    // Update MCU starting pointers
                                    MCU_start = place;
                                    MCU_bit_offset = bit_offset;

                                    edge_cr = 0;
                                }

                                // Regardless, record all the included Cr DC values
                                last_cr_dc = current_cr_dc;
                            }
                        }

                        // Find Huffman decoding for Chrominance (Cr) - AC
                        seen_acs = 0;
                        while (seen_acs < 63)
                        {
                            // Decode the run-length pair byte
                            output = decode_size(place, bit_offset, ac_tree_2);
                            place += output.bytes;
                            bit_offset = output.bits;

                            // Number of zeros skipped is first 4 bits
                            seen_acs += output.value >> 4;
                            
                            // Value category bits from latter 4 are skipped
                            bit_offset += output.value & 0x0f;
                            place += bit_offset / 8;
                            bit_offset %= 8;
                            
                            // 0xf0 = another 0 seen, 0x00 = end of MCU
                            seen_acs += output.value & 0xff ? 1 : 63;
                        }
                    }

                    // If it's inside the crop space, copy outstanding component
                    if (m_y >= min_y_MCU && m_y <= max_y_MCU)
                    {
                        if (m_x >= min_x_MCU && m_x <= max_x_MCU)
                        {
                            if (bit_offset)
                            {
                                mcu_size = place - MCU_start + 1;
                            }
                            else mcu_size = place - MCU_start;

                            cropped_start = cropped->size - 1;
                            if (!cropped_offset) cropped_start++;

                            // Consolidate relevant bits into cropped buffer
                            write_bits(cropped, cropped_offset, MCU_start, 
                                    MCU_bit_offset, mcu_size, bit_offset);
                            cropped_offset += 8 + bit_offset - MCU_bit_offset;
                            cropped_offset %= 8;

                            cropped_end = cropped->size - 1;
                            if (!cropped_offset) cropped_end++;

                            // Escape 0xff values, copy into results buffer
                            for (int i = cropped_start; i < cropped_end; i++)
                            {
                                write_data(results, cropped->data + i, 1);
                                if (cropped->data[i] == 0xff)
                                {
                                    write_data(results, &filler, 1);
                                }
                            }
                            edge_cr = 0;
                        }
                        else edge_cr = 1;
                    }
                    else edge_cr = 1;
                }
            }
        }
    }

    // PROGRESSIVE: Decompress the raw data in layers
    else if (mode == 0xc2)
    {
        // Organize dimensions and trees for easy access
        BYTE dims[color];
        dims[0] = y_dims;
        if (color == 3)
        {
            dims[1] = cb_dims;
            dims[2] = cr_dims;
        }

        Node dc_trees[2];
        dc_trees[0] = dc_tree_1;
        DHT dc_huffmans[2];
        dc_huffmans[0] = dc_huffman_1;
        Node ac_trees[2];
        ac_trees[0] = ac_tree_1;
        if (color == 3)
        {
            dc_trees[1] = dc_tree_2;
            ac_trees[1] = ac_tree_2;
            dc_huffmans[1] = dc_huffman_2;
        }

        size_t dc_starts[2];
        dc_starts[0] = dc_1_start;
        if (color == 3) dc_starts[1] = dc_2_start;

        // Organize DC value trackers for each component
        int current_dcs[color];
        int last_dcs[color];
        for (int i = 0; i < color; i++)
        {
            current_dcs[i] = 0;
            last_dcs[i] = 0;
        }

        BYTE done = 0;
        while (!done)
        {
            if ((sos_start > sos_end) || (num_comps != 1 && sos_end))
            {
                fprintf(stderr, "ERROR: Invalid SOS dimensions\n");
                free(comps);
                free(table_ids);
                destroy_buffer(results);
                destroy_chunk(c);
                destroy_huffman(dc_huffman_1);
                destroy_huffman(ac_huffman_1);
                destroy_tree(dc_tree_1);
                destroy_tree(ac_tree_1);
                if (color == 3)
                {
                    destroy_huffman(dc_huffman_2);
                    destroy_huffman(ac_huffman_2);
                    destroy_tree(dc_tree_2);
                    destroy_tree(ac_tree_2);
                }
                fclose(infile);
                return 13;
            }

            // Write data from SOS to results buffer
            write_data(results, c->type, 2);
            write_data(results, (BYTE *) &(c->size), 2);
            true_l = reverse_bytes(c->size) - 2;
            write_data(results, c->data, true_l);

            // Read in full scan before Huffman decompression
            BYTE eof_check = 1;
            while ((b1 != 0xff || b2 != 0xda) && eof_check)
            {
                b1 = b2;
                eof_check = fread(&b2, sizeof(BYTE), 1, infile);

                // Remove any padding by moving to the next byte
                if (b1 != 0xff || (b2 != 0x00 && b2 != 0xff))
                {
                    write_data(raw_buf, &b2, 1);
                }
            }
            if (!eof_check) break;

            // Set up important variables for iterating
            BYTE edges[3] = {0, 0, 0};
            BYTE comp, t_comp, e_band;

            // Reset dimensions
            MCU_rows = (height0 + 7) / (8 * y_scale);
            MCU_cols = (width0 + 7) / (8 * x_scale);
            min_x_MCU = x / (8 * x_scale);
            min_y_MCU = y / (8 * y_scale);
            max_x_MCU = ((width + 7) / (8 * x_scale)) + min_x_MCU - 1;
            max_y_MCU = ((height + 7) / (8 * y_scale)) + min_y_MCU - 1;

            // If only 1 component, no Chroma Subsampling
            if (num_comps == 1)
            {
                MCU_cols *= x_scales[comps[0]];
                MCU_rows *= y_scales[comps[0]];
                min_x_MCU *= x_scales[comps[0]];
                min_y_MCU *= y_scales[comps[0]];
                max_x_MCU = (width + 7) * x_scales[comps[0]];
                max_x_MCU = (max_x_MCU / (8 * x_scale)) + min_x_MCU - 1;
                max_y_MCU = (height + 7) * y_scales[comps[0]];
                max_y_MCU = (max_y_MCU / (8 * y_scale)) + min_y_MCU - 1;
            }

            // Track first visit vs. refinement for each coefficient
            BYTE visited[color * 8];
            for (int i = 0; i < color * 8; i++)
            {
                // Each bit in the 64 represents one coefficient
                visited[i] = 0;
            }

            // Iterate, decompress, and crop
            for (int m_y = 0; m_y < MCU_rows; m_y++)
            {
                for (int m_x = 0; m_x < MCU_cols; m_x++)
                {
                    // TODO: Reset interval?

                    // Iterate through each component included in this scan
                    for (int comp_ind = 0; comp_ind < num_comps; comp_ind++)
                    {
                        // Update MCU starting pointers & components
                        MCU_start = place;
                        MCU_bit_offset = bit_offset;
                        comp = comps[comp_ind];
                        cs_dims = num_comps == 1 ? 1 : dims[comp];
                        t_comp = comp > 0 ? 1 : 0;

                        for (cs_ind = 0; cs_ind < cs_dims; cs_ind++)
                        {
                            // DC Scan
                            if (!sos_start && !sos_end)
                            {
                                // Find Huffman decoding for this component
                                output = decode_size(place, bit_offset, dc_trees[t_comp]);
                                place += output.bytes;
                                bit_offset = output.bits;

                                // Feed value bits into decode_value, update DC
                                pull_value(decoded_dc, output, place, bit_offset);
                                current_dcs[comp] += decode_value(output.value, decoded_dc) << sa_low;

                                // Scroll past value category bits
                                bit_offset += output.value;
                                place += bit_offset / 8;
                                bit_offset %= 8;

                                // If cropping here, write the new value in advance
                                if (m_y >= min_y_MCU && m_y <= max_y_MCU)
                                {
                                    if (m_x >= min_x_MCU && m_x <= max_x_MCU)
                                    {
                                        if (edges[comp])
                                        {
                                            // Find and format the new value as a difference
                                            new_diff = current_dcs[comp] - last_dcs[comp];
                                            recode_value(new_diff, &re_value);
                                            recode_size(re_value.n, dc_trees[t_comp], &re_size);

                                            // Check for red flag
                                            if (re_size.n == 100)
                                            {
                                                // Insert missing value into Huffman table
                                                add_huffman(re_value.n, dc_huffmans[t_comp]);
                                                destroy_tree(dc_trees[t_comp]);
                                                dc_trees[t_comp] = create_tree(dc_huffmans[t_comp]);
                                                if (rewrite_huffman(results, dc_starts[t_comp], 
                                                                    dc_huffmans[t_comp]))
                                                {
                                                    fprintf(stderr, 
                                                        "ERROR: Incompatible Huffman Table\n");
                                                    free(comps);
                                                    free(table_ids);
                                                    destroy_buffer(results);
                                                    destroy_buffer(raw_buf);
                                                    destroy_buffer(cropped);
                                                    destroy_huffman(dc_huffman_1);
                                                    destroy_huffman(ac_huffman_1);
                                                    destroy_tree(dc_tree_1);
                                                    destroy_tree(ac_tree_1);
                                                    if (color == 3)
                                                    {
                                                        destroy_huffman(dc_huffman_2);
                                                        destroy_huffman(ac_huffman_2);
                                                        destroy_tree(dc_tree_2);
                                                        destroy_tree(ac_tree_2);
                                                    }
                                                    return 15;
                                                }

                                                else
                                                {
                                                    if (!t_comp) dc_starts[1]++;
                                                    recode_size(re_value.n, dc_trees[t_comp], &re_size);
                                                }
                                            }

                                            // Consolidate re_size and re_value into flex
                                            consolidate(re_size, re_value, flex);
                                            edge_offset = (re_size.n + re_value.n) % 8;
                                            mcu_size = (re_size.n + re_value.n + 7) / 8;

                                            cropped_start = cropped->size - 1;
                                            if (!cropped_offset) cropped_start++;

                                            // Fit relevant bits into cropped buffer
                                            write_bits(cropped, cropped_offset, flex, 0, 
                                                    mcu_size, edge_offset);
                                            cropped_offset += 8 + edge_offset;
                                            cropped_offset %= 8;

                                            // Last byte isn't finished unless offset = 0
                                            cropped_end = cropped->size - 1;
                                            if (!cropped_offset) cropped_end++;

                                            // Escape 0xff values, copy into results buffer
                                            for (int i = cropped_start; i < cropped_end; i++)
                                            {
                                                write_data(results, cropped->data + i, 1);
                                                if (cropped->data[i] == 0xff)
                                                {
                                                    write_data(results, &filler, 1);
                                                }
                                            }

                                            // Update MCU starting pointers
                                            MCU_start = place;
                                            MCU_bit_offset = bit_offset;

                                            edges[comp] = 0;
                                        }

                                        // Regardless, record all the included DC values
                                        last_dcs[comp] = current_dcs[comp];
                                    }
                                }

                                // Mark coefficient as seen
                                set_bit(visited + (8 * comp), 0);
                            }

                            // AC Scan
                            else if (sos_start)
                            {
                                seen_acs = sos_start;
                                while (seen_acs < sos_end)
                                {
                                    // Decode the run-length pair byte
                                    output = decode_size(place, bit_offset, ac_trees[t_comp]);
                                    place += output.bytes;
                                    bit_offset = output.bits;

                                    // Number of zeros skipped is the first 4 bits
                                    seen_acs += output.value >> 4;
                                    
                                    // Value category bits from end 4 are skipped
                                    bit_offset += output.value & 0x0f;
                                    place += bit_offset / 8;
                                    bit_offset %= 8;

                                    // AC First Visit: check for End of Band
                                    e_band = output.value >> 4;
                                    if (!e_band && (output.value & 0xf0) && (e_band < 0x0f))
                                    {
                                        // Skipping 2^n + (n bits read) blocks
                                        seen_acs = 64;
                                        m_x += (int) pow(2, e_band);
                                        m_x += read_bits(place, bit_offset, e_band);
                                        if (m_x >= MCU_cols)
                                        {
                                            m_y += m_x % MCU_cols;
                                            m_x %= MCU_cols;
                                        }
                                        if (m_y >= MCU_rows)
                                        {
                                            fprintf(stderr, 
                                                "ERROR: Invalid End Band value\n");
                                            free(comps);
                                            free(table_ids);
                                            destroy_buffer(results);
                                            destroy_buffer(raw_buf);
                                            destroy_buffer(cropped);
                                            destroy_huffman(dc_huffman_1);
                                            destroy_huffman(ac_huffman_1);
                                            destroy_tree(dc_tree_1);
                                            destroy_tree(ac_tree_1);
                                            if (color == 3)
                                            {
                                                destroy_huffman(dc_huffman_2);
                                                destroy_huffman(ac_huffman_2);
                                                destroy_tree(dc_tree_2);
                                                destroy_tree(ac_tree_2);
                                            }
                                            return 16;
                                        }
                                    }
                                    
                                    // 0xf0 = another 0 seen, 0x00 = end of MCU
                                    else seen_acs += output.value & 0xff ? 1 : 63;
                                }

                                // Mark coefficients as seen
                                for (int i = sos_start; i < sos_end; i++)
                                {
                                    set_bit(visited + (8 * comp) + (i / 8), i % 8);
                                }
                            }
                        }

                        // If it's inside the crop space, copy outstanding
                        if (m_y >= min_y_MCU && m_y <= max_y_MCU)
                        {
                            if (m_x >= min_x_MCU && m_x <= max_x_MCU)
                            {
                                mcu_size = place - MCU_start;
                                if (bit_offset) mcu_size++;

                                cropped_start = cropped->size - 1;
                                if (!cropped_offset) cropped_start++;

                                // Fit relevant bits into cropped buffer
                                write_bits(cropped, cropped_offset, MCU_start, 
                                        MCU_bit_offset, mcu_size, bit_offset);
                                cropped_offset += 8 + bit_offset - MCU_bit_offset;
                                cropped_offset %= 8;

                                cropped_end = cropped->size - 1;
                                if (!cropped_offset) cropped_end++;

                                // Escape 0xff values, copy into results buffer
                                for (int i = cropped_start; i < cropped_end; i++)
                                {
                                    write_data(results, cropped->data + i, 1);
                                    if (cropped->data[i] == 0xff)
                                    {
                                        write_data(results, &filler, 1);
                                    }
                                }

                                edges[comp] = 0;
                            } else edges[comp] = 1;
                        } else edges[comp] = 1;
                    }
                }
            }

            // Read next SOS chunk
            if (b1 == 0xff && b2 == 0xda)
            {
                infile -= 2;
                read_chunk(infile, c);
                num_comps = c->num_comps;
                for (int i = 0; i < num_comps; i++)
                {
                    comps[i] = c->comps[i];
                    table_ids[i] = c->table_ids[i];
                }
                sos_start = c->start;
                sos_end = c->end;
                sa_low = c->successive_appr & 0x0f;
            }
            else done = 1;
        }

        // After everything is done here...
        fclose(infile);
        free(comps);
        free(table_ids);
        destroy_chunk(c);
    }

    // Write last byte (if necessary) & EOI to results, destroy cropped buf
    if (!cropped_offset)
    {
        write_data(results, cropped->data + cropped->size - 1, 1);
    }
    BYTE eoi[2] = {0xff, 0xd9};
    write_data(results, eoi, 2);
    destroy_buffer(raw_buf);
    destroy_buffer(cropped);

    // Clean up Huffman Tables & Huffman Trees
    destroy_huffman(dc_huffman_1);
    destroy_huffman(ac_huffman_1);
    destroy_tree(dc_tree_1);
    destroy_tree(ac_tree_1);
    if (color == 3)
    {
        destroy_huffman(dc_huffman_2);
        destroy_huffman(ac_huffman_2);
        destroy_tree(dc_tree_2);
        destroy_tree(ac_tree_2);
    }

    // Write to the output file
    FILE *outfile = fopen(argv[6], "w");
    if (!outfile)
    {
        fprintf(stderr, "ERROR: No outfile found\n");
        destroy_buffer(results);
        fclose(outfile);
        return 17;
    }
    fwrite(results->data, sizeof(BYTE), results->size, outfile);
    printf("\nCrop Complete!\nNew dimensions: %zux%zu\n", width, height);

    // Clean up
    destroy_buffer(results);
    fclose(outfile);
    return 0;
}
