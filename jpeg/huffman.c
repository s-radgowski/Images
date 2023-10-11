#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jpeg.h"

// Additional function definitions for jpeg.h

// Return the right node of n in the current level or NULL
static Node right_shift(const Node n)
{
    if (!n)
    {
        return NULL;
    }

    // If it's a left child, just return the right child
    if (n->parent && n->parent->left == n)
    {
        return n->parent->right;
    }

    // Else, traverse back to find the nearest node to the right
    int count = 0;
    Node current = n;
    while (current->parent && current->parent->right == current)
    {
        current = current->parent;
        count++;
    }

    // If you're at the root, you started at the right
    if (!current->parent) return NULL;

    // Cross to right side of tree
    current = current->parent->right;
    while (count > 0 && current)
    {
        current = current->left;
        count--;
    }

    return current;
}

// Create an empty node
static void create_child(Node parent, char c)
{
    Node n = malloc(sizeof(struct _node));
    assert(n);
    n->value = 0;

    // Add this char to the parent's code
    n->code = calloc(parent->length + 1, sizeof(char));
    assert(n->code);
    for (int i = 0; i < parent->length; i++)
    {
        n->code[i] = parent->code[i];
    }
    n->code[parent->length] = c;
    n->length = parent->length + 1;

    n->left = NULL;
    n->right = NULL;
    n->leaf = 0;
    n->parent = parent;

    // Add this node as its parent's child
    if (c == '0')
    {
        parent->left = n;
    }
    else if (c == '1')
    {
        parent->right = n;
    }
    else
    {
        printf("ERROR: Invalid code digit entered\n");
    }
}

// Populate a Huffman Tree from a DHT chunk
Node create_tree(const DHT d)
{
    // Create root with empty left/right children
    Node root = malloc(sizeof(struct _node));
    assert(root);
    root->value = 0;
    root->code = NULL;
    root->length = 0;
    root->leaf = 0;
    root->parent = NULL;

    // Add placeholders for level 1 children
    create_child(root, '0');
    create_child(root, '1');

    // Iterate through Number of Symbols array
    Node current;
    Node leftmost = root->left;
    size_t seen = 0;
    for (int i = 0; i < 16; i++)
    {
        // Check for valid number of symbols
        if (d->counts[i] > pow(2, i + 1))
        {
            fprintf(stderr, "ERROR: Invalid symbol count\n");
            return NULL;
        }

        // No symbols in this layer
        if (d->counts[i] == 0)
        {
            current = leftmost;
            while (current)
            {
                // Add empty children of current
                create_child(current, '0');
                create_child(current, '1');
                current = right_shift(current);
            }

            leftmost = leftmost->left;
        }

        // Symbols in this layer
        else
        {
            // Add symbols, close those branches
            for (int j = seen; j < seen + d->counts[i]; j++)
            {
                leftmost->value = d->symbols[j];
                leftmost->leaf = 1;
                leftmost = right_shift(leftmost);
            }
            seen += d->counts[i];

            // Add empty children of leftmost
            create_child(leftmost, '0');
            create_child(leftmost, '1');
            current = right_shift(leftmost);
            leftmost = leftmost->left;

            // Add empty children to rest of level
            while (current)
            {
                create_child(current, '0');
                create_child(current, '1');
                current = right_shift(current);
            }
        }
    }

    return root;
}

// Print a branch of a Huffman Tree in 2D space
#define BUFFER (2)
static void print_branch(const Node root, size_t space)
{
    // Base case
    if (!root) return;

    space += BUFFER;
    print_branch(root->right, space);

    // Print the current node after space
    printf("\n");
    for (int i = BUFFER; i < space; i++)
    {
        printf(" ");
    }
    printf("%i\n", root->value);

    print_branch(root->left, space);
}

// Print a full Huffman Tree horizontally in 2D space
void print_tree(const Node root)
{
    print_branch(root, 0);
}

// Recursively destroy the branches of a Huffman Tree
void destroy_tree(Node n)
{
    if (!n) return;

    // Free the left branch
    if (n->left) destroy_tree(n->left);

    // Free the right branch
    if (n->right) destroy_tree(n->right);

    free(n->code);
    free(n);
}

// Decode a series of bits using a Huffman Tree
#define INITIAL_CAPACITY (8)
#define CAPACITY_MULTIPLIER (2)
decoded decode_size(const BYTE *str, BYTE bit, const Node root)
{
    // Set up output with bit and byte offset
    decoded output;
    output.bytes = 0;
    output.bits = bit;
    output.value = 0x00;

    // Keep track of the path that you have traversed
    size_t path_size = 0;
    size_t capacity = INITIAL_CAPACITY;
    char *path = calloc(capacity, sizeof(char));
    assert(path);

    Node current = root;
    while (current && !current->leaf)
    {
        // Scale buffer if necessary
        if (path_size + 1 > capacity)
        {
            capacity *= CAPACITY_MULTIPLIER;
            path = realloc(path, capacity);
        }

        // Traverse the tree
        if (BIT_SET(*(str + output.bytes), output.bits))
        {
            current = current->right;
            path[path_size] = '1';
        }
        else
        {
            current = current->left;
            path[path_size] = '0';
        }

        // Increment the bit/byte offsets and path size
        if (output.bits++ == 7) 
        {
            output.bits = 0;
            output.bytes++;
        }
        path_size++;
    }

    if (current && current->leaf)
    {
        // Check if the codes match as we expect
        if (!strncmp(path, current->code, path_size))
        {
            output.value = current->value;
        }
    }

    free(path);
    return output;
}

// Decode value of n bits from b using bitstream table
int decode_value(size_t n, const BYTE *b)
{
    if (!n) return 0;

    // Determine sign from first bit
    BYTE sign = *b & 0x80;

    // Isolate relevant bits
    size_t adjusted;
    if (n < 8)
    {
        adjusted = *b >> (8 - n);
    }
    else if (n == 8)
    {
        adjusted = (int) *b;
    }
    else
    {
        // n has a max value of 15
        adjusted = *(b + 1) >> (16 - n);
        adjusted += *b << (n - 8);
        adjusted &= 0xffff00ff;
        adjusted += (*b >> (16 - n)) << 8;
    }

    // Find distance from either extreme
    int value;
    if (!sign)
    {
        value = -1 * (int) (pow(2, n) - 1);
        value += adjusted;
    }
    else
    {
        value = adjusted;
    }

    return value;
}

// Feed value bits into decoded_dc byte array (dc) for analysis
void pull_value(BYTE *dc, decoded out, BYTE *place, size_t offset)
{
    dc[0] = *place << offset;
    if (out.value > (8 - offset))
    {
        dc[0] += *(place + 1) >> (8 - offset);
    }

    // If it spills into next output byte
    if (out.value > 8)
    {
        dc[1] = *(place + 1) << offset;
        // Pulls from second byte after place
        if (out.value > (16 - offset))
        {
            dc[1] += *(place + 2) >> (8 - offset);
        }
    }
}

// Recursively searches a branch for value n
static Node search_branch(size_t n, const Node root)
{
    // Avoiding segmentation faults
    if (!root) return NULL;

    // Base case: root is our desired value
    if ((root->value == (BYTE) n) && root->leaf)
    {
        return root;
    }

    Node solution = NULL;

    // Search left branch
    if (root->left)
    {
        solution = search_branch(n, root->left);
        if (solution) return solution;
    }

    // Search right branch
    if (root->right)
    {
        solution = search_branch(n, root->right);
        if (solution) return solution;
    }

    return NULL;
}

// Recode an integer into a bitstring using a Huffman Tree
void recode_size(size_t n, const Node root, recoded *r)
{
    // Find the correct node from the tree
    Node solution = search_branch(n, root);
    if (!solution)
    {
        // n > 16 is impossible, so this is a red flag
        r->n = 100;
        return;
    }

    r->n = solution->length;
    r->b[0] = 0x00;
    r->b[1] = 0x00;

    // Plug each bit into b
    if (solution->length < 9)
    {
        for (int i = 0; i < solution->length; i++)
        {
            if (solution->code[i] == '1')
            {
                r->b[0] += pow(2, 7 - i);
            }
        }
    }

    else
    {
        for (int i = 0; i < 8; i++)
        {
            if (solution->code[i] == '1')
            {
                r->b[0] += pow(2, 7 - i);
            }
        }

        for (int i = 8; i < solution->length; i++)
        {
            if (solution->code[i] == '1')
            {
                r->b[1] += pow(2, 15 - i);
            }
        }
    }
}

// Recode value of n bits from b using bitstream table
void recode_value(int value, recoded *r)
{
    if (!value)
    {
        r->n = 0;
        r->b[0] = 0x00;
        r->b[1] = 0x00;
        return;
    }

    // Determine sign
    BYTE sign = value < 0 ? 0 : 1;

    // Determine nearest power of 2
    r->n = (size_t) ceil(log2(abs(value) + 1));

    // Increment as necessary
    if (sign)
    {
        if (r->n < 8)
        {
            r->b[0] = ((BYTE) value) << (8 - r->n);
        }
        else if (r->n == 8)
        {
            r->b[0] = (BYTE) value;
        }
        else
        {
            int temp = ((value << (16 - r->n)) & 0x0000ff00);
            r->b[0] = (BYTE) (temp >> 8);

            temp = ((value << (16 - r->n)) & 0x000000ff);
            r->b[1] = (BYTE) temp;
        }
    }
    else
    {
        int temp2 = -1 * (int) (pow(2, r->n) - 1);
        if (r->n < 8)
        {
            r->b[0] = (BYTE) ((value - temp2) << (8 - r->n));
        }
        else if (r->n == 8)
        {
            r->b[0] = (BYTE) (value - temp2);
        }
        else
        {
            int temp3 = (value - temp2) << (16 - r->n);
            temp3 &= 0x0000ff00;
            r->b[0] = (BYTE) (temp3 >>  8);

            temp3 = (value - temp2) << (16 - r->n);
            temp3 &= 0x000000ff;
            r->b[1] = (BYTE) temp3;
        }
    }
}

// Consolidate recoded values into one "flex" array
void consolidate(recoded re_size, recoded re_value, BYTE *flex)
{
    flex[0] = re_size.b[0];

    // Recoded size spills into second flex byte
    if (re_size.n > 8)
    {
        // Populate remaining re_size and opening re_value bits
        flex[1] = re_size.b[1];
        if (re_size.n == 16)
        {
            flex[2] = re_value.b[0];
            flex[3] = re_value.b[1];
            return;
        }

        flex[1] &= 0xff << (16 - re_size.n);
        flex[1] += re_value.b[0] >> (re_size.n - 8);

        // Add the rest of the recoded value bits
        if (re_value.n > 16 - re_size.n)
        {
            flex[2] = re_value.b[0] << (16 - re_size.n);
            if (re_value.n > 8)
            {
                flex[2] += re_value.b[1] >> (re_size.n - 8);
                flex[3] = re_value.b[1] << (16 - re_size.n);
            }
        }
    }

    // Recoded size fills just the first flex byte
    else if (re_size.n == 8)
    {
        flex[1] = re_value.b[0];
        flex[2] = re_value.b[1];
    }

    // Extra room at end of first flex byte left over
    else
    {
        flex[0] &= 0xff << (8 - re_size.n);
        flex[0] += re_value.b[0] >> re_size.n;

        if (re_value.n > 8 - re_size.n)
        {
            flex[1] = re_value.b[0] << (8 - re_size.n);
            flex[1] += re_value.b[1] >> re_size.n;
            flex[2] = re_value.b[1] << (8 - re_size.n);
        }
    }
}
