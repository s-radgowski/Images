# IMAGE CROPPER PROGRAM
Version 1.1 - February 9th, 2023
**Author: Shaun Radgowski**

## JPEG Version notes:
- Current version only supports Baseline DCT encoding
- Current version only supports bit depths that are divisible by 8
- Current version only supports color modes 1 & 3 (neither 4 nor 5)

## PNG Version notes:
- Current version only supports Compression Method 0
- Current version only supports Interlace Method 0

## Longer-term fixes
- Add_huffman and search_branch functions could be algo-improved
- But, you can actually implement error checking for its return value...

Currently in the process of adding support for Progressive encoding

## Progressive JPEG problems
- For scans with just 1 component, chroma subsampling is erased