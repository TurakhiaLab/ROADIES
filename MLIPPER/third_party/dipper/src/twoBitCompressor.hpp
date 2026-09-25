/**
 * src/twoBitCompressor.hpp
 *
 * DNA/RNA: 2-bit encoding (A=0, C=1, G=2, T/U=3) for MASH on unaligned sequences.
 * Output (seqLen+31)/32 uint64_ts; 32 bases per word; ambiguous bases as 0 (A).
 *
 * Protein (default): 8-bit encoding per residue (same mapping as
 * fourBitCompressor.hpp); output (seqLen+7)/8 uint64_ts; 8 residues per word.
 *
 * Protein reduced mode: 2-bit residues using 4 groups:
 *   L: LVIMC, A: AGSTP, F: FYW, E: EDNQKRH
 * output (seqLen+31)/32 uint64_ts; 32 residues per word.
 *
 * Protein Murphy8 mode: 3-bit residues (8 states), 21 residues per uint64_t.
 */

#ifndef TWOBITCOMPRESSOR_HPP
#define TWOBITCOMPRESSOR_HPP

#include <string>
#include <iostream>
#include <tbb/parallel_for.h>
#include "kseq.h"

/** Compress seq[0..seqLen-1]; reduced protein uses 2-bit grouped alphabet. */
void twoBitCompressor(
    std::string seq,
    size_t seqLen,
    uint64_t* compressedSeq,
    bool isProtein = false,
    bool useReducedProtein = false,
    bool useMurphy8 = false);

#endif