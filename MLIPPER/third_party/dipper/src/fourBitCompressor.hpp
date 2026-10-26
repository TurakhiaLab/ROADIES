/**
 * src/fourBitCompressor.hpp
 *
 * 4-bit encoding for DNA/RNA (A=0, C=1, G=2, T/U=3; 16 bases per uint64_t) and
 * 8-bit encoding for protein (amino acids 0–19; 8 residues per uint64_t). Used
 * for aligned sequences prior to MSA distance computation. Optional
 * [lowerLimit, upperLimit] slice; caller allocates output buffer.
 */

#ifndef FOURBITCOMPRESSOR_HPP
#define FOURBITCOMPRESSOR_HPP

#include <string>
#include <iostream>
#include <tbb/parallel_for.h>
#include "kseq.h"

/** Compress seq[ lowerLimit..upperLimit (or 0..seqLen-1) ] into compressedSeq.
 *  DNA/RNA: 4 bits/base, (len+15)/16 uint64_ts. Protein: 8 bits/residue, (len+7)/8 uint64_ts.
 *  lowerLimit=0, upperLimit=-1 => full sequence; isProtein=false => DNA/RNA. */
void fourBitCompressor(std::string seq, size_t seqLen, uint64_t* compressedSeq, int lowerLimit=0, int upperLimit=-1, bool isProtein=false);

#endif