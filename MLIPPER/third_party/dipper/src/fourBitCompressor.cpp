/**
 * src/fourBitCompressor.cpp
 *
 * 4-bit encoding for DNA/RNA (16 bases per uint64_t) and 8-bit for protein
 * (8 residues per uint64_t). Used for aligned sequences prior to MSA distance
 * computation. Optional [lowerLimit, upperLimit] slice; ambiguous/unknown
 * encoded as 4 (DNA/RNA) or 20 (protein).
 */

#ifndef FOURBITCOMPRESSOR_HPP
#include "fourBitCompressor.hpp"
#endif

/** Compress seq[ lowerLimit..upperLimit (or 0..seqLen-1) ] into compressedSeq.
 *  DNA/RNA: 4 bits/base, (len+15)/16 uint64_ts. Protein: 8 bits/residue, (len+7)/8 uint64_ts. */
void fourBitCompressor(std::string seq, size_t seqLen, uint64_t* compressedSeq, int lowerLimit, int upperLimit, bool isProtein) {
    /* Resolve effective length and offset for [lowerLimit, upperLimit] slice. */
    if (upperLimit > -1)
        seqLen = upperLimit + 1;
    else
        upperLimit = -1;
    if (lowerLimit > 0)
        seqLen -= lowerLimit;
    else
        lowerLimit = 0;

    size_t compressedSeqLen = isProtein ? (seqLen + 7) / 8 : (seqLen + 15) / 16;

    /* Process each 64-bit word (8 protein residues or 16 DNA/RNA bases). */
    for (size_t i = 0; i < compressedSeqLen; i++) {
        compressedSeq[i] = 0;
        size_t start = isProtein ? 8 * i : 16 * i;
        size_t end = std::min(seqLen, start + (isProtein ? 8 : 16));
        uint64_t shift = 0;

        for (auto j = start; j < end; j++) {
            uint64_t fourBitVal;
            switch (seq[j + lowerLimit]) {
            case 'A': fourBitVal = 0;  break;
            case 'C': fourBitVal = 1;  break;
            case 'G': fourBitVal = 2;  break;
            case 'T':
            case 'U': fourBitVal = 3;  break;
            case 'D': fourBitVal = 4;  break;
            case 'E': fourBitVal = 5;  break;
            case 'F': fourBitVal = 6;  break;
            case 'H': fourBitVal = 7;  break;
            case 'I': fourBitVal = 8;  break;
            case 'K': fourBitVal = 9;  break;
            case 'L': fourBitVal = 10; break;
            case 'M': fourBitVal = 11; break;
            case 'N': fourBitVal = 12; break;
            case 'P': fourBitVal = 13; break;
            case 'Q': fourBitVal = 14; break;
            case 'R': fourBitVal = 15; break;
            case 'S': fourBitVal = 16; break;
            case 'V': fourBitVal = 17; break;
            case 'W': fourBitVal = 18; break;
            case 'Y': fourBitVal = 19; break;
            default:
                fourBitVal = isProtein ? 20 : 4;  /* Ambiguous / unknown. */
                break;
            }
            compressedSeq[i] |= (fourBitVal << shift);
            shift += isProtein ? 8 : 4;
        }
    }
}