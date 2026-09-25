/**
 * src/twoBitCompressor.cpp
 *
 * DNA/RNA: 2-bit encoding for MASH. Protein: 8-bit residue codes aligned with
 * fourBitCompressor.cpp (MSA protein mode).
 */

#ifndef TWOBITCOMPRESSOR_HPP
#include "twoBitCompressor.hpp"
#endif
#include <cctype>

/** Encode one residue for protein MASH (same numeric mapping as fourBitCompressor, protein). */
static uint64_t encodeProteinResidue(char c) {
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T':
    case 'U': return 3;
    case 'D': return 4;
    case 'E': return 5;
    case 'F': return 6;
    case 'H': return 7;
    case 'I': return 8;
    case 'K': return 9;
    case 'L': return 10;
    case 'M': return 11;
    case 'N': return 12;
    case 'P': return 13;
    case 'Q': return 14;
    case 'R': return 15;
    case 'S': return 16;
    case 'V': return 17;
    case 'W': return 18;
    case 'Y': return 19;
    default:
        return 20;
    }
}

static uint64_t encodeProteinReducedResidue(char c) {
    char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    switch (u) {
    case 'L': case 'V': case 'I': case 'M': case 'C': return 0;  // L group
    case 'A': case 'G': case 'S': case 'T': case 'P': return 1;  // A group
    case 'F': case 'Y': case 'W': return 2;                      // F group
    case 'E': case 'D': case 'N': case 'Q': case 'K': case 'R': case 'H':
        return 3;                                                 // E group
    default:
        return 1; // Unknowns fall back to A group.
    }
}

/* Murphy8 alphabet (8 states):
 * L: LVIMC, A: AG, S: ST, F: FYW, E: EDNQ, K: KRH, P: P, X: unknown.
 */
static uint64_t encodeProteinMurphy8Residue(char c) {
    char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    switch (u) {
    case 'L': case 'V': case 'I': case 'M': case 'C': return 0; // L
    case 'A': case 'G': return 1;                                // A
    case 'S': case 'T': return 2;                                // S
    case 'F': case 'Y': case 'W': return 3;                      // F
    case 'E': case 'D': case 'N': case 'Q': return 4;            // E
    case 'K': case 'R': case 'H': return 5;                      // K
    case 'P': return 6;                                           // P
    default: return 7;                                            // X
    }
}

void twoBitCompressor(
    std::string seq,
    size_t seqLen,
    uint64_t* compressedSeq,
    bool isProtein,
    bool useReducedProtein,
    bool useMurphy8) {
    if (isProtein && useReducedProtein) {
        size_t compressedSeqLen = (seqLen + 31) / 32;
        for (size_t i = 0; i < compressedSeqLen; i++) {
            compressedSeq[i] = 0;
            size_t start = 32 * i;
            size_t end = std::min(seqLen, start + 32);
            uint64_t shift = 0;
            for (size_t j = start; j < end; j++) {
                uint64_t v = encodeProteinReducedResidue(seq[j]);
                compressedSeq[i] |= (v << shift);
                shift += 2;
            }
        }
        return;
    }

    if (isProtein && useMurphy8) {
        size_t compressedSeqLen = (seqLen + 20) / 21;
        for (size_t i = 0; i < compressedSeqLen; i++) {
            compressedSeq[i] = 0;
            size_t start = 21 * i;
            size_t end = std::min(seqLen, start + 21);
            uint64_t shift = 0;
            for (size_t j = start; j < end; j++) {
                uint64_t v = encodeProteinMurphy8Residue(seq[j]);
                compressedSeq[i] |= (v << shift);
                shift += 3;
            }
        }
        return;
    }

    if (isProtein) {
        size_t compressedSeqLen = (seqLen + 7) / 8;
        for (size_t i = 0; i < compressedSeqLen; i++) {
            compressedSeq[i] = 0;
            size_t start = 8 * i;
            size_t end = std::min(seqLen, start + 8);
            uint64_t shift = 0;
            for (size_t j = start; j < end; j++) {
                uint64_t v = encodeProteinResidue(seq[j]);
                compressedSeq[i] |= (v << shift);
                shift += 8;
            }
        }
        return;
    }

    size_t compressedSeqLen = (seqLen + 31) / 32;
    for (size_t i = 0; i < compressedSeqLen; i++) {
        compressedSeq[i] = 0;
        size_t start = 32 * i;
        size_t end = std::min(seqLen, start + 32);
        uint64_t shift = 0;
        for (auto j = start; j < end; j++) {
            uint64_t twoBitVal;
            switch (seq[j]) {
            case 'A':
                twoBitVal = 0;
                break;
            case 'C':
                twoBitVal = 1;
                break;
            case 'G':
                twoBitVal = 2;
                break;
            case 'T':
            case 'U':
                twoBitVal = 3;
                break;
            default:
                twoBitVal = 0;
                break;
            }
            compressedSeq[i] |= (twoBitVal << shift);
            shift += 2;
        }
    }
}