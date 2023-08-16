#include <iostream>
#include <boost/program_options.hpp>
#include <time.h>
#include "kseq.h"
#include "zlib.h"

// For parsing the command line values
namespace po = boost::program_options;


bool hasN(std::string str) {
    for (auto c: str) {
        if ((c != 'A') && 
                (c != 'C') &&
                (c != 'G') &&
                (c != 'T') &&
                (c != 'a') &&
                (c != 'c') &&
                (c != 'g') &&
                (c != 't'))
            return true;
    }
    return false;
}

bool upperExceedsThresh(std::string str, size_t h) {
    size_t upperCount = 0;
    for (auto c: str) {
        if ((c == 'A') || 
                (c == 'C') ||
                (c == 'G') ||
                (c == 'T'))
            upperCount++;
    }
    return (upperCount > h);
}

// For reading in the FASTA file
KSEQ_INIT2(, gzFile, gzread)

int main(int argc, char** argv) {

    std::string inFilename, outFilename;
    size_t l, s, e;
    float t;

    // Parse the command line options
    po::options_description desc{"Options"};
    desc.add_options()
    ("input,i", po::value<std::string>(&inFilename)->required(), "Input FASTA file name [REQUIRED].")
    ("output,o", po::value<std::string>(&outFilename)->required(), "Output FASTA file name [REQUIRED].")
    ("len,l", po::value<size_t>(&l)->default_value(1000), "Length of segments.")
    ("start,s", po::value<size_t>(&s)->required(), "start index.")
    ("end,e", po::value<size_t>(&e)->required(), "end index.")
    ("upper-case threshold,t", po::value<float>(&t)->required(), "Lower-case threshold.")
    ("help,h", "Print help messages");

    po::options_description allOptions;
    allOptions.add(desc);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(allOptions).run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << desc << std::endl;
        exit(1);
    }

    std::vector<std::string> record_names;
    std::vector<std::string> record_seqs;
    std::vector<size_t> record_lens;

    size_t total_length = 0;

    // Read input sequence as kseq_t object
    //fprintf(stdout, "Reading input sequence.\n");
    gzFile fp = gzopen(inFilename.c_str(), "r");
    if (!fp) {
        fprintf(stdout, "ERROR: Cannot open file: %s\n", inFilename.c_str());
        exit(1);
    }
    kseq_t *record = kseq_init(fp);
    while (kseq_read(record) >= 0) {
	    record_names.push_back(std::string(record->name.s, record->name.l));
	    record_seqs.push_back(std::string(record->seq.s, record->seq.l));
	    record_lens.push_back(record->seq.l);
        total_length+=record->seq.l;
    }
    gzclose(fp);
    
    size_t num_samples = (e > s) ? e-s+1 : 0;
    if (num_samples < 1) {
	    fprintf(stderr, "Invalid index\n");
	    exit(1);
    }

    fprintf(stdout, "Number of regions: %ld\n", num_samples);
    fprintf(stdout, "ID START: %ld, ID END: %ld\n", s, e);
    fprintf(stdout, "Region length: %ld\n", l);
    fprintf(stdout, "Input file: %s\n", inFilename.c_str());
    fprintf(stdout, "Output file: %s\n", outFilename.c_str());

    size_t threshold = size_t(t*l);

    FILE *f_wr = fopen(outFilename.c_str(), "w");

    size_t count = 0;

    srand(time(0));
    for (size_t j=0; j<num_samples; j++) {
        bool found = false;
        std::string sample;
        while (!found) {
            size_t k = rand() % total_length;
            size_t loc=0, last_idx=0, idx=0, start=0;
            for (auto i: record_lens) {
                idx += i;
                if (k > idx) {
                    loc++;
                }
                else {
                    start = k-last_idx;
                    break;
                }
                last_idx = idx;
            }
            if (start+l > size_t(record_lens[loc])) {
                count++;
                continue;
            }
            std::string frag = record_seqs[loc].substr(start, l);
            if (hasN(frag) || !upperExceedsThresh(frag, threshold)) {
                count++;
                continue;
            }
            fprintf(f_wr, ">gene_%ld\n%s\n", s+j, frag.c_str());
            found = true;
        }
    }

    fclose(f_wr);
    
    fprintf(stdout, "Number of resampling: %ld\n", count);

    return 0;
}

