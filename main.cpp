#include <iostream>
#include <fstream>
#include <zlib.h>
#include <string>
#include "kseq.h"
#include "cmdline.h"
// g++ -std=c++11 main.cpp -lz -o filter

using namespace std;
KSEQ_INIT(gzFile, gzread)

cmdline::parser parameter(int argc, char *argv[]) {
    cmdline::parser opt;
    opt.add<string>("read1", '1', "input read1, compressed or uncompressed, must be Pherd33", true, "");
    opt.add<string>("read2", '2', "input read2", true, "");
    opt.add<string>("out1", '3', "out read1, compressed", true, "");
    opt.add<string>("out2", '4', "out read2, compressed", true, "");
    opt.add<int>("qual", 'q', "mean read quality cut off, default is 30", false, 30);
    opt.add("umi", 'u', "extract UMI sequence, default false");
    // 10U8B3F*R
    opt.add<int>("umiStart", '\0', "start position of UMI sequence, 0 based", false);
    opt.add<int>("umiLength", '\0', "umi length", false);
    opt.add<int>("readStart", '\0', "start position of real sequence, 0 based", false);
    opt.add<int>("readLength", '\0', "real sequence length, it can be ignored when the sequence reach the end", false);
    opt.parse_check(argc, argv);
    return opt;
}

int average_quality(const kseq_t* read, int start, int length) {
    char *quality = read->qual.s;
    int sumQ = 0, end = start + length;
    for(int i = start; i < end; i++)
        sumQ += quality[i];
    // Pherd33
    return sumQ / length - 33;
}

string get_umi(const kseq_t *read1, const kseq_t *read2,
               int start, int length, char connect = ':') {
    string umi, sequence1, sequence2;
    sequence1 = read1->seq.s;
    sequence2 = read2->seq.s;
    umi = sequence1.substr(start, length) +
          connect +
          sequence2.substr(start, length);
    return umi;
}

void write_read(const kseq_t *read, gzFile file,
                string umi, int start, int length) {
    // inspired by chen
    gzsetparams(file, 4, Z_DEFAULT_STRATEGY); // TODO set for level
    gzbuffer(file, 1024*1024);
    gzputc(file, '@'); gzputs(file, read->name.s);
    if(umi.length()) {
        gzputc(file, ' '); gzputs(file, umi.c_str());
    }
    if(read->comment.l) {
        gzputc(file, ' '); gzputs(file, read->comment.s);
    }
    gzputc(file, '\n');
    gzwrite(file, read->seq.s + start, length);
    gzputs(file, "\n+\n");
    gzwrite(file, read->qual.s + start, length);
    gzputc(file, '\n');
}

int main(int argc, char *argv[]) {
    // parameter logic check
    cmdline::parser opt = parameter(argc, argv);
    int cutQ = opt.get<int>("qual");
    bool treat_umi = opt.exist("umi");
    int umi_start = 0, umi_length = 0;
    int seq_start = 0, seq_length = 0;
    if(treat_umi) {
        if(opt.exist("umiStart") && opt.exist("umiLength") && opt.exist("readStart")) {
            umi_start = opt.get<int>("umiStart");
            umi_length = opt.get<int>("umiLength");
            seq_start = opt.get<int>("readStart");
        } else {
            cerr << "Error: --umiStart, --umiLength, --readStart, these 3 are all need for umi treat" << endl;
            return 3;
        }

        if(opt.exist("readLength"))
            seq_length = opt.get<int>("readLength");
    } else {
        if(opt.exist("umiStart") || opt.exist("umiLength") || opt.exist("readStart")) {
            cerr << "Error: it seems that you want do umi treat, but lack of --umi" << endl;
            return 4;
        }
    }
    // file handler
    gzFile fp1, fp2, out1, out2;
    fp1 = gzopen(opt.get<string>("read1").c_str(), "r");
    fp2 = gzopen(opt.get<string>("read2").c_str(), "r");
    out1 = gzopen(opt.get<string>("out1").c_str(), "wb");
    out2 = gzopen(opt.get<string>("out2").c_str(), "wb");
    // initialize read
    kseq_t *reads1 = kseq_init(fp1);
    kseq_t *reads2 = kseq_init(fp2);

    string UMISeq;
    int seq1_end, seq2_end, q1, q2;
    while (kseq_read(reads1) >= 0) {
        kseq_read(reads2);
        // if not exist readLength, use all remain length of reads
        seq1_end = seq_length? seq_length: reads1->seq.l - seq_start;
        seq2_end = seq_length? seq_length: reads2->seq.l - seq_start;
        // todo multi-thread
        q1 = average_quality(reads1, seq_start, seq1_end);
        q2 = average_quality(reads2, seq_start, seq2_end);
        if(q1 >= cutQ && q2 >= cutQ) {
            UMISeq = treat_umi? get_umi(reads1, reads2, umi_start, umi_length, ':'): "";
            // todo multi-thread
            write_read(reads1, out1, UMISeq, seq_start, seq1_end);
            write_read(reads2, out2, UMISeq, seq_start, seq2_end);
        }
    }
    // destroy read
    kseq_destroy(reads1); kseq_destroy(reads2);
    gzclose(fp1); gzclose(fp2);
    gzclose(out1); gzclose(out2);
    return 0;
}
