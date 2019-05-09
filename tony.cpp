#include <iostream>
#include <fstream>
#include <zlib.h>
#include <string>
#include "kseq.h"
#include "cmdline.h"
// g++ -std=c++11 main.cpp -lz -o tony

using namespace std;
KSEQ_INIT(gzFile, gzread)

cmdline::parser parameter(int argc, char *argv[]) {
    cmdline::parser opt;
    opt.add<string>("read1", '1', "Required, input read1, it can be compressed or not.", true);
    opt.add<string>("read2", '2', "Required, input read2.", true);
    opt.add<string>("out1", '3', "Required, out read1, compressed.", true);
    opt.add<string>("out2", '4', "Required, out read2.", true);
    opt.add<int>("qual", 'q', "mean quality threshold for discard reads, default is 30.", false, 30);
    opt.add("umi", '\0', "extract UMI sequence, default is NO.");
    opt.add<int>("umiStart", '\0', "start position(0 based) of UMI sequence, required for UMI.", false);
    opt.add<int>("umiLength", '\0', "UMI sequence length at one single read, required for UMI.", false);
    opt.add<string>("connection", '\0', "the character between readsID and UMI, it can be space(S), colon(C), or underline(U), default is colon",
                 false, "C", cmdline::oneof<string>("S", "C", "U"));
    opt.add("disComment", '\0', "for disable reads's comment, default is NO.");
    opt.add<int>("readStart", '\0', "start position of real sequence, required for UMI.", false);
    opt.add<int>("readLength", '\0', "real sequence length, it can be ignored when the sequence reach the end", false);
    opt.add<int>("level", '\0', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);
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
               int start, int length, char connect = '_') {
    string umi, sequence1, sequence2;
    sequence1 = read1->seq.s;
    sequence2 = read2->seq.s;
    umi = sequence1.substr(start, length) +
          connect +
          sequence2.substr(start, length);
    return umi;
}

void write_read(const kseq_t *read, gzFile file,
                bool add_comment, int level, char p,
                string umi, int start, int length) {
    // cpoy from chen
    gzsetparams(file, level, Z_DEFAULT_STRATEGY);
    gzbuffer(file, 1024*1024);
    gzputc(file, '@'); gzputs(file, read->name.s);
    if(umi.length()) {
        gzputc(file, p); gzputs(file, umi.c_str());
    }
    if(add_comment && read->comment.l) {
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
    int level = opt.get<int>("level");
    bool add_comment = !opt.exist("disComment");
    bool treat_umi = opt.exist("umi");
    int umi_start = treat_umi && opt.exist("umiStart") ? opt.get<int>("umiStart"): 0;
    int umi_length = treat_umi && opt.exist("umiLength")? opt.get<int>("umiLength"): 0;
    string connection = opt.get<string>("connection");
    char prefix;
    if(connection == "S") prefix = ' ';
    if(connection == "C") prefix = ':';
    if(connection == "U") prefix = '_';
    int seq_start = treat_umi && opt.exist("readStart")? opt.get<int>("readStart"): 0; //
    int seq_length = opt.exist("readLength")? opt.get<int>("readLength"): 0; //
    if(treat_umi) {
        if(!opt.exist("umiStart") || !opt.exist("umiLength") || !opt.exist("readStart")) {
            cerr << "Error: --umiStart, --umiLength and --readStart are all need for umi treat" << endl;
            return -1;
        }
    } else {
        if(opt.exist("umiStart") || opt.exist("umiLength") || opt.exist("readStart")) {
            cerr << "Error: it seems that you want do umi treat, but lack of --umi" << endl;
            return -1;
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
            UMISeq = treat_umi? get_umi(reads1, reads2, umi_start, umi_length): "";
            // todo multi-thread
            write_read(reads1, out1, add_comment, level, prefix, UMISeq, seq_start, seq1_end);
            write_read(reads2, out2, add_comment, level, prefix, UMISeq, seq_start, seq2_end);
        }
    }
    // destroy read
    kseq_destroy(reads1); kseq_destroy(reads2);
    gzclose(fp1); gzclose(fp2);
    gzclose(out1); gzclose(out2);
    return 0;
}
