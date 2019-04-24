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

//
int average_quality(const kseq_t *read) {
   char *quality = read->qual.s;
   int sumQ = 0;
   for (; *quality; quality++)
       sumQ += *quality - 33; // Pherd33
   return sumQ / read->seq.l;
}
int average_quality(const kseq_t* read, int start) {
    string quality = read->qual.s;
    quality = quality.substr(start); //
    int sumQ = 0, L = quality.length();
    for(int i = 0; i < L; i++)
        sumQ += quality[i] - 33;
    return sumQ / L;
}
int average_quality(const kseq_t* read, int start, int length) {
    string quality = read->qual.s;
    quality = quality.substr(start, length); //
    int sumQ = 0, L = quality.length();
    for(int i = 0; i < L; i++)
        sumQ += quality[i] - 33;
    return sumQ / L;
}

//
string get_umi(const kseq_t *read1, const kseq_t *read2, int start, int length, char connect = ':') {
    string umi, sequence1, sequence2;
    sequence1 = read1->seq.s;
    sequence2 = read2->seq.s;
    umi = sequence1.substr(start, length) + connect + sequence2.substr(start, length);
    return umi;
}

// 
void write_read(const kseq_t *read, gzFile file) {
    string lines = "@";
    lines.append(read->name.s);
    if(read->comment.l) { // comment.s
        lines.append(" ");
        lines.append(read->comment.s);
    }
    lines.append("\n");
    lines.append(read->seq.s);
    lines.append("\n+\n");
    lines.append(read->qual.s);
    lines.append("\n");
    gzputs(file, lines.c_str()); // IO
}
void write_read(const kseq_t *read, gzFile file, string umi, int start) {
    string lines = "@";
    lines.append(read->name.s);
    lines.append(":" + umi); // umi
    if(read->comment.l) {
        lines.append(" ");
        lines.append(read->comment.s);
    }
    lines.append("\n");
    string sequence = read->seq.s;
    lines.append(sequence.substr(start)); //
    lines.append("\n+\n");
    string quality = read->qual.s;
    lines.append(quality.substr(start)); //
    lines.append("\n");
    gzputs(file, lines.c_str());
}
void write_read(const kseq_t *read, gzFile file, string umi, int start, int length) {
    string lines = "@";
    lines.append(read->name.s);
    lines.append(":" + umi);
    if(read->comment.l) {
        lines.append(" ");
        lines.append(read->comment.s);
    }
    lines.append("\n");
    string sequence = read->seq.s;
    lines.append(sequence.substr(start, length)); //
    lines.append("\n+\n");
    string quality = read->qual.s;
    lines.append(quality.substr(start, length)); //
    lines.append("\n");
    gzputs(file, lines.c_str()); //
}

//
int filter(kseq_t *reads1, kseq_t *reads2, gzFile out1, gzFile out2, int cutQ) {
    int q1 = average_quality(reads1);
    if(q1 < cutQ)
        return 1;
    int q2 = average_quality(reads2);
    if(q2 < cutQ)
        return 2;

    write_read(reads1, out1);
    write_read(reads2, out2);
    return 0;
}
int filter(kseq_t *reads1, kseq_t *reads2, gzFile out1, gzFile out2, int cutQ, int umi_start, int umi_length, int seq_start) {
    int q1 = average_quality(reads1, seq_start);
    if(q1 < cutQ)
        return 1;
    int q2 = average_quality(reads2, seq_start);
    if(q2 < cutQ)
        return 2;
    string umi = get_umi(reads1, reads2, umi_start, umi_length, ':');
    write_read(reads1, out1, umi, seq_start);
    write_read(reads2, out2, umi, seq_start);
    return 0;
}
int filter(kseq_t *reads1, kseq_t *reads2, gzFile out1, gzFile out2, int cutQ, int umi_start, int umi_length, int seq_start, int seq_length) {
    int q1 = average_quality(reads1, seq_start, seq_length);
    if(q1 < cutQ)
        return 1;
    int q2 = average_quality(reads2, seq_start, seq_length);
    if(q2 < cutQ)
        return 2;
    string umi = get_umi(reads1, reads2, umi_start, umi_length, ':');
    write_read(reads1, out1, umi, seq_start, seq_length);
    write_read(reads2, out2, umi, seq_start, seq_length);
}

//
int main(int argc, char *argv[]) {
    // parameter
    cmdline::parser opt = parameter(argc, argv);
    int cutQ = opt.get<int>("qual");
    bool treat_umi = opt.exist("umi");
    int umi_start = -1, umi_length = -1, seq_start = -1, seq_length = 0; //#
    if(treat_umi) {
        if(opt.exist("umiStart") && opt.exist("umiLength") && opt.exist("readStart")) {
            umi_start = opt.get<int>("umiStart");
            umi_length = opt.get<int>("umiLength");
            seq_start = opt.get<int>("readStart");
        } else {
            cerr << "Error: --umiStart, --umiLength, --readStart, these 3 are all need for umi treat" << endl;
            return 3; // exit
        }
        if(opt.exist("readLength"))
            seq_length = opt.get<int>("readLength");
    } else {
        if(opt.exist("umiStart") || opt.exist("umiLength") || opt.exist("readStart")) {
            cerr << "Error: it seems that you want do umi treat, but lack of --umi" << endl;
            return 5;
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

    while (kseq_read(reads1) >= 0) {
        kseq_read(reads2);
        if(treat_umi) {
            if(seq_length) // not 0
                filter(reads1, reads2, out1, out2, cutQ, umi_start, umi_length, seq_start, seq_length);
            else
                filter(reads1, reads2, out1, out2, cutQ, umi_start, umi_length, seq_start);
        } else
            filter(reads1, reads2, out1, out2, cutQ);
    }

    // destroy read
    kseq_destroy(reads1);
    kseq_destroy(reads2);
  
    // close file
    gzclose(fp1); gzclose(fp2);
    gzclose(out1); gzclose(out2);
  
    return 0;
}
