#include <iostream>
#include <fstream>
#include <zlib.h>
#include "kseq.h"
#include "cmdline.h"
// g++ -std=c++11 main.cpp -lz -o filter

using namespace std;
// declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread);

void write_seq(const kseq_t *seq, gzFile file) {
  gzputc(file, '@');
  gzputs(file, seq->name.s);
  // some read has no comment.s
  if (seq->comment.s) {
    gzputc(file, ' ');
    gzputs(file, seq->comment.s);
  }
  gzputc(file, '\n');
  gzputs(file, seq->seq.s);
  gzputs(file, "\n+\n");
  gzputs(file, seq->qual.s);
  gzputc(file, '\n');
}

cmdline::parser parameter(int argc, char *argv[]) {
  cmdline::parser opt;
  opt.add<string>("read1", '1', "input read1, compressed or not", true, "");
  opt.add<string>("read2", '2', "input read2, compressed or not", true, "");
  opt.add<string>("out1", '3', "out read1, compressed", true, "");
  opt.add<string>("out2", '4', "out read2, compressed", true, "");
  opt.add<int>("qual", 'q', "mean reads quality cut off, default is 30", false, 30);
  opt.parse_check(argc, argv);
  return opt;
}

int read_check(const kseq_t *seq) {
  char *quality = seq->qual.s;
  int Q, sumQ = 0;
  for (; *quality; quality++)
    // hard code here
    sumQ += *quality - 33;
  Q = sumQ / seq->seq.l;
  return Q;
}

bool pair_check(const kseq_t *seq1, const kseq_t *seq2, int Q) {
  int q1 = read_check(seq1);
  int q2 = read_check(seq2);
  if (q1 >= Q and q2 >= Q)
    return true;
  else
    return false;
}

int main(int argc, char *argv[]) {
  cmdline::parser opt = parameter(argc, argv);
  gzFile fp1, fp2, out1, out2;
  // open the file handler
  fp1 = gzopen(opt.get<string>("read1").c_str(), "r");
  fp2 = gzopen(opt.get<string>("read2").c_str(), "r");
  out1 = gzopen(opt.get<string>("out1").c_str(), "wb");
  out2 = gzopen(opt.get<string>("out2").c_str(), "wb");
  int qual = opt.get<int>("qual");
  // initialize seq
  kseq_t *seq1 = kseq_init(fp1);
  kseq_t *seq2 = kseq_init(fp2);
  // read sequence
  while (kseq_read(seq1) >= 0) {
      kseq_read(seq2);
      bool qc = pair_check(seq1, seq2, qual);
      if (qc) {
         write_seq(seq1, out1);
         write_seq(seq2, out2);
      }
  }
  // destroy seq
  kseq_destroy(seq1);
  kseq_destroy(seq2);
  // close the file handler
  gzclose(fp1); gzclose(fp2);
  gzclose(out1); gzclose(out2);

  return 0;
}
