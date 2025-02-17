#ifndef UTILS_H
#define UTILS_H

#include <zlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
// #include "edlib.h"

struct Seq
{
    const std::string name;
    std::string seq;
    Seq(const std::string &name, std::string &seq) : name(name), seq(seq) { transform(seq.begin(), seq.end(), seq.begin(), ::toupper); }
    Seq() = default;
};

void read_fa(const char *fa_path, std::vector<Seq> &seqs);

std::string get_rev_comp(const std::string &sequence);

// float cal_edit_identity(const std::string &query, const std::string &ref);

// void get_tems(const std::vector<Seq> &orignal_tems, const Seq &ref_seq, std::vector<Seq> &final_tems); // get the temlate set

#endif