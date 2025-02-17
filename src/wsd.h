#ifndef WSD_H
#define WSD_H

#include <iostream>
#include <queue>
#include <limits> // for std::numeric_limits
#include <algorithm>
#include <iomanip> // for std::fixed and std::setprecision
#include "robin_hood.h"
#include "utils.h"

typedef robin_hood::unordered_map<int, robin_hood::unordered_map<int, int>> offset_dict; // diagonal map to a new hash table with cost-offset pair
typedef robin_hood::unordered_map<int, bool> is_extend_dict;

struct DemInfo
{

    int tem_idx;
    std::string ref_name;
    std::string strand;
    float identity;
    int st_pos;
    int end_pos;

    DemInfo(const int &tem_idx, const std::string &ref_name, const std::string &strand, const float &identity, const int &st_pos, const int &end_pos) : tem_idx(tem_idx), ref_name(ref_name), strand(strand), identity(identity), st_pos(st_pos), end_pos(end_pos) {}

    DemInfo() = default;
};

inline void print_dem_info_lst(const std::vector<DemInfo> &dem_res)
{
    for (const DemInfo &item : dem_res)
    {
        std::cout << item.tem_idx << "\t" << item.st_pos << "\t" << item.end_pos
                  << "\t" << item.identity << "\t" << item.strand << "\n";
    }
}

class WSD
{
public:
    std::vector<Seq> ref_seq_vec;              // repeat reference sequence
    std::vector<std::string> tem_seq_lst;      // a list load repeat repeat units
    std::vector<std::string> tem_seq_name_lst; // a list load repeat repeat units
    std::vector<int> pelenty_score;            // the pelenty score when decompostion is end
    int mismatch;                              // mismatch score
    int g;                                     // gap score
    int adap_max_cost;                         // start heuristic strategy when cost score larger than this value
    int adap_max_dist;                         // extend this wave when the distance between its offset and best offset less than this value
    int batch_size;                            // set size of the batch in parallelization
    int p_max_tem_len = 0;                     // max length of temlate

    WSD(const std::vector<Seq> &ref_seq_vec, const std::vector<Seq> &s_tem_seqs, const int &mismatch, const int &g, const int &adap_max_cost, const int &adap_max_dist, int &batch_size)
        : ref_seq_vec(ref_seq_vec), p_tem_lst_len(0), mismatch(mismatch), g(g), adap_max_cost(adap_max_cost), adap_max_dist(adap_max_dist), batch_size(batch_size)
    {
        tem_seq_lst.resize(s_tem_seqs.size());
        tem_seq_name_lst.resize(s_tem_seqs.size());
        p_tem_lst_len = static_cast<int>(s_tem_seqs.size());
        for (size_t j = 0; j < tem_seq_lst.size(); ++j)
        {
            tem_seq_lst[j] = s_tem_seqs[j].seq;
            tem_seq_name_lst[j] = s_tem_seqs[j].name;
            int tem_len = static_cast<int>(tem_seq_lst[j].length());
            if (tem_len > p_max_tem_len)
                p_max_tem_len = tem_len;
        }
    }

    void waf_for_decompose(const std::string &ref_seq, const bool &p_adap, const std::string &ref_name, std::vector<DemInfo> &dem_res);
    void wsd_adap(const int &s, const std::string &ref_seq, std::vector<offset_dict> &wave_list, int &current_offset, std::vector<is_extend_dict> &lr_dia_lst);
    void waf_extend(const int &s, const std::string &ref_seq, std::vector<offset_dict> &wave_list, int &current_offset, std::vector<is_extend_dict> &lr_dia_lst);
    void waf_next(const int &s, const std::string &ref_seq, std::vector<offset_dict> &wave_list, std::vector<is_extend_dict> &lr_dia_lst);
    void backtrace(const std::string &ref_seq, std::vector<offset_dict> &wave_list, int &pelenty_score, const std::string &ref_name, const int &pos_offset, std::vector<DemInfo> &dem_res);
    void para_decompose(const bool &p_adap, const int &thread_num); // parallel decomsose

private:
    int p_tem_lst_len;

    inline void print_wave_list(const std::vector<offset_dict> &wave_list)
    {
        for (size_t j = 0; j < tem_seq_lst.size(); ++j)
        {
            for (auto &k_v1 : wave_list[j])
            {
                for (auto &k_v2 : k_v1.second)
                {
                    std::cout << j << "," << k_v1.first << "," << k_v2.first << "," << k_v2.second << "\n";
                }
            }
        }
    }

    inline void init_wave_lst(std::vector<offset_dict> &wave_list, const std::string &ref_seq)
    {
        wave_list.resize(tem_seq_lst.size());
        for (size_t j = 0; j < tem_seq_lst.size(); ++j)
        {
            if (tem_seq_lst[j][0] != ref_seq[0])
                wave_list[j][0][0] = 0;
            else
                wave_list[j][0][0] = 1;
        }
    }

    inline void init_lr_dia_lst(std::vector<is_extend_dict> &lr_dia_lst, const std::string &ref_seq)
    {
        lr_dia_lst.resize(tem_seq_lst.size());
        for (size_t j = 0; j < tem_seq_lst.size(); ++j)
        {
            for (int k = -static_cast<int>(tem_seq_lst[j].length()) + 1; k < static_cast<int>(ref_seq.length()); ++k)
                lr_dia_lst[j][k] = 1;
        }
    }

    inline void print_dem_res(const std::vector<std::vector<DemInfo>> &total_dem_res)
    {
        for (size_t i = 0; i < total_dem_res.size(); ++i)
        {
            for (const DemInfo &block : total_dem_res[i])
            {
                std::cout << block.ref_name << "\t"
                          << tem_seq_name_lst[block.tem_idx] << "\t"
                          << block.strand << "\t"
                          << block.st_pos << "\t"
                          << block.end_pos << "\t"
                          << block.end_pos - block.st_pos << "\t"
                          << block.identity << "\n";
            }
        }
    }
};

inline float cal_identity(const int &edit_dist, const int &ref_len)
{
    return 1 - static_cast<float>(edit_dist) / ref_len;
}

inline std::string get_true_tem_name(const std::string &name)
{
    if (name.back() == '-')
    {
        return name.substr(0, name.length() - 1);
    }
    return name;
}

inline std::string get_strand(const int &tem_idx, const int &tem_lst_len)
{
    // std::cout << "Tem_len: " << tem_lst_len << "," << static_cast<int>(tem_lst_len) / 2 << "," << tem_idx << "\n";
    if (tem_idx >= static_cast<int>(tem_lst_len) / 2)
        return "-";
    return "+";
}

#endif