#include "wsd.h"
#include <omp.h>

int batch_size_factor = 3;
void WSD::waf_for_decompose(const std::string &ref_seq, const bool &p_adap, const std::string &ref_name, std::vector<DemInfo> &dem_res)
{

    if (p_adap)
        batch_size = p_max_tem_len * batch_size_factor; // `batch_size_factor` can be revised, default is 3 and it should be fine

    int pos_offset = 0;
    // std::cout << batch_size << "\n";
    std::string sub_ref_seq = ref_seq.substr(0, batch_size);
    // int num = 0;
    while (1)
    {
        int sub_current_offset = 0;
        std::vector<offset_dict> sub_wave_list;
        std::vector<is_extend_dict> sub_lr_dia_lst;
        int sub_pelenty_score = 0;
        int s = 0;
        init_wave_lst(sub_wave_list, sub_ref_seq);
        init_lr_dia_lst(sub_lr_dia_lst, sub_ref_seq);
        while (1)
        {
            waf_extend(s, sub_ref_seq, sub_wave_list, sub_current_offset, sub_lr_dia_lst);
            if (sub_current_offset >= static_cast<int>(sub_ref_seq.length()))
                break;
            if (s > adap_max_cost)
            {
                wsd_adap(s, sub_ref_seq, sub_wave_list, sub_current_offset, sub_lr_dia_lst);
            }

            s++;
            waf_next(s, sub_ref_seq, sub_wave_list, sub_lr_dia_lst);
        }
        backtrace(sub_ref_seq, sub_wave_list, sub_pelenty_score, ref_name, pos_offset, dem_res);
        if (dem_res.back().end_pos >= ref_seq.size())
            break;
        else
        {
            // if (dem_res.back().identity < 0.9)
            //     dem_res.pop_back();
            dem_res.pop_back();
            pos_offset = dem_res.back().end_pos;
            if (pos_offset >= ref_seq.size())
                break;
            sub_ref_seq = ref_seq.substr(pos_offset, batch_size);
        }
        // num++;
    }

    // std::cout << "Num: " << num << "\n";
}

void WSD::wsd_adap(const int &s, const std::string &ref_seq, std::vector<offset_dict> &wave_list, int &current_offset, std::vector<is_extend_dict> &lr_dia_lst)
{
    for (size_t j = 0; j < tem_seq_lst.size(); ++j)
    {
        for (int k = -static_cast<int>(tem_seq_lst[j].length()) + 1; k < static_cast<int>(ref_seq.length()); ++k)
        {
            if (lr_dia_lst[j][k] == 1)
            {
                if (wave_list[j][k].find(s) != wave_list[j][k].end())
                {
                    if (wave_list[j][k][s] == std::min(k + static_cast<int>((tem_seq_lst[j].length())), static_cast<int>(ref_seq.length())))
                    {
                        lr_dia_lst[j][k] = 0;
                        continue;
                    }

                    if (current_offset - wave_list[j][k][s] > adap_max_dist)
                        lr_dia_lst[j][k] = 0;
                }
            }
        }
    }
}

void WSD::waf_extend(const int &s, const std::string &ref_seq, std::vector<offset_dict> &wave_list, int &current_offset, std::vector<is_extend_dict> &lr_dia_lst)
{
    std::queue<std::pair<int, int>> q; // Queue to hold (j, k) for wave extension

    for (size_t j = 0; j < tem_seq_lst.size(); ++j)
    {
        if (current_offset >= static_cast<int>(ref_seq.length()))
        {
            return;
        }

        for (int k = -static_cast<int>(tem_seq_lst[j].length()) + 1; k < static_cast<int>(ref_seq.length()); ++k)
        {
            if (lr_dia_lst[j][k] && wave_list[j][k].find(s) != wave_list[j][k].end())
            {
                int pos_in_ref = wave_list[j][k][s];
                current_offset = std::max(pos_in_ref, current_offset);
                int pos_in_tem = wave_list[j][k][s] - k;

                if (current_offset >= static_cast<int>(ref_seq.length()))
                {
                    return;
                }

                if (pos_in_tem < static_cast<int>(tem_seq_lst[j].length()))
                {
                    while (ref_seq[pos_in_ref] == tem_seq_lst[j][pos_in_tem])
                    {
                        if (pos_in_ref < static_cast<int>(ref_seq.length()) && pos_in_tem < static_cast<int>(tem_seq_lst[j].length()))
                        {
                            wave_list[j][k][s]++;
                            pos_in_ref++;
                            pos_in_tem++;
                        }
                        else
                            break;
                    }
                }

                if (wave_list[j][k][s] == std::min(k + static_cast<int>(tem_seq_lst[j].length()), static_cast<int>(ref_seq.length())))
                {
                    q.push({j, k});
                    lr_dia_lst[j][k] = 0;
                }
            }
        }
    }

    while (!q.empty())
    {
        auto wave_idx = q.front();
        q.pop();
        int extend_i = wave_list[wave_idx.first][wave_idx.second][s];

        if (extend_i >= static_cast<int>(ref_seq.length()))
        {
            return;
        }

        for (size_t j = 0; j < tem_seq_lst.size(); ++j)
        {
            if (current_offset >= static_cast<int>(ref_seq.length()))
            {
                return;
            }

            if (ref_seq[extend_i] == tem_seq_lst[j][0])
            {
                if (wave_list[j].find(extend_i) == wave_list[j].end() || wave_list[j][extend_i].empty())
                {
                    wave_list[j][extend_i][s] = extend_i;
                    int pos_in_ref = wave_list[j][extend_i][s];
                    current_offset = std::max(pos_in_ref, current_offset);
                    int pos_in_tem = wave_list[j][extend_i][s] - extend_i;

                    if (pos_in_tem != static_cast<int>(tem_seq_lst[j].length()) - 1 && current_offset < static_cast<int>(ref_seq.length()))
                    {
                        while (ref_seq[pos_in_ref] == tem_seq_lst[j][pos_in_tem])
                        {
                            if (pos_in_ref < static_cast<int>(ref_seq.length()) && pos_in_tem < static_cast<int>(tem_seq_lst[j].length()))
                            {
                                wave_list[j][extend_i][s]++;
                                pos_in_ref++;
                                pos_in_tem++;
                            }
                            else
                                break;
                        }
                    }

                    if (wave_list[j][extend_i][s] == std::min(extend_i + static_cast<int>(tem_seq_lst[j].length()), static_cast<int>(ref_seq.length())))
                    {
                        q.push({j, extend_i});
                        lr_dia_lst[j][extend_i] = false;
                    }

                    current_offset = std::max(pos_in_ref, current_offset);
                }
            }
            else
            {
                if (wave_list[j][extend_i].empty())
                {
                    wave_list[j][extend_i][s + mismatch] = extend_i + 1;
                }
            }
        }
    }
}

void WSD::waf_next(const int &s, const std::string &ref_seq, std::vector<offset_dict> &wave_list, std::vector<is_extend_dict> &lr_dia_lst)
{
    for (size_t j = 0; j < tem_seq_lst.size(); ++j)
    {
        for (int k = -static_cast<int>(tem_seq_lst[j].length()) + 1; k < static_cast<int>(ref_seq.length()); ++k)
        {
            if (lr_dia_lst[j][k] && wave_list[j][k].find(s) == wave_list[j][k].end())
            {
                // If the wave is not yet added at (j, k, s), handle the next wave extension

                // Initialize mismatch, insertion, and gap values
                int inter_mismatch = -1;
                if (wave_list[j][k].find(s - mismatch) != wave_list[j][k].end())
                {
                    inter_mismatch = wave_list[j][k][s - mismatch]; // Mismatch from previous position
                }

                int insertion = -1;
                if (wave_list[j].find(k + 1) != wave_list[j].end() && wave_list[j][k + 1].find(s - g) != wave_list[j][k + 1].end())
                {
                    insertion = wave_list[j][k + 1][s - g]; // Insertion from the next diagonal
                }

                int gap = -1;
                if (wave_list[j].find(k - 1) != wave_list[j].end() && wave_list[j][k - 1].find(s - g) != wave_list[j][k - 1].end())
                {
                    gap = wave_list[j][k - 1][s - g]; // Gap from the previous diagonal
                }

                // Compute the maximum wave limit
                int max_wave = k + static_cast<int>(tem_seq_lst[j].length());
                if (inter_mismatch != -1 || insertion != -1 || gap != -1)
                {
                    // If there is a mismatch, insertion, or gap, calculate the score
                    int score = std::min(std::max(insertion, std::max(gap + 1, inter_mismatch + 1)), max_wave);
                    wave_list[j][k][s] = score; // Store the score for the current wave position
                }
            }
        }
    }
}

void WSD::backtrace(const std::string &ref_seq, std::vector<offset_dict> &wave_list, int &pelenty_score, const std::string &ref_name, const int &pos_offset, std::vector<DemInfo> &dem_res)
{
    int _j = 0, _k = 0;
    int cost = std::numeric_limits<double>::infinity();
    // Find the initial wave that matches the reference sequence length
    for (size_t j = 0; j < tem_seq_lst.size(); ++j)
    {
        for (int k = -static_cast<int>(tem_seq_lst[j].length()) + 1; k < static_cast<int>(ref_seq.length()); ++k)
        {
            for (auto &s : wave_list[j][k])
            {
                if (s.second == static_cast<int>(ref_seq.length()))
                {
                    if (s.first < cost)
                    {
                        cost = s.first;
                        pelenty_score = s.first;
                        _j = j;
                        _k = k;
                    }
                }
            }
        }
    }
    int pos_in_ref = static_cast<int>(ref_seq.length()) - 1;
    int pos_in_tem = pos_in_ref - _k;
    int idx_tem = _j;
    std::pair<int, int> out_pair = {idx_tem, -1};
    std::vector<std::pair<int, int>> out_pair_list;

    std::vector<int> edit_dist_lst;
    int edit_dist = 0;

    while (pos_in_ref > 0)
    {
        if (pos_in_tem != 0)
        {
            if (wave_list[idx_tem][_k - 1].find(cost - g) != wave_list[idx_tem][_k - 1].end())
            {
                pos_in_ref--;
                cost -= g;
                _k--;
                edit_dist++;
            }
            else if (wave_list[idx_tem][_k + 1].find(cost - g) != wave_list[idx_tem][_k + 1].end())
            {
                pos_in_tem--;
                cost -= g;
                _k++;
                edit_dist++;
            }
            else if (wave_list[idx_tem][_k].find(cost - mismatch) != wave_list[idx_tem][_k].end())
            {
                pos_in_ref--;
                pos_in_tem--;
                cost -= mismatch;
                edit_dist++;
            }
            else
            {
                pos_in_ref -= pos_in_tem;
                pos_in_tem = 0;
            }
        }
        else
        {
            if (wave_list[idx_tem][_k - 1].find(cost - g) != wave_list[idx_tem][_k - 1].end())
            {
                pos_in_ref--;
                cost -= g;
                edit_dist++;
            }
            else
            {
                out_pair = {idx_tem, pos_in_ref};
                out_pair_list.push_back(out_pair);
                edit_dist_lst.push_back(edit_dist);
                edit_dist = 0;

                if (tem_seq_lst[idx_tem][0] != ref_seq[pos_in_ref])
                {
                    cost -= mismatch;
                }

                for (size_t j = 0; j < tem_seq_lst.size(); ++j)
                {
                    _k = (pos_in_ref - 1) - (static_cast<int>(tem_seq_lst[j].length()) - 1);
                    if (wave_list[j][_k].find(cost) != wave_list[j][_k].end() && wave_list[j][_k][cost] == static_cast<int>(tem_seq_lst[j].length() + _k))
                    {
                        idx_tem = j;
                        pos_in_tem = static_cast<int>(tem_seq_lst[j].length()) - 1;
                        pos_in_ref--;
                        out_pair = {idx_tem, -1};
                        break;
                    }
                }
            }
        }
    }

    edit_dist_lst.push_back(edit_dist);
    if ( out_pair_list.empty() )
    {
        int ref_len = static_cast<int>(ref_seq.length());
        std::string strand = get_strand(out_pair.first, p_tem_lst_len);
        dem_res.push_back(DemInfo(out_pair.first, ref_name, strand, cal_identity(edit_dist_lst.back(), ref_len), pos_offset + pos_in_ref, pos_offset + static_cast<int>(ref_seq.length())));
        return;
    }

    // Reverse the output pair list
    std::reverse(out_pair_list.begin(), out_pair_list.end());
    std::reverse(edit_dist_lst.begin(), edit_dist_lst.end());
    // DemInfo tmp;

    std::string strand = get_strand(out_pair.first, p_tem_lst_len);
    dem_res.push_back(DemInfo(out_pair.first, ref_name, strand, cal_identity(edit_dist_lst[0], out_pair_list[0].second), pos_offset, pos_offset + out_pair_list[0].second));

    for (size_t i = 0; i < out_pair_list.size() - 1; ++i)
    {
        std::string strand = get_strand(out_pair_list[i].first, p_tem_lst_len);
        int ref_len = out_pair_list[i + 1].second - out_pair_list[i].second;
        dem_res.push_back(DemInfo(out_pair_list[i].first, ref_name, strand, cal_identity(edit_dist_lst[i], ref_len), pos_offset + out_pair_list[i].second, pos_offset + out_pair_list[i + 1].second));
    }

    int ref_len = static_cast<int>(ref_seq.length()) - out_pair_list.back().second;
    strand = get_strand(out_pair_list.back().first, p_tem_lst_len);
    dem_res.push_back(DemInfo(out_pair_list.back().first, ref_name, strand, cal_identity(edit_dist_lst.back(), ref_len), pos_offset + out_pair_list.back().second, pos_offset + static_cast<int>(ref_seq.length())));

    // Print results
    // std::cout << "TemplateName\tStartPositon\tEndPositon\tRegionLength\tIdentity\n";
    // std::cout << out_pair.first << "\t" << "0" << "\t" << out_pair_list[0].second << "\t" << out_pair_list[0].second << "\t" << cal_identity(edit_dist_lst[0], out_pair_list[0].second) << "\n";
    // for (size_t i = 0; i < out_pair_list.size() - 1; ++i)
    // {
    //     int ref_len = out_pair_list[i + 1].second - out_pair_list[i].second;
    //     std::cout << out_pair_list[i].first << "\t" << out_pair_list[i].second << "\t" << out_pair_list[i + 1].second
    //               << "\t" << out_pair_list[i + 1].second - out_pair_list[i].second << "\t" << cal_identity(edit_dist_lst[i], ref_len) << "\n";
    // }
    // int ref_len = static_cast<int>(ref_seq.length()) - out_pair_list.back().second;
    // std::cout << out_pair_list.back().first << "\t" << out_pair_list.back().second << "\t" << static_cast<int>(ref_seq.length())
    //           << "\t" << static_cast<int>(ref_seq.length()) - out_pair_list.back().second << "\t" << cal_identity(edit_dist_lst.back(), ref_len) << "\n";
}

// int batch_size_factor = 3;
// float batch_overlap_factor = 0.1; // it decides the overlap length between two batch
void WSD::para_decompose(const bool &p_adap, const int &thread_num)
{

    std::vector<std::vector<DemInfo>> total_dem_res(ref_seq_vec.size());
#pragma omp parallel num_threads(thread_num)
    {
#pragma omp for
        for (size_t i = 0; i < ref_seq_vec.size(); ++i)
        {
            std::vector<DemInfo> dem_res;
            WSD::waf_for_decompose(ref_seq_vec[i].seq, p_adap, ref_seq_vec[i].name, dem_res);
            total_dem_res[i] = dem_res;
        }
    }
    print_dem_res(total_dem_res);
}