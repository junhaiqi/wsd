#include "utils.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

void read_fa(const char *fa_path, std::vector<Seq> &seqs)
{
    gzFile fp;
    kseq_t *ks;
    if ((fp = gzopen(fa_path, "r")) == 0)
    {
        printf("%s not found!\n", fa_path);
        exit(0);
    }

    ks = kseq_init(fp);
    while (kseq_read(ks) >= 0)
    {
        std::string seq = ks->seq.s;
        std::string name = ks->name.s;
        seqs.push_back(Seq(name, seq));
    }
}

std::string get_rev_comp(const std::string &sequence)
{
    std::string complement;
    for (int i = sequence.size() - 1; i >= 0; --i)
    {
        switch (sequence[i])
        {
        case 'A':
            complement += 'T';
            break;
        case 'T':
            complement += 'A';
            break;
        case 'C':
            complement += 'G';
            break;
        case 'G':
            complement += 'C';
            break;
        default:
            std::cerr << "Invalid character in sequence: " << sequence[i] << std::endl;
            return "";
        }
    }
    return complement;
}

// float cal_edit_identity(const std::string &query, const std::string &ref)
// {

//     int edist;
//     EdlibAlignResult result =
//         edlibAlign(query.c_str(), query.length(), ref.c_str(), ref.length(),
//                    edlibDefaultAlignConfig());
//     if (result.status == EDLIB_STATUS_OK)
//     {
//         edist = result.editDistance;
//     }
//     edlibFreeAlignResult(result);
//     return 1 - static_cast<float>(edist) / std::min(static_cast<int>(query.length()), static_cast<int>(ref.length()));
// }

// void get_tems(const std::vector<Seq> &orignal_tems, const Seq &ref_seq, std::vector<Seq> &final_tems)
// {
//     std::vector<Seq> rev_tems;
//     for (size_t i = 0; i < orignal_tems.size(); ++i)
//     {
//         std::string rev_name = orignal_tems[i].name + "-";
//         std::string rev_seq = get_rev_comp(orignal_tems[i].seq);
//         rev_tems.push_back( Seq(rev_name, rev_seq) );
//     }

//     float p_ave_iden;
//     for (const Seq & item : orignal_tems)
//         p_ave_iden += cal_edit_identity(ref_seq.seq, item.seq);
//     p_ave_iden /= orignal_tems.size();

//     float n_ave_iden;
//     for (const Seq & item : rev_tems)
//         n_ave_iden += cal_edit_identity(ref_seq.seq, item.seq);
//     n_ave_iden /= orignal_tems.size();
// }