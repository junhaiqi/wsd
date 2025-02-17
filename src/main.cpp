#include "wsd.h"
#include "utils.h"
#include "ketopt.h"
#include "version.h"

void test1()
{
    const char *fa_path1 = "example/test_tems.fa";
    const char *fa_path2 = "example/test_ref.fa";
    std::vector<Seq> tem_seqs;
    std::vector<Seq> ref_seqs;
    std::vector<std::string> tem_lst;
    std::string ref_seq;
    read_fa(fa_path1, tem_seqs);
    read_fa(fa_path2, ref_seqs);
    for (const Seq &item : tem_seqs)
        tem_lst.emplace_back(item.seq);
    int c = tem_seqs.size();
    for (size_t i = 0; i < c; ++i)
    {
        std::string rev = get_rev_comp(tem_seqs[i].seq);
        tem_seqs.emplace_back(tem_seqs[i].name, rev);
    }
    const int mismatch = 1;
    const int g = 1;
    const int adap_max_cost = 10;
    const int adap_max_dist = 100;
    int batch_size = 1000;
    WSD ex1(ref_seqs, tem_seqs, mismatch, g, adap_max_cost, adap_max_dist, batch_size);
    ex1.para_decompose(1, 1);
}

int main(int argc, char *argv[])
{
    ketopt_t o = KETOPT_INIT;
    int32_t c;

    const char *tem_fa_path = "";

    int mismatch = 1;
    int g = 1;
    int adap_max_cost = 10;
    int adap_max_dist = 100;
    int batch_size = 1000;
    bool adap_get_bs = 1;
    int thread_num = 1;

    while ((c = ketopt(&o, argc, argv, 1, "m:t:d:c:b:a:M:G:", 0)) >= 0)
    {
        if (c == 'm')
        {
            if (o.arg == 0)
            {
                fprintf(stderr, "Error: -m requires an argument\n");
                return 1;
            }
            tem_fa_path = o.arg;
        }

        else if (c == 'c')
            adap_max_cost = atoi(o.arg);

        else if (c == 'd')
            adap_max_dist = atoi(o.arg);

        else if (c == 'b')
            batch_size = atoi(o.arg);

        else if (c == 'a')
            adap_get_bs = atoi(o.arg);

        else if (c == 't')
            thread_num = atoi(o.arg);

        else if (c == 'M')
        {
            mismatch = atoi(o.arg);
            if (mismatch <= 0)
            {
                fprintf(stderr, "Error: -M requires an argument larger than 0\n");
                return 1;
            }
        }
        else if (c == 'G')
        {
            g = atoi(o.arg);
            if (g <= 0)
            {
                fprintf(stderr, "Error: -G requires an argument larger than 0\n");
                return 1;
            }
        }
    }

    if (argc - o.ind < 1)
    {
        fprintf(stderr, "Version %d.%d.%d\n", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
        fprintf(stderr, "Usage: %s [Options:] <ref.fa> -m <templates.fa>\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -m STR     Specify the template file path with fasta "
                        "type [required parameters]\n");
        fprintf(stderr, "  -c INT     Specify the cost score threshold for fast wave extension [default = 10]\n");
        fprintf(stderr, "  -d INT     Specify the max offset distance threshold for fast wave extension [default = 100]\n");
        fprintf(stderr, "  -b INT     Specify the batch size for fast wave extension [default = 1000] \n");
        fprintf(stderr, "  -a BOOL    Specify whether (1: yes, 0: no) to adaptively determine the batch size [default = 1]\n");
        fprintf(stderr, "  -t INT     Specify the thread number [default = 1]\n");
        fprintf(stderr, "  -M INT     Specify the mismatch penalty score [default = 1] \n");
        fprintf(stderr, "  -G INT     Specify the gap/insertion penalty score [default = 1]\n");
        return 1;
    }

    if (argv[o.ind] == "")
    {
        fprintf(stderr, "Error: requires a repeat reference file path\n");
        return 1;
    }

    if (tem_fa_path == "")
    {
        fprintf(stderr, "Error: requires a template file path\n");
        return 1;
    }

    std::vector<Seq> tem_seqs;
    std::vector<Seq> ref_seqs;
    std::vector<std::string> tem_lst;
    std::string ref_seq;
    read_fa(tem_fa_path, tem_seqs);
    read_fa(argv[o.ind], ref_seqs);
    for (const Seq &item : tem_seqs)
        tem_lst.emplace_back(item.seq);
    int num_tem = tem_seqs.size();
    for (size_t i = 0; i < num_tem; ++i)
    {
        std::string rev = get_rev_comp(tem_seqs[i].seq);
        tem_seqs.emplace_back(tem_seqs[i].name, rev);
    }

    WSD ex1(ref_seqs, tem_seqs, mismatch, g, adap_max_cost, adap_max_dist, batch_size);
    ex1.para_decompose(adap_get_bs, thread_num);

    return 0;
}