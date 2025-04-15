## Getting Started

```bash
# Requirements
WSD runs Linux and requires gcc 9.3.0+.

# Install WSD
git clone https://github.com/junhaiqi/wsd.git
cd wsd && make -j8  # C++11 required to compile

# Run on test data
./wsd example/test_ref.fa -m example/test_tems.fa

# If you have a fasta file $ref.fa that load tandem repeats (TRs) and a fasta file $tems.fa that load template sequences, you can run this command to obtain decomposition results in a tsv file $out.tsv  
./wsd $ref.fa -m $tems.fa >$out.tsv

# Using 36 threads to decompose the ultra-long tandem repeats (centromeres)  
./wsd $ref.fa -m $tems.fa -t36 -A1 >$out.tsv   

```

## Overview of WSD
WSD is a efficient wavefront-based algorithm for string decomposition problem, When decomposing human centromere regions, WSD is, on average, 2 times faster than StringDecomposer (https://academic.oup.com/bioinformatics/article/36/Supplement_1/i93/5870498) and uses two orders of magnitude less memory. All of the benchmarking data is available at https://doi.org/10.6084/m9.figshare.28601708.

## Table of contents
  * [Usage](#usage)
  * [Output](#output)
  * [Acknowledgments](#acknowledgments)
  * [License](#license)
  * [Cite](#cite)

## Usage
```bash
Version 1.0.1
Usage: ./wsd [Options:] <ref.fa> -m <templates.fa>
Options:
  -m STR     Specify the template file path with fasta type [required parameters]
  -c INT     Specify the cost score threshold for fast wave extension [default = 10]
  -d INT     Specify the max offset distance threshold for fast wave extension [default = 100]
  -b INT     Specify the batch size for fast wave extension [default = 1000]
  -a BOOL    Specify whether (1: yes, 0: no) to adaptively determine the batch size [default = 1]
  -t INT     Specify the thread number [default = 1]
  -A BOOL    Specify whether (1: yes, 0: no) to decompose ultralong tandem repeat assemblies [default = 0]
  -M INT     Specify the mismatch penalty score [default = 1]
  -G INT     Specify the gap/insertion penalty score [default = 1]
```

## Output
The WSD output is a tsv file with columns: TR name/template name/alignment direction ("+" for forward, "-" for reverse complement)/TR start/TR end/length/identity score, like:
```bash
1       0       +       0       183     183     0.994536
1       1       +       183     366     183     0.994536
1       2       +       366     549     183     0.994536
1       3       +       549     732     183     0.994536
1       4       +       732     915     183     0.994536
1       5       +       915     1098    183     0.994536
1       6       +       1098    1281    183     0.994536
```

## Acknowledgments
None.

## License 
MIT License.

## Cite
None.
