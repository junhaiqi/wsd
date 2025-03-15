## Getting Started

```bash
# Requirements
WSD runs Linux and requires gcc 9.3.0+.

# Install WSD
git clone https://github.com/junhaiqi/wsd.git
cd wsd && make -j8  # C++11 required to compile

# Run on test data
./wsd example/test_ref.fa -m example/test_tems.fa

# If you have a fasta file $ref.fa that load tandem repeats and a fasta file $tems.fa that load template sequences, you can run this command to obtain decomposition results in a tsv file $out.tsv  
./wsd $ref.fa -m $tems.fa >$out.tsv  

```

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
