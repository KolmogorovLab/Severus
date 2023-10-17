# Severus

Severus is a tool to call somatic and germline structural variations (SV) in single or multi sample long-read sequencing data. Severus tested for both Pacbio and Oxford Nanopore read data ([see results](## Why use Severus?)).  
Severus builds breakpoint graphs for one or multiple long-read cancer samples to cliuster complex SVs. For more information, [see](docs/README.md)

<p align="center">
  <img src="docs/severus_first.png" alt="Severus overview"/>
</p>



## Installation

Requirements:
* Python3
* networkx
* numpy
* pydot
* graphviz
* pysam

The easiest way to install dependencies is through conda.

```
conda env create --name <env> --file environment.yml
```

## Quick Start

### Somatic SV calling

```
# Single sample somatic SV calling

./severus.py --target-bam phased_tumor.bam --control-bam phased_normal.bam --out-dir severus_out -t 16 --phasing-vcf phased.vcf --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
dot -Tsvg -O severus_out/breakpoint_graph.dot

# Multisample somatic SV calling

./severus.py --target-bam phased_tumor1.bam phased_tumor2.bam --control-bam phased_normal.bam --out-dir severus_out -t 16 --phasing-vcf phased.vcf --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
dot -Tsvg -O severus_out/breakpoint_graph.gv

```
Providing phased bam files and phasing vcf is optional but recommended. For somatic SV calling single control file is supported. 

### Germline SV calling

```
# Single sample SV calling

./severus.py --target-bam phased_tumor.bam --out-dir severus_out -t 16 --phasing-vcf phased.vcf --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
dot -Tsvg -O severus_out/breakpoint_graph.gv

# Multisample SV calling

./severus.py --target-bam phased_tumor1.bam phased_tumor2.bam --out-dir severus_out -t 16 --phasing-vcf phased.vcf --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
dot -Tsvg -O severus_out/breakpoint_graph.gv
```

Providing phased bam files and phasing vcf is optional but recommended.

## Important parameters

### Required

* `--target-bam` path to one or multiple target bam files (must be indexed) 

* `--out-dir` path to output directory

### Optional 

* `--control-bam` path to one or multiple control bam files (must be indexed)

* `--vntr-bed` path to bed file for tandem repeat regions (must be ordered)

* `--phasing-vcf` path to vcf file used for phasing (must for the haplotype specific SV calling)

* `--threads` number of threads [8]
  
* `--min-support` minimum number of reads supporting a breakpoint [3]

* `--TIN-ratio` tumor in normal ratio [0.01]

For all parameters, see [parameter list](docs/README.md)
## Outputs

### breakpoint_graph.gv  

The primary output is the breakpoint graph, like on the example above. Solid edges correspond to the fragments of the reference genome, (L: length C: coverage)
and dashed colored edges correspond to non-reference connections from reads (R: number of support reads). Each breakpoint is defined by its coordinate.

```
# To convert gv format to svg
dot -Tsvg -O severus_out/breakpoint_graph.gv
```

### VCF file

If phased bam and phasing vcf is provided haplotype specific SV calls are reported as `0|1` or `1|0`.

### breakpoints_double.csv

Detailed information for all breakpoints detected in any of the bam files provided.

## Why use Severus?

### Germline benchmarking results using HG002

We compared performance of Severus, [sniffles2](https://github.com/fritzsedlazeck/Sniffles) and [cuteSV](https://github.com/tjiangHIT/cuteSV) in [HG002 GIAB SV benchmark set](https://www.nature.com/articles/s41587-020-0538-8).  

|SV Caller| TP | FP | FN | Recall | F1 score |
|---------|----|----|----|--------|----------|
| Severus |9453| 345| 402| 0.965| 0.962|
| sniffles2|9459| 336| 396| 0.966| 0.963|
| cuteSV |9231| 676| 624| 0.937| 0.934|

### Somatic benchmarking results COLO829

We compared the performance of existing somatic SV callers [nanomonSV](https://github.com/friend1ws/nanomonsv), [SAVANA](https://github.com/cortes-ciriano-lab/savana) and [sniffles2](https://github.com/fritzsedlazeck/Sniffles) in mosaic mode using COLO829 cell line data against multi-platform [Valle-Inclan et al. truthset](https://www.sciencedirect.com/science/article/pii/S2666979X22000726). 

#### Oxford Nanopore

|SV Caller| TP | FN | FP | Recall | F1 score |
|---------|----|----|----|--------|----------|
| Severus | 53 | 15 | 35 | 0.78   | 0.68 | 
| nanomonsv| 44 | 24 | 46 | 0.65 | 0.56 |
| sniffles2| 35 | 33 | 263 | 0.51 | 0.19 | 
| SAVANA | 51 | 17 | 18 | 0.75 | 0.74 |

#### Pacbio

|SV Caller| TP | FP | FN | Recall | F1 score |
|---------|----|----|----|--------|----------|
| Severus | 59 | 9 | 33 | 0.87 | 0.74 | 
| nanomonsv| 47 | 21 | 27 | 0.69 | 0.66 |
| sniffles2| 40 | 28 | 181 | 0.59 | 0.28 |
| SAVANA | 48 | 20 | 24 | 0.70 | 0.69 |

---
### Contact
For advising, bug reporting and requiring help, please contact aysegokce.keskus@nih.gov





