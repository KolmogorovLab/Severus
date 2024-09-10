

# Severus

<p>
<img src="docs/severus_logo_right.png" alt="Severus logo" align="left" style="width:100px;"/>
Severus is a somatic structural variation (SV) caller for long reads (both PacBio and ONT). It is designed for matching tumor/normal analysis,
supports multiple tumor samples, and produces accurate and complete somatic and germline calls. Severus takes
advantage of long-read phasing and uses the breakpoint graph framework to model complex chromosomal rearrangements.

If you use Severus, please cite https://www.medrxiv.org/content/10.1101/2024.03.22.24304756v1
</p>

<br/>

## Contents

* [Installation](#installation)
* [Quick Usage](#quick-usage)
* [Input and Parameters](#inputs-and-parameters)
* [Benchmarking Severus and other SV callers](#benchmarking-severus-and-other-sv-callers)
* [Output Files](#output-files)
* [Overview of the Severus algorithm](#overview-of-the-severus-algorithm)
* [Preparing phased and haplotagged alignments](#preparing-phased-and-haplotagged-alignments)
* [Generating VNTR annotation](#generating-vntr-annotation)
* [Generating PoN file](#generating-pon-file)
* [Additional Fields in the vcf](#additional-fields-in-the-vcf)
* [Breakpoint Graphs](#breakpoint-graphs)

## Installation

The easiest way to install is through conda:

```
conda create -n severus_env severus
conda activate severus_env
severus --help
```

Or alternatively, you can clone the repository and run without installation,
but you'll still need to install the dependencies via conda:

```
git clone https://github.com/KolmogorovLab/Severus
cd Severus
conda env create --name severus_env --file environment.yml
conda activate severus_env
./severus.py
```

## Quick Usage

Single sample SV calling (Tumor-only)

```
severus --target-bam phased_tumor.bam --out-dir severus_out -t 16 --phasing-vcf phased.vcf \
    --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed --PON ./pon/PoN_1000G_hg38.tsv.gz
```
We are currently providing two PoN files for GRCh38 and chm13 which are generated using [1000G data](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.0/). 
If you use any of these please cite [this paper](https://www.biorxiv.org/content/10.1101/2024.04.18.590093v1) . To generate a PoN file please see [below](#generating-pon-file)

Single sample somatic SV calling (Tumor/Normal pair)

```
severus --target-bam phased_tumor.bam --control-bam phased_normal.bam --out-dir severus_out \
    -t 16 --phasing-vcf phased.vcf --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
```

Multi-sample somatic SV calling

```
severus --target-bam phased_tumor1.bam phased_tumor2.bam --control-bam phased_normal.bam \
    --out-dir severus_out -t 16 --phasing-vcf phased.vcf \
    --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
```

Haplotagged (phased) alignment input is highly recommended but not required. See [below](#preparing-phased-and-haplotagged-alignments)
for the detailed instructions on how to prepare haplotagged alignments. 
If using haplotagged bam, the matching phased VCF file should be provided as `--phasing-vcf` option.

`--vntr-bed` argument is optional but highly recommended. VNTR annotations for common references are available
in the `vntrs` folder. See [below](#generating-vntr-annotation) how to generate annotations for a custom reference. 

Without `--control-bam`, only germline variants will be called.

After running, vcf files with somatic and germline(with somatic) calls are available at in the output folder,
along with complexSV clusters and additional information about SVs. See [below][#breakpoint-graphs] for complexSV cluster outputs.


## Inputs and Parameters

### Required

```
--target-bam    path to one or multiple target bam files (e.g. tumor, must be indexed) 
--out-dir       path to output directory
```

### Highly recommended

```
--control-bam     path to the control bam file (e.g. normal, must be indexed)
--vntr-bed        path to bed file for tandem repeat regions (must be ordered)
--phasing-vcf     path to vcf file used for phasing (if using haplotype specific SV calling)
```

### Optional parameters

```
--threads               number of threads [8]
--min-support           minimum number of reads supporting a breakpoint [3]
--vaf-thr               variant allele frequency threshold for SVs
--TIN-ratio             tumor in normal ratio [0.01]
--min-mapq              minimum mapping quality for aligned segment [10]
--max-genomic-len       maximum length of genomic segment to form connected components [2Mb]
--min-sv-size           minimum SV size to be reported [50]
--min-reference-flank   minimum distance between a breakpoint and reference ends [10000]
--write-alignments      write read alignments to file
--bp-cluster-size       maximum distance in bp cluster [50]
--write-collapsed-dup   outputs a bed file with identified collapsed duplication regions
--no-ins-seq            do not output insertion sequences to the vcf file
--resolve-overlaps      resolve overlaps between split alignments.
--between-junction-ins  report unmapped insertions around breakpoints
--single-bp             to add hanging breakpoints 
--output-read-ids       outputs read IDs for support reads
--use-supplementary-tag to use HP tag in supplementary alignments. Need to be added if HiPhase or LongPhase is used for haplotagging.
```

## Benchmarking Severus and other SV callers

For the details of benchmarking and complete results, please check https://www.medrxiv.org/content/10.1101/2024.03.22.24304756v1

### Germline benchmarking results using HG002

First, we verified performance of Severus on a germline SV benchmark. 
We compared Severus, [sniffles2](https://github.com/fritzsedlazeck/Sniffles) and [cuteSV](https://github.com/tjiangHIT/cuteSV) using 
[HG002 GIAB SV benchmark set](https://www.nature.com/articles/s41587-020-0538-8). Comparison was perfromed using [Minda](https://github.com/KolmogorovLab/minda).
The benchamrking was done against the grch38 reference, for consistency with the benchmarks below (confident regions were lifted over).
Severus and Sniffles2 performed similarly, with CuteSV running a bit behind.

|SV Caller | TP   | FN  | FP  | Precision | Recall | F1 score |
|----------|------|-----|-----|-----------|--------|----------|
|Severus   | 9568 |	302 | 216 |	  **0.98** |**0.97**| **0.97** |
|Sniffles2 | 9603 |	277 | 305 |	   0.97    |   0.97 |	  0.97 |
|cuteSV    | 9530 |	465	| 530 |	   0.95    |   0.95 |     0.95 |

### Somatic benchmarking results COLO829

We compared the performance of existing somatic SV callers [nanomonSV](https://github.com/friend1ws/nanomonsv), [SAVANA](https://github.com/cortes-ciriano-lab/savana) 
and [sniffles2](https://github.com/fritzsedlazeck/Sniffles) in mosaic mode using COLO829 cell line data against multi-platform 
[Valle-Inclan et al. truthset](https://www.sciencedirect.com/science/article/pii/S2666979X22000726).  
We compared somatic SVs using [Minda](https://github.com/KolmogorovLab/minda) with 0.1 VAF threshold. Severus had the highest recall and precision on the HiFi dataset,
and highest recall on the ONT dataset, with nanomonsv having highest precision.

|Technology|Caller|TP|FN|FP|Precision|Recall|F1|
|----------|------|--|--|--|--------|-------|--|
|PacBio|Severus|59|10|20|0.75|**0.86**|**0.80**|
|PacBio|SAVANA|56|14|78|0.42|0.80|0.55|
|PacBio|nanomonsv|51|17|15|**0.77**|0.75|0.76|
|PacBio|Sniffles2|49|20|249|0.16|0.71|0.27|
|-----|
|ONT|Severus|54|14|24|0.69|**0.79**|0.74|
|ONT|SAVANA|51|18|14|0.78|0.74|**0.76**|
|ONT|nanomonsv|45|25|11|**0.80**|0.64|0.71|
|ONT|Sniffles2|36|32|262|0.12|0.53|0.20|
|-----|
|Illumina|SvABA|48|20|10|0.83|0.71|0.76|
|Illumina|GRIPSS|55|13|0|**1.00**|**0.81**|**0.89**|
|Illumina|manta|46|22|7|0.87|0.68|0.76|

### Somatic Benchmarking results: Tumor/Normal Cell line pairs

We compared the performance of the somatic SV callers using 5 tumor/normal cell line pairs. Since no ground truth
SV calls are available, we created an ensemble set of SVs supported by 2+ technology and 4+ callers (out off 11) for each dataset.
This assumes that singleton calls are false-positives, and calls supported by multiple tools are more reliable.
Severus consistently had the highest recall and precision against the ensemble SV sets.

<p align="center">
  <img src="docs/bench.png" alt="benchmarking" style="width:90%"/>
</p>


## Output Files

#### VCF file

For each target sample, Severus outputs a VCF file with somatic SV calls.
If the input alignment is haplotagged, haplotype will be reported as HP in INFO. In addition,
Severus outputs a set of all SVs (somatic + germline) for each input sample.
VCF contains additional information about SVs, such as the clustering of complex
variants. Please see the detailed description [below](#additional-fields-in-the-vcf).

#### html plots

Severus outputs a breakpoint graph as interactive plotly graph that describes the derived structure
of tumor haplotypes. See the detailed description [below](#breakpoint-graphs).

#### breakpoint_double.csv

Detailed info about the detected breakpoints for all samples in text format, intended for an advanced user.

#### breakpoint_clusters_list.tsv

A summary file for the complexSVs clusters. 

#### breakpoint_clusters.tsv
Detailed information of the junctions in involved in complex SVs.

## Overview of the Severus algorithm

<p align="center">
  <img src="docs/g1.png" alt="Severus workflow" style="width:90%"/>
</p>

Somatic SVs in cancer are typically more complex compared to germline SVs. For example, breakage-fusion-bridge (BFB) amplifications 
are characterized by multiple foldback inversions and oscillating copy numbers. Current long-read SV algorithms were designed 
for relatively simple germline variation and do not automatically detect complex multi-break rearrangements in cancer.
Severus is designed to provide a comprehensive view of somatic and germline SVs, and implements several algorithmic novelties.
The algorithm has modes for `single`, `two (e.g., tumor/normal)`, and `multiple samples` (e.g., multi-site or time series). 

**Haplotype-specific SVs**. A key advantage of long-read sequencing is the ability to phase variants into longer haplotypes, without the need of any additional data. 
Both germline and somatic variants can be phased, however somatic variants must still be attributed to a single germline haplotype of origin. 
If the normal sample is available, we perform SNP calling and phasing of the normal sample alignment, 
and then haplotag both normal and tumor alignments. In the absence of a normal sample, phasing of tumor alignment data can be performed instead. 

**Mismapped repeats filtering**. Mismapped reads from wrong repeat copies may artificially inflate (or deplete) 
the coverage of the region, and these reads may contain SV signatures, resulting in false-positive calls. 
It is one of the major sources of false positive somatic SV calls, as mismapped read alignment is often unstable and may vary between tumor and normal samples. 
Since the copies of evolutionary divergent repeats are not identical, a common signature of mismapped repeats is an increased sequence divergence. 
Severus filters out mismaped repetitive reads with an increased subtitution rates.

**VNTR synchronization**. Variation inside variable number tandem repeats (VNTRs) is another major source of error in SV analysis. This is because the alignment inside the VNTR regions is often ambiguous. 
Since reads contain errors at random positions, this may affect the scores of the near-optimal alignments in different ways. For each read, Severus transforms the SV signatures 
inside annotated VNTR regions into a uniform coordinate system. 

**Breakpoint graph**. Each SV signature (extracted from an individual read) connects two non-adjacent breakpoints of the reference genome. 
We will denote these regions as `start_pos` and `end_pos`. Each breakpoint has a direction, either `left` or `right`. Thus, an SV signature is defined by a pair of breakpoints: 
`(start_pos, start_dir):(end_pos, end_dir)`. For example, a deletion at position `P` of size `L` corresponds to: `(P, left):(P + L, right)`. 

Given a set of SV signatures, they are clustered based on their start and end breakpoints. For each accepted cluster that passes several quality thresholds, 
we create two nodes on the breakpoint graph (for start and end breakpoints), and connect the nodes with an `adjacency` edge. We then add `genomic` edges
that corrspond to target genome segments using long reads that span multiple consecutive SVs. Genomic edges are also formed for distant breakpoints that
are not connected by long reads, but are in the same (germline) phase and consistent direction. 

As a result, the connected components of the constructed breakpoint graph naturally represent clusters of SVs, for example chains of deletions or translocations. To capture
more complex and potentially overlapping SVs - such as chromoplexy or breakage-fusion-bridge - we perform additional clustering specific to complex SV types.


## Preparing phased and haplotagged alignments

We recommend running Severus with phased and haplotagged tumor and normal alignments. Below is an example
workflow that can be used to produce them, assuming that reads are already aligned using minimap2/pbmm2.
If phasing is difficult or not possible (e.g. haploid genome), unphased alignment can be used as input.

If normal sample is available:
* 1. SNP calling and phasing normal bam. See for [DeepVariant](https://github.com/google/deepvariant) and [margin](https://github.com/UCSC-nanopore-cgl/margin). 

  + Using DeepVariant + Margin

```
# To install DeepVariant and margin using singularity

 singularity pull docker://google/deepvariant:latest
 docker pull kishwars/pepper_deepvariant:latest
 
# To generate phased VCF for ONT R10 data using DeepVariant + Margin

singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ deepvariant_latest.sif run_deepvariant \
--model_type ONT_R104 
--ref ref.fa \
--reads normal.bam \
--output_vcf normal_vcf.vcf \
--num_shards 56

docker run kishwars/pepper_deepvariant:latest \
margin phase normal.bam ref.fa normal_vcf.vcf allParams.haplotag.ont-r104q20.json -t 56 -o output_dir

# To generate phased VCF for PacBio HiFi data using DeepVariant + Margin

singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ deepvariant_latest.sif run_deepvariant \
--model_type PACBIO \
--ref ref.fa \
--reads normal.bam \
--output_vcf normal_vcf.vcf \
--num_shards 56

docker run kishwars/pepper_deepvariant:latest \
margin phase normal.bam ref.fa normal_vcf.vcf allParams.haplotag.pb-hifi.json -t 56 -o output_dir

```

  + Using Clair3

For the complete list of the available models, [see](https://github.com/HKU-BAL/Clair3/tree/main#pre-trained-models). with `--enable_phasing` and `--longphase_for_phasing` [Clair3](https://github.com/HKU-BAL/Clair3) generates unphased and phased vcfs using [LongPhase](https://academic.oup.com/bioinformatics/article/38/7/1816/6519151). 

```
# To install Clair3

 singularity pull docker://google/deepvariant:latest
 docker pull kishwars/pepper_deepvariant:latest
 
# To generate phased VCF using Clair3

INPUT_DIR="[YOUR_INPUT_FOLDER]"       
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      
THREADS="[MAXIMUM_THREADS]"            
MODEL_NAME="[YOUR_MODEL_NAME]"

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    
  --ref_fn=${INPUT_DIR}/ref.fa \       
  --threads=${THREADS} \               
  --platform="ont" \                   
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR} \
  --enable_phasing \
  --longphase_for_phasing

```

  + Using HiPhase for phasing (HiFi only) [see](https://github.com/PacificBiosciences/HiPhase)
    
    
```
# To install HiPhase

    conda create -n hiphase -c bioconda hiphase
    conda install hiphase

# To run HiPhase with DV or Clair3 output
    
    hiphase --bam normal.bam\
        --vcf hifi_vcf.vcf.gz\
        --output-vcf hifi_hiphase.vcf.gz\
        --reference ref.fa --threads 16 --ignore-read-groups
    
```

* 2. Haplotagging normal and tumor bams using phased vcf from step1.

  + Using [whatshap](https://whatshap.readthedocs.io/en/latest/index.html)
    
```
# To install WhatsHap using Conda

 conda create -n whatshap-env whatshap 
 conda activate whatshap-env
 
# To haplotag Normal and Tumor

whatshap haplotag --reference ref.fa phased_vcf.vcf normal.bam -o normal.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=4

whatshap haplotag --reference ref.fa phased_vcf.vcf tumor.bam -o tumor.haplotagged.bam --ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads=4

```

  + Using [longphase](https://github.com/twolinin/longphase)
    
```
# To install longphase

 wget https://github.com/twolinin/longphase/releases/download/v1.6/longphase_linux-x64.tar.xz
 tar -xJf longphase_linux-x64.tar.xz
 
# To haplotag Normal and Tumor

longphase haplotag -r ref.fa -s phased_vcf.vcf -b normal.bam -t 8 -o normal_haplotagged --tagSupplementary --qualityThreshold 10

longphase haplotag -r ref.fa -s phased_vcf.vcf -b tumor.bam -t 8 -o tumor_haplotagged --tagSupplementary --qualityThreshold 10


```

## Generating VNTR annotation

Alignments to highly repetitive regions are often ambigous and leads false positives. Severus clusters SVs inside a single VNTR region to uniform the SV representation for each read.

In the example above, alignment led two deletion inside a vntr. Severus combines two deletions and transform them into a single longer deletion.

VNTR annotation file for the most commonly used genomes are generated using [findTandemRepeats](https://github.com/PacificBiosciences/pbsv/tree/master/annotations) and provided in the [vntrs](vntrs) folder .

To generate VNTR annotation file for a new genome using [findTandemRepeats](https://github.com/PacificBiosciences/pbsv/tree/master/annotations):

```
findTandemRepeats --merge <REF>.fa <REF>.trf.bed

```
## Generating PoN file

To improve the computation time, Severus requires a summary text file generated using the vcf file. To generate a PoN data from a multisample vcf file please see [the script](scripts/pon_from_vcf.py) 

## Additional fields in the vcf

Severus outputs several unique fileds in INFO column.

* `DETAILED_TYPE` : Severus can identify SV types other than standart SV types; INS, DEL, INV, BND. Additional types: 
  + `reciprocal_inv`: a reciprocal inversion with matching (+,+) and (-,-) adjecencies.
  + `dup_inv_segment`: an inversion with a duplicated segment
  + `reciprocal_inv_del`: an inversion with a deletion
  + `foldback`: foldback inversion
  + `BFB_foldback`: Two matching foldback inversions
  + `Reciprocal_tra`: Reciprocal translocation
  + `inv_tra`: An inversion with a translocation
  + `Templated_ins`: templated insertion cycle
  + `Intra_chr_ins`: one chromosome piece inserted into the same chromosome
  

* `INSLEN`: if `--between-junction-ins` added to run severus, it outputs the length of unmapped sequence between two ends.

* `INSSEQ`: if `--between-junction-ins` added to run severus, it outputs the unmapped sequence between two ends.

* `INSIDE_VNTR`: True if the both breakpoints are inside the same VNTR.

* `ALINGED_POS`: Position if the insertion sequence if there is a supporting read.

* `PHASESET_ID` : if phasing vcf is provided, severus outputs phaseset ID from phasing vcf for phased SVs.

* `CLUSTERID` : subgraph ID from breakpoint graph.  
   

## Breakpoint Graphs

<p align="center">
  <img src="docs/graph1.png" alt="Severus overview"/ style="width:90%"/>
</p>

The primary output of Severus is the breakpoint graph generated separately for somatic and germline SVs.
After identifying breakpoint pairs (or a single breakpoint and insertion sequence in unmapped insertions), it filters simple SVs, i.e., INDELS, reciprocal inversions, local duplications.
For the remaining SVs, it adds the adjacency edges (dashed edges) between breakpoint pairs.
The second step is to add genomic segments as it connects breakpoints to the nearest breakpoints based on genomic coordinates and direction: i), (ii) both breakpoints have similar coverage,
and (ii) the length is less than MAX_GENOMIC_DIST (2Mb by default).

Each subgraph with +2 SVs are selected as complex SV clusters and assigned a cluster ID which are also output in vcfs in CLUSTER_ID filed in the INFO column.

The breakpoints in each cluster and the number of support reads are also summarized in breakpoint_cluster.csv. 

Severus outputs breakpoint graphs as interactive plotly graphs in plots folder along with summary files; breakpoint_clusters.csv and breakpoint_clusters_list.tsv.


#### Plotly graphs

Plotly graphs are generated as html files in plots folder. 

Genomic segments are seperated by chromosome and ordered with their respective position in chromosome and represented by green segments. Each segment is labelled with the coverage of the segments (total coverage), length, and haplotype as Bp1 Haplotype|Bp2 Haplotype.

Adjacencies are represented by solid edges and colored according to the direction of the junction i.e., HT(Head-to-Tail or '+-' or DEL-like), TH(Tail-to-Head or '-+' or DUP-like), HH/TT(Head-to-Head or Tail-to-Tail or INV-like) and interchromosomal (translocation like).
Each adjacency is labelled with the supporting sample, number of supporting reads and haplotype.

In the example above, there is a deletion in one of the haplotypes (chr19: 9040865- 9041235),
and there is an insertion of a 73.3kb region in chr10 to chr19 in the other haplotype.

License
-------

Severus is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Severus is developed in Kolmogorov Lab at the National Cancer Institute.

Key contributors:

* Ayse Keskus
* Asher Bryant
* Mikhail Kolmogorov

---
### Contact
For advising, bug reporting and requiring help, please submit an [issue](https://github.com/KolmogorovLab/Severus/issues).
You can also contact the developer: aysegokce.keskus@nih.gov
