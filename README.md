# Severus

A tool to build breakpoint graphs for one or multiple long-read cancer samples. This is work in progress and subject to frequent updates. 

<img width="1373" alt="bp_example" src="https://user-images.githubusercontent.com/2475380/198373853-a606ef10-63f9-4bde-b72f-31a249a38948.png">


## Installation

Requirements:
* Python3
* networkx
* numpy
* pydot
* graphviz

The easiest way to install dependencies is through conda.

## Running

```
./bga.py --target-bam tumor.bam --control-bam normal.bam --out-dir bga_out -t 16 --min-support 5 --max-read-error 0.01
dot -Tsvg -O bga_out/breakpoint_grpah.dot
```
The primary output is the breakpoint graph, like on the example above. Black edges correspond to the fragments of the reference genome,
and dashed colored edges correspond to non-reference connections from reads. Each breakpoint is defined by its coordinate
and sign. Plus sign corresponds to connection to the left of breakpoint (reference coordinates), and minus - to the right.

Providing control bams is optional, but recommended. Ideally, it is a matching normal sample. But this could be any unrelated,
non-cancerous sample as well. This helps to filter out many regions with uncertain alignemnt due to mapping difficulties or reference bias. 
If control bams are provided, the output will prioritize connected components with adjacencies that are unique for "target" bams.

## Important parameters

* `--target-bam` path to one or multiple target bam files (must be indexed)
  
* `--control-bam` path to one or multiple control bam files (must be indexed)
  
* `--min-support` minimum number of reads supporting a breakpoint [5]
  
* `--max-read-error` maximum base alignment error for read [0.1]. Setting to `0.01` is resommended for HiFi.
  
* `--min-mapq` minimum mapping quality for aligned segment [10]

* `--reference-adjacencies` draw reference adjacencies (as dashed black edges)

* `--max-genomic-len` maximum length of genomic segment to form connected components [100000]
