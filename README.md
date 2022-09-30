# Breakpoint graph assembler

A proof-of-consept script to identify breakpoints from long read alignment and build breakpoint graphs.

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
./bga.py -b grcm39/grcm39_m4_c11.bam -o out_breakpoints -t 16 --min-support 3 --max-read-error 0.01 --min-reference-flank 500 --coverage
dot -Tsvg -O out_breakpoints/breakpoint_grpah.dot
```

This will generate a number of files with breakpoint information, and visualization in `breakpoint_graph.dot.svg`.

## Parameters

* `-b` and `-o`: Paths to input bam (must be sorted and indexed) and output difrectory respectively.

* `-t` controlds the number of threads

* `--min-support`: the minimum number of reads that support breakpoint. Should be adjusted depending on the dataset coverage 

* `--max-read-error`: maximum base-level read error. The recommended is `0.1` for ONT and `0.01` for PacBio HiFi

* `--min-reference-flank`: minimum distance from the end of reference sequence for a breakpoint. Could be set to 0 for specific datasets, for example containing HPV.

* `--coverage`: enable computing genomic segments coverage
