MultiSplice
===========

A Robust Method for Transcript Quantification with RNA-seq Data.
The project is actively maintained at: http://www.netlab.uky.edu/p/bioinfo/
Please cite our paper: "A Robust Method for Transcript Quantification with RNA-seq Data."
Journal of Computational Biology, 2013, doi: 10.1089/cmb.2012.0230.

===============

[Yan Huang](http://protocols.netlab.uky.edu/~yan) \(yan at netlab dot uky dot edu\)

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Compilation & Installation](#compilation)
* [Usage](#usage)


* * *

## <a name="introduction"></a> Introduction

The advent of high throughput RNA-seq technology allows
deep sampling of the transcriptome, making it possible to characterize
both the diversity and the abundance of transcript isoforms. Accurate
abundance estimation or transcript quantification of isoforms is critical
for downstream differential analysis (e.g. healthy vs. diseased cells), but
remains a challenging problem for several reasons. First, while various
types of algorithms have been developed for abundance estimation, short
reads often do not uniquely identify the transcript isoforms from which
they were sampled. As a result, the quantification problem may not be
identifiable, i.e. lacks a unique transcript solution even if the read maps
uniquely to the reference genome. In this paper, we develop a general
linear model for transcript quantification that leverages reads spanning
multiple splice junctions to ameliorate identifiability. Second, RNA-seq
reads sampled from the transcriptome exhibit unknown position-specific
and sequence-specific biases. We extend our method to simultaneously
learn bias parameters during transcript quantification to improve accuracy.
Third, transcript quantification is often provided with a candidate
set of isoforms, not all of which are likely to be significantly expressed
in a given tissue type or condition. By resolving the linear system with
LASSO our approach can infer an accurate set of dominantly expressed
transcripts while existing methods tend to assign positive expression to
every candidate isoform. Using simulated RNA-seq datasets, our method
demonstrated better quantification accuracy than existing methods. The
application of our method on real data experimentally demonstrated
that transcript quantification is effective for differential analysis of
transcriptomes.

## <a name="compilation"></a> Compilation & Installation

To compile MultiSplice, simply run
   
    make

To install, simply put the rsem directory in your environment's PATH
variable.

### Prerequisites

C++ and Matlab are required to be installed. 


## <a name="usage"></a> Usage

Please refer: http://www.netlab.uky.edu/p/bioinfo/MultiSplice/userguide
for detailed usage of the software and examples.
