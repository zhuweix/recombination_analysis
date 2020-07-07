# Analysis of the non-allelic homologous recombination events using a k-mer counting approach

This package aims to predict non-allelic homologous recombination (NAHR) events in ORFs and other genetic elements between genomes.

## Prerequisite

matplotlib

numpy

pandas

bitarray

## How to identify NAHR events

For each query sequence in query strain, we identified sequence variations ratlive to the orthologue sequence in the reference strain. We speculate that NAHR between the orthologue and the paralogue(s) in the referenence strain lead to sequence variations that shared between the query sequence and the paralogous sequence(s) rather than the orthologous sequence. On the contrary, sequence mutations between the orthologue and the query gene do not lead to shared variations between the query and the orthologue(s).

We identified the shared sequence variations using k-mers. For each query sequence and all the reference sequences homologous to the query sequence, we performed k-mer counting. K-mers which are shared with only between the query and the paralog, but not with the corresponding ortholog in reference, indicate shared sequence variations, therefore, indicating NAHR events.


## What is the input of the analysis

We need the sequences of query genes, and the sequences of reference genes, and the table of the homologous pairs between query genes and reference genes. The orthologue of the query gene in the reference genes must have the same name with the query gene. Our package can omit the requirement for orthologue for NAHR identification. However, the orthologue is required to draw the dotplot illustration of NAHR events.

## Gerenal Model for NAHR identification

<details>
  <summary>Click for details</summary>
<p>
Each query sequence is separated into minimal recombination regions (including non-recombined regions as the special case for minimal recombination regions) and the homologous flanking regions.
 </p>
<p> 
  Minimal recombination regions
  
  Minimal recombination regions are regions containing the k-mers shared only with the query sequence and the recombined paralogues in NAHR (recombined k-mers).
  
  1. In order to minimize the NAHR events and the associated recombined paralogues, the adjacent minimal recombination regions will be combined if they share a paralogue, i.e. all the adjacent minimal recombination regions are associated with different paralogues.
  
  2. The length of the minimal recombination regions satisfying 1. is minimized by trimming the minimal recombination regions to the first and last recombined k-mers within the regions. 
  
  Non-recombined regions
  
  Non-recombined regions is a operationaly a special form of minimal recombination regions which are identified using k-mers shared only between the query sequence and the orthologue sequence. The boundaries of the non-recombined regions are also set up with the first and last k-mer shared only between the query and the orthoogue. In practice, our package do not specially label the non-recombined regions. They are the minimal recombination regions in which the orthologue gene leads to "NAHR".
  
  Homologous flanking regions
  
  Homologous flanking regions are regions between the minimal recombination regions and the non-recombined regions. By definition, all the shared k-mers between the query sequence and the reference sequences are shared between the query sequence in the query strain, the orthologue sequence and the recombined paralogues. Therefore, this region is shared between the orthologue and the recombined paralogue for NAHR.
  
  The homologous flanking regions can also occur between two minimal recombination regions, which indicate an NAHR event between the two (sets of) parlogues associated with the two regions. Similarly, all the shared k-mers within these regions are shared between the query, and the two (sets of) paralogues.
  
</p>
  
  
  
</details>

## Bi-directional k-mer scanning to identify recombination regions

<details>
  <summary>Click for details</summary>
  The minimal recombination regions and the homologous flanking regions are identified using bi-directional k-mer scanning.
  <ol>
  <li> K-mer counting in reference sequences </li>
  
  For a query sequence with N reference sequences, we count the k-mers in 1XN binrary arrays. Each array records the existance of one kmer in reference sequence 1, 2, 3, ..., N. We call the {k-mer: binary array of the associated reference sequences} the reference dictionary.
  
  <li> Forward scanning </li>
  
  We scan the k-mer from the 5' end of the query sequence. 
  <ol>
  <li> We find the first k-mer shared between the query and the references, which is the start location of the first minimal recombination region. </li>
 <li> We iniate the associated refereence sequences of the first minimal recombination region with the reference sequences accounting for the first k-mer to the reference dictionary </li>
  <li> Proceed to the next k-mer. Skip it if it is not shared with the reference. </li>
  <li> For the next shared k-mer, perform intersection of the associated refereence sequences with the next shared k-mer and those with the recombination region.</li>
  <li>
  If the intersection is not empty, update the associated refereence sequences of the first minimal recombination with the intersection set. Repeat the iii-v/vi until the end of the query. </li>
  <li>
 If the intersection is empty, this indicates that non of the reference sequences contains all the shared k-mers. This location is the start of the next minimal recombination region, and we initiated the the associated refereence sequences of the next minimal recombination region with the reference sequences accounting for the first k-mer to the reference dictionary. Repeat iii-v/vi until the end of the query.</li>
  </ol>
  After the forward scanning, we obtained the start locations of all the minimal recombination regions (as well as the non-recombined regions as the special case).
  
  
  <li> Backward scanning </li>
  
  We performed similar scanning from the 3' end, and obtained the end locations of  all the minimal recombination regions (as well as the non-recombined regions as the special case).
  
  <li> Combination of forward and backward scanning to identify minimal recombination regions, non-recombined regions and the homologous flanking regions. </li>
  
  The forward scaning identifies recombination regions separated with the start locations of the minimal recombination regions, therefore, they are constituted with two parts: The left part is the minimal recombination region, and the right part is the homologous flanking region. Similarly, backward scanning identifies recombination regions separated with the end locations of the minimal recombination regions, therefore, they are constituted with two parts: The left part is the homologous flanking region and the right part is the minimal recombination region. 
  
  When we intersect the two sets of recombination regions, the forward and backward recombination regions identifying the same NAHR event will associated with the same (or at least one shared) reference genes, and the overlap of the two regions define the minimal recombination region; In contrast, the forward and backward recombination regions identifying the recombined region and the adjacent non-recombined regions are associated with different reference genes, thus, the overlap of the two regions define the homologous flanking region.
  </ol>

  
  
</details>
  
## Integrate Method for two-step NAHR prediction between strains

We used a two-step idenfication procedure in the integrate method. We first performed the NAHR identification using the references sequences in the provided homology pair file using the k-mer_size_first; then, we performed a second round of NAHR identification with the recombined paralogs identified in the first step using the second k-mer_size_second to refine the NAHR identification. Note that k-mer_size_first > k-mer_size_second, otherwise, it is no different from one-step identifcation using k-mer_size_second.

The integrate method:

```
orfeome_comparison(
  query_fn: str, 
  ref_fn: str, 
  homo_pair_fn: str,  
  window_first: int,
  window_second: int, 
  assign_fig_dir='', 
  compare_fig_dir='', 
  result_prefix='out',
  fig_prefix='Out',
  comp_top=3, skip_novel=True)
```                 

Arguments:

- query_fn (str): fasta file with the query sequences.
- ref_fn (str): fasta file with reference sequences.
- homo_pair_fn (str):  name of the tsv file for homologous gene pairs.
- window_fist (int): k-mer size of the first round NAHR identification.
- window_second (int): k-mer size of the second round NHAR identification.
- assign_fig_dir (str): dir name of the NAHR assignment figures. If empty, no figures will be generated.
- compare_fig_dir (str): dir name of the dotplot comparison for NAHR. If empty, no figures will be generated.
- resulf_prefix (str): prefix of the NAHR identification files.
- fig_prefix (str): prefix of the figures.
- com_top (int): Top N recombined paralogues to be illustrated for dotplot comparison.
- skip_novel=True: Skip ORFs with no orthologue in reference. Dotplot comparison REQUIRES orthologue, the comparison figure will SKIP genes without orthologous even if skip_novel=False.

