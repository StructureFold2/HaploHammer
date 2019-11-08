# HaploHammer

Tools for generating and comparing whole haplotypes from re-sequenced data. This readme is incomplete as of right now
and the enclosed manual.pdf should serve as a reference in the interm.

**Dependencies**
+ [Python 2.7.X](https://www.python.org/)
+ [BioPython](https://biopython.org/)
+ [vcftools](https://vcftools.github.io/index.html)

HaploHammer is a suite of Python scripts designed to quickly apply a large number of <.vcf>
files back to their respective genome, yielding all inferred changes to whole transcripts or particular
transcript features. These features may then be sorted into a whole haplotype matrix, cataloging the full
range of haplotypes present among all surveyed <.vcf>s, along with an accompanying ‘master’ <.fasta>
file containing an entry for every haplotype referenced by the matrix. Thus, it effectively aims to
summarize the variation among re-sequenced lines and organize that information into a useful and
halotype aware format.

**Contents**
1. vcf_manager.py
2. haplo_hammer.py
3. batch_translate.py
4. allele_matrix.py
5. matrix_expand.py
6. allele_extract.py
7. annotation_parse.py

**Input data**
This tool matches variation from vcf files back to the reference genome they were generated
against, thus you must have the same version of the reference genome that was used to generate the vcf
files, i.e. what the variants were called from. There is no known limit on the number of vcfs that can be
incorporated into a single haplotype matrix. It is possible to use an updated GFF3/GTF file as long as it
still references the same build of the genome, thus giving more or updated gene models without
perturbing the genomic coordinates. The pipeline was built around the IRRI 3000 genomes project and
the 1001 Arabidopsis genomes project, where each re-sequenced line is both diploid and contained
within it’s own vcf file; support for vcf files detailing the changes from multiple lines may be added in
the future.


