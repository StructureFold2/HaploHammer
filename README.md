#HaploHammer

Tools for generating and comparing whole haplotypes from re-sequenced data. This readme is incomplete as of right now
and the enclosed manual.pdf should serve as a reference in the interm.

**Dependencies**
+ [Python 2.7.X](https://www.python.org/)
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

