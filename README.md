# variant_annotation v0.1
Prototype tool for annotating allele variants, developed as part of the Tempus challenge.

The annotation required are:
1. Type of variation, and their effect from ExAC API. If there are multiple effects, the most deleterious is selected.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from ExAC API.

