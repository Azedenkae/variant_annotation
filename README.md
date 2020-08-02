# variant_annotation v0.1
Prototype tool for annotating allele variants, developed as part of the Tempus challenge.

### Default annotations:
1. Type of variation, and their effect from ExAC API. If there are multiple effects, the most deleterious is selected.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from ExAC API.

# Dependencies:
* Python libraries: argparse, os, requests, json

# Installation:
* variant_annotation does not require installation.
* argparse, os, and json are default components of the Python 3 library (https://docs.python.org/3/library/).
* requests can be installed following these instructions: https://requests.readthedocs.io/en/master/user/install/.

# Running the tool:
    variant_annotation.py -i [input]
Input must be in standard vcf format.

### Optional parameters:
* -o [output], specify output file name in tsv format
* -def, add definition subheaders
* -dump, use this argument for straight data dump

Only use -dump if straight data dump is desired (publish all INFO and FORMAT annotations existing in the vcf file).

### Recommended command:
    variant_annotation.py -i [input] -def
