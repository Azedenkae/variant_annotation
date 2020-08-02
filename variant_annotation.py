##################################################################
### Prototype variant annotation tool designed to extract data ###
### for allelic variants as specified by the Tempus challenge. ###
##################################################################

import argparse # For parsing arguments.
import os.path # For parsing arguments.
import requests # For API.
import json # For API.

# Define ext_check as a variable type. The purpose is to check that the input has a .vcf extension. Does not allow
# tool to progress unless input is valid.
def ext_check(input_file):
    base, ext = os.path.splitext(input_file)
    if '.vcf' not in ext.lower():
        raise argparse.ArgumentTypeError('File must have a vcf extension.')
    return input_file

parser = argparse.ArgumentParser()

# Allow user to define (compulsory) input and (optional) output files. Also allows the user to add ID definitions as
# sub-headers, and what type of output they want.
parser.add_argument("-i", required=True, type=ext_check, help="Specify variant call input file in vcf format.")
parser.add_argument("-o", required=False, default="", type=str, help="Specify output file name in tsv format.")
parser.add_argument("-dump", required=False, action="store_true", help="Use this argument for straight data dump.")
parser.add_argument("-def", required=False, action="store_true", help="Add definition sub-headers.")

args = vars(parser.parse_args())
vcf_infile = args["i"]
tsv_outfile = args["o"]
data_dump = args["dump"]
id_def_include = args["def"]

# If -o not specified, outfile use infile name.
if tsv_outfile == "":
    tsv_outfile = vcf_infile.split(".")[0]

# Adds the .tsv extension, especially useful if the user forgot to specify. Unlike the input file, there is no check for
# if the user decided to add any other extension, in case that is intentional... for whatever reason.
if ".tsv" not in tsv_outfile[-4:]:
    tsv_outfile += ".tsv"

# List of consequences, as taken from the Gemini database:
# https://gemini.readthedocs.io/en/latest/content/database_schema.html, which lists the same consequences as returned
# after querying the ExAC browser. Consequence IDs have been manually swapped to VEP IDs where necessary to be in line
# with the IDs returned from the ExAC querying. Consequences are listed in order of highest to lowest severity. "unknown"
# is added at the end because not all variants have consequences listed in the ExAC browser.
consequences = ["exon_loss_variant", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost",
                "stop_gained", "stop_lost", "initiator_codon_variant", "rare_amino_acid_variant", "chromosomal_deletion",
                "missense_variant", "inframe_insertion", "inframe_deletion", "coding_sequence_variant",
                "disruptive_inframe_deletion", "disruptive_inframe_insertion", "5_prime_UTR_truncation + exon_loss_variant",
                "3_prime_UTR_truncation + exon_loss_variant", "splice_region_variant", "mature_miRNA_variant",
                "regulatory_region_variant", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
                "TFBS_ablation", "TFBS_amplification", "stop_retained_variant", "synonymous_variant", "5_prime_UTR_variant",
                "3_prime_UTR_variant", "intron_variant", "coding_sequence_variant", "upstream_gene_variant", "downstream_gene_variant",
                "intergenic_variant", "intragenic_variant", "gene_variant", "transcript_variant", "exon_variant",
                "5_prime_UTR_premature_start_codon_gain_variant", "start_retained_variant", "conserved_intron_variant",
                "nc_transcript_variant", "NMD_transcript_variant", "incomplete_terminal_codon_variant",
                "non_coding_transcript_exon_variant", "transcript_ablation", "transcript_amplification",
                "feature_elongation", "feature_truncation","unknown"]

# Prepping for later data processing.
vcf_dict = {}
info_abbreviations = []
format_abbreviations = []
info_definition = {}
format_definition = {}
var_counter = 0 # Count number of variants.
for each_line in open(vcf_infile):
    # Identify info and format IDs, and put them into lists to be used later on. We could have also specified the list
    # to save space/time, but having the tool determine the IDs available is probably better to somewhat future-proof
    # the tool in case we want to add/remove/change/update the INFO and FORMATs available. The -dump option should work
    # regardless of what changes.
    if each_line.startswith('#'):
        if "INFO=<ID=" in each_line:
            info_abbreviations.append(each_line.split(",")[0][11:]) # Grabs INFO IDs.
            info_definition[each_line.split(",")[0][11:]] = each_line.split("\"")[1] # Grabs INFO ID definitions.
        elif "FORMAT=<ID=" in each_line:
            format_abbreviations.append(each_line.split(",")[0][13:]) # Grabs FORMAT IDs.
            format_definition[each_line.split(",")[0][13:]] = each_line.split("\"")[1] # Grabs FORMAT ID definitions.
    else:
        # Creates keys from CHROM and POS ID.
        each_line_split = each_line.strip().split('\t')
        vcf_key = '%s_%s' % (each_line_split[0], each_line_split[1])
        vcf_value = each_line_split[2:]
        # Store INFO values into a dict.
        vcf_info_str_split = each_line_split[7].split(';')
        vcf_info_dict = {}
        for each_info in vcf_info_str_split:
            each_info_split = each_info.split('=')
            vcf_info_dict[each_info_split[0]] = each_info_split[1]
        # Replace the string of INFO values with the INFO dict.
        vcf_value[5] = vcf_info_dict
        vcf_dict[vcf_key] = vcf_value
        var_counter += 1

# Adds custom headers and its definition to INFO dictionary.
var_effect = "Effect"
info_definition[var_effect] = "Most significant (negative) genotypic consequence of variant."
var_ref_sup_id = "Variant/reference support."
info_definition[var_ref_sup_id] = "Percentage of reads supporting the variant versus those supporting reference reads."
exac_header = "Allele frequency (ExAC)."
info_definition[exac_header] = "Allele frequency of variant from ExAC API."

# Final lists to iterate over. Default is specified as per the Tempus challenge.
final_infos = ["TYPE",var_effect,"DP","DPB","DPRA",var_ref_sup_id,exac_header]
final_formats = []
if data_dump == True:
    final_infos = info_abbreviations
    final_formats = ["FORMAT","normal","vaf5"]

# Definitions for first four headers.
contig_def = "The name of the sequence (typically a chromosome) on which the variation is being called."
pos_def = "The 1-based position of the variation on the given sequence."
ref_def = "The reference base (or bases in the case of an indel) at the given position on the given reference sequence."
alt_allele_def = "The list of alternative alleles at this position."

variant_annotation = "Contig\tPosition\tReference\tAlternate Allele"
definition_subheadings = ""

# Add definition sub-headers for headings, if the option is included.
if id_def_include == True:
    definition_subheadings += "%s\t%s\t%s\t%s" % (contig_def, pos_def, ref_def, alt_allele_def)

# Generates headers for tsv file.
for header in final_infos:
    variant_annotation += "\t" + header
    if id_def_include == True:
        definition_subheadings += "\t" + info_definition[header]
        # Just for formatting consistency. We would probably generally want the (end) user to see consistency in the formatting.
        if not info_definition[header][-1].endswith("."):
            definition_subheadings += "."
for header in final_formats:
    variant_annotation += "\t" + header
    if header == "FORMAT":
        definition_subheadings += "\t" + "An extensible list of fields for describing the samples: "
        for definition in format_definition:
            definition_subheadings += "%s:%s, " % (definition, format_definition[definition])
        definition_subheadings = definition_subheadings[0:-2] + "\tFORMAT values for reference allele\tFORMAT values for variant"
variant_annotation += "\n" + definition_subheadings + "\n"

# Adds data for each variant.
completion_counter = 0
csq_severity = 0
for variant in vcf_dict:
    variant_split = variant.split("_")
    alternate_allele = vcf_dict[variant][2].split(",")
    alternate_counter = 0
    for alternate in alternate_allele:
        first_four_data = (variant_split[0], variant_split[1], vcf_dict[variant][1], alternate)
        variant_annotation += "%s\t%s\t%s\t%s" % first_four_data
        for abbrev in final_infos:
            variant_annotation += "\t"
            if abbrev in vcf_dict[variant][5]:
                if "," in vcf_dict[variant][5][abbrev]:
                    variant_annotation += str(vcf_dict[variant][5][abbrev].split(",")[alternate_counter])
                else:
                    variant_annotation += str(vcf_dict[variant][5][abbrev])

            # Querying ExAC for consequence information.
            elif abbrev == var_effect:
                exac_response = requests.get("http://exac.hms.harvard.edu/rest/variant/ordered_csqs/" + "%s-%s-%s-%s" % first_four_data)
                exac_data = json.loads(exac_response.content)
                # Data returned can exist as null (None in Python), an empty list, or a list of consequences. This if/else
                # statement is necessary to check for all three of the possibilities. If no data is found from the querying,
                # returns "unknown" rather than "na".
                if exac_data == None:
                    variant_annotation += "unknown"
                else:
                    exac_data = [item for item in exac_data if item] # Removes empty strings from list.
                    if not len(exac_data) == 0:
                        csq_severity = consequences.index("unknown") # Resets rank to lowest severity.
                        for csq in exac_data:
                            if consequences.index(csq) < csq_severity:
                                csq_severity = consequences.index(csq)
                        variant_annotation += str(consequences[csq_severity])
                    else:
                        variant_annotation += "unknown"

            # Concatenates proportion of paired reads supporting allele versus reference.
            elif abbrev == var_ref_sup_id:
                variant_annotation += "%s/%s" % (vcf_dict[variant][5]["PAIRED"], vcf_dict[variant][5]["PAIREDR"])

            # Querying ExAC for allele frequency information.
            elif abbrev == exac_header:
                exac_response = requests.get("http://exac.hms.harvard.edu/rest/variant/variant/" + "%s-%s-%s-%s" % first_four_data)
                exac_data = json.loads(exac_response.content)
                if "allele_freq" in exac_data:
                    variant_annotation += str(exac_data["allele_freq"])
                else:
                    variant_annotation += "na"
            else:
                # Just in case we ever get data where a specific bit of info is missing.
                variant_annotation += "na"
        if data_dump == True:
            variant_annotation += "\t%s\t%s\t%s" % (vcf_dict[variant][6], vcf_dict[variant][7], vcf_dict[variant][8])
        variant_annotation += "\n"
        alternate_counter += 1
    completion_counter += 1
    # So we know what the progress is.
    print("%s of %s variants processed." % (completion_counter, var_counter), end="\r")

# Publish data in output tsv file.
out_file = open(tsv_outfile,'w')
out_file.write(variant_annotation)
out_file.close()
