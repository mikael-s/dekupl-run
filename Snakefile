#!/bin/python3

#######################################################################
# The MIT License
#
# Copyright (c) 2017, Jérôme Audoux (jerome.audoux@inserm.fr)
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files 
# (the “Software”), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be 
# included in all copies or substantial portions of the Software.
#
# The Software is provided “as is”, without warranty of any kind, 
# express or implied, including but not limited to the warranties of 
# merchantability, fitness for a particular purpose and 
# noninfringement. In no event shall the authors or copyright holders
# be liable for any claim, damages or other liability, whether in an 
# action of contract, tort or otherwise, arising from, out of or in 
# connection with the software or the use or other dealings in the 
# Software. 
#######################################################################

import os
import gzip
from snakemake.utils import R

__author__ = "Jérôme Audoux (jerome.audoux@inserm.fr)"


def getbasename(fileName):
    fileName = os.path.basename(fileName)
    *name, extension, compression = fileName.split(os.path.extsep)
    baseName = '.'.join(name)
    return(baseName)


configfile: "config.json"

# COMMON VARIABLES
SAMPLE_NAMES    = [i['name'] for i in config["samples"]]
CONDITION_COL   = "condition"
CONDITION_A     = config['diff_analysis']['condition']['A']
CONDITION_B     = config['diff_analysis']['condition']['B']
PVALUE_MAX      = config['diff_analysis']['pvalue_threshold']
LOG2FC_MIN      = config['diff_analysis']['log2fc_threshold']
MIN_REC         = config['dekupl_counter']['min_recurrence']
MIN_REC_AB      = config['dekupl_counter']['min_recurrence_abundance']
LIB_TYPE        = config['lib_type']    if 'lib_type'     in config else "rf"
R1_SUFFIX       = config['r1_suffix']   if 'r1_suffix'    in config else "_1.fastq.gz"
R2_SUFFIX       = config['r2_suffix']   if 'r2_suffix'    in config else "_2.fastq.gz"
CHUNK_SIZE      = config['chunk_size']  if 'chunk_size'   in config else 1000000
TMP_DIR         = config['tmp_dir']     if 'tmp_dir'      in config else os.getcwd()
KMER_LENGTH     = config['kmer_length'] if 'kmer_length'  in config else 31
DIFF_METHOD     = config['diff_method'] if 'diff_method'  in config else 'DESeq2'
OUTPUT_DIR      = config['output_dir']
FASTQ_DIR       = config['fastq_dir']

# DIRECTORIES
BIN_DIR         = "bin"
TMP_DIR         = temp(TMP_DIR + "/dekupl_tmp")
GENE_EXP_DIR    = OUTPUT_DIR + "/gene_expression"
COUNTS_DIR      = OUTPUT_DIR + "/kmer_counts"
KMER_DE_DIR     = OUTPUT_DIR + "/" + CONDITION_A + "_vs_" + CONDITION_B + "_kmer_counts"
METADATA_DIR    = OUTPUT_DIR + "/metadata"
REFERENCE_DIR   = OUTPUT_DIR + "/references"
LOGS            = OUTPUT_DIR + "/Logs"

# FILES
RAW_COUNTS                  = COUNTS_DIR    + "/raw-counts.tsv.gz"
NORMALIZATION_FACTORS       = COUNTS_DIR  + "/normalization_factors.tsv"
DIFF_COUNTS                 = KMER_DE_DIR   + "/diff-counts.tsv.gz"
PVALUE_ALL                  = KMER_DE_DIR   + "/raw_pvals.txt.gz"
MERGED_DIFF_COUNTS          = KMER_DE_DIR   + "/merged-diff-counts.tsv.gz"
ASSEMBLIES_FASTA            = KMER_DE_DIR   + "/merged-diff-counts.fa.gz"
ASSEMBLIES_BAM              = KMER_DE_DIR   + "/merged-diff-counts.bam"
SAMPLE_CONDITIONS           = METADATA_DIR  + "/sample_conditions.tsv"
SAMPLE_CONDITIONS_FULL      = METADATA_DIR  + "/sample_conditions_full.tsv"
DEGS                        = GENE_EXP_DIR  + "/" + CONDITION_A + "vs" + CONDITION_B + "-DEGs.tsv"
CHECKING_PLOTS              = KMER_DE_DIR   + "/checking_plots.pdf"
DIST_MATRIX                 = GENE_EXP_DIR  + "/clustering_of_samples.pdf"
NORMALIZED_COUNTS           = GENE_EXP_DIR  + "/normalized_counts.tsv"
PCA_DESIGN                  = GENE_EXP_DIR  + "/pca_design.tsv"

# binaries
REVCOMP         = BIN_DIR + "/revCompFastq.pl"
DEKUPL_COUNTER  = BIN_DIR + "/dekupl-counter"
DIFF_FILTER     = BIN_DIR + "/diffFilter.pl"
TTEST_FILTER    = BIN_DIR + "/TtestFilter"
JOIN_COUNTS     = BIN_DIR + "/joinCounts"
MERGE_COUNTS    = BIN_DIR + "/mergeCounts.pl"
MERGE_TAGS      = BIN_DIR + "/mergeTags"
COMPUTE_NF      = BIN_DIR + "/compute_norm_factors.R"
JELLYFISH       = "jellyfish"
JELLYFISH_COUNT = JELLYFISH + " count"
JELLYFISH_DUMP  = JELLYFISH + " dump"
PIGZ            = "pigz"
ZCAT            = "gunzip -c"
SORT            = "sort"

# SET MEMORY/THREAD USAGE FOR EACH RULE
MAX_MEM_JELLYFISH = 8000
MAX_MEM_SORT      = 3000

MAX_CPU           = 1000 
MAX_CPU_JELLYFISH = 10
MAX_CPU_SORT      = 10

# GET THE METHOD USED FOR DETECT DE KMERS
if DIFF_METHOD == "DESeq2":
    TEST_DIFF_SCRIPT   = BIN_DIR + "/DESeq2_diff_method.R"
elif DIFF_METHOD == "Ttest":
    TEST_DIFF_SCRIPT   = BIN_DIR + "/Ttest_diff_method.R"
else:
    sys.exit("Invalid value for 'diff_method', possible choices are: 'DESeq2' and 'Ttest'")

# VERIFY LIB_TYPE VALUE
if LIB_TYPE not in ['rf', 'fr', 'unstranded']:
    sys.exit("Invalid value for 'lib_type', possible choices are: 'rf', 'rf' and 'unstranded'")

rule all:
  input: MERGED_DIFF_COUNTS


###############################################################################
#
# SOFTWARE INSTALLATION
#
rule compile_joinCounts:
  output: JOIN_COUNTS
  run:
    shell("cd share/joinCounts && make")
    shell("ln -s -f ../share/joinCounts/joinCounts bin/")

rule compile_mergeTags:
  output: MERGE_TAGS
  input: "share/mergeTags/mergeTags.c"
  run:
    shell("cd share/mergeTags && make")
    shell("ln -s -f ../share/mergeTags/mergeTags bin/")

rule compile_TtestFilter:
  input: "share/TtestFilter/TtestFilter.c"
  output: TTEST_FILTER
  run:
    shell("cd share/TtestFilter && make")
    shell("ln -s -f ../share/TtestFilter/TtestFilter bin/")

###############################################################################
#
# UTILS
# Creates :
#   1. A tabulated file with the sample names and conditions
#   2. A tabulated file with the sample names and normalization factors
#   3. A tabulated file with the sample names, condition and normalization factors

rule sample_conditions:
  output: SAMPLE_CONDITIONS
  run:
    with open(output[0], "w") as f:
      f.write("\t".join(["sample",CONDITION_COL]) + "\n")
      for sample in config["samples"]:
        f.write("\t".join([sample["name"],sample[CONDITION_COL]]) + "\n")

rule compute_normalization_factors:
  input:
    raw_counts = RAW_COUNTS
  output: 
    nf      = NORMALIZATION_FACTORS,
    tmp_dir = temp(TMP_DIR + "/NF")
  log: LOGS + "/compute_norm_factors.log"
  script: COMPUTE_NF

rule sample_conditions_full:
  output:
    SAMPLE_CONDITIONS_FULL
  input:
    sample_conditions     = SAMPLE_CONDITIONS,
    normalization_factors  = NORMALIZATION_FACTORS
  shell: "join --header {input.sample_conditions} {input.normalization_factors} > {output}"


###############################################################################
#
# STEP 2: KMER COUNTS
#         Compiple DEkupl counter and count k-mers on all the samples
#
rule jellyfish_count:
  input:
    r1 = FASTQ_DIR + "/{sample}" + R1_SUFFIX,
    r2 = FASTQ_DIR + "/{sample}" + R2_SUFFIX
  output: COUNTS_DIR + "/{sample}.jf"
  log:
    exec_time = LOGS + "/{sample}_jellyfishRawCounts_exec_time.log"
  threads: MAX_CPU_JELLYFISH
  resources: ram = MAX_MEM_SORT
  run:
    options = "-L 2 -m {KMER_LENGTH} -s 10000 -t {threads} -o {output} -F 2"

    shell("echo -e \"******\" >{log.exec_time}")
    shell("echo -e \"start of rule jellyfish_count (raw counts) : $(date)\n\" >>{log.exec_time}")

    if LIB_TYPE == "rf":
      options += " <({ZCAT} {input.r1} | {REVCOMP}) <({ZCAT} {input.r2})"
      shell("echo -e \"R1 is rev comp\n\" >>{log.exec_time}")
    elif LIB_TYPE == "fr":
      options += " <({ZCAT} {input.r1}) <({ZCAT} {input.r2} | {REVCOMP})"
      shell("echo -e \"R2 is rev comp\n\" >>{log.exec_time}")
    elif LIB_TYPE == "unstranded":
      options += " -C <({ZCAT} {input.r1}) <({ZCAT} {input.r2})"
    else:
      sys.exit('Unknown library type')

    shell("{JELLYFISH_COUNT} " + options)
    shell("echo -e \"\nend of rule jellyfish_count : $(date)\n\" >>{log.exec_time}")
    shell("echo -e \"******\" >>{log.exec_time}")

rule jellyfish_dump:
  input: COUNTS_DIR + "/{sample}.jf"
  output: COUNTS_DIR + "/{sample}.txt.gz"
  threads: MAX_CPU_SORT
  resources: ram = MAX_MEM_SORT
  log :
    exec_time = LOGS + "/{sample}_jellyfishDumpRawCounts_exec_time.log"
  shell: """

         echo -e \"******\" >{log.exec_time}
         echo -e \"start of rule jellyfish_dump : $(date)\n\" >>{log.exec_time}

         {JELLYFISH_DUMP} -c {input} | {SORT} -k 1 -S {resources.ram}G --parallel {threads}| pigz -p {threads} -c > {output}

         echo -e \"\nend of rule jellyfish_dump : $(date)\n\" >>{log.exec_time}
         echo -e \"******\" >>{log.exec_time}

         """

rule join_counts:
  input:
    fastq_files = expand("{counts_dir}/{sample}.txt.gz",counts_dir=COUNTS_DIR,sample=SAMPLE_NAMES),
    binary = JOIN_COUNTS
  params:
    sample_names = "\t".join(SAMPLE_NAMES)
  output: RAW_COUNTS
  log :
    exec_time = LOGS + "/joinRawCounts_exec_time.log"
  run:
    shell("echo 'tag\t{params.sample_names}' | gzip -c > {output}")
    shell("""

           echo -e \"******\" >{log.exec_time}
           echo -e \"start of rule join_counts : $(date)\n\" >>{log.exec_time}

           {JOIN_COUNTS} -r {MIN_REC} -a {MIN_REC_AB} \
          {input.fastq_files} | gzip -c >> {output}

          echo -e \"\nend of rule dekupl_counter : $(date)\n\" >>{log.exec_time}
          echo -e \"******\" >>{log.exec_time}

          """)

###############################################################################
#
# STEP 4: SELECT DIFFERENTIALLY EXPRESSED K-MERS
#         Apply a T-test on all new k-mers to select only those that are
#         differentially expressed.
#
rule test_diff_counts:
  input:
    counts = RAW_COUNTS, # this is the difference with the transcriptome version, we take the raw counts here, not the transcriptome-filtered ones
    sample_conditions = SAMPLE_CONDITIONS_FULL,
    binary = TTEST_FILTER
  output:
    diff_counts = DIFF_COUNTS,
    pvalue_all  = PVALUE_ALL,
    tmp_dir     = TMP_DIR + "/test_diff"
    #tmp_dir     = temp(TMP_DIR + "/test_diff")
  params:
    conditionA  = CONDITION_A,
    conditionB  = CONDITION_B,
    pvalue_threshold = PVALUE_MAX,
    log2fc_threshold = LOG2FC_MIN,
    chunk_size = CHUNK_SIZE,
  threads: MAX_CPU
  log: LOGS + "/test_diff_counts.logs"
  script: TEST_DIFF_SCRIPT

rule merge_tags:
  input:
    counts = DIFF_COUNTS,
    binary = MERGE_TAGS
  output:
    MERGED_DIFF_COUNTS
  log:
    exec_time = LOGS + "/merge_tags_exec_time.log"
  run:
    options = "-k {KMER_LENGTH}"

    if LIB_TYPE == "unstranded":
      options += " -n"

    shell("echo -e \"******\" >{log.exec_time}")
    shell("echo -e \"start of merge_tags : $(date)\n\" >>{log.exec_time}")
    shell("{MERGE_TAGS} " + options + " {input.counts} | gzip -c > {output}")
    shell("echo -e \"\nend of merge_tags : $(date)\n\" >>{log.exec_time}")
    shell("echo -e \"******\" >>{log.exec_time}")
