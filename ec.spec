#### common ####
# input_for_database: filtered subreads in fasta or fastq format
input_for_database pacbio.filtered_subreads.fasta

# min_len_for_query: the subreads longer than or equal to this value will be corrected
min_len_for_query 500

#if you don't know the estimated genome size, give a large number
estimated_genome_size 50000
#if you don't know the estimated depth of coverage, give 0
estimated_depth 100

# ca_path: where Celera Assembler exist in
# If have it and want to run a full assembly with it,
# set this variable
# ca_path /home/imai/wgs/Linux-amd64/bin/

# the number of processes used by all vs. all alignment
# = 'partition' (in single node mode)
# = 'pre_partition' * 'partition' (in many node mode)
pre_partition 2
partition 12

#### many node mode (advanced) ####

#sge: options for all the SGE jobs (used by ezez4qsub_vx1.pl)
#sge -soft -l ljob,lmem,sjob
#queue_req: additional options for all the SGE jobs (used by ezez4qsub_vx1.pl and ezez4makefile_v4.pl)
#queue_req -l s_vmem=4G -l mem_req=4
#longestXx_queue_req: if valid, displaces queue_req (used by ezez4qsub_vx1.pl)
#longestXx_queue_req -l s_vmem=64G -l mem_req=64
#BLAST_RREQ: additional options for SGE jobs of all vs. all alignment (used by ezez4qsub.pl and ezez4makefile_v4.pl)
#BLAST_RREQ -pe def_slot 4
#ec_rreq: options for error correction (used by ezez4makefile_v4.pl)
#ec_rreq -l s_vmem=4G -l mem_req=4

#### common (advanced) ####

# used by blastn
word_size 18
evalue 1e-50
num_threads 1
max_target_seqs 100

#valid_voters 11

#trim: both ends of each alignment by blastn will be trimmed 'trim' bases to detect chimeric reads
trim 42

# if not 0, use only one subread per one molecule
use_one_subread 1

# direct_vote & copy_blastdb are used by ezez4makefile_v4.pl
direct_vote 0
# skip writing the blast results once to disk before selecting
# voters (default 1), or write to the disk to allow multiple use
# of the blast results for different number of voters (0)
copy_blastdb 1
# copy the blastdb to $TMP, presumably set by gridengine to local node,
# and use it during execution. This could reduce NFS access per job, but may
# increase the total transfer if each node is running multiple jobs.

