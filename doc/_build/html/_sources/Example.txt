=====================================
Example: assemble phage genome
=====================================

After installing Sprai, assemble a small genome.

For example, we showed how to assemble the phage genome.

Prepare data
=======================
Go to `pacbiotoca wiki <http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=PacBioToCA>`_ and download phage data ( \http://www.cbcb.umd.edu/software/PBcR/data/sampleData.tar.gz ).

::

   mkdir tmp
   cd tmp
   wget http://www.cbcb.umd.edu/software/PBcR/data/sampleData.tar.gz
   tar xvzf sampleData.tar.gz

Convert fasta to fastq. You can use a fa2fq.pl script in Sprai.

(Sprai version 0.9.9.13 or newer can be fed both fastq and fasta format. So you can skip converting fasta to fastq.)

::

   fa2fq.pl sampleData/pacbio.filtered_subreads.fasta > pacbio.filtered_subreads.fq

Prepare Sprai
==========================================

::

   mkdir tmp2
   cd tmp2
   ln -s ../pacbio.filtered_subreads.fq all.fq
   cp /your/sprai/directory/asm.spec .
   cp /your/sprai/directory/ec.spec .

Open ec.spec and change values.
::

   #>- params -<#
   input_fastq all.fq
   estimated_genome_size 50000
   estimated_depth 100
   partition 12
   evalue 1e-50
   trim 42
   ca_path /path/to/your/wgs/Linux-amd64/bin/
   word_size 18
   #>- params -<#

*input_fastq* is your input file name.

*estimated_genome_size* is the number of nucleotides of your target.
If you do not know it, set large number. For example, set 1e+12.

*estimated_depth* is the depth of coverage of input_fastq of your target.
If you do not know it, set 0.

*partition* is the number of processors Sprai uses.

*evalue* is used by blastn.

*trim* is the number of nucleotides Sprai cut from both sides of alignments.

*ca_path* is the path to your wgs-assembler (Celera Assembler) installed.

*word_size* is used by blastn.

Correct errors & assemble
==============================
::

   ezez_vx1.pl ec.spec pbasm.spec > log.txt 2>&1 &

ezez_vx1.pl outputs into a *result_yyyymmdd_hhmmss* directory. The improved reads will be in a *c01.fin.idfq.gz* file and contigs will be in a *CA/9-terminator/asm.ctg.fasta* file.

We explain temporary files in a *tmp* directory.
Sprai gives ID to each read in *all.fq* and outputs *c00.idfq.gz*.  
Sprai partitions *all.fq* to *c00_00xx.fa*.
*c00.nin*, *c00.nhr* and *c00.nsq* are output by makeblastdb.
Sprai aligns *c00_00xx.fa* to *c00.nxx* database by using blastn, correct errors and output *c00_00xx.dfq.gz*.
dfq format contains the IDs of reads, bases, aligned depths and quality values.
(Current Sprai outputs dummy quality values.)
Sprai cuts low depth regions of each read, deletes'-' characters to convert to FASTQ format and outputs corrected reads in *c01.fin.idfq.gz*.
Sprai extracts longest 20X reads of the *estimated_genome_size* from *c01.fin.idfq.gz*.
And feed them to Celera Assembler.
Celera Assembler outputs files into *CA* directory.

If you only correct errors and don't assemble, do

::

   ezez_vx1.pl ec.spec -ec_only > log.txt 2>&1 &

or

::

   ezez_vx1.pl ec.spec > log.txt 2>&1 &

After error correction, if you want to assemble corrected reads using Celera Assembler, do

::

   ca_ikki_v5.pl pbasm.spec estimated_genome_size \
     -d directory in which fin.idfq.gzs exist \
     -ca_path /path/to/your/wgs/Linux-amd64/bin \
     -sprai_path the path to get_top_20x_fa.pl installed 



Find contigs
===================
You will find contigs in a *CA/9-terminator/asm.ctg.fasta* file.
Files in *CA* are produced by the wgs-assembler (Celera Assembler).
Read `the wgs-assembler document <http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=Main_Page>`_ for details. 

Get N50 value
================
If you install Statistics::Descriptive module to a directory Celera Assembler can see, Celera Assembler outputs N50 value and so on in a CA/do_*_c01.fin.top20x.log file.
You will find lines in the file like below:
::

   [Contigs]
   TotalContigsInScaffolds         1         
   TotalBasesInScaffolds           49848     
   TotalVarRecords                 0         
   MeanContigLength                49848
   MinContigLength                 49848     
   MaxContigLength                 49848     
   N25ContigBases                  49848
   N50ContigBases                  49848
   N75ContigBases                  49848

These are the statistics of contigs of the assembly.

Notes
==========================================
If you would like to use more than about 1000 processors, we recommend to use a *pre_partition* parameter in the ec.spec file like:

::

   # ec.spec
   pre_partition 4
   partition 300

And

::

   ezez4qsub_vx1.pl ec.spec pbasm.spec > log 2>&1 &

In this example, Sprai will use 4*300 = 1200 processors.

