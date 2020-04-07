=======
FAQ
=======
.. | Q. How many files does sprai create in each cycle?
.. | A. Given P partition in ec.spec, the number of jobs in each cycle is 1 (by makeblastdb)+1 (by partition_fa.pl)+P (by blastn, bfmt72s, nss2v and myrealigner)+ceil(P/CURRENT) (by dfq->fq)+1 (by cat) = P+ceil(P/CURRENT)+3, where CURRENT is 1 in ezez_v3.pl and ezez4make.pl and 12 in ezez4qsub_v4.pl. The number of files created in each cycle is {3 (by makeblastdb)+P (by partition_fa.pl)+P (by blastn, bfmt72s, nss2v and myrealigner)+ceil(P/CURRENT) (by dfq->fq)+1 (by cat)}+2*(# of jobs)(.e and .o files)+2*P(if you specify -pe thet blastn,... will use, and otherwise 0) = {2*P+ceil(P/CURRENT)+4}+2*(P+ceil(P/CURRENT)+3)+(2*P (if -pe) or 0 (otherwise)

Q. How much sequencing coverage depth would be required?
---------------------------------------------------------
A. It depends. If you can afford to sequence 100x of the estimated genome size, that would be the best.
When the genome is large (e.g., 400Mbp or larger), you may reduce the depth to 30x or 40x.
If you want to sequence less than 20x, you might consider hybrid assembly in which you combine Illumina reads and PacBio reads.
Higher sequencing depth than 200x does not usually improve the assembled result.
Note that the sequencing depth itself is not a good predictor of the assembly result.
What matters is the sequencing depth of reads longer than repetitive elements in a target genome.
If reads are short and only few of them are longer than some threshold, namely 2kbp (for bacteria), even 100x data would be beated by 20x data with longer reads.
We strongly recommend the use of the size selection protocol using Blue Pippin, which is highly recommended by PacBio.

Q. How large genomes can Sprai assemble?
--------------------------------------------------------
A. The largest genome we tested is 2Gbp, but there is no theoretical limit for Sprai (although you need more computational resources as the target genome gets larger).
Error-correction process is highly scalable, so most probably, Celera assembler would be the bottleneck.

Q. Can Sprai assemble both PacBio CLRs and other reads such as Illumina short reads?
-------------------------------------------------------------------------------------
A. No. Sprai is designed to take only PacBio long reads.
If you want to combine both reads, you may first correct errors in PacBio long reads by Sprai, and then you can use whatever genome assembler that accepts both error-corrected long reads and Illumina reads.

Q. Why is the number of the output nucleotides by Sprai smaller than the number of the nucleotides in the input file?
----------------------------------------------------------------------------------------------------------------------
A. Sprai trims nucleotides that will not be useful for subsequent analysis such as de novo assembly or mapping to a reference genome. Such nucleotides include low coverage regions (valid_depth paramter; default = 4), regions of too short reads (valid read length: default = 500 bp), and both ends of reads (default = 42 bp).
The exact ratio of such nucleotides depends on input data, but however, it usually falls between 10% to 50% in our experience.
