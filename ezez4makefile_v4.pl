#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);

# Environmental variables (deprecated):
#   MAKEDB_RREQ    options for qsub in creating databases (blast)
#   BLAST_RREQ     options for qsub in searching databases (blast)
#   BLAST_OPT      options for blast in searching databases
#   ASSEMBLY_RREQ  options for qsub in assembly (celera)

my $valid_depth=4;
my $valid_read_length=500;

my $now = `date +%Y%m%d_%H%M%S`;
chomp $now;
my $submit=0;
my $jobs=0;
my $opt_help;

GetOptions(
  "now=s"=>\$now, # string
  "jobs" => \$jobs,
  "submit" => \$submit,
  "h" => \$opt_help
);

my %params;

my @msgs = (
  "USAGE: <this> <ec.spec> <asm.spec> or <this> <ec.spec>",
  "[-jobs: create job files for a grid engine]",
  "[-submit: actually submit the job to a grid engine]",
  "[-now string: a suffix of some directories]",
  "[-h: show this message]",
);

if(@ARGV == 0 || @ARGV >= 3 || $opt_help){
  my $msg = join "\n\t", @msgs;
  die "$msg\n";
  #die "USAGE: <this> <ec.spec> <asm.spec>\n\t[-jobs: create job files for a grid engine]\n\t[-submit: actually submit the job to a grid engine]\n\t[-now string: a suffix of some directories]\n";
}

my $spec = "/dev/null";
if(defined($ARGV[1])){
  $spec = $ARGV[1];
}

{
  my $ec_spec = $ARGV[0];
  open my $fh,"<",$ec_spec or die "cannot open $ec_spec :$!\n";
  while($_ = <$fh>){
    next if($_ =~ /^\s+$/);
    next if($_ =~ /^\s*#/);
    chomp;
    my @line = split /\s+/,$_;
    for(my $i=0; $i<@line; ++$i){
      if($line[$i] =~ /^\s*#/){
        @line = @line[0..$i-1];
        last;
      }
    }
    $params{$line[0]}=join(" ",@line[1..@line-1]);
  }
  close $fh;
}

my @input_libs;
my $nlibs;
my $pre_partition=1;
my $partition=12;
my $word_size="";
my $evalue=1e-50;
my $num_threads=1;
my $max_target_seqs=100;
my $valid_voters=11;
my $max_read_length=128000;
my @voters=();
my $direct_vote = 1;
my $copy_blastdb = 0;
my $trim=42;
my $estimated_genome_size=0;
my $estimated_depth = 0;
my @gsizes=();
my $ca_path="/opt/wgs-8.1/Linux-amd64/bin/";
my $blast_path="/usr/bin/";
my $sprai_path="/usr/lib/sprai/";
my $param_makedb_rreq = ""; # not so much resource 
my $param_blast_rreq = ""; # desired to have a high CPU slot for latency optimization
my $param_ec_rreq = ""; # not so much resource 
my $param_blast_opt = "";
my $param_gettop_rreq = ""; # requires relatively high memory
my $param_assembly_rreq = "";
my $filter_same_lib = "| cat";
my $min_len_for_query=1;
my $max_len_for_query=1000000000000000;
my $use_one_subread=0;

my $command="";

if(defined($params{input_for_database})){
  # new multi-lib mode
  # use complementary set of data for each library to avoid correction of
  # chimeric data as authentic data
  my $libs;
  $libs = $params{input_for_database};
  @input_libs = split(/\s+/, $libs);
  $nlibs = @input_libs;
  printf STDERR "number of libraries: $nlibs\n";
  printf STDERR "libs: @input_libs\n";
}else{
  die "specify input_for_database in ec.spec\n";
}

if(defined($params{estimated_genome_size})){
  $estimated_genome_size = $params{estimated_genome_size};
  @gsizes = split(/\s+/, $estimated_genome_size);
}
else{
  die "specify estimated_genome_size in ec.spec\n";
}
if(defined($params{estimated_depth})){
  $estimated_depth=$params{estimated_depth};
}

foreach my $eg (@gsizes){
  if($eg <= 0){die "estimated_genome_size must be > 0\n";}
}

if(defined($params{pre_partition})){
  $pre_partition = $params{pre_partition};
}
if(defined($params{partition})){
  $partition = $params{partition};
}
if(defined($params{word_size})){
  $word_size = $params{word_size};
}
if(defined($params{evalue})){
  $evalue = $params{evalue};
}
if(defined($params{num_threads})){
  $num_threads = $params{num_threads};
}
if(defined($params{max_target_seqs})){
  $max_target_seqs = $params{max_target_seqs};
}
if(defined($params{valid_voters})){
  $valid_voters = $params{valid_voters};
  @voters = split(/\s+/, $valid_voters);
  undef($valid_voters);
}
else{
  if($estimated_depth == 0) {
    my $n_base = 0;
    foreach my $input_fastq(@input_libs){
     $n_base += -s $input_fastq;
    }
    $n_base /= 2;
    $estimated_depth=$n_base/$estimated_genome_size;
  }
  $valid_voters = int(0.8*$estimated_depth);
  $valid_voters = ($valid_voters < 11) ? 11 : $valid_voters;
  $valid_voters = ($valid_voters > 30) ? 30 : $valid_voters;
  @voters = ($valid_voters);
  undef($valid_voters);
}
if(defined($params{direct_vote})){
  $direct_vote = $params{direct_vote};
}
if(defined($params{copy_blastdb})){
  $copy_blastdb = $params{copy_blastdb};
}
if(defined($params{filter_same_lib})){
  if($params{filter_same_lib} =~ /t/i){
    $filter_same_lib = "| perl -e 'while(<>){\@a=split;next if(\$\$a[0] ne \$\$a[3] && \$\$a[0]=~s/_.*//&& \$\$a[3]=~s/_.*//&& \$\$a[0] eq \$\$a[3]);print;}'"
  }
}
if(defined($params{max_read_length})){
  $max_read_length = $params{max_read_length};
}
if(defined($params{trim})){
  $trim = $params{trim};
}
if(defined($params{use_one_subread})){
  $use_one_subread = $params{use_one_subread};
}
if(defined($params{ca_path})){
  $ca_path = $params{ca_path};
}
if(defined($params{blast_path})){
  $blast_path = $params{blast_path};
}
if(defined($params{sprai_path})){
  $sprai_path = $params{sprai_path};
}
if(exists $ENV{'MAKEDB_RREQ'}) {
  $param_makedb_rreq = $ENV{'MAKEDB_RREQ'};
  print STDERR "WARNING: environmental variable MAKEDB_RREQ is deprecated. Please set queue_req in ec.spec file.\n";
}
if(defined($params{queue_req})) {
  $param_assembly_rreq = $params{queue_req};
  $param_blast_rreq = $params{queue_req};
  $param_makedb_rreq = $params{queue_req};
  $param_ec_rreq = $params{queue_req};
}
if(defined($params{BLAST_RREQ})) {
  $param_blast_rreq = $params{BLAST_RREQ};
}
if(exists $ENV{'BLAST_OPT'}) {
  $param_blast_opt = $ENV{'BLAST_OPT'};
  print STDERR "WARNING: environmental variable BLAST_OPT is deprecated. Please set BLAST_OPT in ec.spec file.\n";
}
if(defined($params{blast_opt})) {
  $param_blast_opt = $params{blast_opt};
}
if($param_blast_opt =~ /num_threads/) {
  print STDERR "ERROR: you cannot specify -num_threads in blast_opt. Please use num_threads in ec.spec file.\n";
  exit 2;
}
if($param_blast_opt =~ /word_size/) {
  print STDERR "ERROR: you cannot specify -word_size in blast_opt. Please use word_size in ec.spec file.\n"+
  exit 2;
}
if(defined($params{makedb_rreq})) {
  $param_makedb_rreq = $params{makedb_rreq};
}
if(defined($params{ec_rreq})) {
  $param_ec_rreq = $params{ec_rreq};
}
if(defined($params{gettop_rreq})) {
  $param_gettop_rreq = $params{gettop_rreq};
}else{
  $param_gettop_rreq = $param_makedb_rreq;
}
if(exists $ENV{'ASSEMBLY_RREQ'}) {
  $param_assembly_rreq = $ENV{'ASSEMBLY_RREQ'};
  print STDERR "WARNING: environmental variable ASSEMBLY_RREQ is deprecated. Please set assembly_rreq in ec.spec file.\n";
}
if(defined($params{assembly_rreq})) {
  $param_assembly_rreq = $params{assembly_rreq};
}
if(defined($params{min_len_for_query})){
  $min_len_for_query = $params{min_len_for_query};
}
if(defined($params{max_len_for_query})){
  $max_len_for_query = $params{max_len_for_query};
}

printf STDERR ("-- params --\n");
printf STDERR ("input_fastq %s\n", join(" ", @input_libs));
printf STDERR ("estimated_genome_size %s\n",join(" ", @gsizes));
printf STDERR ("pre_partition %s\n",$pre_partition);
printf STDERR ("partition %s (total %d)\n",$partition, $pre_partition * $partition);
$partition *= $pre_partition;
printf STDERR ("word_size %d\n",$word_size);
printf STDERR ("evalue %g\n",$evalue);
printf STDERR ("num_threads %d\n",$num_threads);
if($min_len_for_query){
  printf STDERR ("min_len_for_query %d\n",$min_len_for_query);
}
if($max_len_for_query){
  printf STDERR ("max_len_for_query %d\n",$max_len_for_query);
}
if($max_target_seqs){
  printf STDERR ("max_target_seqs %d\n",$max_target_seqs);
}
printf STDERR ("valid_voters %s\n",join(" ", @voters));
printf STDERR ("trim %d\n",$trim);
printf STDERR ("use_one_subread %d\n",$use_one_subread);
printf STDERR ("ca_path %s\n",$ca_path);
printf STDERR ("blast_path %s\n",$blast_path);
printf STDERR ("sprai_path %s\n",$sprai_path);
printf STDERR ("BLAST_RREQ %s\n", $param_blast_rreq);
printf STDERR ("blast_opt %s\n", $param_blast_opt);
printf STDERR ("makedb_rreq %s\n", $param_makedb_rreq);
printf STDERR ("gettop_rreq %s\n", $param_gettop_rreq);
printf STDERR ("assembly_rreq %s\n", $param_assembly_rreq);
printf STDERR ("ec_rreq %s\n", $param_ec_rreq);
printf STDERR ("direct_vote %d\n", $direct_vote);
printf STDERR ("copy_blastdb %d\n", $copy_blastdb);
printf STDERR ("-- params --\n");

my $log="./ezez4makefile.log";

if(length($blast_path) > 0 && $blast_path !~ /\/$/){
  $blast_path .= "/";
}
if(length($sprai_path) > 0 && $sprai_path !~ /\/$/){
  $sprai_path .= "/";
}


my $date="";
my $message="";

my $bindir=$sprai_path;
my $path2blast=$blast_path;

if(!-d $bindir){
  die "$bindir does not exist\n"
}
if(!-d $path2blast){
  die "$path2blast does not exist\n"
}

my $pwd = `pwd`;
chomp $pwd;

my $outdir = "result";
if(!-d $outdir){
  mkdir "$outdir" or die "cannot mkdir $outdir: $!\n";
}

my @allscripts=();
my $script;
my @alloutputfiles=();
my @allerrfiles=();
my $errfile;

my $makefile = "$pwd/Makefile";
open my $mh, ">", $makefile or die "cannot open $makefile: $!\n";
my $prefix="c01.fin";
my $rreq = "";
if(defined $param_assembly_rreq) {
  $rreq = "SGE_RREQ=\"$param_assembly_rreq\" ";
}
my @fintargets = ();
foreach my $egsize (@gsizes){
  foreach my $vv (@voters){
    my $ps = "$prefix.v$vv.20x$egsize";
    my $target = "$outdir/$ps/9-terminator/$ps.scf.fasta";
    push(@fintargets, $target);
  }
}
printf $mh "all: %s\n\n", join(" ", @fintargets);
foreach my $egsize (@gsizes){
  foreach my $vv (@voters){
    my $ps = "$prefix.v$vv.20x$egsize";
    my $target = "$outdir/$ps/9-terminator/$ps.scf.fasta";
    printf $mh "$outdir/$ps/9-terminator/$ps.scf.fasta: $outdir/$ps.frg\n";
    printf $mh "\t${rreq}$ca_path/runCA -dir $outdir/`basename \$< .frg` -p `basename \$< .frg` -s $spec \$< \n";
  }
}


my $libi;
my @current_dep;
my $tmp_dfq;
my $tmp_blastdb = sprintf("tmp_blastdb");
if(!-d "$outdir/$tmp_blastdb"){
    mkdir "$outdir/$tmp_blastdb" or die "cannot mkdir $outdir/$tmp_blastdb: $!\n";
} # this is not per library directory; thus outside the loop;
my $aggr_db ="$outdir/$tmp_blastdb/all_libs";
my $nalf ="$aggr_db.nal";
for($libi = 0; $libi < $nlibs; $libi++){
# $libi is a safe string that represents an integer; safe to embed directly in a string
  my $iprefix = "$outdir/c_l$libi";
  my $idfqgz = "$iprefix.idfq.gz";
  my $qidfqgz = "$iprefix.q.idfq.gz";
  # longer scope required
  {
    # construct a idfq file from a fastq file

    my $PG= "$bindir/fq2idfq.pl";
    my $PG0="$bindir/fa2fq.pl";
    my $PG2="$bindir/fqfilt.pl";
    my $PG3="$bindir/dumbbell_filter.pl";
    if(!$use_one_subread){
      $PG3="cat";
    }
    my $command = "STMP=`mktemp -d --tmpdir sprai-XXXXXX` ; set -o pipefail; $PG0 \$\< | $PG3 - | $PG - --prefix $iprefix | gzip -c -1 > \$\$STMP/c_l${libi}.idfq.gz && cp -p \$\$STMP/c_l${libi}.idfq.gz \$@ && openssl dgst \$\$STMP/c_l${libi}.idfq.gz \$@; status=\$\$?; rm -fr \$\$STMP;exit \$\$status";

    $rreq = "";
    if(defined $param_makedb_rreq) {
      $rreq = "SGE_RREQ=\"$param_makedb_rreq\" ";
    }
    printf $mh "$idfqgz: %s\n", $input_libs[$libi];
    printf $mh "\t${rreq}$command\n\n";

    if($min_len_for_query > 0){
      $command = "STMP=`mktemp -d --tmpdir sprai-XXXXXX`; set -o pipefail; gzip -d -c \$\< | $PG2 - $min_len_for_query -max_len $max_len_for_query | gzip -c -1 > \$\$STMP/c_l${libi}.q.idfq.gz && cp -p \$\$STMP/c_l${libi}.q.idfq.gz \$@ && openssl dgst \$\$STMP/c_l${libi}.q.idfq.gz \$@; status=\$\$?; rm -fr \$\$STMP; exit \$\$status";
      printf $mh "$qidfqgz: %s\n", $idfqgz;
      printf $mh "\t${rreq}$command\n\n";
    }else{
      $qidfqgz=$idfqgz;
    }
  }

  my $tmp_idfq;
  my $dfq2fa = "$bindir/dfq2fq_v2.pl -f -";
  my $addlibn = "sed -e 's/^>/>l${libi}_s/'";
  #construct blast database from idfq file
  {
    $rreq = "";
    if(defined $param_makedb_rreq) {
      $rreq = "SGE_RREQ=\"$param_makedb_rreq\" ";
    }
    my $PG = "$path2blast/makeblastdb";
    $command = "gzip -d -c $idfqgz | $dfq2fa | $addlibn | $PG -in - -dbtype nucl -out $outdir/$tmp_blastdb/c_l$libi -title c_l$libi ";

    my $makeblastdone = "$outdir/c${libi}_imakeblastdb.done";

    printf $mh "$makeblastdone: $idfqgz\n";
    printf $mh "\t${rreq}$command\n";
    printf $mh "\ttouch \$\@\n\n";
  }

  if($libi==0){  # %.nal depends on makeblastdb.done of other libs
# construct aggregate complementary blast database
# don't need real data. just write to .nal file 
# write a file with the two lines stating TITLE and DBLIST
#TITLE libs_excludinglib1
#DBLIST lib2 lib3 
     my @clibs=();
     my @current_dep=();
     my $k;
     for($k=0; $k < $nlibs; $k++){
       push(@clibs, "c_l$k");
       push(@current_dep, "$outdir/c${k}_imakeblastdb.done");
     }
     printf $mh "$nalf: %s\n", join(" ", @current_dep);
     printf $mh "\techo TITLE all_libs > \$\@\n";
     printf $mh "\techo DBLIST %s >> \$\@\n\n", join(" ", @clibs);
  }

  my $tmp_fa;
  {
    $rreq = "";
    if(defined $param_makedb_rreq) {
      $rreq = "SGE_RREQ=\"$param_makedb_rreq\" ";
    }
    $tmp_fa = "tmp_fa_l$libi";
    if(!-d "$outdir/$tmp_fa"){
      mkdir "$outdir/$tmp_fa" or die "cannot mkdir $outdir/$tmp_fa: $!\n";
    }
    my $PG2 = "$bindir/partition_fa.pl";
    $command = "${rreq}gzip -d -c $qidfqgz | $dfq2fa | $addlibn | $PG2 - $partition -p $outdir/$tmp_fa/c_l$libi";

    my $ipartdone = "$outdir/c_l${libi}_partition_fa.done";

    printf $mh "$ipartdone: $qidfqgz\n";
    printf $mh "\t$command\n";
    printf $mh "\ttouch \$\@\n\n";
  }

  { # main error correction job
    $rreq = "";
    if(defined $param_blast_rreq && $param_blast_rreq ne '') {
      $rreq = "SGE_RREQ=\"$param_blast_rreq\" ";
    }
    my $tmp_evalue;
    $tmp_evalue = $evalue;
    $tmp_dfq = "tmp_dfq_l$libi";
    if(!-d "$outdir/$tmp_dfq"){
      mkdir "$outdir/$tmp_dfq" or die "cannot mkdir $outdir/$tmp_dfq: $!\n";
    }
    my $BLAST_CMD = "$path2blast/blastn -dbsize 1 -num_threads $num_threads";
    if($word_size){
      $BLAST_CMD .= " -word_size $word_size";
    }
    if($max_target_seqs){
      $BLAST_CMD .= " -max_target_seqs $max_target_seqs";
    }
    if(defined $param_blast_opt) {
      $BLAST_CMD .= " $param_blast_opt";
    }
    my $outfmt = "7 qseqid sstart send sacc qstart qend bitscore evalue pident qseq sseq";
    $BLAST_CMD .= " -evalue $tmp_evalue -outfmt '$outfmt'";
    my $PG2=sprintf("$bindir/bfmt72s -c $trim -u -i",$libi);
    if($copy_blastdb){ # copy the blastdb per each job to reduce the number of NFS access; this might increase transer if many jobs are run on the same node 
      $BLAST_CMD .= " -db \$\$STMP/$aggr_db";
    }else{
      $BLAST_CMD .= " -db $aggr_db";
    }
    if($direct_vote){
      foreach my $valid_voters(@voters){
        my $PG3="$bindir/nss2v_v3 -r $max_read_length -v $valid_voters";
        my $PG4="$bindir/myrealigner -f -B $valid_voters -b 3 -d 0.5 -l $max_read_length";
        if(!-d "$outdir/$tmp_dfq.v$valid_voters"){
          mkdir "$outdir/$tmp_dfq.v$valid_voters" or die "cannot mkdir $outdir/$tmp_dfq.v$valid_voters: $!\n";
        }
        for(my $i=0; $i < $partition; ++$i){
          my $inputfile = sprintf("$outdir/$tmp_fa/c_l${libi}_%04d.fa", $i);
          printf $mh "$outdir/$tmp_dfq.v$valid_voters/c_l${libi}_%04d.v$valid_voters.dfq.gz: $outdir/c_l${libi}_partition_fa.done $nalf\n", $i;
          printf $mh "\t${rreq} STMP=`mktemp -d --tmpdir sprai-XXXXXX`; set -o pipefail; mkdir -p \$\$STMP/$outdir $outdir/$tmp_dfq.v$valid_voters &&";
          if($copy_blastdb){ printf $mh " cp -rp $outdir/$tmp_blastdb \$\$STMP/$outdir &&";}
          printf $mh " cp -p $inputfile \$\$STMP/input && openssl dgst $inputfile \$\$STMP/input &&";
          printf $mh " $BLAST_CMD -query \$\$STMP/input $filter_same_lib | $PG2 - | $PG3 - | $PG4 - | gzip -c -1 > \$\$STMP/output.tmp && cp -p \$\$STMP/output.tmp \$\@ ; status=\$\$?; rm -fr \$\$STMP; exit \$\$status\n\n";
        }
      }
    }else{ # not direct_vote
      for(my $i=0; $i < $partition; ++$i){
        my $outputfile = sprintf("$outdir/$tmp_dfq/c_l${libi}_%04d.blast.sam.gz",$i);
        my $inputfile = sprintf("$outdir/$tmp_fa/c_l${libi}_%04d.fa", $i);
        printf $mh "$outputfile : $outdir/c_l${libi}_partition_fa.done $nalf\n";
        printf $mh "\t${rreq}";
        printf $mh " STMP=`mktemp -d --tmpdir sprai-XXXXXX`; set -o pipefail; mkdir -p \$\$STMP/$outdir &&";
        if($copy_blastdb){ printf $mh " cp -rp $outdir/$tmp_blastdb \$\$STMP/$outdir &&";}
        printf $mh " cp -p $inputfile \$\$STMP/input && openssl dgst $inputfile \$\$STMP/input &&";
        printf $mh " $BLAST_CMD -query \$\$STMP/input $filter_same_lib | $PG2 - | gzip -c -1 > \$\$STMP/output.tmp && cp -p \$\$STMP/output.tmp \$\@ ; status=\$\$?; rm -fr \$\$STMP; exit \$\$status\n\n", $i;
      }

      $rreq = ""; # smaller resource request
      if(defined $param_makedb_rreq) {
        $rreq = "SGE_RREQ=\"$param_makedb_rreq\" ";
      }
      foreach my $valid_voters(@voters){
        my $PG3="$bindir/nss2v_v3 -r $max_read_length -v $valid_voters";
        my $PG4="$bindir/myrealigner -f -B $valid_voters -b 3 -d 0.5";
        if(!-d "$outdir/$tmp_dfq.v$valid_voters"){
          mkdir "$outdir/$tmp_dfq.v$valid_voters" or die "cannot mkdir $outdir/$tmp_dfq.v$valid_voters: $!\n";
        }
        # instead of using pattern match write every pattern for tge_make
        for(my $i=0; $i < $partition; ++$i){
          printf $mh "$outdir/$tmp_dfq.v$valid_voters/c_l${libi}_%04d.v$valid_voters.dfq.gz: $outdir/$tmp_dfq/c_l${libi}_%04d.blast.sam.gz\n", $i, $i;
          printf $mh "\t${rreq}gzip -dc \$< | $PG3 - | $PG4 - | gzip -c -1 > \$\@.tmp && mv \$\@.tmp \$\@\n\n";
        }
      }
    }
  }

  {
    $rreq = "";
    if(defined $param_makedb_rreq) {
      $rreq = "SGE_RREQ=\"$param_makedb_rreq\" ";
    }
    foreach my $valid_voters(@voters){
      $tmp_idfq = "tmp_idfq_l$libi";
      if(!-d "$outdir/$tmp_idfq.v$valid_voters"){
        mkdir "$outdir/$tmp_idfq.v$valid_voters" or die "cannot mkdir $outdir/$tmp_idfq: $!\n";
      }
      for(my $i=0; $i < $partition; ++$i){
        printf $mh "$outdir/$tmp_idfq.v$valid_voters/c_l${libi}.refined_%04d.v$valid_voters.idfq.gz: $outdir/$tmp_dfq.v$valid_voters/c_l${libi}_%04d.v$valid_voters.dfq.gz\n", $i, $i;
        printf $mh "\t${rreq}gzip -d -c \$< | $bindir/dfq2fq_v2.pl - | gzip -c -1 > \$\@.tmp && mv \$\@.tmp \$\@\n\n";
      }
    }
  }

  foreach my $valid_voters(@voters){
    my $outputfile = "$outdir/c_l${libi}.1.v$valid_voters.idfq.gz";
    my @part_refined_dfqgz = (); # initialize each time
    for(my $i=0; $i < $partition; ++$i){
      my $pidfqfiles = sprintf("$outdir/$tmp_idfq.v$valid_voters/c_l${libi}.refined_%04d.v$valid_voters.idfq.gz", $i);
      push(@part_refined_dfqgz, $pidfqfiles);
    }
    my $files=join(" ", @part_refined_dfqgz);

    printf $mh "$outputfile : $files\n";
    printf $mh "\tcat \$^ > \$\@.tmp && mv \$\@.tmp \$\@\n\n";
  }

  {
    $rreq = "";
    if(defined $param_makedb_rreq) {
      $rreq = "SGE_RREQ=\"$param_makedb_rreq\" ";
    }
    $tmp_dfq = "tmp_dfq_l$libi";
    foreach my $valid_voters(@voters){
      my $tmp_finish = "tmp_finish_l$libi.v$valid_voters";
      for(my $i=0; $i < $partition; ++$i){
        if(!-d "$outdir/$tmp_finish"){
          mkdir "$outdir/$tmp_finish" or die "cannot mkdir $outdir/$tmp_finish: $!\n";
        }
        printf $mh "$outdir/$tmp_finish/c_l${libi}.refined_%04d.fin.v$valid_voters.idfq.gz: $outdir/$tmp_dfq.v$valid_voters/c_l${libi}_%04d.v$valid_voters.dfq.gz\n", $i, $i;
        printf $mh "\t${rreq}gzip -d -c \$< | $bindir/dfq2fq_v2.pl - -finish -valid_depth $valid_depth -valid_read_length $valid_read_length | gzip -c -1 > \$\@.tmp && mv \$\@.tmp \$\@\n\n";
      }
    }
  }

  foreach my $valid_voters(@voters){
    my $tmp_finish = "tmp_finish_l$libi.v$valid_voters";
    my @part_finidfqgz = ();
    for(my $i = 0; $i < $partition; ++$i){
      my $outputfile = sprintf("$outdir/$tmp_finish/c_l${libi}.refined_%04d.fin.v$valid_voters.idfq.gz",$i);
      push(@part_finidfqgz,$outputfile);
    }
    my $files=join(" ", @part_finidfqgz);
    my $outputfile = "$outdir/c_l${libi}.1.fin.v$valid_voters.idfq.gz";

    printf $mh "\n";
    printf $mh "$outputfile: $files\n";
    printf $mh "\tcat \$^ > \$\@.tmp && mv \$\@.tmp \$\@\n\n";
  }
}

my @cdep=();
foreach my $valid_voters(@voters){
  for($libi = 0; $libi < $nlibs; $libi ++){
    push(@cdep, "$outdir/c_l${libi}.1.fin.v$valid_voters.idfq.gz");
  }
}
my $ecdone = "ec.done";
printf $mh "$ecdone: %s\n", join(" ", @cdep);
printf $mh "\ttouch \$\@;\n\n";

my $ec_only = "ec_only: $ecdone";
printf $mh "$ec_only\n\n";

my $bashcommand="";
my $suffix = "top20x";


foreach my $estimated_genome_size (@gsizes){
  foreach my $valid_voters (@voters){
    my $combined = "$outdir/$prefix.v$valid_voters.20x$estimated_genome_size";
    {
      $rreq = "";
      if(defined $param_gettop_rreq) {
        $rreq = "SGE_RREQ=\"$param_gettop_rreq\" ";
      }
      my $PG1 = "$bindir/get_top_20x_fa.pl";
      my @vdep = ();
      for($libi = 0; $libi < $nlibs; $libi ++){
        push(@vdep, "$outdir/c_l${libi}.1.fin.v$valid_voters.idfq.gz");
      }
      printf $mh "%s.fq: %s\n", $combined, join(" ", @vdep);
      printf $mh "\t${rreq}STMP=`mktemp -d --tmpdir sprai-XXXXXX`;gzip -d -c \$^ > \$\$STMP/ec.idfq && $PG1 \$\$STMP/ec.idfq -l -g $estimated_genome_size -q > \$\@.tmp && mv \$\@.tmp \$\@; status=\$\$?; rm -fr \$\$STMP; exit \$\$status \n\n";

      printf $mh "%s.frg: %s.fq\n", $combined, $combined;
      printf $mh "\t$ca_path/fastqToCA -libraryname foo -technology pacbio-corrected -reads \$<  > \$\@\n\n";
    }
  }
}
my @frgs = ();
foreach my $egsize (@gsizes){
  foreach my $vv (@voters){
    my $frg = "$outdir/$prefix.v$vv.20x$egsize.frg";
    push(@frgs, $frg);
  }
}
printf $mh ".PHONY: prepfrg\n";
printf $mh "prepfrg: %s\n\n", join(" ", @frgs);

close $mh;


if($submit || $jobs){
  my $headstr ="#!/bin/sh\n#\$ -S /bin/sh\n#\$ -cwd\n";
# construct jobs for grid engines
# each job just calls a make command to construct certain targets
# $partition should be run as array job and accessed as $SGE_TASK_ID
# STR0=`expr $SGE_TASK_ID - 1`
# PART=`printf '%04d' $STR0`
  my $partgen = q(STR0=`expr $SGE_TASK_ID - 1`
PART=`printf '%04d' $STR0`
);
  my $jobdir = "jobs";
  if(!-d $jobdir){
    make_path($jobdir) or die "cannot make_path $jobdir: $!\n";
  }
  make_path('log');
  my $rreq = $param_makedb_rreq;
  my @idfqjobids = (); # store the jobids for each library

  make_path('log/idfqgz');
  for($libi = 0; $libi < $nlibs; $libi++){
    my $iprefix = "$outdir/c_l$libi";
    my $idfqgz = "$iprefix.idfq.gz";
    my $jf = "$jobdir/idfq$libi.job";
    open my $jh, ">", $jf or die "cannot open $jf: $!\n";
      print $jh $headstr;
      printf $jh "#\$ %s\n", $rreq;
      print $jh "make -f $makefile $idfqgz\n"; 
    close $jh;
    if($submit){
      my $retv=`qsub -o log/idfqgz -e log/idfqgz $jf`;
      if($retv =~ /Your job (\d+)/){
        push(@idfqjobids, $1);
        print "$1: qsub -o log/idfqgz -e log/idfqgz $jf\n";
      }else{
        die "job submission error: qsub $jf said $retv";
      }
    }
  }

  my @jobids = (); # store the jobids for each library
  make_path('log/makeblastdb');
  for($libi = 0; $libi < $nlibs; $libi++){
    my $jf = "$jobdir/blastdb_partition$libi.job";
    open my $jh, ">", $jf or die "cannot open $jf: $!\n";
      print $jh $headstr;
      printf $jh "#\$ %s\n", $rreq;
      print $jh "make -f $makefile $outdir/c${libi}_imakeblastdb.done\n"; 
    close $jh;
    if($submit){
      my $retv=`qsub -o log/makeblastdb -e log/makeblastdb -hold_jid $idfqjobids[$libi] $jf`;
      if($retv =~ /Your job (\d+)/){
        push(@jobids, $1);
        print "$1: qsub -o log/makeblastdb -e log/makeblastdb -hold_jid $idfqjobids[$libi] $jf\n";
      }else{
        die "job submission error: qsub $jf said $retv";
      }
    }
  }
  # prepare the target database consisting of all the libraries by writing a single file
  # pseudo dependency to synchronize
  my $jf = "$jobdir/blastdb_all.job";
  {
    open my $jh, ">", $jf or die "cannot open $jf: $!\n";
    print $jh $headstr;
    printf $jh "#\$ %s\n", $rreq;
    print $jh "make -f $makefile $nalf\n";
    close $jh;
  }
  my $nalfjid;
  if($submit){
    my $depends = join(",",@jobids);
    my $retv=`qsub -o log/makeblastdb -e log/makeblastdb -hold_jid $depends $jf`;
    if($retv =~ /Your job (\d+)/){
      $nalfjid = $1;
      print "$1: qsub -o log/makeblastdb -e log/makeblastdb -hold_jid $depends $jf\n";
    }else{
      die "job submission error: qsub $jf said $retv";
    }
  }

  make_path('log/partition');
  my $preq = sprintf " -t 1-%d ", $partition;
  for($libi = 0; $libi < $nlibs; $libi++){
    my $ipjid;
    my $jf = "$jobdir/input_partition_l$libi.job";
    {
      open my $jh, ">", $jf or die "cannot open $jf: $!\n";
      print $jh $headstr;
      printf $jh "#\$ %s\n", $rreq;
      print $jh "make -f $makefile $outdir/c_l${libi}_partition_fa.done\n";
      print $jh "sleep 5s\n"; 
      close $jh;
    }
    if($submit){
      my $retv=`qsub -o log/partition -e log/partition -hold_jid $idfqjobids[$libi] $jf`;
      if($retv =~ /Your job (\d+)/){
        $ipjid = $1;
        print "$1: qsub -o log/partition -e log/partition -hold_jid $idfqjobids[$libi] $jf\n";
      }else{
        die "job submission error: qsub $jf said $retv";
      }
    }

    my $tmp_dfq = "tmp_dfq_l$libi";
    make_path('log/blast');
    #blast individual partitions and convert to sam like format
    $jf = "$jobdir/blast_exec$libi.job";
    {
      open my $jh, ">", $jf or die "cannot open $jf: $!\n";
      print $jh $headstr;
      printf $jh "#\$ %s\n", $param_blast_rreq;
      printf $jh "#\$ %s\n", $preq;
      print $jh $partgen;
      if($direct_vote){
        # This code do not separate the jobs for different valid_voters, which is apparently inefficient.
        # This is justified because for $direct_vote mode, only a single valid_voters are expected.
        # Not reusing the blast output for the multiple valid_voters is by itself inefficient, and
        # parallelization merit for multiple valid_voters is not implemented.
        foreach my $valid_voters(@voters){
          print $jh "make -f $makefile $outdir/$tmp_dfq.v$valid_voters/c_l${libi}_\$PART.v$valid_voters.dfq.gz\n";
        }
      }else{
        print $jh "make -f $makefile $outdir/$tmp_dfq/c_l${libi}_\$PART.blast.sam.gz\n";
      }
      close $jh;
    }
    if($submit){
      my $retv=`qsub -o log/blast -e log/blast -hold_jid $nalfjid,$ipjid $jf`;
      if($retv =~ /Your job-array (\d+)/){
        print "$1: qsub -o log/blast -e log/blast -hold_jid $nalfjid,$ipjid $jf\n";
        $jobids[$libi] = $1;
      }else{
        die "job submission error: qsub $jf said $retv";
      }
    }
  }
  $rreq = $param_ec_rreq;
  my %jid_vv;
  make_path('log/part_ec');
  foreach my $valid_voters(@voters){
    my @t_jid_lib = ();
    for($libi = 0; $libi < $nlibs; $libi++){
      #error correct individual partitions
      my $jf="$jobdir/part_ec$libi.v$valid_voters.job";
      my $tmp_dfq = "tmp_dfq_l$libi";
      my $tmp_finish = "tmp_finish_l$libi.v$valid_voters";
      {
        open my $jh, ">", $jf or die "cannot open $jf: $!\n";
        print $jh $headstr;
        printf $jh "#\$ %s\n", $rreq;
        printf $jh "#\$ %s\n", $preq;
        print $jh $partgen;
        print $jh "make -f $makefile $outdir/$tmp_finish/c_l${libi}.refined_\$PART.fin.v$valid_voters.idfq.gz\n";
        close $jh;
      }
      if($submit){
        my $depend = $jobids[$libi];
        my $retv=`qsub -o log/part_ec -e log/part_ec -hold_jid_ad $depend $jf`;
        if($retv =~ /Your job-array (\d+)/){
          $t_jid_lib[$libi] = $1;
          print "$1: qsub -o log/part_ec -e log/part_ec -hold_jid_ad $depend $jf\n";
        }else{
          die "job submission error: qsub $jf said $retv";
        }
      }
      #finalize ec proc (merge error corrected partitions)
      # consecutive job with array dependency
      $jf="$jobdir/merge_gc$libi.v$valid_voters.job";
      {
        open my $jh, ">", $jf or die "cannot open $jf: $!\n";
        print $jh $headstr;
        printf $jh "#\$ %s\n", $rreq;
        print $jh "make -f $makefile $outdir/c_l${libi}.1.fin.v$valid_voters.idfq.gz\n";
        close $jh;
      }
      if($submit){
        my $depend = $t_jid_lib[$libi];
        my $retv=`qsub -o log -e log -hold_jid $depend $jf`;
        if($retv =~ /Your job (\d+)/){
          $jid_vv{$valid_voters} = $1;
          print "$1: qsub -o log -e log -hold_jid $depend $jf\n";
        }else{
          die "job submission error: qsub $jf said $retv";
        }
      }
    }
  }
  $rreq = $param_gettop_rreq;
  foreach my $estimated_genome_size (@gsizes){
    foreach my $valid_voters (@voters){
      #prepCA (select top 20 x based on estimated genome size)
      my $jf="$jobdir/prepCA.v$valid_voters.e$estimated_genome_size.job";
      my $ps = "$prefix.v$valid_voters.20x$estimated_genome_size";
      {
        open my $jh, ">", $jf or die "cannot open $jf: $!\n";
        print $jh $headstr;
        printf $jh "#\$ %s\n", $rreq;
        print $jh "make -f $makefile $outdir/$ps.frg\n";
        close $jh;
      }
      my $t_jid;
      if($submit){
        my $depend = $jid_vv{$valid_voters};
        my $retv=`qsub -o log -e log -hold_jid $depend $jf`;
        if($retv =~ /Your job (\d+)/){
          $t_jid = $1;
          print "$1: qsub -o log -e log -hold_jid $depend $jf\n";
        }else{
          die "job submission error: qsub $jf said $retv";
        }
      }
      #execCA
      $jf="$jobdir/execCA.v$valid_voters.e$estimated_genome_size.job";
      {
        open my $jh, ">", $jf or die "cannot open $jf: $!\n";
        print $jh $headstr;
        printf $jh "#\$ %s\n", $rreq;
        print $jh "make -f $makefile $outdir/$ps/9-terminator/$ps.scf.fasta\n";
        close $jh;
      }
      if($submit){
        my $depend = $t_jid;
        my $retv=`qsub -o log -e log -hold_jid $depend $jf`;
        if($retv =~ /Your job (\d+)/){
          $t_jid = $1;
          print "$1: qsub -o log -e log -hold_jid $depend $jf\n";
        }else{
          die "job submission error: qsub $jf said $retv";
        }
      }
    }
  }
}

