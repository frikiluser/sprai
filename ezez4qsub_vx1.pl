#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util;

my $DEBUG=0;
my $DEVEL=0;

my $CURRENT=12;
my $valid_depth=4;
my $valid_read_length=500;
my $scriptdir = "sgescript";
my $data_dir_prefix="data";
my $preprefix="c";
my $ppp="p";

my $confident_depth = $valid_depth;
my $confident_length_coefficient = 0.75;

my $opt_dryrun;
my $opt_ec_only;
my $opt_hgc;

my $now = `date +%Y%m%d_%H%M%S`;
chomp $now;

my %original_time_stamps = ();
my @modified_file_names = ();

GetOptions(
  "n" => \$opt_dryrun,
  "devel" => \$DEVEL,
  "debug"=>\$DEBUG,
  "now=s"=>\$now,
  "ec_only"=>\$opt_ec_only,
  "hgc"=>\$opt_hgc
);

if($DEBUG){
  die "debug mode is not supported yet, sorry.\nDon't use -debug option.\n";
}

my %params;

my @msgs=(
  'USAGE: <this> <ec.spec> <asm.spec>',
  "or: <this> <ec.spec> -ec_only",
  #'[-debug: outputs intermediate files (not implemented)]',
  '[-n: outputs qsub scripts and does NOT qsub]',
  '[-now yyyymmdd_hhmmss: use a XXX_yyyymmdd_hhmmss directories, detect unfinished jobs and restart at the appropriate stage.]',
  "[-ec_only: does error correction and does NOT assemble]",
);

if(@ARGV == 0 || @ARGV > 2){
  my $msg = join "\n\t",@msgs;
  printf STDERR ("%s\n",$msg);
  exit(1);
}
if(@ARGV == 1 && !$opt_ec_only){
  printf STDERR ("WARNING: %s\n", "-ec_only was added");
  $opt_ec_only = 1;
}

my $pwd = `pwd`;
chomp $pwd;

my $asm_spec="";
if(@ARGV == 2){
  $asm_spec = $ARGV[1];
  $asm_spec =~ s/^\s+//;
  if(!-e $asm_spec){
    die "$asm_spec does not exist.\n";
  }
  if($asm_spec =~ /^\//){
    # real path; do nothing
  }
  else{
    $asm_spec = "$pwd/$asm_spec";
  }
}

{
  my $ec_spec = $ARGV[0];
  open my $fh,"<",$ec_spec or die "cannot open $ec_spec :$!\n";
  while($_ = <$fh>){
    next if($_ =~ /^\s+$/);
    next if($_ =~ /^\s*#/);
    chomp;
    my @line = split /\s+/,$_;
  #  if(@line < 2){
  #    die "strange line in $ec_spec\n$_\n";
  #  }
    for(my $i=0; $i<@line; ++$i){
      if($line[$i] =~ /^\s*#/){
        @line = @line[0..$i-1];
        last;
      }
    }
    $params{$line[0]}=join(" ",@line[1..@line-1]);
    if($params{$line[0]}){
  #    printf("%s %s#\n",$line[0],$params{$line[0]});
    }
  }
  close $fh;
}

my $input_for_database;
my $from=0;
my $to=1;
my $pre_partition=1;
my $partition="";
my $word_size=11;
my $evalue=1e-50;
#my $second_evalue="";
my $valid_voters="";
my $trim="";
my $estimated_genome_size="";
my $ca_path="";
my $blast_path="/usr/bin";
my $sprai_path="/usr/lib/sprai";
my $queue_req="";
my $longestXx_queue_req="";
my $blast_rreq="";
my $blasr_path="";
my $blasr_opt="";
my $num_threads=1;
my $max_target_seqs=100;
my $min_len_for_query=1;
my $max_len_for_query=1000000000000000;
my $sge="";
my $estimated_depth=0;
my $use_one_subread=0;

if(defined($params{input_for_database})){
  $input_for_database = $params{input_for_database};
  if(!-e $input_for_database){
    die "$input_for_database does not exist.\n";
  }
}
else{
  die "specify input_for_database in ec.spec\n";
}

if(defined($params{estimated_genome_size})){
  $estimated_genome_size = $params{estimated_genome_size};
}
else{
  die "specify estimated_genome_size in ec.spec\n";
}
if($estimated_genome_size<=0){
  die "estimated_genome_size must be > 0\n";
}
if(defined($params{estimated_depth})){
  $estimated_depth = $params{estimated_depth};
}
else{
  die "specify estimated_depth in ec.spec\n";
}

if(defined($params{from})){
  $from = $params{from};
}
if(defined($params{to})){
  $to = $params{to};
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
#if(defined($params{second_evalue})){
#  $second_evalue = $params{second_evalue};
#}
if(defined($params{valid_voters})){
  $valid_voters = $params{valid_voters};
}
else{
  $valid_voters = int(0.8*($estimated_depth+0.0));
  if($valid_voters > 30){
    $valid_voters =  30;
  }
  #my $n_base = -s $input_for_database;
  #$n_base /= 2;
  #$valid_voters = int(0.8*$n_base/$estimated_genome_size);
  $valid_voters = ($valid_voters < 11) ? 11 : $valid_voters;
}
if(defined($params{trim})){
  $trim = $params{trim};
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
if(defined($params{queue_req})){
  $queue_req = $params{queue_req};
}
if(defined($params{longestXx_queue_req})){
  $longestXx_queue_req = $params{longestXx_queue_req};
}
if(defined($params{BLAST_RREQ})){
  $blast_rreq = $params{BLAST_RREQ};
}
if(defined($params{blasr_path})){
  $blasr_path = $params{blasr_path};
}
if(defined($params{blasr_opt})){
  $blasr_opt = $params{blasr_opt};
}
if(defined($params{min_len_for_query})){
  $min_len_for_query = $params{min_len_for_query};
}
if(defined($params{max_len_for_query})){
  $max_len_for_query = $params{max_len_for_query};
}
if(defined($params{sge})){
  $sge = $params{sge};
}

if(defined($params{max_target_seqs})){
  $max_target_seqs = $params{max_target_seqs};
}

if(defined($params{use_one_subread})){
  $use_one_subread = $params{use_one_subread};
}

printf STDERR ("#>- params -<#\n");
printf STDERR ("input_for_database %s\n",$input_for_database);
printf STDERR ("estimated_genome_size %g\n",$estimated_genome_size);
printf STDERR ("estimated_depth %d\n",$estimated_depth);
#printf STDERR ("from %s\n",$from);
#printf STDERR ("to %s\n",$to);
printf STDERR ("pre_partition %s\n",$pre_partition);
printf STDERR ("partition %s\n",$partition);
if($word_size){
  printf STDERR ("word_size %d\n",$word_size);
}
printf STDERR ("evalue %g\n",$evalue);
printf STDERR ("num_threads %d\n",$num_threads);
printf STDERR ("valid_voters %s\n",$valid_voters);
printf STDERR ("trim %d\n",$trim);
printf STDERR ("use_one_subread %d\n",$use_one_subread);
if($ca_path){
  printf STDERR ("ca_path %s\n",$ca_path);
}
printf STDERR ("blast_path %s\n",$blast_path);
printf STDERR ("sprai_path %s\n",$sprai_path);
if($queue_req){
  printf STDERR ("queue_req %s\n",$queue_req);
}
if($longestXx_queue_req){
  printf STDERR ("longestXx_queue_req %s\n",$longestXx_queue_req);
}
if($blast_rreq){
  printf STDERR ("blast_rreq %s\n",$blast_rreq);
}
if($DEVEL){
  printf STDERR ("development_mode %s\n","true");
}
if($DEBUG){
  printf STDERR ("debug_mode %s\n","true");
}
if($blasr_path){
  printf STDERR ("blasr_path %s\n",$blasr_path);
  if($blasr_opt){
    printf STDERR ("blasr_opt %s\n",$blasr_opt);
  }
}
if($min_len_for_query){
  printf STDERR ("min_len_for_query %d\n",$min_len_for_query);
}
if($max_len_for_query){
  printf STDERR ("max_len_for_query %d\n",$max_len_for_query);
}
if($sge){
  printf STDERR ("sge %s\n",$sge);
}
if($max_target_seqs){
  printf STDERR ("max_target_seqs %d\n",$max_target_seqs);
}
printf STDERR ("#>- params -<#\n");

if(length($blast_path) > 0 && $blast_path !~ /\/$/){
  $blast_path .= "/";
}
if(length($sprai_path) > 0 && $sprai_path !~ /\/$/){
  $sprai_path .= "/";
}

my $outfmt = "7 qseqid sstart send sacc qstart qend bitscore evalue pident qseq sseq";
#my $outfmt = "7 qseqid qlen qstart qend sacc slen sstart send bitscore evalue pident qseq sseq";

my $command="";

my $date="";
my $message="";

my $bindir=$sprai_path;
my $path2blast=$blast_path;

if($DEVEL){
  $bindir = `pwd`;
  chomp $bindir;
}

if(!-d $bindir){
  die "$bindir does not exist\n"
}
if(!-d $path2blast && !-d $blasr_path){
  die "specify aligner dir. $path2blast or $blasr_path does not exist\n"
}
if(!-e "$bindir/bfmt72s"){
  die "$bindir/bfmt72s does not exist in $bindir\ninstall sprai programs in $bindir\n"
}
if(!-e "$path2blast/blastn"){
  die "$path2blast/blastn does not exist in $path2blast\n"
}
if(!$opt_ec_only && !-e $ca_path){
  die "ca_path $ca_path does not exist.\n";
}


$scriptdir = "$pwd/${scriptdir}_$now";
if(!-d $scriptdir){
  `mkdir $scriptdir`;
}

my @pre_array_jobs=();
my @array_jobs=();
my @post_array_jobs=();

my @do_qsub_preaj=();
my @do_qsub_aj=();
my @do_qsub_postaj=();

my $script;
my $outputfile;
my $errfile;

my @list=(0..$partition-1);

my $datadir="$pwd/${data_dir_prefix}_$now";
if(!-d $datadir){
  mkdir $datadir or die "cannot mkdir $datadir: $!\n";
}
my $logdir = sprintf("$pwd/log_$now");
if(!-d $logdir){
  mkdir $logdir or die "cannot mkdir $logdir: $!\n";
}

my $orig_idfq = sprintf("$datadir/$preprefix%02d.idfq.gz",0);

my $db_idfq_gz = sprintf("$datadir/$preprefix%02d.db.idfq.gz",0);
my $qu_idfq_gz = sprintf("$datadir/$preprefix%02d.qu.idfq.gz",0);

my $pipefail = "set -o pipefail;";

if($from == 0)
{
  # fq -> idfq (& id2n)
  my $PG = "$bindir/fq2idfq.pl";
  my $PG2 = "$bindir/fqfilt.pl";
  my $dumbbell_filter = "$bindir/dumbbell_filter.pl";
  if(!$use_one_subread){
    $dumbbell_filter = "cat";
  }

  if($input_for_database !~ /^\//){
    $input_for_database = sprintf("$pwd/$input_for_database");
  }

  my $command;
  my $do_qsub=0;

  my $f_do=1;
  my @parents;
  my $parent = $input_for_database;
  push @parents,$parent;
  my $child=$db_idfq_gz;
  $f_do = &do_or_not(\@parents,\$child);

  my $PG0 = "$bindir/fa2fq.pl";
  $command = sprintf("$pipefail cat $parent | $PG0 - | $dumbbell_filter - | $PG - -flag --prefix $datadir/$preprefix%02d | gzip -c -1 > $child.tmp && mv $child.tmp $child",0);
  if($f_do){
    $do_qsub=1;
  }
  else{
    $command = sprintf("#%s",$command);
  }

  $parent = $db_idfq_gz;
  @parents = ();
  push @parents, $parent;
  $child = $qu_idfq_gz;
  $f_do = &do_or_not(\@parents,\$child);

  my $com2;
  if($min_len_for_query > 1){
    $com2 = sprintf("$pipefail gzip -d -c $parent | $PG2 - $min_len_for_query -max_len $max_len_for_query | gzip -c -1 > $child.tmp && mv $child.tmp $child",0);
  }
  else{
    $com2 = sprintf("ln -s $parent $child");
  }
  if($f_do){
    $do_qsub=1;
  }
  else{
    $com2 = sprintf("#%s",$com2);
  }

#  $command .= " & wait";

  $script = sprintf("$scriptdir/fq2idfq.sh");
  open my $fh, ">", $script or die $!;
  printf $fh ("#!/bin/bash\n");
  printf $fh ("#\$ -S /bin/bash\n");
  printf $fh ("#\$ -cwd\n");
  printf $fh ("#\$ -V\n");
  printf $fh ("#\$ -N $preprefix%02d_fq2idfq_$now\n",0);
  printf $fh ("#\$ -o $logdir\n");
  printf $fh ("#\$ -e $logdir\n");
  printf $fh ("time ($command)\n");
  printf $fh ("time ($com2)\n");
  close $fh;
  push(@pre_array_jobs, $script);
  push(@do_qsub_preaj, $do_qsub);
}

if($pre_partition and $pre_partition < 1){
  printf STDERR ("WARNING: given pre_partition %d was changed to 1\n",$pre_partition);
  $pre_partition = 1;
}

my @makeblastdb_holdjids=();
my @prepartition_holdjids=();

for(my $index=$from; $index<$to; ++$index){
  if($index>0){
    @pre_array_jobs=();
    @array_jobs=();
    @post_array_jobs=();
    @do_qsub_preaj=();
    @do_qsub_aj=();
    @do_qsub_postaj=();
  }

  my $prepadir = sprintf("$datadir/$preprefix%02d",$index);
  if(!-d $prepadir){
    mkdir $prepadir or die "cannot mkdir $prepadir: $!\n";
  }
  my $script_sub_dir = sprintf("$scriptdir/$preprefix%02d",$index);
  if(!-d $script_sub_dir){
    mkdir $script_sub_dir or die "cannot mkdir $script_sub_dir: $!\n";
  }
  # XXX
  my $errdir = sprintf("$logdir/$preprefix%02d/err/",$index);
  my $doutdir = sprintf("$logdir/$preprefix%02d/dout/",$index);
  {
    my $tmp = sprintf("$logdir/$preprefix%02d/",$index);
    if(!-d $tmp){
      mkdir $tmp or die "cannot mkdir $tmp: $!\n";
    }
  }
  if(!-d $errdir){
    mkdir $errdir or die "cannot mkdir $errdir: $!\n";
  }
  if(!-d $doutdir){
    mkdir $doutdir or die "cannot mkdir $doutdir: $!\n";
  }

  my @ec_holdjids=();

  # makeblastdb
  {
    my $PG;
    my $PG2;
    my $do_qsub=0;
    if($blasr_path){
      $PG2= "$bindir/dfq2fq_v2.pl";
      $command = sprintf("gzip -d -c $db_idfq_gz | $PG2 -f - > $datadir/$preprefix%02d.fasta ",$index);
      $command .= sprintf(" && $blasr_path/sawriter $datadir/$preprefix%02d.fasta ",$index);
    }
    else{
      my $f_do=1;
      my $parent = $db_idfq_gz;
      my @parents;
      push @parents,$parent;
      my $child = sprintf("$datadir/makeblastdb_%02d.done",$index);
      $f_do = &do_or_not(\@parents,\$child);

      $PG = "$path2blast/makeblastdb";
      $PG2= "$bindir/dfq2fq_v2.pl";
      $command = sprintf("$pipefail gzip -d -c $parent | $PG2 -f - | $PG -in - -dbtype nucl -out $datadir/$preprefix%02d -title $preprefix%02d 1>$child.tmp && mv $child.tmp $child",$index,$index);
      if($f_do){
        $do_qsub=1;
      }
      else{
        $command = sprintf("#%s",$command);
      }

    }

    $script = sprintf("$script_sub_dir/$preprefix%02d_makeblastdb.sh",$index);
    push(@pre_array_jobs,$script);
    push(@do_qsub_preaj,$do_qsub);
    open my $fh, ">", $script or die $!;
    printf $fh ("#!/bin/bash\n");
    printf $fh ("#\$ -S /bin/bash\n");
    printf $fh ("#\$ -cwd\n");
    printf $fh ("#\$ -V\n");
    my $jobname = sprintf("$preprefix%02d_makeblastdb_$now",$index);
    printf $fh ("#\$ -N $jobname\n");
    push(@ec_holdjids,$jobname);
    if($queue_req){
      printf $fh ("#\$ $queue_req\n");
    }
    printf $fh ("#\$ -o $logdir\n");
    printf $fh ("#\$ -e $logdir\n");
    if($index > $from){
      my $holdlist=join(',',@makeblastdb_holdjids);
      printf $fh ("#\$ -hold_jid $holdlist\n");
    }
    elsif($index == 0){
      printf $fh ("#\$ -hold_jid $preprefix%02d_fq2idfq_$now\n",0);
    }
    else{
      printf STDERR ("strange index. %d\n", $index);
    }
    printf $fh ("time ($command)\n");
    close $fh;
  }

  # prepartitioning
  {
    my $PG1= "$bindir/dfq2fq_v2.pl";
    my $PG2="$bindir/partition_fa.pl -1origin";
    my $do_qsub=0;

    my $f_do=1;
    my $parent = $qu_idfq_gz;
    my @parents;
    push @parents,$parent;
    my $child = sprintf("$prepadir/pre_partition_fa_%02d.done",$index);
    $f_do = &do_or_not(\@parents,\$child);
    $command = sprintf("gzip -d -c $parent | $PG2 -q - %d -p $prepadir/$preprefix%02d 1>$child", $pre_partition,$index);
    #$command = sprintf("gzip -d -c $parent | $PG1 -f - | $PG2 - %d -p $prepadir/$preprefix%02d 1>$child", $pre_partition,$index);
    if($f_do){
      $do_qsub=1;
    }
    else{
      $command = sprintf("#%s",$command);
    }

    $script = sprintf("$script_sub_dir/$preprefix%02d_pre_partition_fa.sh",$index);
    push(@pre_array_jobs,$script);
    push(@do_qsub_preaj,$do_qsub);
    open my $fh, ">", $script or die $!;
    printf $fh ("#!/bin/bash\n");
    printf $fh ("#\$ -S /bin/bash\n");
    printf $fh ("#\$ -cwd\n");
    printf $fh ("#\$ -V\n");
    printf $fh ("#\$ -N $preprefix%02d_pre_partition_fa_$now\n",$index);
    if($queue_req){
      printf $fh ("#\$ $queue_req\n");
    }
    printf $fh ("#\$ -o $logdir\n");
    printf $fh ("#\$ -e $logdir\n");
    if($index > $from){
      my $holdlist=join(',',@prepartition_holdjids);
      printf $fh ("#\$ -hold_jid $holdlist\n");
    }
    elsif($index == 0){
      printf $fh ("#\$ -hold_jid $preprefix%02d_fq2idfq_$now\n",0);
    }
    printf $fh ("time ($command)\n");
    close $fh;
  }


  my @second_cat_pnp_finish_holdjids = ();

  for(my $pp = 1,my $orig_datadir = $datadir; $pp<=$pre_partition; ++$pp,$datadir = $orig_datadir){
    my $ppdir = sprintf("$datadir/$preprefix%02d/$pp/",$index);
    if(!-d $ppdir){
      mkdir $ppdir or die "cannot mkdir $ppdir: $!\n";
    }
    my $maindir= $ppdir;
    if(!-d $maindir){
      mkdir $maindir or die "cannot mkdir $maindir: $!\n";
    }

    # sub partitioning
    {
      my $PG="$bindir/partition_fa.pl";
      my $do_qsub=0;

      my $f_do=1;
      my $parent = sprintf("$prepadir/$preprefix%02d_%04d.fq",$index,$pp);
      #my $parent = sprintf("$prepadir/$preprefix%02d_%04d.fa",$index,$pp);
      my $dummy_parent_1 = sprintf("$prepadir/pre_partition_fa_%02d.done",$index);
      my @parents;
      push @parents,$parent;
      push @parents,$dummy_parent_1;
      my $child = sprintf("$maindir/partition_fa_%02d.done",$index);
      $f_do = &do_or_not(\@parents,\$child);

      $command = sprintf("cat $parent | $PG -q - %d -p $maindir/$preprefix%02d 1>$child",scalar(@list),$index);
      #$command = sprintf("cat $parent | $PG - %d -p $maindir/$preprefix%02d 1>$child",scalar(@list),$index);
      if($f_do){
        $do_qsub=1;
      }
      else{
        $command = sprintf("#%s",$command);
      }
      $script = sprintf("$script_sub_dir/$preprefix%02d$ppp%02d_partition_fa.sh",$index,$pp);
      push(@pre_array_jobs,$script);
      push(@do_qsub_preaj,$do_qsub);
      open my $fh, ">", $script or die $!;
      printf $fh ("#!/bin/bash\n");
      printf $fh ("#\$ -S /bin/bash\n");
      printf $fh ("#\$ -cwd\n");
      printf $fh ("#\$ -V\n");
      #my $jobname = sprintf("$preprefix%02d$ppp%02d_partition_fa_$now",$index,$pp);
      my $jobname = sprintf("$preprefix%02d_partition_fa_$now",$index);
      printf $fh ("#\$ -N $jobname\n");
      #printf $fh ("#\$ -N $preprefix%02d$ppp%02d_partition_fa_$now\n",$index,$pp);
      if($pp==1){
        push(@ec_holdjids,$jobname);
      }
      if($queue_req){
        printf $fh ("#\$ $queue_req\n");
      }
      printf $fh ("#\$ -o $doutdir\n");
      printf $fh ("#\$ -e $errdir\n");
      printf $fh ("#\$ -hold_jid $preprefix%02d_pre_partition_fa_$now\n",$index);
      printf $fh ("time ($command)\n");
      close $fh;
    }
  }
#  print STDERR "pre_array_jobs printed\n";
  
  # array jobs
  my $ppdir = sprintf("$datadir/$preprefix%02d/\${SGE_TASK_ID}/",$index);
  my $maindir= $ppdir;
  my $mdprefix = sprintf("$datadir/$preprefix%02d/",$index);
  my @dqa_offsets=();
  my $dqa_offset=0;
  push @dqa_offsets,$dqa_offset;
  my $task_l=$pre_partition;
  for(my $i=0; $i<@list; ++$i){
    for(my $j=0; $j<$task_l; ++$j){
      push @do_qsub_aj,0;
    }
  }
  {
    for(my $i=0; $i<@list; ++$i){
      my $PG="$path2blast/blastn -dbsize 1 -num_threads \${nthreads}";
      if($word_size){
        $PG .= " -word_size $word_size";
      }
      if($max_target_seqs){
        $PG .= " -max_target_seqs $max_target_seqs";
      }
      #my $PG1="$bindir/overlap_finder";
      my $PG2;
      {
        my $t_trim;
        #if($index == 0){
          $t_trim = $trim;
        #}
        #else{
        #  $t_trim = 0;
        #}
        $PG2=sprintf("$bindir/bfmt72s -c %d -u -i",$t_trim);
      }
      #my $PG3="$bindir/nss2v_v3 -v $valid_voters -q";
      my $PG3="$bindir/nss2v_v3 -v $valid_voters";
      my $PG4="$bindir/myrealigner -f -B $valid_voters -b 3 -d 0.5 -l 131072";
      #my $PG3="$bindir/nss2v_v3 -q";
      #my $PG3="$bindir/nss2v_v3 -q -s";
      #my $PG3="$bindir/nss2v_v3 -v $valid_voters";
      #my $distinguishbale = 0.5;
      #my $PG4="$bindir/realigner_v2 -f -b $valid_depth -d $distinguishbale -q";
      #my $PG4="$bindir/realigner_v2 -f -b $valid_depth -d $distinguishbale -q -g";
      #if($index+1 == $to){
      #  $PG4 = "$PG4 -k";# mask option
      #}
      #my $PG4="$bindir/realigner_v2 -f -B $valid_voters -b 3 -d 0.5";
      my $do_qsub=0;
      my @othercom = ();
      if($DEBUG){# TODO
        die "sorry. debug mode is under construction.\n";
        #$PG4="$bindir/myrealigner -q -B $valid_voters -b 3";
        #my $PG5="$bindir/myrealigner -f -n -q -B $valid_voters -b 3";
        #$command = sprintf("cat $maindir/$preprefix%02d_%04d.fa | $PG -db $maindir/$preprefix%02d -query - -evalue $evalue -outfmt '$outfmt' | tee $datadir/$preprefix%02d_%04d.blastn | $PG2 - | tee $datadir/$preprefix%02d_%04d.nss | $PG3 - | tee $datadir/$preprefix%02d_%04d.vertical | $PG4 - | tee $datadir/$preprefix%02d_%04d.realigned | $PG5 - | gzip -c -1 > $maindir/$outputfile", $index,$list[$i],$index,$index,$list[$i],$index,$list[$i],$index,$list[$i],$index,$list[$i]);
      }
      else{
        $outputfile = sprintf("$maindir/$preprefix%02d_%04d.dfq.gz",$index,$list[$i]);
        if($blasr_path){
          die "sorry. under construction.\nPlease use blastn, not blasr.\n";
#          my $PG6 = "$bindir/m52bfmt7";
#          $command = sprintf("$blasr_path/blasr $maindir/$preprefix%02d_%04d.fa $datadir/$preprefix%02d.fasta -sa $datadir/$preprefix%02d.fasta.sa -m 5 $blasr_opt -out $maindir/tmp%02d_%04d.m5 ",$index,$list[$i],0,0,$index,$list[$i]);
#          $command .= sprintf(" && cat $maindir/tmp%02d_%04d.m5 | $PG6 - |  $PG2 - | $PG3 - | $PG4 - | gzip -c -1 > $outputfile", $index,$list[$i]);
          }
        else{
          my $f_do=1;
          my $input = sprintf("$maindir/$preprefix%02d_%04d.fq",$index,$list[$i]);
          #my $input = sprintf("$maindir/$preprefix%02d_%04d.fa",$index,$list[$i]);
          $command = sprintf("$pipefail cat $input | $bindir/fq2fa.pl - | $PG -db $datadir/$preprefix%02d -query - -evalue $evalue -outfmt '$outfmt' | $PG2 - | $PG3 - | $PG4 - | gzip -c -1 > $outputfile.tmp && mv $outputfile.tmp $outputfile", $index);
          #$command = sprintf("cat $input | $bindir/fq2fa.pl - | $PG -db $datadir/$preprefix%02d -query - -evalue $evalue -outfmt '$outfmt' | $PG2 - -f $input | $PG3 - | $PG4 - | gzip -c -1 > $outputfile.tmp && mv $outputfile.tmp $outputfile", $index);
          push @othercom, "gzip -d -t $outputfile";
          push @othercom, "while [ \$? -ne 0 ]";
          push @othercom, "do";
          push @othercom, "  rm $outputfile";
          push @othercom, "  time ($command)";
          push @othercom, "  gzip -d -t $outputfile";
          push @othercom, "done";
          #push @othercom, "echo test";
          for(my $j=1; $j<=$task_l; ++$j){
            my $mdextended = sprintf("$mdprefix/%d",$j);
            my $parent = sprintf("$mdextended/$preprefix%02d_%04d.fq",$index,$list[$i]);
            my $dummy_parent_1 = sprintf("$datadir/makeblastdb_%02d.done",$index);
            my $dummy_parent_2 = sprintf("$mdextended/partition_fa_%02d.done",$index);
            my @parents=();
            push @parents,$parent;
            push @parents,$dummy_parent_1;
            push @parents,$dummy_parent_2;
            my $child = sprintf("$mdextended/$preprefix%02d_%04d.dfq.gz",$index,$list[$i]);
            $f_do = &do_or_not(\@parents,\$child);
            my $idx = $task_l*$i+($j-1);
#            printf STDERR ("%d\n",$idx);
            if($f_do){
              $do_qsub_aj[$idx] = 1;
            }
            else{
              #$do_qsub_aj[$idx] = 0;
              #$command = sprintf("#%s",$command);
            }
          }
        }
      }

      $script = sprintf("$script_sub_dir/$preprefix%02d_ec_%04d.sh",$index,$list[$i]);
      push(@array_jobs,$script);
      #push(@do_qsub_aj,$do_qsub);
      open my $fh, ">", $script or die $!;

      printf $fh ("#!/bin/bash\n");
      printf $fh ("#\$ -S /bin/bash\n");
      printf $fh ("#\$ -cwd\n");
      printf $fh ("#\$ -V\n");
      #printf $fh ("#\$ -N $preprefix%02d_ec_%04d_$now\n",$index,$list[$i]);
      printf $fh ("#\$ -N $preprefix%02d_ec_$now\n",$index);
      if($queue_req){
        printf $fh ("#\$ $queue_req\n");
      }
      if($blast_rreq){
        printf $fh ("#\$ $blast_rreq\n");
      }
      printf $fh ("#\$ -o $doutdir\n");
      printf $fh ("#\$ -e $errdir\n");

      my $holdlist = join(",",@ec_holdjids);
      printf $fh ("#\$ -hold_jid $holdlist\n");

      if($opt_hgc){
        printf $fh ("hostname=`hostname`\n");
        printf $fh ("nthreads=1\n");
        printf $fh ("if [[ \${hostname} =~ ^ncl ]]; then\n");
        printf $fh ("  nthreads=1\n");
        printf $fh ("else\n");
        printf $fh ("  nthreads=%d\n",$num_threads);
        printf $fh ("fi\n");
      }
      else{
        printf $fh ("nthreads=%d\n",$num_threads);
      }
      printf $fh ("time ($command)\n");
      for(my $i=0; $i<@othercom; ++$i){
        printf $fh ("%s\n",$othercom[$i]);
      }
      close $fh;
    }
  }
  $dqa_offset+=@list*$task_l;
  push @dqa_offsets,$dqa_offset;

  # finish
  {
    my $loop;
    my $diff4dqa;
    {
      my $tmp = (@list+0.0)/($CURRENT+0.0);
=pod
      if($tmp == int($tmp)){
        $loop = int($tmp);
      }
      else{
        $loop = int($tmp+1.0);
      }
      $tmp = &round_up($tmp);
      if($tmp != $loop){
        die "strange $tmp $loop\n";
      }
=cut
      $loop = &round_up($tmp);
    }
    for(my $i=0; $i<$loop; ++$i){
      for(my $j=0; $j<$task_l; ++$j){
        push @do_qsub_aj,0;
      }
    }
    $diff4dqa = $loop*$task_l;
    for(my $k=1; $k<=$task_l; ++$k)
    {
      my @tmp;
      for(my $i=0; $i<@list; ++$i){
        #$tmp[$i]=sprintf("$preprefix%02d_ec_%04d_$now",$index,$list[$i]);
      }
      $tmp[0]=sprintf("$preprefix%02d_ec_$now",$index);
      my $holdlist=join(',',@tmp);

      my $mdextended = sprintf("$mdprefix/%d",$k);
      for(my $i=0; $i<@list; $i+=$CURRENT){
        $outputfile = sprintf("$preprefix%02d_pnp_finish_%04d.idfq.gz",$index,$list[$i]);
        my $child = sprintf("$mdextended/$outputfile");
        $outputfile = sprintf("$maindir/$outputfile");
        my $do_qsub=0;

        my $f_do=1;
        my $parent;
        my @parents;
        for(my $j=0; $j<$CURRENT && $i+$j<@list; ++$j){
          $parent = sprintf("$mdextended/$preprefix%02d_%04d.dfq.gz",$index,$list[$i+$j]);
          push @parents,$parent;
        }
        $f_do = &do_or_not(\@parents,\$child);

        {
          my @ta;
          for(my $j=0; $j<$CURRENT && $i+$j<@list; ++$j){
            my $p = sprintf("$maindir/$preprefix%02d_%04d.dfq.gz",$index,$list[$i+$j]);
            push @ta,$p;
          }
          my $tmp = join " ", @ta;
          if($index+1 == $to){
            $command = sprintf("$pipefail gzip -d -c $tmp | $bindir/dfq2fq_v2.pl - --finish --valid_depth $valid_depth -valid_read_length $valid_read_length | gzip -c -1 > $outputfile.tmp && mv $outputfile.tmp $outputfile");
          }
          else{
            $command = sprintf("$pipefail gzip -d -c $tmp | $bindir/dfq2fq_v2.pl - | gzip -c -1 > $outputfile.tmp && mv $outputfile.tmp $outputfile");
          }
        }

        my $idx = $dqa_offset+&round_up((@list+0.0)/($CURRENT+0.0))*($k-1)+$i/$CURRENT;
#            printf STDERR ("%d\n",$idx);
        if($f_do){
          $do_qsub_aj[$idx] = 1;
        }
        else{
          #$do_qsub_aj[$idx] = 0;
          #$command = sprintf("#%s",$command);
        }
        if($k==1){# write once
          $script = sprintf("$script_sub_dir/$preprefix%02d_pnp_finish_%04d.sh",$index,$list[$i]);
          push(@array_jobs,$script);
#            push(@do_qsub_aj,$do_qsub);
          open my $fh, ">", $script or die $!;

          printf $fh ("#!/bin/bash\n");
          printf $fh ("#\$ -S /bin/bash\n");
          printf $fh ("#\$ -cwd\n");
          printf $fh ("#\$ -V\n");
          #printf $fh ("#\$ -N $preprefix%02d_pnp_finish_%04d_$now\n",$index,$list[$i]);
          printf $fh ("#\$ -N $preprefix%02d_pnp_finish_$now\n",$index);
          if($queue_req){
            printf $fh ("#\$ $queue_req\n");
          }
          printf $fh ("#\$ -o $doutdir\n");
          printf $fh ("#\$ -e $errdir\n");
          printf $fh ("#\$ -hold_jid $holdlist\n");
          printf $fh ("time ($command)\n");
          close $fh;
        }
      }
    }
=pod
    {
      my $loop;
      my $tmp = (@list+0.0)/($CURRENT+0.0);
      if($tmp == int($tmp)){
        $loop = int($tmp);
      }
      else{
        $loop = int($tmp+1.0);
      }
      $dqa_offset += $loop*$task_l;
    }
=cut
    $dqa_offset += $diff4dqa;
    push @dqa_offsets,$dqa_offset;

    for(my $j=0; $j<$task_l; ++$j){
      push @do_qsub_aj,0;
    }
    for(my $k=1; $k<=$task_l; ++$k)
    {
      my $mdextended = sprintf("$mdprefix/%d",$k);
      my @tmp;
      for(my $i=0; $i<@list; $i+=$CURRENT){
        #$tmp[$i/$CURRENT]=sprintf("$preprefix%02d_pnp_finish_%04d_$now",$index,$list[$i]);
      }
      $tmp[0]=sprintf("$preprefix%02d_pnp_finish_$now",$index);
      my $holdlist=join(',',@tmp);

      my $do_qsub=0;
      my $f_do=1;
      my $parent;
      my @parents = ();

      my $files="";
      for(my $i=0; $i<@list; $i+=$CURRENT){
        my $input = sprintf("$preprefix%02d_pnp_finish_%04d.idfq.gz",$index,$list[$i]);
        my $parent = sprintf("$mdextended/$input");
        push @parents,$parent;
        $input = "$maindir/$input";
        $files .= $input;
        $files .= " ";
      }
      chop $files;

      $outputfile = sprintf("$prepadir/$preprefix%02d$ppp\${SGE_TASK_ID}.pnp.fin.idfq.gz",$index+1);
      my $child = sprintf("$prepadir/$preprefix%02d$ppp$k.pnp.fin.idfq.gz",$index+1);
      $f_do = &do_or_not(\@parents,\$child);

      $command = sprintf("cat $files > $outputfile");

      my $idx = $dqa_offset+($k-1);
#      printf STDERR ("pnpfin %s\n",$f_do);
      if($f_do){
        $do_qsub_aj[$idx] = 1;
      }
      else{
        #$do_qsub_aj[$idx] = 0;
        #$command = sprintf("#%s",$command);
      }

      if($k==1){
        $script = sprintf("$script_sub_dir/$preprefix%02d_cat_pnp_finish.sh",$index);
        push(@array_jobs,$script);
#          push(@do_qsub_aj,$do_qsub);
        open my $fh, ">", $script or die $!;

        printf $fh ("#!/bin/bash\n");
        printf $fh ("#\$ -S /bin/bash\n");
        printf $fh ("#\$ -cwd\n");
        printf $fh ("#\$ -V\n");
        my $jobname = sprintf("$preprefix%02d_cat_pnp_finish_$now",$index);
        push(@second_cat_pnp_finish_holdjids,$jobname);
        printf $fh ("#\$ -N $jobname\n");
        if($queue_req){
          printf $fh ("#\$ $queue_req\n");
        }
        printf $fh ("#\$ -o $doutdir\n");
        printf $fh ("#\$ -e $errdir\n");
        printf $fh ("#\$ -hold_jid $holdlist\n");
        printf $fh ("time ($command)\n");

        close $fh;
      }
    }
    $dqa_offset+=$task_l;
    push @dqa_offsets,$dqa_offset;
  }
#  print STDERR "array_jobs printed\n";

  # second cat
  my @get_topXx_holdjids=();
  my $second_catted_file;
  {
    my $holdlist=join(',',@second_cat_pnp_finish_holdjids);

    my $do_qsub = 0;
    my $f_do=1;
    my @parents=();

    my $files="";
    for(my $pp=1; $pp<=$pre_partition; ++$pp){
      my $parent = sprintf("$prepadir/$preprefix%02d$ppp%d.pnp.fin.idfq.gz",$index+1,$pp);
      push @parents,$parent;
      $files .= $parent;
      $files .= " ";
    }

    $outputfile = sprintf("$datadir/$preprefix%02d.fin.idfq.gz",$index+1);
    my $child = $outputfile;

    $f_do = &do_or_not(\@parents,\$child);

    $second_catted_file = "$outputfile";
    $command = sprintf("cat $files > $outputfile.tmp && mv $outputfile.tmp $outputfile");
    if($f_do){
      $do_qsub=1;
    }
    else{
      $command = sprintf("#%s",$command);
    }

    $script = sprintf("$script_sub_dir/$preprefix%02d_second_cat_pnp_finish.sh",$index);
    push(@post_array_jobs,$script);
    push(@do_qsub_postaj,$do_qsub);
    open my $fh, ">", $script or die $!;

    printf $fh ("#!/bin/bash\n");
    printf $fh ("#\$ -S /bin/bash\n");
    printf $fh ("#\$ -cwd\n");
    printf $fh ("#\$ -V\n");
    my $jobname = sprintf("$preprefix%02d_second_cat_pnp_finish_$now",$index);
    push(@get_topXx_holdjids,$jobname);
    printf $fh ("#\$ -N $jobname\n");
    if($queue_req){
      printf $fh ("#\$ $queue_req\n");
    }
    printf $fh ("#\$ -o $logdir\n");
    printf $fh ("#\$ -e $logdir\n");
    printf $fh ("#\$ -hold_jid $holdlist\n");
    printf $fh ("time ($command)\n");

    close $fh;
  }

  my @fastqToCA_holdjids=();
  my @prep_holdjids=@get_topXx_holdjids;
  my $fastqToCA_input;
  my $longestXx;
  #if(!$opt_ec_only)
  if($index+1 == $to)
  {
    my $holdlist=join(',',@get_topXx_holdjids);
    my $f_do=1;
    my $do_qsub=0;
    my @parents=();

    $outputfile = sprintf("$datadir/$preprefix%02d.fin.topXx.fq",$index+1);
    my $child = $outputfile;
    $fastqToCA_input = "$child";
    $longestXx = "$child";

    my $PG1 = "$bindir/get_top_20x_fa.pl";

    my $parent = $second_catted_file;
    push @parents,$parent;
    $f_do = &do_or_not(\@parents,\$child);

    my $uuid = $now;
    #my $uuid = `uuidgen`;
    chomp $uuid;
    $command = sprintf("$pipefail gzip -d -c $second_catted_file > $datadir/$uuid.tmp && $PG1 $datadir/$uuid.tmp -l -g $estimated_genome_size -q -c 20 > $child.tmp && mv $child.tmp $child && rm $datadir/$uuid.tmp");
    #$command = sprintf("gzip -d -c $second_catted_file | $PG1 - -g $estimated_genome_size -q -c 20 > $child.tmp && mv $child.tmp $child");
    if($f_do){
      $do_qsub=1;
    }
    else{
      $command = sprintf("#%s",$command);
    }

    $script = sprintf("$script_sub_dir/$preprefix%02d_get_topXx.sh",$index);
    push(@post_array_jobs,$script);
    push(@do_qsub_postaj,$do_qsub);
    open my $fh, ">", $script or die $!;

    printf $fh ("#!/bin/bash\n");
    printf $fh ("#\$ -S /bin/bash\n");
    printf $fh ("#\$ -cwd\n");
    printf $fh ("#\$ -V\n");
    my $jobname = sprintf("$preprefix%02d_get_topXx_$now",$index);
    push @fastqToCA_holdjids,$jobname;
    push @prep_holdjids,$jobname;
    printf $fh ("#\$ -N $jobname\n");
    if($longestXx_queue_req){
      printf $fh ("#\$ $longestXx_queue_req\n");
    }
    elsif($queue_req){
      printf $fh ("#\$ $queue_req\n");
    }
    printf $fh ("#\$ -o $logdir\n");
    printf $fh ("#\$ -e $logdir\n");
    printf $fh ("#\$ -hold_jid $holdlist\n");
    printf $fh ("time ($command)\n");

    close $fh;
  }

  # prep data
  $qu_idfq_gz = sprintf("$datadir/$preprefix%02d.qu.idfq.gz",$index+1);
  $db_idfq_gz = sprintf("$datadir/$preprefix%02d.db.idfq.gz",$index+1);
  if($index+1 != $to)
  {
=pod
    @makeblastdb_holdjids = ();
    @prepartition_holdjids = ();
    my $jobname = sprintf("$preprefix%02d_data_prep_$now",$index);
    push(@makeblastdb_holdjids,$jobname);
    push(@prepartition_holdjids,$jobname);

    my $f_do=1;
    my $do_qsub=0;
    my @parents;
    my $parent = $second_catted_file;
    #my $parent = $longestXx;
    push @parents,$parent;
    my $child = $qu_idfq_gz;
    $f_do = &do_or_not(\@parents,\$child);

    my $PG = sprintf("$bindir/fq2idfq.pl");
    my $com1 = sprintf("ln -s $parent $child");
    #my $com1 = sprintf("cat $parent | $PG - --prefix $datadir/qu%02d | gzip -c -1 > $child.tmp && mv $child.tmp $child",$index+1);
    if($f_do){
      $do_qsub=1;
    }
    else{
      $com1 = sprintf("#%s",$com1);
    }

    @parents=();
    my $parent_1 = $second_catted_file;
    my $parent_2 = sprintf("$datadir/$preprefix%02d.db.idfq.gz",0);
    push @parents,$parent_1;
    push @parents,$parent_2;
    #push @parents,$longestXx;
    #push @parents,$input_for_database;
    $child = $db_idfq_gz;
    $f_do = &do_or_not(\@parents,\$child);
    my $com2;
    {
      $com2 = "ln -s $parent_2 $child";
      #my $tmp = join " ",@parents;
      #$com2 = sprintf("cat $tmp | $PG - --prefix $datadir/db%02d | gzip -c -1 > $child.tmp && mv $child.tmp $child",$index+1);
    }
    if($f_do){
      $do_qsub=1;
    }
    else{
      $com2 = sprintf("#%s",$com2);
    }

    my $script = sprintf("$script_sub_dir/$preprefix%02d_data_prep.sh",$index);
    push(@post_array_jobs,$script);
    push(@do_qsub_postaj,$do_qsub);
    open my $fh, ">", $script or die $!;

    printf $fh ("#!/bin/bash\n");
    printf $fh ("#\$ -S /bin/bash\n");
    printf $fh ("#\$ -cwd\n");
    printf $fh ("#\$ -V\n");
    printf $fh ("#\$ -N $jobname\n");
    if($queue_req){
      printf $fh ("#\$ $queue_req\n");
    }
    printf $fh ("#\$ -o $datadir\n");
    printf $fh ("#\$ -e $datadir\n");
    my $holdlist=join(',',@prep_holdjids);
    printf $fh ("#\$ -hold_jid $holdlist\n");
    printf $fh ("$com1 & $com2 & wait\n");

    close $fh;
=cut
  }

  my @runCA_holdjids=();
  my $runCA_input;
  if($index+1 == $to)
  {
    my $holdlist=join(',',@fastqToCA_holdjids);
    $outputfile = sprintf("$datadir/$preprefix%02d.fin.topXx.frg",$index+1);
    $runCA_input = "$outputfile";

    my $f_do=1;
    my $do_qsub=0;
    my $parent = $fastqToCA_input;
    my @parents;
    push @parents,$parent;
    my $child = $outputfile;
    $f_do = &do_or_not(\@parents,\$child);

    my $PG1 = "$ca_path/fastqToCA";
    $command = sprintf("$PG1 -libraryname foo -technology pacbio-corrected -reads $parent > $child");

    if($f_do){
      $do_qsub = 1;
    }
    else{
      $command = sprintf("#%s",$command);
    }

    $script = sprintf("$script_sub_dir/$preprefix%02d_fastqToCA.sh",$index);
    if(!$opt_ec_only){
      push(@post_array_jobs,$script);
      push(@do_qsub_postaj,$do_qsub);

      open my $fh, ">", $script or die $!;

      printf $fh ("#!/bin/bash\n");
      printf $fh ("#\$ -S /bin/bash\n");
      printf $fh ("#\$ -cwd\n");
      printf $fh ("#\$ -V\n");
      my $jobname = sprintf("$preprefix%02d_fastqToCA_$now",$index);
      push @runCA_holdjids,$jobname;
      printf $fh ("#\$ -N $jobname\n");
      if($queue_req){
        printf $fh ("#\$ $queue_req\n");
      }
      printf $fh ("#\$ -o $logdir\n");
      printf $fh ("#\$ -e $logdir\n");
      printf $fh ("#\$ -hold_jid $holdlist\n");
      printf $fh ("time ($command)\n");

      close $fh;
    }
  }
  if($index+1 == $to)
  {
    my $holdlist=join(',',@runCA_holdjids);

    my $outdir = sprintf("$datadir/CA_c%02d",$index+1);
    if(-d $outdir){
      $outdir = "${outdir}_$now";
    }
    if(!-d $outdir){
      mkdir "$outdir" or die "cannot mkdir $outdir: $!\n";
    }

    my $PG1 = "$ca_path/runCA";
    my $foo = `date +%Y%m%d_%H%M%S`;
    chomp $foo;

    my $f_do=1;
    my $parent = $runCA_input;
    my @parents;
    push @parents,$parent;
    #my $child = sprintf("$datadir/do_${now}_c%02d.fin.topXx.log",$index+1);
    my $child = sprintf("$outdir/9-terminator/asm_$foo.ctg.fasta");
    $f_do = &do_or_not(\@parents,\$child);
    my $do_qsub=0;
    if($f_do){
      $do_qsub=1;
    }
    else{
    }

    $command = sprintf("$PG1 -dir $outdir -p asm_$foo -s $asm_spec $parent");
    #$command = sprintf("$PG1 -dir $outdir -p asm_$foo -s $asm_spec $parent | tee -a $child");

    $script = sprintf("$script_sub_dir/$preprefix%02d_runCA.sh",$index);
    if(!$opt_ec_only){
      push(@post_array_jobs,$script);
      push(@do_qsub_postaj,$do_qsub);

      open my $fh, ">", $script or die $!;

      printf $fh ("#!/bin/bash\n");
      printf $fh ("#\$ -S /bin/bash\n");
      printf $fh ("#\$ -cwd\n");
      printf $fh ("#\$ -V\n");
      my $jobname = sprintf("$preprefix%02d_runCA_$now",$index);
      printf $fh ("#\$ -N $jobname\n");
      if($queue_req){
        printf $fh ("#\$ $queue_req\n");
      }
      printf $fh ("#\$ -o $logdir\n");
      printf $fh ("#\$ -e $logdir\n");
      printf $fh ("#\$ -hold_jid $holdlist\n");
      printf $fh ("time ($command)\n");

      close $fh;
    }
  }

#  print STDERR "post_array_jobs printed\n";

  {
    my $qsub = `which qsub`;
    chomp $qsub;
    if($sge){
      $qsub = sprintf("%s %s",$qsub,$sge);
    }

    if($index+1 != $to){
      my $script = sprintf("$script_sub_dir/$preprefix%02d_next_qsub.sh",$index);
      push(@post_array_jobs,$script);
      push(@do_qsub_postaj,1);

      open my $fh, ">", $script or die $!;
      printf $fh ("#!/bin/bash\n");
      printf $fh ("#\$ -S /bin/bash\n");
      printf $fh ("#\$ -cwd\n");
      printf $fh ("#\$ -V\n");
      my $jobname = sprintf("$preprefix%02d_next_qsub_$now",$index);
      printf $fh ("#\$ -N $jobname\n");
      printf $fh ("#\$ -o $logdir\n");
      printf $fh ("#\$ -e $logdir\n");

      my $prep = sprintf("$preprefix%02d_data_prep_$now",$index);
      my $holdlist = $prep;
      printf $fh ("#\$ -hold_jid $holdlist\n");

      my $next_script_sub_dir = sprintf("$scriptdir/$preprefix%02d",$index+1);
      my $next_qsub_script = sprintf("$next_script_sub_dir/$preprefix%02d_qsub.sh",$index+1);
      printf $fh ("$qsub $next_qsub_script\n");
      close $fh;
    }


    my $qsub_script = sprintf("$script_sub_dir/$preprefix%02d_qsub.sh",$index);
    open my $fh, ">", $qsub_script or die $!;

    printf $fh ("#!/bin/bash\n");
    printf $fh ("#\$ -S /bin/bash\n");
    printf $fh ("#\$ -cwd\n");
    printf $fh ("#\$ -V\n");
    my $jobname = sprintf("$preprefix%02d_qsub_$now",$index);
    printf $fh ("#\$ -N $jobname\n");
    printf $fh ("#\$ -o $logdir\n");
    printf $fh ("#\$ -e $logdir\n");

    my $task_f=1;
    my $task_l=$pre_partition;
    if($task_l<1){
      die "task_l must be >= 1\n";
    }

#    if($#pre_array_jobs != $#do_qsub_preaj){
#      die "#pre_array_jobs != #do_qsub_preaj\n";
#    }
#    if($#array_jobs != $#do_qsub_aj){
#      die "#array_jobs != #do_qsub_aj\n";
#    }
#    if($#post_array_jobs != $#do_qsub_postaj){
#      die "#post_array_jobs != #do_qsub_postaj\n";
#    }
    for(my $i=0; $i<@pre_array_jobs; ++$i){
      my $tmp = sprintf("$qsub $pre_array_jobs[$i]");
      if($do_qsub_preaj[$i]){
        printf $fh ("$tmp\n");
      }
      else{
        printf $fh ("#$tmp\n");
      }
    }

=pod
    my $foofoo = @list*$task_l;
    for(my $i=0; $i<@do_qsub_aj; ++$i){
      printf STDERR ("%d ",$do_qsub_aj[$i]);
      if($i % $task_l == $task_l-1){
        if($i >= $foofoo){
          printf STDERR (" post\n");
        }
        else{
          printf STDERR (" pre\n");
        }
      }
      if($i % (@list*$task_l) == (@list*$task_l)-1){
        printf STDERR ("\n");
      }
    }
    for(my $i=0; $i<@dqa_offsets; ++$i){
      printf STDERR ("# %d\n",$dqa_offsets[$i]);
    }
=cut
    for(my $i=0; $i<@array_jobs; ++$i){
      my $c=$do_qsub_aj[$task_l*$i];
      $task_f=1;
      my $j;
      for($j=0; $j<$task_l; ++$j){
        if($do_qsub_aj[$task_l*$i+$j] == $c){
        }
        else{
          my $tmp;
          $tmp = sprintf("$qsub -t $task_f-$j:1 $array_jobs[$i]");
          if($c){
            printf $fh ("$tmp\n");
          }
          else{
            printf $fh ("#$tmp\n");
          }
          $c = $do_qsub_aj[$task_l*$i+$j];
          $task_f=$j+1;
        }
      }
      my $tmp;
      $tmp = sprintf("$qsub -t $task_f-$j:1 $array_jobs[$i]");
      if($c){
        printf $fh ("$tmp\n");
      }
      else{
        printf $fh ("#$tmp\n");
      }
    }
    for(my $i=0; $i<@post_array_jobs; ++$i){
      my $tmp = sprintf("$qsub $post_array_jobs[$i]");
      if($do_qsub_postaj[$i]){
        printf $fh ("$tmp\n");
      }
      else{
        printf $fh ("#$tmp\n");
      }
    }
    close $fh;

    if(!$opt_dryrun && $index == $from){
      `qsub $qsub_script`;
    }
  }
}

# runCA
#my $ca_command="ca_ikki_v3.pl ";
#my $ikki_log="$pwd/ca_ikki_v3.log";
#$ca_command .= "-d $datadir -from $from -to $to $asm_spec $estimated_genome_size -ca_path $ca_path > $ikki_log 2>&1";
#if($opt_dryrun){
#  printf("%s\n",$ca_command);
#}
#else{
#  `$ca_command`;
#}

{
  open my $fh, ">>$logdir/$now.tss" or die "cannot open $now.tss: $!\n";

  #foreach my $key (keys %original_time_stamps){
  foreach my $key (@modified_file_names){
    printf $fh ("%s\t%d\n",$key,$original_time_stamps{$key});
  }

  close $fh;
}

sub do_or_not($$){
  my $a = shift;
  my $b = shift;

  my @parents = @$a;
  my $child = $$b;

  my $f_do=0;

#  printf STDERR ("c %s\n",$child);
#  printf STDERR ("p %s\n",$parents[0]);
  if(!-e $child){
    $f_do=1;
  }
  else{
    for(my $i=0; $i<@parents; ++$i){
      if(!-e $parents[$i]){
        $f_do=1;
        my @c_fs=stat $child;
        push @modified_file_names,$child;
        $original_time_stamps{"$child"} = $c_fs[9];
        `touch $child`;
        last;
      }
      else{
        my @p_fs=stat $parents[$i];
        my @c_fs=stat $child;
        my $p_bd = $p_fs[9];
        my $c_bd = $c_fs[9];
        if($p_bd > $c_bd){
          $f_do=1;
          push @modified_file_names,$child;
          $original_time_stamps{"$child"} = $c_fs[9];
          `touch $child`;
          last;
        }
        else{
          #$f_do=0;
        }
      }
    }
=pod
    if(!$f_do && $child =~ /\.dfq\.gz$/){
      my $ret = system("gzip -d -c -t $child 2> /dev/null");
      if($ret){
        printf STDERR "broken $child\n";
        $f_do=1;
        `touch $child`;
      }
    }
=cut
  }
  return $f_do;
}

sub round_up($){
  my $val = shift;
  my $tmp = Scalar::Util::looks_like_number($val);
  if(!$tmp){
    die "$val does not look line number\n";
  }
  my $ret;
  if($val == int($val)){
    $ret = int($val);
  }
  else{
    $ret = int($val+1.0);
  }
  return $ret;
}
