#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $DEVEL;

my $from=0;
my $to=1;
my $fastqdir="./";
my $ca_path="";
my $out_dir="CA";
#my $tmp_dir=$out_dir;
my $sprai_path="";
my $coverage=20;
my $raw_fastq="";

my @msg=(
"USAGE: <this> <asm.spec> estimated_genome_size",
#"[-from integer]",
#"[-to integer ]",
"[-d directory in which fin.idfq.gzs exist (default: $fastqdir)]",
"[-ca_path /path/to/your/wgs/Linux-amd64/bin (default: $ca_path)]",
#"[-tmp_dir temporary directory (default: $tmp_dir)]",
"[-out_dir output directory (default: $out_dir)]",
"[-sprai_path the path to get_top_20x_fa.pl installed (default: $sprai_path)]",
"[-coverage int : use longer than N(coverage) reads for assembly (default: $coverage)]",
"",
"[-raw_fastq in.fq : use all reads in in.fq (default: off)]",
);

GetOptions(
  'from=i' => \$from,
  'to=i' => \$to,
  'd=s' => \$fastqdir,
  'devel' => \$DEVEL,
  'ca_path=s'=>\$ca_path,
#  'tmp_dir=s'=>\$tmp_dir,
  'out_dir=s'=>\$out_dir,
  'coverage=i'=>\$coverage,
  'raw_fastq=s'=>\$raw_fastq,
  'sprai_path=s'=>\$sprai_path
  );

if(@ARGV != 2){
  my $tmp = join "\n\t",@msg;
  die "$tmp\n";
  #die "USAGE: <this> <asm.spec> estimated_genome_size [-from integer -to integer]\n\t[-d directory in which fin.fq.gzs exist]\n\t[-ca_path /path/to/your/wgs/Linux-amd64/bin]\n";
}

my $spec = $ARGV[0];
printf STDERR ("%s is given\n",$spec);
my $estimated_genome_size = $ARGV[1];
if($estimated_genome_size <= 0){
  die "estimated_genome_size must be > 0\n";
}

printf STDERR ("#>- params -<#\n");
printf STDERR ("spec\t%s\n",$spec);
printf STDERR ("estimated_genome_size\t%s\n",$estimated_genome_size);
#printf STDERR ("from\t%s\n",$from);
#printf STDERR ("to\t%s\n",$to);
printf STDERR ("fastq_dir\t%s\n",$fastqdir);
#printf STDERR ("tmp_dir\t%s\n",$tmp_dir);
printf STDERR ("out_dir\t%s\n",$out_dir);
printf STDERR ("sprai_path\t%s\n",$sprai_path);

if($DEVEL){
  printf STDERR ("development mode\t%s\n","true");
}
if($ca_path){
  printf STDERR ("ca_path\t%s\n",$ca_path);
}
if($coverage){
  printf STDERR ("coverage\t%s\n",$coverage);
}
if($raw_fastq){
  printf STDERR ("raw_fastq\t%s\n",$raw_fastq);
}
printf STDERR ("#>- params -<#\n");
#exit;

my $PWD=`pwd`;
chomp $PWD;

my $now = `date +%Y%m%d_%H%M%S`;
chomp $now;

if($out_dir !~ /^\//){
  $out_dir = "$PWD/$out_dir";
}

if($out_dir =~ /\/$/){
  chop $out_dir;
}
my @out_dirs=();
my @p2=();
for(my $i=$from+1; $i<=$to; ++$i){
  my $tmp;
  my $now_used=0;
  if($to-$from > 1){
    $tmp = sprintf("%s_%02d_%s",$out_dir,$i,$now);
  }
  else{
    $tmp = sprintf("%s",$out_dir);
  }
  if(-d $tmp){
    my $now = `date +%Y%m%d_%H%M%S`;
    chomp $now;
    $tmp = sprintf("%s_%02d_%s",$out_dir,$i,$now);
    if($to-$from>1){
      redo;
    }
    $now_used=1;
  }
  mkdir "$tmp" or die "cannot mkdir $tmp: $!\n";
  $out_dirs[$i] = $tmp;
  my @foo = split /\//,$out_dirs[$i];
  if($to-$from>1 || $now_used){
    $p2[$i] = $foo[$#foo];
  }
  else{
    my $now = `date +%Y%m%d_%H%M%S`;
    chomp $now;
    $p2[$i] = sprintf("%s_%02d_%s",$foo[$#foo],$from+1,$now);
  }
}

=pod
if($tmp_dir !~ /^\//){
  $tmp_dir = "$PWD/$tmp_dir";
}
if(-d "$tmp_dir"){
  my $tmp = "${tmp_dir}_$now";
  mkdir "$tmp" or die "cannot mkdir $tmp: $!\n";
}
else{
  my $tmp = "${tmp_dir}";
  mkdir "$tmp" or die "cannot mkdir $tmp: $!\n";
}
=cut

#$tmp_dir = $out_dir;

my $bashcommand="";

my $suffix = "top20x";


if(!$raw_fastq){
  #printf STDERR ("start idfq2fq\n");
  for(my $i=$from+1; $i<=$to; ++$i){
    my $prefix=sprintf("c%02d.fin",$i);
    my $PG1 = "get_top_20x_fa.pl";
    if($sprai_path){
      $PG1 = "$sprai_path/$PG1";
    }
    my $uuid = $now;
    #my $uuid = `uuidgen`;
    chomp $uuid;
    $bashcommand .= "gzip -d -c $fastqdir/$prefix.idfq.gz > $PWD/$uuid.tmp && $PG1 $PWD/$uuid.tmp -l -c $coverage -g $estimated_genome_size -q > $out_dirs[$i]/$prefix.$suffix.fq && rm $PWD/$uuid.tmp &\n";
  }

  `
    $bashcommand
    wait
  `;
  #printf STDERR ("done idfq2fq\n");
}
else{
  if($to-$from != 1){
    printf STDERR ("strange 'from' and 'to'.\n");
    printf STDERR ("from - to must be one if you specify -raw_fastq.\n");
    exit(1);
  }
  if($raw_fastq !~ /^\//){
    $raw_fastq = "$PWD/$raw_fastq";
  }
  my $prefix=sprintf("c%02d.fin",$from+1);
  `ln -s $raw_fastq $out_dirs[$from+1]/$prefix.$suffix.fq`;
}

#printf STDERR ("start fastqToCA\n");
for(my $i=$from+1; $i<=$to; ++$i){
  my $prefix=sprintf("c%02d.fin",$i);
  my $fastqToCA = "fastqToCA";
  if($ca_path){
    $fastqToCA = "$ca_path/$fastqToCA";
  }
  `$fastqToCA -libraryname foo -technology pacbio-corrected -reads $out_dirs[$i]/$prefix.$suffix.fq > $out_dirs[$i]/$prefix.$suffix.frg`;
}
#printf STDERR ("done fastqToCA\n");

#printf STDERR ("start CA (stopAfter=unitigger)\n");

my $now_used=0;

my $runCA="runCA";
if($ca_path){
  $runCA = "$ca_path/$runCA";
}

=pod
for(my $i=$from+1; $i<=$to; ++$i){
  my $prefix=sprintf("c%02d.fin",$i);
  `$runCA stopAfter=unitigger -dir $out_dirs[$i] -p asm_$p2[$i] -s $spec $out_dirs[$i]/$prefix.$suffix.frg 2>> $out_dirs[$i]/do_$prefix.$suffix.$now.log`;
}
#printf STDERR ("done CA (stopAfter=unitigger)\n");

#printf STDERR ("start CA (the rest)\n");
$bashcommand="";
for(my $i=$from+1; $i<=$to; ++$i){
  my $prefix=sprintf("c%02d.fin",$i);
  $bashcommand .= "($runCA -dir $out_dirs[$i] -p asm_$p2[$i] -s $spec $out_dirs[$i]/$prefix.$suffix.frg 2>> $out_dirs[$i]/do_$prefix.$suffix.$now.log) &\n";
}
`
  $bashcommand
`;
=cut
$bashcommand="";
for(my $i=$from+1; $i<=$to; ++$i){
  my $prefix=sprintf("c%02d.fin",$i);
  $bashcommand .= "$runCA -dir $out_dirs[$i] -p asm_$p2[$i] -s $spec $out_dirs[$i]/$prefix.$suffix.frg 2>> $out_dirs[$i]/do_$prefix.$suffix.$now.log &\n";
}
`
  $bashcommand
  wait
`;
