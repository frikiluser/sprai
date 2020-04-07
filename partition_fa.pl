#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
#use Compress::Zlib;

#my $opt_n;
#GetOptions("n" => \$opt_n);

my $opt_prefix=`date +%Y%m%d%H%M%S`;
chomp $opt_prefix;
my $opt_gz;
my $opt_1origin=0;
my $opt_fq=0;

#print $opt_prefix,"\n";
#exit 1;
GetOptions(
  "p=s" => \$opt_prefix,
  "q" => \$opt_fq,
  "1origin"=>\$opt_1origin
);

my $ext="fa";
if($opt_fq){
  $ext = "fq";
}

#GetOptions("p=s" => \$opt_prefix, "z" => \$opt_gz);
$opt_gz=0;

if(@ARGV != 2){
  printf STDERR ("USAGE: <this> <in.fa> <int (number of partitions)>\n");
  printf STDERR ("\t[-p <prefix>: prefix to the outputs. default is yymmddhhmmss]\n");
  printf STDERR ("\t[-1origin: for file naming\n");
  printf STDERR ("\t[-q : regard input as in fastq]\n");
#  printf STDERR ("\t[-z: output gz compressed file instead of uncompressed fasta file]\n");
  #printf STDERR ("\t[-n (name: output only names , not sequences)\n");
  exit 1;
}

#printf("%s\n", $ARGV[0]);
#printf("%s\n", $ARGV[1]);
#exit 1;

my $in_fa = $ARGV[0];
chomp $in_fa;
#my $IN_FA;
#if($in_fa eq '-'){
#  open $IN_FA,"<",STDIN or die "cannot open $in_fa : $!\n";
#}
#else{
  open my $IN_FA,"<".$in_fa or die "cannot open $in_fa : $!\n";
#}

my $partitions = $ARGV[1];
if($partitions !~ /^\+?\d+$/){
  printf STDERR ("nth must be an integer\n");
  exit 1;
}

my @out;
for(my $i=0; $i<$partitions; ++$i){
  if($opt_gz){
#    $out[$i] = gzopen(sprintf("%s_%04d.fa.gz",$opt_prefix,$i),"wb") or die "cannot open output file:$!\n";
  }
  else{
    if($opt_1origin){
      open $out[$i],">",sprintf("%s_%04d.$ext",$opt_prefix,$i+1) or die "cannot open output file:$!\n";
    }
    else{
      open $out[$i],">",sprintf("%s_%04d.$ext",$opt_prefix,$i) or die "cannot open output file:$!\n";
    }
  }
}

my $counter=-1;

my $name = <$IN_FA>;

while(!eof($IN_FA)){
  ++$counter;
  chomp $name;
  my $tmp_name = $name;
  my $bases="";
  my $opts="";
  my $qvs="";
  if(!$opt_fq){
    while(!eof($IN_FA) && (($name = <$IN_FA>) !~ /^>/)){
      chomp $name;#bases
      $bases .= $name;
    }
  }
  else{
    while(!eof($IN_FA) && (($name = <$IN_FA>) !~ /^\+/)){
      chomp $name;#bases
      $bases .= $name;
    }
    chomp $name;
    $opts = $name;
    while(!eof($IN_FA) && (($name = <$IN_FA>) !~ /^\@/)){
      chomp $name;#qvs
      $qvs .= $name;
    }
    if(length($bases) != length($qvs)){
      printf STDERR ("|bases| != |qvs|: %d %d\n",length($bases),length($qvs));
      exit(1);
    }
  }
  #next if ($counter % $mod != $nth);
  my $fh = $out[$counter % $partitions];
  if($opt_gz){
#    $fh->gzwrite(sprintf("%s\n",$tmp_name));
#    $fh->gzwrite(sprintf("%s\n",$bases));
  }
  else{
    printf $fh ("%s\n",$tmp_name);
    printf $fh ("%s\n",$bases);
    if($opt_fq){
      printf $fh ("%s\n",$opts);
      printf $fh ("%s\n",$qvs);
    }
  }

  #printf $out[$counter % $partitions] ("%s\n",$tmp_name);
  #printf $out[$counter % $partitions] ("%s\n",$bases);
  #my $num=70;
  #my $loop = length($bases);
  #for(my $stt=0; $stt < $loop; $stt += $num){
  #  printf("%s\n",substr($bases,$stt,$num));
  #}
}

close $IN_FA;
for(my $i=0; $i<$partitions; ++$i){
  if($opt_gz){
#    $out[$i]->gzclose;
  }
  else{
    close $out[$i];
  }
}

