#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $pacbio;
my $min_len = 1;
my $max_len = 1000000000000;

GetOptions('p'=>\$pacbio,"max_len=i"=>\$max_len);

if(@ARGV != 2){
  die "USAGE: <this> <in.fq> min_len\n";
}

my $id_head_character="\@";

if($pacbio){
  $id_head_character=">";
}
else{
}

$min_len = $ARGV[1];

my $line = <>;
while(!eof){
  chomp $line;
  my $name = $line;
  #print $line,"\n";
  
  my $bases="";
  $line =<>;
  while($line !~ /^\+/){
    chomp $line;
    $bases .= $line;
    $line = <>;
  }
  #print $bases,"\n";

  my $qvs="";
  $line =<>;# qvs
  while($line !~ /^$id_head_character/ || length($qvs) < length($bases)){
    chomp $line;
    $qvs.=$line;
    # do nothing
    if(eof){
      last;
    }
    $line = <>;
  }
  if(length($bases)>=$min_len && length($bases)<=$max_len){
    print $name,"\n";
    print $bases,"\n";
    print "+\n";
    print $qvs,"\n";
  }
}
