#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $pacbio;

GetOptions('p'=>\$pacbio);

if(@ARGV != 1){
  die "USAGE: <this> <in.fq>\n\t[-p: against '\@'->'>' input]\n";
}

my $id_head_character="\@";

if($pacbio){
  $id_head_character=">";
}
else{
}

my $line = <>;
while(!eof){
  chomp $line;
  my $result = $line =~ s/^$id_head_character/>/;
  if(!$result){
    if($id_head_character eq "\@"){
      $id_head_character = ">";
      redo;
    }
    elsif($id_head_character eq ">"){
      die "1. strange input $result\n$line\n";
    }
  }

  #$line =~ s/^\@/>/;
  print $line,"\n";
  
  my $bases="";
  $line =<>;
  while($line !~ /^\+/){
    chomp $line;
    $bases .= $line;
    $line = <>;
  }
  print $bases,"\n";

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
}
