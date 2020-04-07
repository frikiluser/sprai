#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $opt_largerthan;
my $opt_shortname=1;
my $prefix=`date +%Y%m%d%H%M%S`;
chomp $prefix;
my $dbsuffix = ".id2n";
my $output_is_fasta=0;
my $opt_flag_fq=0;

my $valid_read_length=1;

GetOptions(
  'l' => \$opt_largerthan,
  'prefix=s' => \$prefix,
  'valid_read_length=i'=>\$valid_read_length,
  'flag'=>\$opt_flag_fq,
  'output_is_fasta'=>\$output_is_fasta
);

my @msg = (
  "USAGE: <this> <in.fq> > out.idfq",
  "[--prefix prefix]",
  "[-l : '\@'->'>']",
  "[-output_is_fasta]",
  "[-flag : discard quality values and use quality lines for iterative error correction]",
);

if(@ARGV != 1){
  my $msg = join "\n\t",@msg;
  die "$msg\n";
}

my $id_head_character="\@";

if($opt_largerthan){
  $id_head_character=">";
}
else{
}

my $fh;
if($opt_shortname){
  open $fh, ">", $prefix.$dbsuffix or die "cannot open $prefix.$dbsuffix:$!\n";
}

my $counter=0;
my $printed_line=0;

my $line = <>;
++$counter;
my $result;
while(!eof){
  chomp $line;
  $result = $line =~ s/^$id_head_character/>/;
  if(!$result){
    if($id_head_character eq "\@"){
      $id_head_character = ">";
      redo;
    }
    elsif($id_head_character eq ">"){
      die "1. strange input $result\n$line\n";
    }
  }
  my($idfqline,$nameline,$baseline,$optionalline,$qvline);
  if($opt_shortname){
    #printf $fh ("%d\t%s\n",$counter,substr($line,1));
    $idfqline = sprintf("%d\t%s",$counter,substr($line,1));
    if($output_is_fasta){
      $line = sprintf(">%d",$counter);
    }
    else{
      $line = sprintf("\@%d",$counter);
    }
  }
  $nameline = $line;
  #print $line,"\n";# name
  #++$printed_line;
  
  my $bases="";
  if(eof){
    last;
  }
  $line =<>;
  my $line_c = 1;
  chomp $line;
  while(1){# read bases
    $bases .= $line;
    if(eof){
      last;
    }
    $line = <>;
    chomp $line;
    if($line =~ /^\+/){
      chomp $line;
      $optionalline = $line;
      last;
    }
    else{
      ++$line_c;
    }
  }
  $baseline = $bases;

  if(eof){
    last;
  }
  #print "+\n";
  #++$printed_line;
  my $qvs="";
  for(my $i = 0; $i<$line_c; ++$i){# # of lines of bases and qvs must be =
    if(eof){
      last;
    }
    $line = <>;
    chomp $line;
    $qvs .= $line;
  }
  $qvline = $qvs;
  if(length($baseline) != length($qvline)){
    printf STDERR "strange input: # of bases and qvs differs\n$baseline\n$qvline\nthis record was skipped.\n";
  }
  #print $qvs,"\n";
  #++$printed_line;

  elsif(length($bases)>=$valid_read_length){
    print $fh $idfqline,"\n";
    print $nameline,"\n";
    ++$printed_line;
    print $baseline,"\n";
    ++$printed_line;
    if(!$output_is_fasta){
      print $optionalline,"\n";
      ++$printed_line;
      if($opt_flag_fq){
        my $len = length($qvline);
        $qvline = "";
        for(my $i=0; $i<$len; ++$i){
          $qvline .= "2";
        }
      }
      print $qvline,"\n";
      ++$printed_line;
    }
  }

  if(!eof){
    $line = <>;# next unit's name
  }
  ++$counter;
}

if($opt_shortname){
  close $fh;
}
#if($printed_line % 4 != 0){
#  printf STDERR "WARNING: the input fq file may be broken\n";
#}
