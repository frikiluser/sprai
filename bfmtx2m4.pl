#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $opt_shortname;
#GetOptions('n'=>\$opt_shortname);

if(@ARGV != 2){
  die "USAGE: <this> <m4.pre> <all_norm.fa>\n";
}

my $m4pre = $ARGV[0];
my $all_norm = $ARGV[1];

#print $m4pre,"\n";
#print $all_norm,"\n";
#exit;

my %name2len=();

open my $fa_fh, "<", $all_norm or die "cannot open $all_norm: $!\n";

my $counter=0;

my $name = <$fa_fh>;
chomp $name;
$name =~ s/^>//;

my $bases = "";
my $qval = "";
while(1){
  while(my $buf=<$fa_fh>){
    chomp $buf;
    if($buf =~ /^>/){
      $name2len{$name} = length($bases);

      $name = $buf;
      $bases= "";
      $qval = "";
      $name =~ s/^>//;
      last;
    }
    else{
      $bases .= $buf;
    }
  }
  if(eof){
    last;
  }
}

$name2len{$name} = length($bases);
close $fa_fh;

open my $pre_fh, "<", $m4pre or die "cannot open $m4pre: $!\n";
while(my $line=<$pre_fh>){
  chomp $line;
  if($line =~ /^#/){
    next;
  }
  # qseqid sacc bitscore pident qstart qend sstrand sstart send
  my @tmp = split /\s+/,$line;
  if(@tmp != 9){
    die "strange format: $line\n";
  }
  my ($qseqid,$sacc,$bitscore,$pident,$qstart,$qend,$sstrand,$sstart,$send) = @tmp;
  $bitscore = (-1.0)*$bitscore*5.0/2.0;
  # qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart tend tseqlength mapqv ncells clusterScore probscore numSigClusters
  # qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart tend tseqlength mapqv
  if($sstrand eq "plus"){
    $sstrand = 0;
  }
  elsif($sstrand eq "minus"){
    $sstrand = 1;
  }
  else{
    die "souteigai: $sstrand\n";
  }
  printf("%s %s %d %f %s %d %d %d %s %d %d %d %d\n",$qseqid, $sacc, $bitscore, $pident, "0", $qstart-1, $qend, $name2len{$qseqid}, $sstrand, $sstart-1, $send, $name2len{$sacc}, 254);
}

close $pre_fh;

