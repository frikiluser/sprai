#!/usr/bin/env perl

use strict;
use Getopt::Long;
# use Data::Dump qw(dump);

my $debug = 0;
my $flag_force = 0;
# Minimum match ratio to consider
my $param_min_ratio = 0.94;
my $param_min_len_ratio = 0.98;
# This length may not have a good match to a larger contig.
my $param_max_hang_len = 100;

GetOptions(
  "force"             => \$flag_force,
  "max_hang_len=i"    => \$param_max_hang_len,
  "min_align_ratio=f" => \$param_min_ratio,
  "debug"             => \$debug
);

my $input_fasta_file_name = shift;

unless(defined $input_fasta_file_name) {
  print STDERR "Usage: check_redundancy.pl <input FASTA (assembly)>\n";
  exit 0;
}

my $command_line = "blastn -task blastn -subject $input_fasta_file_name -query $input_fasta_file_name -evalue 1e-100 -outfmt 7";
my @results;
print STDERR "\$ $command_line\n" if($debug > 0);
@results = `$command_line`;
print @results if($debug > 0);

my %id_2_len;
my $current_seq;
my @ids;
for(@results) {
  chomp; chop if(/\r$/);
  next if(m/^#/);
  my ($query_id, $subj_id, $ident_percent, $align_len, $mismatches, $gap_opens, $qstart, $qend, $sstart, $send, $evalue, $bit_score) = split(/\t/);
  if($query_id eq $subj_id && $current_seq ne $query_id) {
    $current_seq = $query_id;
    $id_2_len{$current_seq} = $align_len;
    push(@ids, $current_seq);
  }
}

my %id_2_redundant_arr;
for(@results) {
  chomp; chop if(/\r$/);
  next if(m/^#/);
  # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
  my ($query_id, $subj_id, $ident_percent, $align_len, $mismatches, $gap_opens, $qstart, $qend, $sstart, $send, $evalue, $bit_score) = split(/\t/);
  next if($query_id eq $subj_id);
  next if($ident_percent < $param_min_ratio * 100.0);
  next if($align_len < $id_2_len{$query_id} * $param_min_len_ratio);
  push(@{$id_2_redundant_arr{$query_id}}, "\[$qstart-$qend\]($id_2_len{$query_id})=($ident_percent\%)=>$subj_id\[$sstart-$send\]($id_2_len{$subj_id})");
}

print "#ID\tunique\tcomment\n";
for my $id (@ids) {
  print "$id\t";
  if(exists $id_2_redundant_arr{$id}) {
    print "redundant\t";
    print join(' ', @{$id_2_redundant_arr{$id}});
  } else {
    print "unique\t";
  }
  print "\n";
}

# dump(%id_2_redundant_arr);

