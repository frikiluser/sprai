#!/usr/bin/env perl

use strict;
use Getopt::Long;
# use Data::Dump qw(dump);

my $debug = 0;
my $flag_force = 0;
# Maximum Read Length. We assume that no reads are longer than this (bp).
my $param_max_read_len = 30000;
# We assume that a circular chromosome must yield a contig whose ends
# overlap longer than this length (bp) and more identical than this match ratio.
my $param_min_overlap_len = 1500;
my $param_min_overlap_match_ratio = 0.97;
# When the match ratio of the overlap is lower than this, show a warning. (This aims at not missing potential overlaps)
my $param_alert_overlap_match_ratio = 0.90;
# Circular chromosomes (a contig whose ends overlap) may have a region of
# bad sequence quality at the both ends up to this length.
my $param_max_overlap_hang_len = 100;

GetOptions(
  "force"               => \$flag_force,
  "max_read_len=i"      => \$param_max_read_len,
  "min_overlap_len=i"   => \$param_min_overlap_len,
  "min_overlap_ratio=f" => \$param_min_overlap_match_ratio,
  "alt_overlap_ratio=f" => \$param_alert_overlap_match_ratio,
  "max_hang_len=i"      => \$param_max_overlap_hang_len,
  "debug"               => \$debug
);

my $input_fasta_file_name = shift;
my $temporary_directory_name = shift;

unless(defined $input_fasta_file_name && defined $temporary_directory_name) {
  print STDERR "Usage: check_circurarity.pl <input FASTA (assembly)> <temporary dir>\n";
  exit 0;
}

if(-e $temporary_directory_name) {
  if($flag_force) {
    print STDERR "WARNING: '$temporary_directory_name' already exists, but you gave --force option.\n";
    print STDERR "         We remove the entire directory and start from scratch.\n";
    `rm -rf "$temporary_directory_name"`;
    if($?) {
      print STDERR "ERROR: $temporary_directory_name could not be removed.\n";
      exit 1;
    }
  } else {
    print STDERR "ERROR: '$temporary_directory_name' already exists. We abort for safety.\n";
    print STDERR "       If you wish to remove the directory (and the things in it), and start from scratch,\n";
    print STDERR "       you can give '--force' to do that, but use it carefully.\n";
    exit 1;
  }
}

mkdir $temporary_directory_name or die "ERROR: We could not create '$temporary_directory_name'.";

sub split_fasta_into_separate_fasta
{
  my ($fasta, $dir) = @_;
  my @array_of_sequence_name_and_its_length = ();
  open my $fh, "<", $input_fasta_file_name or die "ERROR: We could not open '$input_fasta_file_name' for input";
  my $ofh;
  my $current_sequence_name = undef;
  while(<$fh>) {
    chomp; chop if(/\r$/);
    if(m|^>(\S+)|) {
      my $sequence_name = $1;
      if(defined $current_sequence_name) { close $ofh; }
      $current_sequence_name = $sequence_name;
      $current_sequence_name =~ s|[^\w\d_ -]||g;
      my $output_file_name = "$dir/$current_sequence_name.fa";
      push(@array_of_sequence_name_and_its_length, {seq_name => $current_sequence_name, len => 0, file_name => $output_file_name});
      open $ofh, ">", $output_file_name or die "ERROR: Could not open '$output_file_name' for output";
      print $ofh ">$current_sequence_name\n";
    } else {
      s|[^ACGT]||ig;
      my $l = length($_);
      if(0 < $l) {
        $array_of_sequence_name_and_its_length[-1]->{len} += $l;
        print $ofh "$_\n";
      }
    }
  }
  if(defined $current_sequence_name) { close $ofh; }
  close $fh;
  return @array_of_sequence_name_and_its_length;
}

sub output_with_folding
{
  my ($fh, $seq) = @_;
  my $l = length($seq);
  my $t = 70;
  for(my $i = 0; $i < $l; ) {
    my $tl = $l - $i; if($t < $tl) { $tl = $t; }
    print $fh substr($seq, $i, $tl), "\n";
    $i += $tl;
  }
}

my @sequence_objects = split_fasta_into_separate_fasta($input_fasta_file_name, $temporary_directory_name);
for my $seq_obj (@sequence_objects) {
  my $command_line = "blastn -task megablast -subject $seq_obj->{file_name} -query $seq_obj->{file_name} -evalue 1e-100 -outfmt 7";
  my @results;
  unless($seq_obj->{len} < $param_max_read_len * 2) {
    my $end_query_pos = $param_max_read_len - 1;
    my $start_subj_pos = $seq_obj->{len} - $param_max_read_len + 1;
    $command_line .= " -subject_loc $start_subj_pos-$seq_obj->{len} -query_loc 1-$end_query_pos";
  }
  print STDERR "\$ $command_line\n" if($debug > 0);
  @results = `$command_line`;
  print @results if($debug > 0);
  my $is_circular = 0;
  my $preserve_pos_from_here_0origin;
  my $preserve_pos_to_here_0origin;
  my $comment;
  for(@results) {
    chomp; chop if(/\r$/);
    next if(m/^#/);
    # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    my ($query_id, $subj_id, $ident_percent, $align_len, $mismatches, $gap_opens, $qstart, $qend, $sstart, $send, $evalue, $bit_score) = split(/\t/);
    next if($qstart == $sstart && $qend == $send);
    next if($sstart >=$send);
    next if($align_len < $param_min_overlap_len);
    next if($ident_percent < $param_alert_overlap_match_ratio * 100.0);
    print STDERR "  len=$seq_obj->{len} [$qstart-$qend] => [$sstart-$send]\n" if($debug);
    if($ident_percent < $param_min_overlap_match_ratio * 100.0) {
      unless($is_circular) {
        $comment = "Possibly circular. [$qstart-$qend] matches [$sstart-$send] with $ident_percent\% identity.";
      }
    } else {
      if($qstart <= $param_max_overlap_hang_len && $seq_obj->{len} - $param_max_overlap_hang_len <= $send) {
        $comment = "[$qstart-$qend] matches [$sstart-$send] with $ident_percent\% identity.";
        $is_circular = 1;
        $preserve_pos_from_here_0origin = $qstart - 1;
        $preserve_pos_to_here_0origin   = $sstart - 1;
      }
    }
  }
  print "$seq_obj->{seq_name}\t" . ($is_circular ? "circular" : "linear") . "\t$comment\n";
  if($is_circular) {
    my $fname = $seq_obj->{file_name} . ".cut.fa";
    my $fname2 = $seq_obj->{file_name} . ".cut_halfrot.fa";
    open my $ifh, "<", $seq_obj->{file_name} or die "ERROR: Could not open '$seq_obj->{file_name}'.";
    open my $ofh, ">", $fname or die "ERROR: Could not open '$fname'";
    open my $ofh2, ">", $fname2 or die "ERROR: Could not open '$fname2'";
    my $header = <$ifh>;
    print $ofh $header;
    print $ofh2 $header;
    my @lines;
    while(<$ifh>) {
      chomp; chop if(/\r$/);
      push(@lines, $_);
    }
    my $seq = join('', @lines);
    my $cut_seq = substr($seq, $preserve_pos_from_here_0origin, $preserve_pos_to_here_0origin - $preserve_pos_from_here_0origin + 1);
    output_with_folding($ofh, $cut_seq);
    {
      my $l = length($cut_seq);
      my $first_half_len = int($l / 2);
      output_with_folding($ofh2, substr($cut_seq, $first_half_len) . substr($cut_seq, 0, $first_half_len));
    }
    close $ifh;
    close $ofh;
    close $ofh2;
  }
}

