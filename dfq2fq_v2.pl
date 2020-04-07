#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $VALIDDEPTH=1;
my $VALIDLENGTH=1;
my $fasta=0;
my $finish=0;
my $opt_dfq_check=0;
my $opt_pnp=0;
my $opt_list=0;
my $opt_nlist=0;
my $opt_orig_depth=0;

my $valid_voter=11;
my $trim=42;

my $confident_depth_coefficient=0.0;
my $confident_length_coefficient=0.0;

my $confident_depth=0;
my $confident_length=0;

GetOptions(
  'valid_depth=i' => \$VALIDDEPTH,
  'valid_read_length=i' => \$VALIDLENGTH,
  'f'=>\$fasta,
  'finish'=>\$finish,
  'check'=>\$opt_dfq_check,
  'pnp'=>\$opt_pnp,
  'list'=>\$opt_list,
  'valid_voter=i'=>\$valid_voter,
  'trim=i'=>\$trim,
  'cdc=f'=>\$confident_depth_coefficient,
  'clc=f'=>\$confident_length_coefficient,
  'cd=i'=>\$confident_depth,
  'cl=i'=>\$confident_length,
  'nlist'=>\$opt_nlist,
  'orig_depth'=>\$opt_orig_depth,
);

my $cdc = $confident_depth_coefficient;
my $clc = $confident_length_coefficient;
my $cd = $confident_depth;
my $cl = $confident_length;

if(@ARGV != 1){
  die "USAGE: <this> <in.dfq>\n\t[-f (outputs in fasta)\n\t --valid_depth int\n\t --valid_read_length int\n\t --finish (chops low depth (<valid_depth) regions)\n\t --check (outputs not broken dfq records and discard the rest ('broken' was defined in this code))\n\t --pnp (outputs confidently corrected reads only ('confident' was defined in this code))\n\t --list (outputs confidently corrected read names)\n\t --nlist (outputs NOT confidently corrected read names)]\n";
}
if($opt_pnp && $opt_list){
  die "pnp and list options are incompatible\n";
}
if(($opt_list || $opt_pnp) && $opt_nlist){
  die "nlist option is incompatible with list or pnp options\n";
}
if($cdc && ($cdc < 0.0 or $cdc > 1.0)){
  die "must be 0.0 <= cdc <= 1.0\n";
}
if($clc && ($clc < 0.0 or $clc > 1.0)){
  die "must be 0.0 <= clc <= 1.0\n";
}

if($cd && $cd <1){
  die "must be 0 < cd\n";
}
if($cl && $cl <1){
  die "must be 0 < cl\n";
}
if($cdc && $cd){
  die "cdc and cd are incompatible\n";
}
if($clc && $cl){
  die "clc and cl are incompatible\n";
}

if($VALIDDEPTH <= 0){
  $VALIDDEPTH = 1;
}
if($VALIDLENGTH <= 0){
  $VALIDLENGTH = 1;
}

my $counter=0;
my $printed_line=0;

my $line=<>;
++$counter;
my $result;

while(!eof){
  chomp $line;
  $result = $line =~ s/^\@//;
  if(!$result){
    die "1. strange input\n";
  }
  my $chr = $line;
  ++$printed_line;

  my $consensus = "";
  $line=<>;
  my $line_c = 1;
  chomp $line;
  while(1){# read bases
    $consensus.=$line;
    if(eof){
      die "2. strange input\n";
    }
    $line=<>;
    chomp $line;
    if($line =~ /^\+/){
      last;
    }
    else{
      ++$line_c;
    }
  }
  ++$printed_line;
  ++$printed_line;

  my $depths="";
  my $orig_depths="";
  if(!$fasta){
    chomp $line;
    ($depths,$orig_depths) = (split /\t/, $line)[1..2];
    if(!defined($orig_depths)){
      $orig_depths = "";
    }
  }
  my $qvs = "";
  for(my $i=0; $i<$line_c; ++$i){# # of lines of bases and qvs must be =
    if(eof){
      last;
    }
    $line=<>;
    chomp $line;
    $qvs.=$line;
  }
  ++$printed_line;

  my $strange=0;
  if($consensus =~ /[^acgtnACGTN\-]/){
    $strange = 1;
  }

  if($strange){
    # discard this record.
    printf STDERR ("strange record found\n");
    printf STDERR ("\@%s\n",$chr);
    printf STDERR ("%s\n",$consensus);
    printf STDERR ("+\t%s\t%s\n",$depths,$orig_depths);
    printf STDERR ("%s\n",$qvs);
  }
  elsif(!$fasta){
    if($opt_orig_depth){
      printf("%s\n",$orig_depths);
    }
    elsif($finish){
      &print_bases(\$chr,\$consensus,\$depths,\$qvs);
    }
    elsif($opt_dfq_check){
      if(length($consensus) == length($depths) && length($depths) == length($qvs)){
        printf("\@%s\n",$chr);
        printf("%s\n",$consensus);
        printf("+\t%s\n",$depths);
        printf("%s\n",$qvs);
      }
      else{
        die "broken record\n";
      }
    }
    elsif($opt_pnp || $opt_list || $opt_nlist){
      &flush4pnp(\$chr,\$consensus,\$qvs,\$depths);
    }
    else{
      &flush(\$chr,\$consensus,\$qvs,\$depths);
    }
  }
  else{
    if(length($consensus) >= $VALIDLENGTH){
      # writes in fasta
      $consensus =~ s/\-//g;
      printf(">%s\n",$chr);
      printf("%s\n",$consensus);
    }
  }

  if(!eof){
    $line = <>;# next unit's name
  }
  ++$counter;
}

if(!eof){
  printf STDERR "strange fastq\n";
}
if($printed_line%4 !=0){
  printf STDERR ("strange fastq: the number of lines is $printed_line, $printed_line % 4 = %d\n",$printed_line%4);
}

sub print_bases($$$$){
  my ($chr,$consensus,$depths,$qvs)=@_;
  my $bases="";
  my $qvs_tbp="";
  my $part=0;
  my $loop = length($$depths);
  for(my $i=0; $i<$loop; ++$i){
    if(ord(substr($$depths,$i,1))-33 >= $VALIDDEPTH){
      $bases .= substr($$consensus,$i,1);
      $qvs_tbp .= substr($$qvs,$i,1);
    }
    else{
      &flush($chr,\$bases,\$qvs_tbp,\$part);
      $bases="";
      $qvs_tbp="";
    }
  }
  &flush($chr,\$bases,\$qvs_tbp,\$part);
  $bases="";
  $qvs_tbp="";
}

sub flush($$$$){
  my($chr,$bases,$qvs,$part) = @_;
  if(length($$bases) != length($$qvs)){
    die "FATAL: bases and qvs have different lengths\n";
  }
  #$$bases =~ s/\-//g;
  {
    my $tmp_bases="";
    my $tmp_qvs="";
    my $tmp_depths="";
    my $loop = length($$bases);
    for(my $i=0; $i<$loop; ++$i){
      my $base = substr($$bases,$i,1);
      my $qv = substr($$qvs,$i,1);
      my $depth;
      if(!$finish){
        $depth = substr($$part,$i,1);#depth
      }
      if($base eq '-' or $base eq '~'){
        next;
      }
      $tmp_bases .= $base;
      $tmp_qvs .= $qv;
      if(!$finish){
        $tmp_depths .= $depth;
      }
    }
    $$bases = $tmp_bases;
    $$qvs = $tmp_qvs;
    if(!$finish){
      $$part = $tmp_depths;
    }
  }

  $$bases =~ tr/a-z/A-Z/;
  if(length($$bases)>= $VALIDLENGTH){
    if($finish){
      printf("\@%s/%06d/%d\n",$$chr,$$part++,length($$bases));
      printf("%s\n",$$bases);
      printf("+\n");
      printf("%s\n",$$qvs);
    }
    else{
      printf("\@%s\n",$$chr);
      printf("%s\n",$$bases);
      printf("+\t%s\n",$$part);# part <- depths
      printf("%s\n",$$qvs);
    }
  }
}

sub flush4pnp($$$$){
  my($chr,$bases,$qvs,$depths) = @_;
  if(length($$bases) != length($$qvs)){
    die "FATAL: bases and qvs have different lengths\n";
  }
  
  {
    my $confident_bases=0;
    for(my $i=0; $i<length($$depths); ++$i){
      my $depth = substr($$depths,$i,1);
      my $d_passed = 0;
      if($cdc && ord($depth)-33 >= int($cdc*$valid_voter)){
        $d_passed = 1;
      }
      elsif($cd && ord($depth)-33 >= $cd){
        $d_passed = 1;
      }
      if(!$cdc && !$cd){
        die "specify cdc or cd\n";
      }

      if($d_passed){
        ++$confident_bases;
      }
    }

    my $l_passed=0;
    if($clc && $confident_bases >= int($clc*(length($$bases)-2*$trim))){
      $l_passed = 1;
    }
    elsif($cl && $confident_bases >= $cl){
      $l_passed = 1;
    }
    if(!$clc && !$cl){
      die "specify clc or cl\n";
    }

    if($l_passed){
      if($opt_pnp){
        $$bases =~ tr/a-z/A-Z/;
        printf("\@%s\n",$$chr);
        printf("%s\n",$$bases);
        printf("+\t%s\n",$$depths);
        printf("%s\n",$$qvs);
      }
      else{
        if($opt_list){
          printf("%s\n",$$chr);
        }
      }
    }
    elsif($opt_nlist){
      printf("%s\n",$$chr);
    }
  }
}
