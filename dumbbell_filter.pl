#!/usr/bin/perl
use strict;
use warnings;

my @msgs=(
#"USAGE: <this> <in.fasta> [<in2.fasta> ...]"
"USAGE: <this> <in.fq> [<in2.fq> ...]"
);

my $prev_name_prefix="";

=pod
my $name = <>;
chomp $name;
$name =~ s/^>//;
{
  my @tmp = split /\//,$name;
  $prev_name_prefix = join "/",@tmp[0..$#tmp-1];
}
my $bases = "";
=cut


my @names=();
my @reads=();
my @opts=();
my @qvs=();

=pod
while(1){
  while(my $buf=<>){
    chomp $buf;
    if($buf =~ /^>/){
      push @names,$name;
      push @reads,$bases;

      $name = $buf;
      $name =~ s/^>//;
      {
        my @tmp = split /\//,$name;
        my $p = join "/",@tmp[0..$#tmp-1];
        if($prev_name_prefix ne $p){
          &flush(\@names,\@reads);
          @names=();
          @reads=();
          $prev_name_prefix = $p;
        }
      }

      $bases= "";
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
push @names,$name;
push @reads,$bases;
&flush(\@names,\@reads);

sub flush($$){
  my $n = shift;
  my $r = shift;
  my @names = @$n;
  my @reads = @$r;
  my $longest = length($reads[0]);
  my $i_longest=0;
  for(my $i=1; $i<@reads; ++$i){
    if(length($reads[$i]) > $longest){
      $longest = length($reads[$i]);
      $i_longest = $i;
    }
  }
  printf(">%s\n",$names[$i_longest]);
  printf("%s\n",$reads[$i_longest]);
}
=cut

my $id_head_character = "\@";

my $name;
my $bases;
my $opt;
my $qv;

my $line = <>;
while(!eof){
  chomp $line;
  $line =~ s/^\@//;
  $name = $line;
  
  $bases="";
  $line =<>;
  while($line !~ /^\+/){
    chomp $line;
    $bases .= $line;
    $line = <>;
  }
  chomp $line;
  $opt = $line;

  $qv="";
  $line =<>;# qv
  while($line !~ /^$id_head_character/ || length($qv) < length($bases)){
    chomp $line;
    $qv .= $line;
    if(eof){
      last;
    }
    $line = <>;
  }
  {
    # name must be in @xxxx/yyyy/z1_z2 format
    my @tmp = split /\//,$name;
    my $p = join "/",@tmp[0..$#tmp-1];
    if(!$prev_name_prefix){
      $prev_name_prefix = $p;
    }
    if((@tmp != 3 && @names>0) || $prev_name_prefix ne $p){
      &flush(\@names,\@reads,\@opts,\@qvs);
      @names=();
      @reads=();
      @opts=();
      @qvs=();
      $prev_name_prefix = $p;
    }
    push @names,$name;
    push @reads,$bases;
    push @opts,$opt;
    push @qvs,$qv;
  }
}
&flush(\@names,\@reads,\@opts,\@qvs);

sub flush($$$$){
  my $n = shift;
  my $r = shift;
  my $o = shift;
  my $q = shift;

  my @names = @$n;
  my @reads = @$r;
  my @opts = @$o;
  my @qvs = @$q;

  my $longest = length($reads[0]);
  my $i_adequate=0;
  # use longest subread
  if(@reads > 1){
    for(my $i=1; $i<@reads; ++$i){
      if(length($reads[$i]) > $longest){
        $longest = length($reads[$i]);
        $i_adequate = $i;
      }
    }
  }
=pod
  # discard last one (because it may have low quality)
  # if # of subreads > 2, then use 1st or 2nd read.
  # else use 1st read
  if(@reads>2){
    if(length($reads[1]) > $longest){
      $longest = length($reads[1]);
      $i_adequate = 1; 
    }
  }
=cut
  printf("\@%s\n",$names[$i_adequate]);
  printf("%s\n",$reads[$i_adequate]);
  printf("%s\n",$opts[$i_adequate]);
  printf("%s\n",$qvs[$i_adequate]);
}
