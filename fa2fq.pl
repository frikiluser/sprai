#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $opt_help=0;
my $opt_shortname=0;
my $QV="5";

GetOptions(
  'q=s'=>\$QV,
  'n'=>\$opt_shortname,
  'help'=>\$opt_help
);
if(length($QV) != 1){
  printf STDERR ("-q <one letter>\n");
  printf STDERR ("you gave %s (length %d>1)\n",$QV,length($QV));
  exit(1);
}

my @msgs=(
  "USAGE: <this> <in.fasta>",
  "[-q=s specify quality value (default: '5')]",
  "[-n convert each read name to 1,2,3,...]",
  "[-h show this message]"
);

if($opt_help){
  my $msg = join("\n\t",@msgs);
  printf STDERR ("%s\n",$msg);
  exit(0);
}

my $counter=0;

my $heading;
my $name;

{
  my $l1 = <>;
  chomp $l1;
  $heading = substr($l1,0,1);
  my $l2 = <>;
  chomp $l2;
  if(eof){
    $l1 =~ s/^$heading/\@/;
    if($opt_shortname){
      printf("\@%d\n",++$counter);
    }
    else{
      printf("%s\n",$l1);
    }
    my $qvs = "";
    for(my $i=0; $i<length($l2); ++$i){
      $qvs .= $QV;
    }
    printf("%s\n",$l2);
    printf("+\n");
    printf("%s\n",$qvs);
    exit 0;
  }

  my $l3 = <>;
  chomp $l3;
  while(substr($l3,0,1) ne "+" && substr($l3,0,1) ne $heading){
    $l2 .= $l3;
    if(eof){
      $l1 =~ s/^$heading/\@/;
      if($opt_shortname){
        printf("\@%d\n",++$counter);
      }
      else{
        printf("%s\n",$l1);
      }
      my $qvs = "";
      for(my $i=0; $i<length($l2); ++$i){
        $qvs .= $QV;
      }
      printf("%s\n",$l2);
      printf("+\n");
      printf("%s\n",$qvs);
      exit 0;
    }
    $l3 = <>;
    chomp $l3;
  }
  if(substr($l3,0,1) eq $heading){
    # fasta
    my $qvs = "";
    for(my $i=0; $i<length($l2); ++$i){
      $qvs .= $QV;
    }
    $l1 =~ s/^$heading/\@/;
    if($opt_shortname){
      printf("\@%d\n",++$counter);
    }
    else{
      printf("%s\n",$l1);
    }
    printf("%s\n",$l2);
    printf("+\n");
    printf("%s\n",$qvs);
    $name = $l3;
  }
  else{
    # fastq
    my $l4 = <>;
    chomp $l4;
    my $qvs = $l4;
    while($l4 = <>){
      chomp $l4;
      if(substr($l4,0,1) eq $heading){
        last;
      }
      $qvs .= $l4;
    }
    $l1 =~ s/^$heading/\@/;
    if($opt_shortname){
      printf("\@%d\n",++$counter);
    }
    else{
      printf("%s\n",$l1);
    }
    printf("%s\n",$l2);
    printf("%s\n",$l3);
    printf("%s\n",$qvs);
    if(eof){
      exit 0;
    }
    while(1){
      $name = $l4;
      my $bases="";
      my $opts="";
      $qvs="";
      $bases = <>;
      chomp $bases;
      $opts = <>;
      chomp $opts;
      while(substr($opts,0,1) ne "+"){
        $bases .= $opts;
        $opts = <>;
        chomp $opts;
      }
      while($l4 = <>){
        chomp $l4;
        if(substr($l4,0,1) eq $heading){
          last;
        }
        $qvs .= $l4;
      }
      $name =~ s/^$heading/\@/;
      printf("%s\n",$name);
      printf("%s\n",$bases);
      printf("%s\n",$opts);
      printf("%s\n",$qvs);
      if(eof){
        exit 0;
      }
      $name = $l4;
    }
  }
}

# fasta

#my $name = <>;
#chomp $name;
$name =~ s/^$heading/\@/;
++$counter;
my $bases = "";
my $qval = "";
while(1){
  while(my $buf=<>){
    chomp $buf;
    if($buf =~ /^$heading/){
      for(my $i=0; $i < length($bases); ++$i){
        $qval .= $QV;
      }
      if($opt_shortname){
        printf("\@%d\n",$counter);
      }
      else{
        printf("%s\n",$name);
      }
      printf("%s\n",$bases);
      printf("+\n");
      printf("%s\n",$qval);
      $name = $buf;
      $bases= "";
      $qval = "";
      $name =~ s/^$heading/\@/;
      ++$counter;
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
for(my $i=0; $i < length($bases); ++$i){
  $qval .= $QV;
}
if($opt_shortname){
  printf("\@%d\n",$counter);
}
else{
  printf("%s\n",$name);
}
printf("%s\n",$bases);
printf("+\n");
printf("%s\n",$qval);
