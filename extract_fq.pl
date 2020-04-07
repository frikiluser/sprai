#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $opt_makelist=0;
my $id_head_character="\@";
my $opt_include="";

GetOptions('makelist'=>\$opt_makelist,'include=s'=>\$opt_include);

my $err_msg = "USAGE:\t<this> <all.fq> <exclude.list> > <the.rest.fq>\n\t\thow to make 'exclude.list':\n\t\t\t<this> -makelist <exclude.fq> > <exclude.list>\n\t<this> <all=db.sorted_by_id.fq> -include <target_ids.sorted_by_id.list> > <target.fq>\n\t\thow to make 'target_ids.sorted_by_id.list':\n\t\t\tdfq2fq.pl -list <in.idfq> | sort -n > <target_ids.sorted_by_id.list>";

if($opt_makelist){
  if(@ARGV != 1){
    die "$err_msg\n";
  }
  my $alldb = $ARGV[0];
  open my $dbh,"<".$alldb or die "cannot open $alldb: $!\n";

  my $line = <$dbh>;
  while(!eof){
    chomp $line;
    $line =~ s/^$id_head_character//;
    my $name = $line;
  
    my $bases="";
    $line =<$dbh>;
    while($line !~ /^\+/){
      chomp $line;
      $bases .= $line;
      $line = <$dbh>;
    }
    chomp $line;
    my $depths = $line;

    my $qvs="";
    $line =<$dbh>;# qvs
    while($line !~ /^$id_head_character/ || length($qvs) < length($bases)){
      chomp $line;
      $qvs.=$line;
      if(eof){
        last;
      }
      $line = <$dbh>;
    }
    print "$name\n";
  }
  close $dbh;
  exit 0;
}
else{
  if(@ARGV != 1){
    die "$err_msg\n";
  }
}

my $alldb = $ARGV[0];
my $list = $opt_include;


if(!$opt_include){
  my %allfq=();
  my %list=();# names of fastq records you would like to exclude from all.fq

  open my $lh,"<".$list or die "cannot open $list: $!\n";
  while(my $el=<$lh>){
    chomp $el;
    $list{$el} = 1;
  }
  close $lh;

  open my $dbh,"<".$alldb or die "cannot open $alldb: $!\n";

  my $line = <$dbh>;
  while(!eof){
    chomp $line;
    $line =~ s/^$id_head_character//;
    my $name = $line;
    
    my $bases="";
    $line =<$dbh>;
    while($line !~ /^\+/){
      chomp $line;
      $bases .= $line;
      $line = <$dbh>;
    }
    chomp $line;
    my $depths = $line;

    my $qvs="";
    $line =<$dbh>;# qvs
    while($line !~ /^$id_head_character/ || length($qvs) < length($bases)){
      chomp $line;
      $qvs.=$line;
      if(eof){
        last;
      }
      $line = <$dbh>;
    }
    if(!defined($list{$name})){
      print "\@$name\n$bases\n$depths\n$qvs\n";
    }
  }

  close $dbh;
}
else{# for include.list
  open my $dbfh, "<".$alldb or die "cannot open $alldb : $!\n";
  open my $listfh, "<".$list or die "cannot open $list : $!\n";
  while(my $target=<$listfh>){
    chomp $target;
    my $printed=0;
    while(my $name=<$dbfh>){
      chomp $name;
      $name =~ s/^\@//;
      my $bases = <$dbfh>;
      chomp $bases;
      my $options = <$dbfh>;
      chomp $options;
      my $qvs = <$dbfh>;
      chomp $qvs;
      if($name eq $target){
        printf("\@%s\n",$target);
        printf("%s\n",$bases);
        printf("%s\n",$options);
        printf("%s\n",$qvs);
        $printed=1;
        last;
      }
      else{
        next;
      }
    }
    if(!$printed){
      die "ERROR: NOT hit. \@$target\n";
    }
  }
  close $dbfh;
  close $listfh;
}
