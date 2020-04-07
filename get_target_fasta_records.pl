#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $sort_by_length=0;;

GetOptions('s'=>\$sort_by_length);

my $error_message ="USAGE: <this> <all.fa> <part.target_names (generated by\n\t\$ cut -f 3 part.sam > part.target_names)>\n\t[-s: sort by length]";

if(@ARGV != 2){
  die "$error_message\n";
}

my $in_fa=$ARGV[0];
my $in_names=$ARGV[1];

open my $fh, "<", $in_fa or die "cannot open $in_fa: $!\n";

my $name = <$fh>;
chomp $name;

my $bases = "";

my %reads;

while(1){
  while(my $buf=<$fh>){
    chomp $buf;
    if($buf =~ /^>/){
      $name =~ s/^>//;
      # confirm $name was not added
      if(exists($reads{$name})){
        printf STDERR ("WARNING: the record %s conflicted.\n",$name);
      }

#      $bases =~ s/^N+//i;
#      $bases =~ s/N+$//i;
      $reads{$name} = $bases;

      $name = $buf;
      $bases= "";
      last;
    }
    else{
      $bases .= $buf;
    }
  }
  if(eof){
#    $bases =~ s/^N+//i;
#    $bases =~ s/N+$//i;
    $reads{$name} = $bases;
    last;
  }
}
close $fh;

my @names = keys %reads;
my %exists_in_targets;
if($sort_by_length){
  @names = sort { length($reads{$b}) <=> length($reads{$a}) } keys %reads;
}
for(my $i=0; $i<@names; ++$i){
  $exists_in_targets{$names[$i]} = 0;
}

open my $targets_fh, "<", $in_names or die "cannot open $in_names: $!\n";

while(my $target_name = <$targets_fh>){
  chomp $target_name;
  if(exists($reads{$target_name})){
    $exists_in_targets{$target_name} = 1;
#    printf(">%s\n",$target_name);
#    printf("%s\n",$reads{$target_name});
  }
  else{
    printf STDERR ("WARNING: $target_name does not exist in $in_fa\n");
  }
  if(!$sort_by_length){
    printf(">%s\n",$target_name);
    printf("%s\n",$reads{$target_name});
  }
}

close $targets_fh;

if(!$sort_by_length){
  exit;
}
else{
  for(my $i=0; $i<@names; ++$i){
    if($exists_in_targets{$names[$i]}){
      printf(">%s\n",$names[$i]);
      printf("%s\n",$reads{$names[$i]});
    }
  }
}