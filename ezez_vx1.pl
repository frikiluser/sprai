#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $DEVEL;
my $DEBUG;

my $opt_dryrun;
my $opt_ec_only;
my $opt_foobar;
my $sprai_path="/usr/lib/sprai";
my $blast_path="";

my $now = `date +%Y%m%d_%H%M%S`;
chomp $now;

my %original_time_stamps = ();
my @modified_file_names = ();

GetOptions(
  "n" => \$opt_dryrun,
  "devel" => \$DEVEL,
  "debug"=>\$DEBUG,
  "ec_only"=>\$opt_ec_only,
  "sprai_path=s"=>\$sprai_path,
  "now=s"=>\$now,
  "foobar"=>\$opt_foobar
);

my %params;

my @msgs = (
  "USAGE: <this> <ec.spec> <asm.spec>",
  "or: <this> <ec.spec> -ec_only",
  "[-n: only shows parameters in ec.spec and exit.]",
  "[-ec_only: does error correction and does NOT assemble]",
  #'[-debug: outputs intermediate files (not implemented)]',
  '[-now yyyymmdd_hhmmss: use a result_yyyymmdd_hhmmss directory, detect unfinished jobs and restart at the appropriate stage.]',
);

if(@ARGV == 0 || @ARGV > 2){
  my $msg = join "\n\t",@msgs;
  printf STDERR ("%s\n",$msg);
  exit(1);
}
if(@ARGV == 1 && !$opt_ec_only){
  printf STDERR ("WARNING: %s\n", "-ec_only was added");
  $opt_ec_only = 1;
}
#print "@ARGV\n";

my $spec="";
if(@ARGV == 2){
  $spec=$ARGV[1];
}

{
  my $ec_spec = $ARGV[0];
  open my $fh,"<",$ec_spec or die "cannot open $ec_spec :$!\n";
  while($_ = <$fh>){
    next if($_ =~ /^\s+$/);
    next if($_ =~ /^\s*#/);
    chomp;
    my @line = split /\s+/,$_;
  #  if(@line < 2){
  #    die "strange line in $ec_spec\n$_\n";
  #  }
    for(my $i=0; $i<@line; ++$i){
      if($line[$i] =~ /^\s*#/){
        @line = @line[0..$i-1];
        last;
      }
    }
    $params{$line[0]}=join(" ",@line[1..@line-1]);
    if($params{$line[0]}){
  #    printf("%s %s#\n",$line[0],$params{$line[0]});
    }
  }
  close $fh;
}

my $input_for_database;
my $from=0;
my $to=1;
my $partition=12;
my $evalue=1e-50;
my $num_threads=1;
my $max_target_seqs=100;
my $valid_voters=11;
my $trim=42;
my $estimated_genome_size=0;
my $ca_path="";
my $word_size=0;
my $min_len_for_query=1;
my $max_len_for_query=1000000000000000;
my $estimated_depth=0;
my $use_one_subread=0;

my $blasr_path="";
my $blasr_opt="";

if(defined($params{input_for_database})){
  $input_for_database = $params{input_for_database};
  if(!-e $input_for_database){
    die "$input_for_database does not exist.\n";
  }
}
else{
  die "specify input_for_database in ec.spec\n";
}

if(defined($params{estimated_genome_size})){
  $estimated_genome_size = $params{estimated_genome_size};
}
else{
  die "specify estimated_genome_size in ec.spec\n";
}
if($estimated_genome_size<=0){
  die "estimated_genome_size must be > 0\n";
}
if(defined($params{estimated_depth})){
  $estimated_depth = $params{estimated_depth};
}
else{
  die "specify estimated_depth in ec.spec\n";
}

if(defined($params{from})){
  $from = $params{from};
}
if(defined($params{to})){
  $to = $params{to};
}
if(defined($params{partition})){
  $partition = $params{partition};
}
if(defined($params{evalue})){
  $evalue = $params{evalue};
}
if(defined($params{num_threads})){
  $num_threads = $params{num_threads};
}
if(defined($params{valid_voters})){
  $valid_voters = $params{valid_voters};
}
else{
#  $valid_voters = 11;

  $valid_voters = int(0.8*($estimated_depth+0.0));
  if($valid_voters > 30){
    $valid_voters = 30;
  }
  #my $n_base = -s $input_for_database;
  #$n_base /= 2;
  #$valid_voters = int(0.8*$n_base/$estimated_genome_size);
  $valid_voters = ($valid_voters < 11) ? 11 : $valid_voters;
}
if(defined($params{trim})){
  $trim = $params{trim};
}
if(defined($params{ca_path})){
  $ca_path = $params{ca_path};
}
if(defined($params{sprai_path})){
  $sprai_path = $params{sprai_path};
  if(!-e "$sprai_path/nss2v_v3" || !-e "$sprai_path/fq2idfq.pl"){
    die "there is no $sprai_path/nss2v_v3 or $sprai_path/fq2idfq.pl.\nset right sprai_path in ec.spec\n";
  }
}
if(defined($params{blast_path})){
  $blast_path = $params{blast_path};
}
if(defined($params{word_size})){
  $word_size = $params{word_size};
}
if(defined($params{blasr_path})){
  $blasr_path = $params{blasr_path};
}
if(defined($params{blasr_opt})){
  $blasr_opt = $params{blasr_opt};
}
if(defined($params{min_len_for_query})){
  $min_len_for_query = $params{min_len_for_query};
}
if(defined($params{max_len_for_query})){
  $max_len_for_query = $params{max_len_for_query};
}
if(defined($params{max_target_seqs})){
  $max_target_seqs = $params{max_target_seqs};
}
if(defined($params{use_one_subread})){
  $use_one_subread = $params{use_one_subread};
}

printf STDERR ("#>- params -<#\n");
printf STDERR ("input_for_database %s\n",$input_for_database);
printf STDERR ("estimated_genome_size %g\n",$estimated_genome_size);
printf STDERR ("estimated_depth %d\n",$estimated_depth);
#printf STDERR ("from %s\n",$from);
#printf STDERR ("to %s\n",$to);
printf STDERR ("partition %s\n",$partition);
printf STDERR ("evalue %g\n",$evalue);
printf STDERR ("num_threads %d\n",$num_threads);
printf STDERR ("valid_voters %s\n",$valid_voters);
printf STDERR ("trim %d\n",$trim);
printf STDERR ("use_one_subread %d\n",$use_one_subread);
if($ca_path){
  printf STDERR ("ca_path %s\n",$ca_path);
}
if($sprai_path){
  printf STDERR ("sprai_path %s\n",$sprai_path);
}
if($blast_path){
  printf STDERR ("blast_path %s\n",$blast_path);
}
if($word_size){
  printf STDERR ("word_size %d\n",$word_size);
}
if($blasr_path){
  printf STDERR ("blasr_path %s\n",$blasr_path);
  if($blasr_opt){
    printf STDERR ("blasr_opt %s\n",$blasr_opt);
  }
}
if($min_len_for_query){
  printf STDERR ("min_len_for_query %d\n",$min_len_for_query);
}
if($max_len_for_query){
  printf STDERR ("max_len_for_query %d\n",$max_len_for_query);
}
if($max_target_seqs){
  printf STDERR ("max_target_seqs %d\n",$max_target_seqs);
}
printf STDERR ("#>- params -<#\n");

if(!$opt_ec_only && !-e $ca_path){
  die "ca_path $ca_path does not exist.\n";
}

if(!$opt_ec_only && !-e $spec){
  die "$spec does not exist.\n";
}

if($opt_dryrun){
  exit;
}

my $preprefix="c";
my $CURRENT=$partition;

my $outfmt = "7 qseqid sstart send sacc qstart qend bitscore evalue pident qseq sseq";

my $index=0;
my $command="";

my $valid_depth = 4;
my $valid_read_length = 500;
my $confident_depth = $valid_depth;
my $confident_length_coefficient = 0.75;

my $date="";
my $message="";

my $pipefail = "set -o pipefail;";

# mkdirs
my $PWD = `pwd`;
chomp $PWD;

my $result_dir="$PWD/result_$now/";
if(!-d $result_dir){
  mkdir $result_dir or die "cannot mkdir $result_dir: $!";
}
my $tmp_dir="$result_dir/tmp/";
if(!-d $tmp_dir){
  mkdir $tmp_dir or die "cannot mkdir $tmp_dir: $!\n";
}

# ec

my $orig_idfq = sprintf("$tmp_dir/$preprefix%02d.idfq.gz",0);

my $db_idfq_gz = sprintf("$tmp_dir/$preprefix%02d.db.idfq.gz",0);
my $qu_idfq_gz = sprintf("$tmp_dir/$preprefix%02d.qu.idfq.gz",0);

if($from == 0)
{
  # fq -> idfq (& id2n)
  my $PG = "fq2idfq.pl";
  if($sprai_path){
    $PG = "$sprai_path/$PG";
  }
  my $PG2 = "fqfilt.pl";
  if($sprai_path){
    $PG2 = "$sprai_path/$PG2";
  }
  my $dumbbell_filter = "dumbbell_filter.pl";
  if($sprai_path){
    $dumbbell_filter = "$sprai_path/$dumbbell_filter";
  }
  if(!$use_one_subread){
    $dumbbell_filter = "cat";
  }
  my $com1;
  my $com2;

  {
    my $f_do=1;
    my @parents;
    push @parents,$input_for_database;
    my $child=$db_idfq_gz;
    $f_do = &do_or_not(\@parents,\$child);
#    printf STDERR ("%d\n",$f_do);

    my $parent=$input_for_database;
    my $PG0 = "fa2fq.pl";
    if($sprai_path){
      $PG0 = "$sprai_path/$PG0";
    }
    if($f_do){
      $com1 = sprintf("$pipefail time cat $parent | $PG0 - | $dumbbell_filter - | $PG - --prefix $tmp_dir/db%02d | gzip -c -1 > $child.tmp && mv $child.tmp $child",$from);
      `$com1`;
    }
    else{
      $com1 = "sleep 0";
    }
  }

  {
    my $f_do=1;
    my @parents;
    push @parents,$db_idfq_gz;
    my $child=$qu_idfq_gz;
    $f_do = &do_or_not(\@parents,\$child);
#    printf STDERR ("%d\n",$f_do);

    my $parent=$db_idfq_gz;
    if($f_do){
      if($min_len_for_query > 1){
        $com2 = sprintf("$pipefail time gzip -d -c $parent | $PG2 - $min_len_for_query -max_len $max_len_for_query | gzip -c -1 > $child.tmp && mv $child.tmp $child",$from);
      }
      else{
        $com2 = sprintf("ln -s $parent $child");
      }
      `$com2`;
    }
    else{
      $com2 = "sleep 0";
    }
  }
}

for(my $index=$from; $index<$to; ++$index){
  {
    my $com1;
    my $com2;
    my $PG1 = "dfq2fq_v2.pl";
    if($sprai_path){
      $PG1 = "$sprai_path/$PG1";
    }
    if($blasr_path){
      die "under construction. sorry\n";
      #$com1 = sprintf("gzip -d -c $qu_idfq_gz | $PG1 -f - > $tmp_dir/$preprefix%02d.fasta ",$index);
      #$com1 .= sprintf(" && $blasr_path/sawriter $tmp_dir/$preprefix%02d.fasta ",$index);
    }
    else{
      my $dummy_target = sprintf("$tmp_dir/makeblastdb_%02d.done",$index);
      my @parents;
      push @parents,$db_idfq_gz;
      my $child = $dummy_target;
      my $f_do=1;
      $f_do = &do_or_not(\@parents,\$child);
      my $parent = $parents[0];
      my $MAKEBLASTDB = "makeblastdb";
      if($blast_path){
        $MAKEBLASTDB = "$blast_path/$MAKEBLASTDB";
      }
      if($f_do){
        $com1 = sprintf("$pipefail time gzip -d -c $parent | $PG1 -f - | $MAKEBLASTDB -in - -dbtype nucl -out $tmp_dir/$preprefix%02d -title $preprefix%02d 1>$child.tmp && mv $child.tmp $child",$index,$index);
      }
      else{
        $com1 = "sleep 0";
      }
    }

    my $PG2="partition_fa.pl";
    if($sprai_path){
      $PG2 = "$sprai_path/$PG2";
    }
    {
      my $dummy_target = sprintf("$tmp_dir/partition_fa_%02d.done",$index);
      my @parents;
      push @parents,$qu_idfq_gz;
      my $child = $dummy_target;
      my $f_do=1;
      $f_do = &do_or_not(\@parents,\$child);
      my $parent = $parents[0];
      if($f_do){
        $com2 = sprintf("$pipefail time gzip -d -c $parent | $PG1 -f - | $PG2 - $partition -p $tmp_dir/$preprefix%02d 1>$child.tmp && mv $child.tmp $child", $index);
      }
      else{
        $com2 = "sleep 0";
      }
    }

    `$com1 & $com2 & wait`;
  }

  {
    $command="";
    for(my $i=0; $i<$partition; ++$i){
      my $t_trim = $trim;
      if($index>0){
        $t_trim = 0;
      }
      my $PG2=sprintf("bfmt72s -c %d -u -i", $t_trim);
      if($sprai_path){
        $PG2 = "$sprai_path/$PG2";
      }
      my $PG3="nss2v_v3 -v $valid_voters";
      if($sprai_path){
        $PG3 = "$sprai_path/$PG3";
      }
      my $PG4="myrealigner -f -B $valid_voters -b 3 -d 0.5 -l 131072";
      if($sprai_path){
        $PG4 = "$sprai_path/$PG4";
      }
      my $BLASTN = "blastn -dbsize 1 -num_threads $num_threads";
      if($blast_path){
        $BLASTN = "$blast_path/$BLASTN";
      }
      if($word_size){
        $BLASTN .= " -word_size $word_size";
      }
      if($max_target_seqs){
        $BLASTN .= " -max_target_seqs $max_target_seqs";
      }
      if($DEBUG){
        die "sorry. under construction\n";
=pod
        $PG4="realigner_v2 -B $valid_voters -b 3";
        if($DEVEL){
          $PG4 = "./".$PG4;
        }
        my $PG5="realigner_v2 -f -n -B $valid_voters -b 3";
        if($DEVEL){
          $PG5 = "./".$PG5;
        }
        $command.= sprintf("cat $preprefix%02d_%04d.fa | $BLASTN -db $preprefix%02d -query - -evalue $evalue -outfmt '$outfmt' | tee $preprefix%02d_%04d.blastn | $PG2 - | tee $preprefix%02d_%04d.nss | $PG3 - | tee $preprefix%02d_%04d.vertical | $PG4 - | tee $preprefix%02d_%04d.realigned | $PG5 - |  gzip -c -1 > $preprefix%02d_%04d.dfq.gz & ", $index,$i,$index,$index,$i,$index,$i,$index,$i,$index,$i,$index,$i);
=cut
      }
      else{
        if($blasr_path){
          die "under construction. sorry 2\n";
          #$command .= sprintf("$blasr_path/blasr $tmp_dir/$preprefix%02d_%04d.fa $tmp_dir/$preprefix%02d.fasta -sa $tmp_dir/$preprefix%02d.fasta.sa -m 5 $blasr_opt -out $tmp_dir/tmp%02d_%04d.m5 ",$index,$i,0,0,$index,$i);
          #$command .= sprintf(" && cat $tmp_dir/tmp%02d_%04d.m5 | m52bfmt7 - |  $PG2 - | $PG3 - | $PG4 - | gzip -c -1 > $tmp_dir/$preprefix%02d_%04d.dfq.gz & ", $index,$i,$index,$i);
        }
        else{
          my $target = sprintf("$tmp_dir/$preprefix%02d_%04d.dfq.gz",$index,$i);
          my $parent = sprintf("$tmp_dir/$preprefix%02d_%04d.fa",$index,$i);
          my $dummy_parent_1 = sprintf("$tmp_dir/makeblastdb_%02d.done",$index);
          my $dummy_parent_2 = sprintf("$tmp_dir/partition_fa_%02d.done",$index);
          my @parents;
          push @parents,$parent;
          push @parents,$dummy_parent_1;
          push @parents,$dummy_parent_2;
          my $child = $target;
          my $f_do=1;
          $f_do = &do_or_not(\@parents,\$child);
          if($f_do){
            $command .= sprintf("$pipefail time cat $parent | $BLASTN -db $tmp_dir/$preprefix%02d -query - -evalue $evalue -outfmt '$outfmt' | $PG2 - | $PG3 - | $PG4 - | gzip -c -1 > $child.tmp && mv $child.tmp $child & ", $index);
          }
          else{
            $command .= "sleep 0 & ";
          }
        }
      }
    }
    $command .= " wait ";

    `$command`;
  }

  #finish
  my $finished_file = sprintf("$result_dir/$preprefix%02d.fin.idfq.gz",$index+1);
  {
    my @n_parents=();

    my $current = $CURRENT;
    for(my $i=0,$command=""; $i<$partition; $i+=$current,$command=""){
      my $right = $i+$current;
      for(my $j=$i; $j<$right && $j<$partition; ++$j){
        my $PG1 = "dfq2fq_v2.pl";
        if($sprai_path){
          $PG1 = "$sprai_path/$PG1";
        }
        my $parent = sprintf("$tmp_dir/$preprefix%02d_%04d.dfq.gz",$index,$j);
        my $child = sprintf("$tmp_dir/$preprefix%02d.refined_%04d.fin.idfq.gz",$index,$j);
        push @n_parents, $child;

        my @parents;
        push @parents,$parent;
        my $f_do=1;
        $f_do = &do_or_not(\@parents,\$child);
        if($f_do){
          $command .= sprintf("$pipefail time gzip -d -c $parent | $PG1 - -finish -valid_depth $valid_depth -valid_read_length $valid_read_length | gzip -c -1 > $child.tmp && mv $child.tmp $child & ");
        }
        else{
          $command .= "sleep 0 & ";
        }
      }
      $command .= " wait ";

      `$command`;
    }

    {
      my $child = $finished_file;
      my $f_do=1;
      my @parents = @n_parents;
      $f_do = &do_or_not(\@parents,\$child);
      if($f_do){
        if(-e $child){
          `rm $child`;
        }
        {
          my $current = $CURRENT;
          for(my $i=0; $i<$partition; $i+=$current){
            my $files="";
            for(my $j=$i; $j<$i+$current && $j<$partition; ++$j){
              $files .=  sprintf(" $parents[$j]");
            }
            $command="";
            $command .= sprintf("cat $files >> $child.tmp && mv $child.tmp $child");
            `$command`;
          }
        }
      }
      else{
      }
    }
  }

  # prep for next iter
  my $topXx = sprintf("$result_dir/$preprefix%02d.fin.longestXx.fq",$index+1);
  {
    my $PG0 = "get_top_20x_fa.pl";
    if($sprai_path){
      $PG0 = "$sprai_path/$PG0";
    }
    my $com;

    my $parent = $finished_file;
    my @parents;
    push @parents,$parent;
    my $child = $topXx;
    my $f_do=1;
    $f_do = &do_or_not(\@parents,\$child);
    if($f_do){
      my $uuid = $now;
      #my $uuid = `uuidgen`;
      chomp $uuid;
      $com = "time gzip -d -c $parent > $result_dir/$uuid.tmp && $PG0 $result_dir/$uuid.tmp -l -g $estimated_genome_size -q -c 20 > $child.tmp && mv $child.tmp $child && rm $result_dir/$uuid.tmp";
    }
    else{
      $com = "sleep 0";
    }

    `$com`;
=pod
    # fq -> idfq (& id2n)
    my $PG = "fq2idfq.pl";
    if($sprai_path){
      $PG = "$sprai_path/$PG";
    }
    my $input_fastq_for_query = $topXx;
    $qu_idfq_gz = sprintf("$tmp_dir/$preprefix%02d.qu.idfq.gz",$index+1);
    $db_idfq_gz = sprintf("$tmp_dir/$preprefix%02d.db.idfq.gz",$index+1);

    $parent = $input_fastq_for_query;
    @parents = ();
    push @parents,$parent;
    $child = $qu_idfq_gz;
    $f_do = &do_or_not(\@parents,\$child);

    my $com1;
    if($f_do){
      $com1 = sprintf("(cat $parent | $PG - --prefix $tmp_dir/qu%02d | gzip -c -1 > $child.tmp && mv $child.tmp $child)",$index+1);
    }
    else{
      $com1 = "sleep 0";
    }

    my $parent_1 = $input_fastq_for_query;
    my $parent_2 = $input_for_database;
    @parents = ();
    push @parents,$parent_1;
    push @parents,$parent_2;
    $child = $db_idfq_gz;
    $f_do = &do_or_not(\@parents,\$child);

    my $com2;
    if($f_do){
      # don't 'cat database query | fq2idfq.pl - ... '
      # do    'cat query database | fq2idfq.pl - ... '
      $com2 = sprintf("(cat $parent_1 $parent_2 | $PG - --prefix $tmp_dir/db%02d | gzip -c -1 > $child.tmp && mv $child.tmp $child)",$index+1);
    }
    else{
      $com2 = "sleep 0";
    }

    `$com1 & $com2 & wait`;
=cut

  }
}

my $ikki="ca_ikki_v5.pl";
my $ikki_log="$PWD/ca_ikki.log";

my $ca_dir="$result_dir/CA/";

my $ca_command_2="$ikki ";
if($sprai_path){
  $ca_command_2 = "$sprai_path/$ca_command_2";
}
if($ca_path){
  $ca_command_2 .= "-ca_path $ca_path ";
}
if($sprai_path){
  $ca_command_2 .= "-sprai_path $sprai_path ";
}
my $t_from = $to-1;
$ca_command_2 .= "-from $t_from -to $to $spec $estimated_genome_size -d $result_dir -out_dir $ca_dir > $ikki_log 2>&1";
#$ca_command_2 .= "-from $from -to $to $spec $estimated_genome_size -d $result_dir -out_dir $ca_dir > $ikki_log 2>&1";

if(!$opt_ec_only && !$opt_foobar){
  my $finished_file = sprintf("$result_dir/$preprefix%02d.fin.idfq.gz",$to);
  my @parents;
  push @parents,$finished_file;
  my $child = $ikki_log;
  my $f_do=1;
  $f_do = &do_or_not(\@parents,\$child);
  if($f_do){
    `$ca_command_2`;
  }
  else{
  }
}

{
  open my $fh, ">>$tmp_dir/$now.tss" or die "cannot open $now.tss: $!\n";

  foreach my $key (@modified_file_names){
    printf $fh ("%s\t%d\n",$key,$original_time_stamps{"$key"});
  }

  close $fh;
}

# test code
if($opt_foobar){
  my $outfmt_2 = "7 qseqid sstart send sacc qstart qend bitscore evalue pident";
  my $unique_limit=-1;
  {
    `gzip -d -c $tmp_dir/c00_0000.dfq.gz | dfq2fq_v2.pl --orig_depth - | count_chars.pl - > $tmp_dir/tmp.count`;
    my $count=0;
    open my $fh,"<","$tmp_dir/tmp.count" or die "cannot open $tmp_dir/tmp.count: $!\n";
    while(my $line = <$fh>){
      ++$count;
      chomp $line;
      if($count == 3){
        $line =~ /(\d+)\s+\d+/;
        $unique_limit = $1;
        close $fh;
        last;
      }
    }
  }
  if(!defined($unique_limit) || $unique_limit<0){
    die "strange unique_limit: $unique_limit\n";
  }

  my $PG1 = "get_top_20x_fa.pl";
  if($sprai_path){
    $PG1 = "$sprai_path/$PG1";
  }
  `(gzip -d -c $result_dir/c01.fin.idfq.gz | $PG1 - -g $estimated_genome_size -q | fq2fa.pl - > $tmp_dir/c01.fin.longestXx.fa) `;

  {
    my $command="";
    $command .= "fa2fq.pl $tmp_dir/c01.fin.longestXx.fa | fq2idfq.pl - | fq2fa.pl - > $tmp_dir/d01.fin.longestXx.fa";
    `$command`;
  }

  {
    my $command1="";
    $command1.= sprintf("makeblastdb -in $tmp_dir/d01.fin.longestXx.fa -dbtype nucl -out $tmp_dir/d01 -title d01 ");
    my $command2="";
    $command2.= sprintf("partition_fa.pl $tmp_dir/d01.fin.longestXx.fa $partition -p $tmp_dir/d01");
    `$command1 & $command2 & wait`;
  }

  {
    my $command="";
    for(my $i=0; $i<$partition; ++$i){
      my $BLASTN = "blastn -dbsize 1 -num_threads $num_threads";
      if($word_size){
        $BLASTN .= " -word_size $word_size";
      }
      $command.= sprintf("cat $tmp_dir/d01_%04d.fa | $BLASTN -db $tmp_dir/d01 -query - -evalue $evalue -outfmt '$outfmt_2' | gzip -c -1 > $tmp_dir/d01_%04d.blastn.gz & ",$i,$i);
    }
    `$command wait`;
  }
}


sub do_or_not($$){
  my $a = shift;
  my $b = shift;

  my @parents = @$a;
  my $child = $$b;

  my $f_do=0;

#  printf STDERR ("%s\n",$child);
#  printf STDERR ("%s\n",$parents[0]);
  if(!-e $child){
    $f_do=1;
  }
  else{
    for(my $i=0; $i<@parents; ++$i){
      if(!-e $parents[$i]){
        $f_do=1;
        my @c_fs=stat $child;
        push @modified_file_names,$child;
        $original_time_stamps{"$child"} = $c_fs[9];
        `touch $child`;
        last;
      }
      else{
        my @p_fs=stat $parents[$i];
        my @c_fs=stat $child;
        my $p_bd = $p_fs[9];
        my $c_bd = $c_fs[9];
        if($p_bd > $c_bd){
          $f_do=1;
          push @modified_file_names,$child;
          $original_time_stamps{"$child"} = $c_fs[9];
          `touch $child`;
          last;
        }
        else{
          #$f_do=0;
        }
      }
    }
  }
  return $f_do;
}
