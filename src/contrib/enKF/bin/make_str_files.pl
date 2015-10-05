#!/usr/bin/perl
use strict;

my $assdir = $ARGV[0];
my $strfile = $ARGV[1];
my $nrens = $ARGV[2];

unless(open(STR0,"<$assdir/../$strfile")) {die "No such file: $assdir/../$strfile\n" };

my $in_para = 0; my $in_title = 0; my $in_name = 0;
while(<STR0>){
   chomp;
   my $line = $_;
   $in_title = 1 if ($line =~ /^\$title.*/);
   $in_title = 0 if ($line =~ /^\$end.*/);
   $in_para = 1 if ($line =~ /^\$para.*/);
   $in_para = 0 if ($line =~ /^\$end.*/);
   $in_name = 1 if ($line =~ /^\$name.*/);
   $in_name = 0 if ($line =~ /^\$end.*/);

   my ($itanf) = $line =~ /itanf\s*=\s*(\d+)/ if ($in_para == 1);
   my ($itend) = $line =~ /itend\s*=\s*(\d+)/ if ($in_para == 1);
   my ($itend) = $line =~ /itend\s*=\s*(\d+)/ if ($in_para == 1);
   
   #my $itanf = $line[2] if (($in_para == 1) and ($line[0] =~ /itanf/));
   #my $itend = $line[2] if (($in_para == 1) and ($line[0] =~ /itend/));
   #print "@line\n";
   print "$in_title $in_para $in_name\n";
}

close(STR0);
