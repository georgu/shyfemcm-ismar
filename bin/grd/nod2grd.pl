#!/usr/bin/perl -w
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# writes node list to GRD format
#
# Usage: nod2grd.pl basin.grd node_list.txt
#
#----------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/lib/perl","$ENV{HOME}/shyfem/lib/perl");

use grd;
use strict;

my $grid = new grd;
my $basin = $ARGV[0];
my $nodes = $ARGV[1];
Usage() unless $basin;
Usage() unless $nodes;

my @nodes = ();
my $n = 0;

$grid->readgrd($basin);

open(FILE,"<$nodes") || die "cannot open file: $nodes\n";

while(<FILE>) {

  chomp;
  next unless $_;
  my @f = split;
  my $nnode = $f[0];
  last if $nnode == 0;

  my $node = $grid->get_node($nnode);
  my $x = $node->{x};
  my $y = $node->{y};
  print "1 $nnode 0 $x $y\n";

  push(@nodes,$nnode);
}

my $nn = @nodes;

print "\n";
print "3 1 0 $nn\n";

my $i = 0;
foreach my $nnode (@nodes) {
  $i++;
  print " $nnode";
  print "\n" if $i%8 == 0;
}
print "\n";

#----------------------------------------------------

sub Usage {

  die "Usage: nod2grd.pl basin.grd node_list.txt\n";
}
