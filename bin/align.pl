#!/usr/bin/perl -w
#
# aligns items separated by white space
#
# Usage: align.pl file
#
#--------------------------------------

use strict;

@::length = ();
my @lines = ();

while(<>) {

  chomp;

  push(@lines,$_);

  my @f = split;

  make_length(@f);

}

my $line = join(" ",@::length);
#print "length: $line\n";

foreach $line (@lines) {
  
  my @f = split(/\s+/,$line);

  my $i=0;
  foreach my $item (@f) {
    next unless $item;		# get rid of leading white space
    my $l = length($item);
    my $ll = $::length[$i];
    #print STDERR "in loop: $i $l $ll\n";
    my $dl = $ll - $l;
    if( $dl > 0 ) {
      while( $dl-- ) {
        print " ";
      }
    }
    print "$item ";
    $i++;
  }
  print "\n";
}

#--------------------------------------

sub make_length {

  my $i=0;
  foreach my $item (@_) {
    my $l = length($item);
    my $ll = $::length[$i];
    $ll = 0 unless $ll;
    $::length[$i] = $l if $l > $ll;
    $i++;
  }

}

#--------------------------------------

