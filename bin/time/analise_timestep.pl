#!/usr/bin/perl -sw
#
# analyses minimum timestep used
#
# use "-min=secs" to print timesteps smaller then secs
#
#--------------------------------------------------

use strict;

$::min = 0 unless $::min;

my $timemin = 100000.;

while(<>) {

  chomp;

  my @f = split;

  my $time = $f[0];
  my $sync = $f[3];
  my $timestep = $f[4];

  if( $::min and $sync == 0 and $timestep <= $::min ) {
    print "$_\n";
  }

  if( $sync == 0 and $timestep <= $timemin ) {
    $timemin = $timestep;
  }

}

print "minimum timestep used: $timemin\n";

#--------------------------------------------------

