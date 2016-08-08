#!/usr/bin/perl
# Change a sea level file from ispra format to enkf, removing tide.

# Ancona 13.506014       43.624636
# Otranto 18.497192       40.2
# Venezia 12.425681       45.418781

use lib ("$ENV{SHYFEMDIR}/femlib/perl","$ENV{HOME}/shyfem/femlib/perl");
use date;
use strict;


# input values
my $date0 = $ARGV[0];
my $file_ispra = $ARGV[1];
my $file_astro = $ARGV[2];
my $x = $ARGV[3];
my $y = $ARGV[4];
my $stdv = $ARGV[5];

Usage() if ( not $stdv );

# output file
my $filout = $file_ispra.'.enkf';

# Init date at time 0
my $datep = new date;
$datep->init_date($date0);

# reads ispra obs file and astro file
%::zo = (); %::astro = ();
&read_ispra($file_ispra);
&read_astro($file_astro);

# writes the output file
open(FOUT,">$filout");

my @keys_s = sort {$a <=> $b} (keys %::zo);

foreach my $k (@keys_s){
  
  my $val = $::zo{$k} - $::astro{$k};
  print FOUT "$k level $x $y 0 $val $stdv\n";

}
close(FOUT);

######################################

###
sub Usage(){
 print "\nUsage: ispra2enkf_level.pl date0 ispra-file astro-file xcoord ycoord stdv\n\n";
 print "date0 = date of time 0 (yyyymmdd)\n";
 print "ispra-file = sea level file in ispra format\n";
 print "astro-file = annual tide file in Tappy format\n";
 print "xcoord = x coordinate of the station (same reference of the grid)\n";
 print "ycoord = y coordinate of the station (same reference of the grid)\n";
 print "stdv = measurement error (in metres)\n\n";
 exit;
}

####
sub read_ispra(){
use strict;
my $filin = $_[0];

unless(open(FO,"<$filin")) {die 'no such file'}

# skip first lines
<FO>; <FO>; <FO>; <FO>; <FO>;

while(<FO>){
  chomp;
  my ($date,$hour,$lev) = split(/;/);
  my $yyyy=substr $date,0,4;
  my $mm=substr $date,4,2;
  my $dd=substr $date,6,2;
  my ($HH,$MM) = split(/:/,$hour);

  $lev =~ s/,/./;
  $lev = $lev + 0.;     # To remove the final ^M

  my $it = $datep->convert_to_it($yyyy,$mm,$dd,$HH,$MM,0);

  $::zo{$it} = $lev;
}
close(ZO);
}

####
sub read_astro(){
use strict;
my $filin = $_[0];

unless(open(ZA,"<$filin")) {die 'no such file'}

while(<ZA>){
  chomp;
  my ($itt,$lev,$date,$hour) = split;
  my $yyyy=substr $date,0,4;
  my $mm=substr $date,5,2;
  my $dd=substr $date,8,2;
  my ($HH,$MM,$SS) = split(/:/,$hour);

  my $it = $datep->convert_to_it($yyyy,$mm,$dd,$HH,$MM,$SS);

  $::astro{$it} = $lev;
}
close(ZA);
}


