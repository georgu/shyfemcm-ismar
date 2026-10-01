#!/usr/bin/perl -w -s
#
# plots a square around an item to better identify it
#
# now supports -node
# should also support -elem and -line
#
#-------------------------------------------------------------------

use lib ("$ENV{SHYFEMDIR}/lib/perl","$ENV{HOME}/shyfem/lib/perl");

use grd;
use strict;

$::n = 0;
$::l = 0;

$::help = 0 unless $::help;
$::h = 0 unless $::h;

$::node = 0 unless $::node;
$::elem = 0 unless $::elem;
$::line = 0 unless $::line;
$::nsquare = 1 unless $::nsquare;  # how many squares to plot (increasing size)
$::zoom = 50 unless $::zoom;

if( $::help or $::h ) {
  print STDERR "Usage: grid_zoom.pl {-node|-elem|-line}=n [options] grid\n";
  print STDERR "  -nsquare     how many squares to plot for the zoom\n";
  print STDERR "  -zoom        smallest square is 1/zoom of whole grid\n";
  die          "  n is item number\n"
}

#-------------------------------------------------------------------

my $grid = new grd;
my $file = $ARGV[0];
$grid->readgrd($file);                          #FEM grid

$grid->is_latlon();
my ($xmin,$ymin,$xmax,$ymax) = $grid->get_xy_minmax();
my $dx = $xmax - $xmin;
my $dy = $ymax - $ymin;
my $dxy = $dx;
$dxy = $dy if $dx < $dy;

if( $::node ) {
  my $nitem = $grid->get_node($::node);
  my $x = $nitem->{x};
  my $y = $nitem->{y};
  print STDERR "$::node $x $y\n";
  my $dd = $dxy / $::zoom;
  for( my $i=0; $i<$::nsquare; $i++ ) {
    write_square($x,$y,$dd);
    $dd = 2 * $dd;
  }
} else {
  die "can handle only -node yet...\n";
}

#-------------------------------------------------------------------

sub write_square
{
  my ($x,$y,$dd) = @_;

  my $x0 = $x - $dd;
  my $y0 = $y - $dd;
  my $x1 = $x + $dd;
  my $y1 = $y + $dd;

  my $list = "";
  ++$::n; print "1 $::n 3 $x0 $y0\n"; $list .= " $::n";
  ++$::n; print "1 $::n 3 $x1 $y0\n"; $list .= " $::n";
  ++$::n; print "1 $::n 3 $x1 $y1\n"; $list .= " $::n";
  ++$::n; print "1 $::n 3 $x0 $y1\n"; $list .= " $::n";
  my $nfirst = $::n - 3; $list .= " $nfirst";

  $::l++;
  print "3 $::l 3 5 $list\n";
  
}

#-------------------------------------------------------------------

