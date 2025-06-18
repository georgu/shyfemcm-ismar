#!/usr/bin/perl -w
#
# computes bias, rmse, correlation coefficient, and averages of two timeseries
#
# interpolates values of one timeseries to the other
#
#-------------------------------------------------------

use strict;
use lib ("$ENV{SHYFEMDIR}/lib/perl","$ENV{HOME}/shyfemcm/lib/perl");
use date;
use utils;

my $obs = shift;
my $sim = shift;

$::debug = 0;

unless( $sim ) {
  die "Usage: ./rms.pl file1.txt file2.txt\n";
}

#print "$obs $sim\n";

my ($tobs,$sobs) = read_ts($obs);
my ($tsim,$ssim) = read_ts($sim);

($tobs,$sobs,$tsim,$ssim) = choose_values($tobs,$sobs,$tsim,$ssim);

my $nn = @$tobs;

#print "==============================\n";
#for( my $i=0; $i<$nn; $i++ ) {
#  print "$tobs->[$i] $tsim->[$i] $sobs->[$i] $ssim->[$i]\n";
#}
#print "==============================\n";

my $obsaver = compute_aver($sobs);
my $simaver = compute_aver($ssim);

my $n = 0;
my $acum = 0.;
my $bias = 0.;
my $aver = 0.;
my $square = 0.;
my $sstot = 0.;

my $xy = 0.;
my $x = 0.;
my $y = 0.;
my $x2 = 0.;
my $y2 = 0.;

my $j = 0;

while( 1 ) {

  my $t1 = shift(@$tobs);
  my $s1 = shift(@$sobs);
  my $t2 = shift(@$tsim);
  my $s2 = shift(@$ssim);

  last unless defined $t1;

  die "time error: $t1 $t2\n" if $t1 ne $t2;
  $j++;
  print STDERR "processing: $j $t1 $t2 $s1 $s2\n" if $::debug;

  $xy += $s1*$s2;
  $x += $s1;
  $y += $s2;
  $x2 += $s1*$s1;
  $y2 += $s2*$s2;

  $aver += $s1;
  $square += $s1 * $s1;
  #my $err = $s1 - $obsaver;
  my $err = $s1 - $simaver;
  my $err2 = $err*$err;
  $sstot += $err2;

  my $diff = $s2 - $s1;		# sim - obs
  my $diff2 = $diff*$diff;
  $acum += $diff2;
  $bias += $diff;

  $n++;

  #print "$n  $t1  $t2   $s1  $s2\n";
}

my $rxy1 = $n*$xy - $x*$y;
my $rxy2 = sqrt( $n*$x2 - $x*$x );
my $rxy3 = sqrt( $n*$y2 - $y*$y );
my $rxy = $rxy1 / ( $rxy2 * $rxy3 );

$bias /= $n;
$acum /= $n;
my $rms = sqrt( $acum );

$rms = round($rms);
$bias = round($bias);
$rxy = round($rxy);
$obsaver = round($obsaver);
$simaver = round($simaver);

print "$bias   $rms   $rxy   $obsaver  $simaver\n";

#-------------------------------------------------------

sub read_ts {

  # reads timeseries (only one column apart the time column)

  my $file = shift;

  my @t = ();
  my @s = ();

  open(FILE,"<$file") || die "cannot open file $file\n";

  while(<FILE>) {
    next if /^\s*#/;
    my @f = split;
    last if scalar @f == 1;
    push(@t,$f[0]);
    push(@s,$f[1]);
  }

  close(FILE);

  return (\@t,\@s);
}

sub compute_aver {

  # computes average of timeseries

  my $ra = shift;

  my $acum = 0;
  foreach my $entry (@$ra) {
    $acum += $entry;
  }

  my $n = @$ra;
  $acum /= $n if $n;
  return $acum;
}

#-------------------------------------------------------

sub choose_values {

  # interpolates one timeseries (the shorter one) to the time of the other one

  my ($t1,$s1,$t2,$s2) = @_;

  my $debug = $::debug;

  my $n1 = @$t1;
  my $n2 = @$t2;

  my ($n,$tstat,$sstat,$naux,$taux,$saux);

  if( $n1 < $n2 ) {
    $n = $n1;
    $tstat = $t1;
    $sstat = $s1;
    $naux = $n2;
    $taux = $t2;
    $saux = $s2;
  } else {
    $n = $n2;
    $tstat = $t2;
    $sstat = $s2;
    $naux = $n1;
    $taux = $t1;
    $saux = $s1;
  }

  print STDERR "choosing: initial size: $n1  $n2  $n\n" if $debug;

  my $date = new date;

  my @tnew = ();
  my @snew = ();
  my @tstat = ();
  my @sstat = ();

  my $j0 = 0;
  for( my $i=0; $i<$n; $i++ ) {
    my $tline = $tstat->[$i];
    my $ss = $sstat->[$i];
    my $t = $date->unformat_abs($tline);
    last unless defined($t);
    for( my $j=$j0+1; $j<$naux; $j++ ) {
      my $ta1 = $date->unformat_abs($taux->[$j-1]);
      my $ta2 = $date->unformat_abs($taux->[$j]);
      if( $ta1 <= $t and $ta2 >= $t ) {
        print STDERR "$i $j  $taux->[$j-1] $tline $taux->[$j]\n" if $debug;
        die "time error\n" if $taux->[$j-1] gt $tline or $tline gt $taux->[$j];
        my $r = ($t-$ta1)/($ta2-$ta1);
        my $s = $saux->[$j-1] + ($saux->[$j]-$saux->[$j-1])*$r;
        push(@tnew,$tline);
        push(@snew,$s);
        push(@tstat,$tline);
        push(@sstat,$ss);
	$j0 = $j;
	last;
      }
    }
  }

  my $nn = @tstat;
  print STDERR "choosing: final size: $nn\n" if $debug;

  if( $n1 < $n2 ) {
    return(\@tstat,\@sstat,\@tnew,\@snew);
  } else {
    return(\@tnew,\@snew,\@tstat,\@sstat);
  }
  
}

#-------------------------------------------------------

