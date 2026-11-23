#!/usr/bin/perl -sw
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2026  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# reads revision log of file and returns lines after compare_date
# 
# checks if dates are in increasing order
#
#----------------------------------------------

use strict;

$::debug = 0;
$::exit_on_error = 1;
$::write_on_error = 0;
$::has_copyright = 0;
$::has_other_copyright = 0;
$::out_of_header = 0;

#$::check_rev_log = 0;
$::check_only = 0 unless $::check_only;
$::compare_date = "none" unless $::compare_date;

my $file = $ARGV[0];

my $aux = $file;
$aux =~ s/.*\.\.\///;				#file starts with ../
$aux =~ s/\.\///;				#file starts with ./
my $ext = $aux;
$ext =~ s/^.*\.//;				#pop everything before .
$ext = "NONE" if $ext eq $aux;			#no extension
$ext = "HIDDEN" if $aux =~ /^\./;		#file starts with .

if( $::debug ) {
  print STDERR "--------------------------\n";
  print STDERR "file: $file\n";
  print STDERR "extension: $ext\n";
  print STDERR "compare date: $::compare_date\n";
  print STDERR "check_only: $::check_only\n";
  print STDERR "--------------------------\n";
}

if( $ext eq "f90" ) {
  treat_f90();
} elsif( $ext eq "f" ) {
  treat_f90();
} elsif( $ext eq "sh" ) {
  treat_for_copyright();
} elsif( $ext eq "pl" ) {
  treat_for_copyright();
} elsif( $ext eq "txt" ) {
  treat_for_copyright();
} elsif( $ext eq "make" ) {
  treat_for_copyright();
} elsif( $ext eq "o" ) {
  ignore();
} elsif( $ext eq "NONE" ) {
  if( $file eq "Makefile" ) {
    treat_for_copyright();
  } elsif( $file eq "README" ) {
    treat_for_copyright();
  } else {
    if( is_script() ) {
      treat_for_copyright();
    } else {
      treat_for_copyright();
      exit 17;				#no extension, no script - warning
    }
  }
} elsif( $ext eq "HIDDEN" ) {
  ignore();
} else {
  if( $::write_on_error ) {
    print STDERR "  *** cannot yet handle extension $ext\n";
  }
  if( $::exit_on_error ) {
    exit 3
  }
}

if( $::debug ) {
  print STDERR "using compare date: $::compare_date\n";
}

#----------------------------------------------

sub treat_f90
{
  my $in_revision_log = 0;
  my $has_revision_log = 0;
  my $revs = 0;
  my $date_old = "";
  my $error = 0;

  while(<>) {
    chomp;
    check_header();
    last if $::out_of_header;		# we are out of header
    check_copyright();
    if( /^[!cC]\s+revision log/ ) {
      $in_revision_log = 1;
      $has_revision_log = 1;
      <>;				# skip next empty line
    } elsif( $in_revision_log ) {
      last if /^\s*$/;
      last if /^[!cC]\s*$/;
      last if /^[!cC]\s*\*\*\*/;
      next if /^[!cC]\s*\.\.\./;
      my @f = split;
      my $date = $f[1];
      $date = invert_date($date,".");
      if( $::debug ) {
        print STDERR "line = $_\n";
        print STDERR "date = $date\n";
      }
      if( $date_old gt $date ) {
        if( $::write_on_error ) {
          print STDERR "  *** dates are not in order\n";
        }
	print STDERR "  *** date error: $date_old  $date\n";
	$error++;
      }
      if( $date ge $::compare_date ) {
	$revs++;
        print STDERR "  $_\n" unless $::check_only;
      }
      $date_old = $date;
    }
  }

  if( $error ) {
    if( $::write_on_error ) {
      print STDERR "  *** some dates were out of order\n";
    }
    exit 7;
  }
  if( $::has_copyright == 0 and $::has_other_copyright == 0 ) {
    if( $::write_on_error ) {
      print STDERR "  *** no copyright found\n";
    }
    exit 9;
  }
  if( $has_revision_log == 0 ) {
    if( $::write_on_error ) {
      print STDERR "  *** no revision log found\n";
    }
    exit 13;
  }
  if( $revs == 0 ) {
    if( $::write_on_error ) {
      print STDERR "  *** no update to revision log found\n";
    }
    exit 5;
  }
  if( $::has_other_copyright == 1 ) {
    exit 19;
  }

  exit 0;
}

#----------------------------------------------

sub treat_for_copyright
{

  while(<>) {
    chomp;
    check_header();
    last if $::out_of_header;		# we are out of header
    check_copyright();
  }
  if( $::has_other_copyright == 1 ) {
    exit 19;
  } elsif( $::has_copyright == 0 ) {
    if( $::write_on_error ) {
      print STDERR "  *** no copyright found\n";
    }
    exit 9;
  }
}

#----------------------------------------------

sub is_script
{
  $_ = <>;
  if( /^\#!\/bin\/bash/ ) {
    return 1;
  } elsif( /^\#!\/bin\/sh/ ) {
    return 1;
  } else {
    return 0;
  }
}

#----------------------------------------------

sub ignore
{
}

#----------------------------------------------

sub invert_date
{
  my ($date,$sep) = @_;

  #my @f = split(/$sep/,$date);
  my @f = split(/\./,$date);
  #print "$date $sep @f\n";

  return $f[2] . '-' . $f[1] . '-' . $f[0];
}

#----------------------------------------------

sub check_copyright
{
  if( /Copyright \(C\)/ ) {
    #print STDERR "  info: Copyright found\n";
    $::has_copyright = 1;
  } elsif( /COPYRIGHT/i ) {
    $::has_other_copyright = 1;
  }
}

#----------------------------------------------

sub check_header
{
  return if /^\s*$/;
  return if /^\s*[!cC]/;
  return if /^\s*[#]/;

  $::out_of_header = 1;
}

#----------------------------------------------

