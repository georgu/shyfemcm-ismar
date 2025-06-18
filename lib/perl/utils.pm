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
# general utilities
#
# Usage:
#
# #!/usr/bin/env perl
#
# use lib ("$ENV{SHYFEMDIR}/lib/perl","$ENV{HOME}/shyfemcm/lib/perl");
# use utils;
# 
#------------------------------------------------------------------------

use strict;

#------------------------------------------------------------------------

sub example_sort_hash_table
{

    my %planets;

    foreach my $distance (sort {$a <=> $b} values %planets) {
        say $distance;
    }

    foreach my $name (sort { $planets{$a} <=> $planets{$b} } keys %planets) {
        printf "%-8s %s\n", $name, $planets{$name};
    }

    foreach my $name (sort { $planets{$a} <=> $planets{$b} or $a cmp $b } 
							keys %planets) {
        printf "%-8s %s\n", $name, $planets{$name};
    }
}

sub round {

  my ($val,$ipos) = @_;

  $ipos = 2 unless defined $ipos;

  my $result= sprintf "%.${ipos}f", "${val}";

  return $result;
}

sub print_to_var($$) {
   $_[0] .= $_[1];
}

#------------------------------------------------------------------------
1;
#------------------------------------------------------------------------

