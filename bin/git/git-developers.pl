#!/usr/bin/perl
#
# returns list of developers
#
#---------------------------------------------

use lib ("$ENV{SHYFEMDIR}/lib/perl","$ENV{HOME}/shyfem/lib/perl");
use lib ("$ENV{SHYFEMDIR}/var/copyright");

use utils_develop

make_dev_names();

foreach my $key (sort keys %::dev_names) {
  my $name = $::dev_names{$key};
  print "  $key   $name\n";
}

#---------------------------------------------

