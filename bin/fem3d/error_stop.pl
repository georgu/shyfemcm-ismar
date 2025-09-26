#!/usr/bin/perl
#
#----------------------------------------

my $error = 0;

while(<>) {

  if( /^\s*stop\s*\'error\s+stop/ ) {
    $error++;
    print "$_";
  }

}

print "  error stop statments found: $error\n";

#----------------------------------------

