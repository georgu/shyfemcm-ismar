#!/usr/bin/perl -w
#
# parses subroutine param_init to write namelist
#
#-------------------------------------------------------------

use strict;

my $in_routine = 0;
my $section = "";
my ($var,$val,$info);
#my (@var,@val,@info);

open(FILE,">weutro.nml");
print FILE "\n";

while(<>)
{
  chomp;
  next unless $_;
  next if /^\s*$/;

  if( /^\s+subroutine param_init/ ) { $in_routine = 1; }
  if( /^\s+end/ ) { $in_routine = 0; }
  next unless $in_routine;

  next if /^!.*=/;	#commented parameter values
  next if /^!---/;	#separator
  next if /^!!/;	#real comment

  if( /^!\s+(.*)/ ) {
    handle_section($section);
    $section = $1;
    $section =~ s/\s+\(.*//;
    $section =~ s/subroutine\s*//;
    $section =~ s/ /_/g;
    push(@::sections,$section);
    next;
  }
  next unless $section;

  if( /^\s+(\S+)\s*=\s*(\S+)\s*$/ ) {
    $var = $1;
    $val = $2;
    $info = "";
  } elsif( /^\s+(\S+)\s*=\s*(\S+)\s*!(.*)$/ ) {
    $var = $1;
    $val = $2;
    $info = $3;
  } else {
    die("cannot parse: $_\n");
  };
  unless( is_number($val) ) {
    #die "not a regular value: $val in line: $_\n";
    print STDERR " *** not a regular value: $val in line: $_\n";
    push(@::later,"$_\n");
  }
  push(@::var,$var);
  push(@::val,$val);
  push(@::info,$info);

  #print "$var | $val | $info\n";
  #print "$_\n";
}

handle_section($section);

write_fortran();

print "\n";
print "sections found:\n";
foreach $section (@::sections) {
  print "   $section\n";
}
print "\n";

print "later lines found:\n";
foreach my $line (@::later) {
  print "   $line";
}
print "\n";

#-------------------------------------------------------------

sub handle_section
{
  my $section = shift;

  return unless $section;

  my $n = @::var;
  my $nml = "";
  my $decl = "";
  my $vars = "";

  print FILE "!-----------------------------------------------\n";
  print FILE "! $section\n";
  print FILE "!\n";
  for(my $i=0; $i<$n; $i++ ) {
    my $var = insert_space($::var[$i],10);
    print FILE "! $var     $::info[$i]\n";
    $decl .= "$::var[$i],";
  }
  $decl =~ s/,$//;
  $vars = $decl;
  $decl = split_line($decl);
  $nml = $decl;
  $nml =~ s/\(\d\)//g;

  print FILE "!\n";
  print FILE "!-----------------------------------------------\n";

  print FILE "\n";
  print FILE " &nml_$section\n";
  for(my $i=0; $i<$n; $i++ ) {
    my $var = insert_space($::var[$i],10);
    print FILE "  $var  =   $::val[$i]\n";
  }
  print FILE " /\n";
  print FILE "\n";

  my $line = "\tnamelist /nml_$section/ $nml\n";
  print STDERR "$line";
  push(@::namelist,$line);
  $line = "\tread(iu,nml=nml_$section,iostat=ios)\n";
  push(@::read,$line);
  $line = "\twrite(6,*) $decl\n";
  #push(@::write,$line);
  push(@::vars,$vars);

  #print STDERR "\treal $decl\n";
  @::var = ();
  @::val = ();
  @::info = ();
}

#-------------------------------------------------------------

sub split_line
{
  my $string = shift;

  my $nmax = 8;
  my @strings = split(/,/,$string);

  my $newstring = "";
  my $i = 0;

  foreach my $s (@strings) {
    $i++;
    if( $i%$nmax == 0 ) {
      $newstring .= " &\n     &\t\t";
    }
    $newstring .= "$s,";
  }
  $newstring =~ s/,$//;

  return $newstring;
}

#-------------------------------------------------------------

sub write_block
{
  foreach my $l (@_) {
    print FORTRAN "$l";
  }
}

#-------------------------------------------------------------

sub is_number
{
  my $s = shift;

  return 0 if $s =~ /[a-zA-Z_]/;

  return 1;
}

#-------------------------------------------------------------

sub write_fortran
{

my $sep = "*" x 68;
$sep = "!" . $sep;

open(FORTRAN,">read_nml.f90");

print FORTRAN "\n";
print FORTRAN "$sep\n";
print FORTRAN "\n";

print FORTRAN "\tsubroutine read_weutro_param(param_file)\n";
print FORTRAN "\n";
print FORTRAN "\timplicit none\n";
print FORTRAN "\n";
print FORTRAN "\tcharacter*(*) param_file\n";
print FORTRAN "\n";
print FORTRAN "\tinclude 'weutro.h'\n";
print FORTRAN "\n";
print FORTRAN "\treal iavpar\n";
print FORTRAN "\tcommon /iavpar/ iavpar\n";
print FORTRAN "\tsave /iavpar/\n";
print FORTRAN "\n";

write_block(@::namelist);

print FORTRAN "\n";
print FORTRAN "\tinteger iu,ios\n";
print FORTRAN "\n";
print FORTRAN "\tiu = 1\n";
print FORTRAN "\topen(iu,file=param_file,status='old')\n";
print FORTRAN "\n";

write_block(@::read);

print FORTRAN "\n";
print FORTRAN "\tclose(iu)\n";
print FORTRAN "\n";
write_block(@::later);
print FORTRAN "\n";
print FORTRAN "\tend\n";
print FORTRAN "\n";

print FORTRAN "$sep\n";
print FORTRAN "\n";

print FORTRAN "\tsubroutine write_weutro_param(iunit)\n";
print FORTRAN "\n";
print FORTRAN "\timplicit none\n";
print FORTRAN "\n";
print FORTRAN "\tinteger iu,iunit\n";
print FORTRAN "\tinclude 'weutro.h'\n";
print FORTRAN "\n";
print FORTRAN "\treal iavpar\n";
print FORTRAN "\tcommon /iavpar/ iavpar\n";
print FORTRAN "\tsave /iavpar/\n";
print FORTRAN "\n";
print FORTRAN "\tiu = iunit\n";
print FORTRAN "\tif( iu == 0 ) iu = 6 \n";
print FORTRAN "\n";

my $n = @::sections;
  for(my $i=0; $i<$n; $i++ ) {
  my $section = $::sections[$i];
    print FORTRAN "\twrite(iu,*) 'section = $section'\n";
    my @vars = split(/,/,$::vars[$i]);
    foreach my $v (@vars) {
      #print FORTRAN "! $v\n";
      $v = insert_space($v,10);
      print FORTRAN "\t   write(iu,*) '   $v = ',$v\n";
    }
  }

print FORTRAN "\n";
print FORTRAN "\tend\n";
print FORTRAN "\n";

print FORTRAN "$sep\n";

print FORTRAN "\n";
print FORTRAN "\tprogram read_nml\n";
print FORTRAN "\timplicit none\n";
print FORTRAN "\tcall read_weutro_param('weutro.nml')\n";
print FORTRAN "\tcall write_weutro_param(0)\n";
print FORTRAN "\tend\n";
print FORTRAN "\n";

print FORTRAN "$sep\n";
print FORTRAN "\n";
}

sub fstatement {
  print FORTRAN "\t$_\n";
}

sub fcomment {
  print FORTRAN "!$_";
}

sub fempty {
  print FORTRAN "\n";
}

sub insert_space
{
  my ($text,$n) = @_;

  my $t = length($text);
  return $text if( $t >= $n );

  my $missing = $n - $t;
  my $space = " " x $missing;

  return $text . $space;
}

