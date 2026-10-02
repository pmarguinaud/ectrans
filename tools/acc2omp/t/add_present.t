#!/usr/bin/env perl
# Tests of add_present.pl; run with: prove tools/acc2omp/t
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Copy qw(copy);
use FindBin qw($Bin);

my $tool = "$Bin/../add_present.pl";
my $tmp = tempdir(CLEANUP => 1);

sub slurp { my $f = shift; open(my $fh, '<', $f) or die "$f: $!"; local $/; return <$fh>; }

for my $case (['', 'add_present.ref'], ['--conditional', 'add_present_conditional.ref'])
{
  my ($opt, $ref) = @$case;
  my $f = "$tmp/add_present.F90";
  copy("$Bin/add_present.F90", $f) or die;
  system("$^X $tool --write $opt $f > /dev/null");
  is($?, 0, "add_present $opt: exit status 0");
  is(&slurp($f), &slurp("$Bin/$ref"), "add_present $opt: output file");
  my @n1 = split(/\n/, &slurp("$Bin/add_present.F90"));
  my @n2 = split(/\n/, &slurp($f));
  is(scalar(@n2), scalar(@n1), "add_present $opt: same number of lines");
  my $before = &slurp($f);
  system("$^X $tool --write $opt $f > /dev/null");
  is(&slurp($f), $before, "add_present $opt: second run changes nothing");
}

done_testing();
