#!/usr/bin/env perl
#
# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

=head1 NAME

unpair.pl - remove hand-written ACCGPU/OMPGPU pairs, keeping line numbers

=head1 SYNOPSIS

  unpair.pl --list file.F90 ...
  unpair.pl [--dry-run] file.F90:line ...

=head1 DESCRIPTION

A pair is an "#ifdef ACCGPU" block and an "#ifdef OMPGPU" block next to each
other (in either order, empty lines allowed between them). Once a pair is
removed, the OpenACC directives are plain source lines, and acc2omp.pl
translates them for the OpenMP build.

--list prints the blocks of each file, one per line:

  file:line: pair|acc-only|omp-only, N OpenACC directives: removable
  file:line: ...: NOT removable: reason

where line is the line of the first "#ifdef". A pair or block that contains
code (lines other than directives, comments and empty lines) is not
removable.

file.F90:line removes the pair or block that starts at that line, in place.
The number of lines of the file never changes:

=over

=item R1

the "#ifdef" and "#endif" lines, the OMPGPU half and the comments of the
OMPGPU half become empty lines;

=item R2

each OpenACC directive moves to the first line of its OpenMP counterpart
(same kind, in order: compute construct, DATA, HOST_DATA, END DATA, UPDATE,
ENTER DATA, EXIT DATA), so that the translated directive is on the line of
the hand-written one, and amdflang gives its kernel the same name. A
directive without counterpart (WAIT, LOOP SEQ, ...) goes on the next free
line. A compute construct that cannot be put on the line of its counterpart
is an error. Comments of the OpenACC half are kept if their line is free.

=back

An "#ifdef ACCGPU" block alone (acc-only) loses its "#ifdef" and "#endif"
lines. An "#ifdef OMPGPU" block alone (omp-only) is removable only if it is
empty or only has END TARGET TEAMS DISTRIBUTE PARALLEL DO directives (optional
in Fortran). A pair whose OMPGPU half has a directive with no OpenACC
counterpart (other than END TARGET TEAMS DISTRIBUTE PARALLEL DO) is not
removable.

--dry-run prints what would be done, without changing the files.

=cut

use strict;
use warnings;
use Getopt::Long;

# Kind of a directive line (text after the sentinel), for the OpenACC / OpenMP matching
sub accKind
{
  my $t = uc(shift);
  return 'END-COMPUTE' if ($t =~ m/^END\s+PARALLEL\s+LOOP\b/o);
  return 'COMPUTE'     if ($t =~ m/^PARALLEL\s+LOOP\b/o);
  return 'END-DATA'    if ($t =~ m/^END\s+(HOST_)?DATA\b/o);
  return 'HOST_DATA'   if ($t =~ m/^HOST_DATA\b/o);
  return 'DATA'        if ($t =~ m/^DATA\b/o);
  return 'UPDATE'      if ($t =~ m/^UPDATE\b/o);
  return 'ENTER'       if ($t =~ m/^ENTER\s+DATA\b/o);
  return 'EXIT'        if ($t =~ m/^EXIT\s+DATA\b/o);
  return 'OTHER';
}

sub ompKind
{
  my $t = uc(shift);
  return 'END-COMPUTE' if ($t =~ m/^END\s+TARGET\s+TEAMS\b/o);
  return 'COMPUTE'     if ($t =~ m/^TARGET\s+TEAMS\b/o);
  return 'END-DATA'    if ($t =~ m/^END\s+TARGET\s+DATA\b/o);
  return 'HOST_DATA'   if ($t =~ m/^TARGET\s+DATA\b.*\bUSE_DEVICE_(ADDR|PTR)\b/o);
  return 'DATA'        if ($t =~ m/^TARGET\s+DATA\b/o);
  return 'UPDATE'      if ($t =~ m/^TARGET\s+UPDATE\b/o);
  return 'ENTER'       if ($t =~ m/^TARGET\s+ENTER\s+DATA\b/o);
  return 'EXIT'        if ($t =~ m/^TARGET\s+EXIT\s+DATA\b/o);
  return 'OTHER';
}

# The "#ifdef ACCGPU|OMPGPU" block starting at line index $i: (kind, first, last) or ()
sub block
{
  my ($lines, $i) = @_;
  return () unless ($lines->[$i] =~ m/^\s*#\s*ifdef\s+(ACCGPU|OMPGPU)\s*$/o);
  my ($kind, $depth) = ($1, 1);
  for my $j ($i + 1 .. $#$lines)
  {
    $depth++ if ($lines->[$j] =~ m/^\s*#\s*if/o);
    $depth-- if ($lines->[$j] =~ m/^\s*#\s*endif/o);
    return ($kind, $i, $j) unless $depth;
  }
  die "unterminated #ifdef at line " . ($i + 1) . "\n";
}

# Content of the lines first+1 .. last-1 of a block:
# (directives [first line, last line, text], comment lines, reason it cannot be removed or undef)
sub content
{
  my ($lines, $first, $last, $sentinel) = @_;
  my (@dir, @cmt, $why);
  my $i = $first + 1;
  while ($i < $last)
  {
    my $s = $lines->[$i];
    if ($s =~ m/^\s*!\$$sentinel\b\s*(.*)$/i)
    {
      my ($j, $text) = ($i, $1);
      while ($lines->[$j] =~ m/&\s*(!.*)?$/o && $j + 1 < $last)
      {
        $j++;
      }
      push @dir, [$i, $j, $text];
      $i = $j + 1;
      next;
    }
    if ($s =~ m/^\s*!\$/o)
    {
      $why //= "directive of another kind at line " . ($i + 1);
    }
    elsif ($s =~ m/^\s*#/o)
    {
      $why //= "cpp line at line " . ($i + 1);
    }
    elsif ($s =~ m/^\s*!/o)
    {
      push @cmt, $i;
    }
    elsif ($s !~ m/^\s*$/o)
    {
      $why //= "code at line " . ($i + 1);
    }
    $i++;
  }
  return (\@dir, \@cmt, $why);
}

# All blocks of a file: [{type, first, last (whole region), acc => [...], omp => [...]}]
sub blocks
{
  my $lines = shift;
  my @b;
  my $i = 0;
  while ($i < @$lines)
  {
    my ($k1, $s1, $e1) = &block($lines, $i);
    unless ($k1)
    {
      $i++;
      next;
    }
    my $j = $e1 + 1;
    $j++ while ($j < @$lines && $lines->[$j] =~ m/^\s*$/o);
    my ($k2, $s2, $e2) = ($j < @$lines) ? &block($lines, $j) : ();
    if ($k2 && $k2 ne $k1)
    {
      push @b, {type => 'pair', first => $s1, last => $e2,
                acc => ($k1 eq 'ACCGPU') ? [$s1, $e1] : [$s2, $e2],
                omp => ($k1 eq 'ACCGPU') ? [$s2, $e2] : [$s1, $e1]};
      $i = $e2 + 1;
    }
    else
    {
      push @b, {type => ($k1 eq 'ACCGPU') ? 'acc-only' : 'omp-only', first => $s1, last => $e1,
                ($k1 eq 'ACCGPU' ? 'acc' : 'omp') => [$s1, $e1]};
      $i = $e1 + 1;
    }
  }
  return @b;
}

# Compute the new lines of a block; return (hash line index => new text, number of OpenACC
# directives, list of notes) or die with the reason it cannot be removed
sub plan
{
  my ($lines, $b) = @_;
  my ($A, $cA, $wA) = $b->{acc} ? &content($lines, @{ $b->{acc} }, 'acc') : ([], [], undef);
  my ($O, $cO, $wO) = $b->{omp} ? &content($lines, @{ $b->{omp} }, 'omp') : ([], [], undef);
  die "$wA\n" if $wA;
  die "$wO\n" if $wO;

  my %new = map { ($_ => '') } ($b->{first} .. $b->{last});
  my @notes;

  if ($b->{type} eq 'omp-only')
  {
    for (@$O)
    {
      die "OpenMP directive at line " . ($_->[0] + 1) . " has no OpenACC counterpart\n"
        unless (&ompKind($_->[2]) eq 'END-COMPUTE');
    }
    return (\%new, 0, \@notes);
  }

  if ($b->{type} eq 'acc-only')
  {
    for my $d (@$A)
    {
      $new{$_} = $lines->[$_] for ($d->[0] .. $d->[1]);
    }
    $new{$_} = $lines->[$_] for @$cA;
    return (\%new, scalar(@$A), \@notes);
  }

  # pair: match the OpenACC directives with the OpenMP ones, by kind, in order
  my @okind = map { &ompKind($_->[2]) } @$O;
  my @used = (0) x @$O;
  my $next = 0;         # first OpenMP directive not yet passed
  my $free = $b->{first};
  my %taken;
  my @cmt = @$cA;                      # comments of the OpenACC half, in line order
  for my $d (@$A)
  {
    my ($l0, $l1) = @$d;
    my $n = $l1 - $l0 + 1;
    my @before;                         # comments just before this directive move with it
    push(@before, shift(@cmt)) while (@cmt && $cmt[0] < $l0);
    my $kind = &accKind($d->[2]);
    my $target;
    if ($kind ne 'OTHER')
    {
      for my $q ($next .. $#$O)
      {
        if ($okind[$q] eq $kind)
        {
          ($target, $used[$q], $next) = ($O->[$q][0], 1, $q + 1);
          last;
        }
      }
    }
    my $p = (defined($target) && $target >= $free) ? $target : $free;
    if ($kind eq 'COMPUTE' && (! defined($target) || $p != $target))
    {
      die "compute construct at line " . ($l0 + 1) . " cannot be put on the line of its OpenMP counterpart\n";
    }
    die "no room for the directive at line " . ($l0 + 1) . "\n" if ($p + $n - 1 > $b->{last});
    die "no room for the comments before the directive at line " . ($l0 + 1) . "\n" if ($p - @before < $free);
    $new{$p - @before + $_} = $lines->[$before[$_]] for (0 .. $#before);
    $new{$p + $_} = $lines->[$l0 + $_] for (0 .. $n - 1);
    $taken{$p + $_} = 1 for (0 .. $n - 1);
    push @notes, sprintf("line %d -> %d%s", $l0 + 1, $p + 1, $kind eq 'COMPUTE' ? ' (kernel)' : '') if ($p != $l0);
    $free = $p + $n;
  }
  for my $q (0 .. $#$O)
  {
    next if ($used[$q] || $okind[$q] eq 'END-COMPUTE');
    die "OpenMP directive at line " . ($O->[$q][0] + 1) . " has no OpenACC counterpart\n";
  }
  # comments after the last directive: on the next free lines
  die "no room for the comments at the end of the OpenACC half\n" if ($free + @cmt - 1 > $b->{last});
  $new{$free + $_} = $lines->[$cmt[$_]] for (0 .. $#cmt);
  return (\%new, scalar(@$A), \@notes);
}

my ($list, $dry);
GetOptions('list' => \$list, 'dry-run' => \$dry)
  or die "usage: unpair.pl --list file.F90 ...\n       unpair.pl [--dry-run] file.F90:line ...\n";

sub readFile
{
  my $file = shift;
  open(my $fh, '<', $file) or die "$file: $!\n";
  my @lines = <$fh>;
  close($fh);
  my $eol = (@lines && $lines[-1] !~ m/\n$/o) ? 0 : 1;
  chomp(@lines);
  return (\@lines, $eol);
}

if ($list)
{
  my ($nrem, $ndir) = (0, 0);
  for my $file (@ARGV)
  {
    my ($lines) = &readFile($file);
    for my $b (&blocks($lines))
    {
      my ($new, $n, $notes);
      eval { ($new, $n, $notes) = &plan($lines, $b) };
      my $nacc = $b->{acc} ? scalar(@{ (&content($lines, @{ $b->{acc} }, 'acc'))[0] }) : 0;
      if ($@)
      {
        chomp(my $why = $@);
        printf("%s:%d: %s, %d OpenACC directives: NOT removable: %s\n", $file, $b->{first} + 1, $b->{type}, $nacc, $why);
      }
      else
      {
        printf("%s:%d: %s, %d OpenACC directives: removable\n", $file, $b->{first} + 1, $b->{type}, $n);
        $nrem++; $ndir += $n;
      }
    }
  }
  print "# removable: $nrem blocks, $ndir OpenACC directives\n";
  exit(0);
}

die "usage: unpair.pl [--dry-run] file.F90:line ...\n" unless @ARGV;
my %todo;
for (@ARGV)
{
  m/^(.+):(\d+)$/o or die "$_: expected file:line\n";
  push @{ $todo{$1} }, $2 - 1;
}
my $err = 0;
for my $file (sort keys %todo)
{
  my ($lines, $eol) = &readFile($file);
  my %start = map { ($_->{first} => $_) } &blocks($lines);
  my @new = @$lines;
  for my $l (@{ $todo{$file} })
  {
    my $b = $start{$l};
    unless ($b)
    {
      printf STDERR ("%s:%d: error: no ACCGPU/OMPGPU block starts here\n", $file, $l + 1);
      $err++;
      next;
    }
    my ($plan, $n, $notes);
    eval { ($plan, $n, $notes) = &plan($lines, $b) };
    if ($@)
    {
      chomp(my $why = $@);
      printf STDERR ("%s:%d: error: %s not removable: %s\n", $file, $l + 1, $b->{type}, $why);
      $err++;
      next;
    }
    $new[$_] = $plan->{$_} for keys %$plan;
    printf("%s:%d: %s removed, %d OpenACC directives%s\n", $file, $l + 1, $b->{type}, $n,
           @$notes ? '; ' . join(', ', @$notes) : '');
  }
  next if ($dry || $err);
  die "$file: line count changed\n" if (@new != @$lines);
  open(my $fh, '>', "$file.tmp") or die "$file.tmp: $!\n";
  print $fh join("\n", @new), ($eol && @new ? "\n" : '');
  close($fh) or die "$file.tmp: $!\n";
  rename("$file.tmp", $file) or die "$file: $!\n";
}
exit($err ? 1 : 0);
