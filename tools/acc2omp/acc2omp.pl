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

acc2omp.pl - OpenACC to OpenMP Target converter for the ectrans GPU sources

=head1 SYNOPSIS

  acc2omp.pl --check [--pairs] [--summary] file.F90 ...

=head1 DESCRIPTION

Check mode: verify that every OpenACC directive of the given free-form
Fortran files belongs to the subset that can be translated to OpenMP
Target. Directives inside a cpp block whose condition names ACCGPU or
OMPGPU are hand-written OpenACC/OpenMP pairs: they are not checked (they
will be copied unchanged). With --pairs, they are checked too: this
shows what must change before the OMPGPU half of a pair can be removed.

Each directive outside the subset is reported as

  file:line: error: message: directive text

and the exit status is 1 if there is at least one error.

The subset:

  DATA            PRESENT COPYIN COPYOUT COPY CREATE IF ASYNC
  END DATA
  ENTER DATA      COPYIN CREATE IF ASYNC
  EXIT DATA       DELETE COPYOUT IF ASYNC
  UPDATE          HOST DEVICE IF ASYNC
  HOST_DATA       USE_DEVICE
  END HOST_DATA
  PARALLEL LOOP   COLLAPSE PRIVATE FIRSTPRIVATE DEFAULT(NONE) TILE ASYNC
                  PRESENT COPYIN COPYOUT COPY CREATE
                  (TILE and COLLAPSE are exclusive)
  END PARALLEL LOOP
  LOOP SEQ
  WAIT [(n)]

=cut

use strict;
use warnings;
use Getopt::Long;

# Directive names, longest first so that "PARALLEL LOOP" wins over "PARALLEL"
my @DIRECTIVES = (
  'END PARALLEL LOOP', 'END KERNELS LOOP', 'END SERIAL LOOP',
  'PARALLEL LOOP', 'KERNELS LOOP', 'SERIAL LOOP',
  'END HOST_DATA', 'END PARALLEL', 'END KERNELS', 'END SERIAL', 'END ATOMIC', 'END DATA',
  'ENTER DATA', 'EXIT DATA',
  'HOST_DATA', 'PARALLEL', 'KERNELS', 'SERIAL', 'DATA', 'LOOP', 'UPDATE', 'WAIT',
  'ROUTINE', 'DECLARE', 'ATOMIC', 'CACHE', 'INIT', 'SHUTDOWN', 'SET',
);

my @DATA_CLAUSES = qw(PRESENT COPYIN COPYOUT COPY CREATE);

# Allowed clauses for each directive of the subset
my %SUBSET = (
  'DATA'              => [@DATA_CLAUSES, qw(IF ASYNC)],
  'END DATA'          => [],
  'ENTER DATA'        => [qw(COPYIN CREATE IF ASYNC)],
  'EXIT DATA'         => [qw(DELETE COPYOUT IF ASYNC)],
  'UPDATE'            => [qw(HOST DEVICE IF ASYNC)],
  'HOST_DATA'         => [qw(USE_DEVICE)],
  'END HOST_DATA'     => [],
  'PARALLEL LOOP'     => [@DATA_CLAUSES, qw(COLLAPSE PRIVATE FIRSTPRIVATE DEFAULT TILE ASYNC)],
  'END PARALLEL LOOP' => [],
  'LOOP'              => [qw(SEQ)],
  'WAIT'              => [],
);

# Read a file; return a list of logical directives:
#   { line => first line, last => last line, text => joined text, skip => reason or undef }
# and a list of errors found while reading (continuations broken by cpp, ...)
sub readDirectives
{
  my $file = shift;
  open(my $fh, '<', $file) or die "$file: $!\n";
  my @lines = <$fh>;
  close($fh);

  my (@dir, @err);
  my @cpp;   # stack of booleans: true if the block is a hand-written ACCGPU/OMPGPU pair
  my $cur;   # directive being continued

  for my $i (0 .. $#lines)
  {
    my $ln = $i + 1;
    my $s = $lines[$i];
    $s =~ s/\r?\n$//;

    if ($s =~ m/^\s*#\s*(\w+)(.*)/)
    {
      my ($kw, $rest) = (lc($1), $2);
      if ($cur)
      {
        push @err, [$ln, "cpp line inside a directive continuation", $cur->{text}];
        push @dir, $cur; $cur = undef;
      }
      if ($kw =~ m/^if/)
      {
        push @cpp, ($rest =~ m/\b(ACCGPU|OMPGPU)\b/) ? 1 : 0;
      }
      elsif ($kw eq 'endif')
      {
        pop @cpp;
      }
      next;
    }

    if ($s =~ m/^\s*!\$acc(&?)(.*)$/io)
    {
      my ($amp, $body) = ($1, $2);
      $body =~ s/!.*$//o;   # trailing comment
      if (! $cur && $amp)
      {
        push @err, [$ln, "continuation line without a directive", $body];
        next;
      }
      if ($cur)
      {
        $cur->{text} .= ' ' . $body;
        $cur->{last} = $ln;
      }
      else
      {
        $cur = {line => $ln, last => $ln, text => $body,
                skip => (grep { $_ } @cpp) ? 'pair' : undef};
      }
      if ($cur->{text} =~ s/\s*&\s*$//o)
      {
        next;   # continued on the next line
      }
      push @dir, $cur; $cur = undef;
      next;
    }

    if ($cur)
    {
      push @err, [$cur->{line}, "directive continuation not followed by !\$ACC", $cur->{text}];
      push @dir, $cur; $cur = undef;
    }
  }
  push @dir, $cur if $cur;

  for (@dir)
  {
    $_->{text} =~ s/\s+/ /go;
    $_->{text} =~ s/^ | $//go;
  }

  return (\@dir, \@err);
}

# Parse a directive text: return (name, [[clause, args or undef], ...]) or die with a message
sub parseDirective
{
  my $text = shift;

  my $name;
  for my $d (@DIRECTIVES)
  {
    (my $re = $d) =~ s/ /\\s+/g;
    if ($text =~ s/^($re)\b\s*//i)
    {
      $name = $d;
      last;
    }
  }
  die "unknown directive\n" unless $name;

  my @clauses;
  if ($name eq 'WAIT' && $text =~ s/^\(([^()]*)\)\s*//o)
  {
    push @clauses, ['(wait-arg)', $1];
  }
  while (length($text))
  {
    $text =~ s/^[\s,]+//o;
    last unless length($text);
    $text =~ s/^(\w+)\s*//o or die "cannot parse clauses at '$text'\n";
    my $clause = uc($1);
    my $args;
    if ($text =~ m/^\(/o)
    {
      my $depth = 0;
      my $i;
      for ($i = 0; $i < length($text); $i++)
      {
        my $c = substr($text, $i, 1);
        $depth++ if ($c eq '(');
        $depth-- if ($c eq ')');
        last if ($depth == 0);
      }
      die "unbalanced parentheses\n" if ($depth != 0);
      $args = substr($text, 1, $i - 1);
      $text = substr($text, $i + 1);
    }
    push @clauses, [$clause, $args];
  }

  return ($name, \@clauses);
}

# Check a parsed directive against the subset; return a list of messages
sub checkDirective
{
  my ($name, $clauses) = @_;

  return ("directive $name is not in the subset") unless exists $SUBSET{$name};

  my %allowed = map { ($_ => 1) } @{ $SUBSET{$name} };
  my %seen;
  my @msg;
  for (@$clauses)
  {
    my ($c, $args) = @$_;
    next if ($c eq '(wait-arg)');
    $seen{$c}++;
    if (! $allowed{$c})
    {
      push @msg, "clause $c is not allowed on $name";
    }
    elsif ($c eq 'DEFAULT' && uc($args // '') !~ m/^\s*NONE\s*$/)
    {
      push @msg, "only DEFAULT(NONE) is allowed";
    }
    elsif ($c ne 'SEQ' && $c ne 'ASYNC' && ! defined($args))
    {
      push @msg, "clause $c needs an argument";
    }
  }
  push @msg, "TILE and COLLAPSE on the same directive" if ($seen{TILE} && $seen{COLLAPSE});
  push @msg, "LOOP without SEQ" if ($name eq 'LOOP' && ! $seen{SEQ});
  push @msg, "UPDATE without HOST or DEVICE" if ($name eq 'UPDATE' && ! ($seen{HOST} || $seen{DEVICE}));
  return @msg;
}

my ($check, $pairs, $summary);
GetOptions('check' => \$check, 'pairs' => \$pairs, 'summary' => \$summary)
  or die "usage: acc2omp.pl --check [--pairs] [--summary] file.F90 ...\n";
die "usage: acc2omp.pl --check [--pairs] [--summary] file.F90 ...\n" unless ($check && @ARGV);

my $nerr = 0;
my (%count, %skipped);
for my $file (@ARGV)
{
  my ($dir, $err) = &readDirectives($file);
  my @msg = @$err;   # [line, message, directive text]
  for my $d (@$dir)
  {
    $skipped{$file}++ if $d->{skip};
    next if ($d->{skip} && ! $pairs);
    my ($name, $clauses);
    eval { ($name, $clauses) = &parseDirective($d->{text}) };
    if ($@)
    {
      chomp(my $msg = $@);
      push @msg, [$d->{line}, $msg, $d->{text}];
      next;
    }
    $count{$name}++;
    push @msg, [$d->{line}, $_, $d->{text}] for (&checkDirective($name, $clauses));
  }
  for (sort { $a->[0] <=> $b->[0] } @msg)
  {
    my ($ln, $msg, $text) = @$_;
    $text =~ s/\s+/ /go; $text =~ s/^ //o;
    print "$file:$ln: error: $msg: !\$ACC $text\n";
    $nerr++;
  }
}

if ($summary)
{
  my $nskip = 0; $nskip += $_ for values %skipped;
  print "# directives checked, by name:\n";
  printf("#   %-20s %5d\n", $_, $count{$_}) for (sort keys %count);
  printf("# directives in hand-written ACCGPU/OMPGPU pairs (%s): %d in %d files\n", $pairs ? "checked" : "not checked",
         $nskip, scalar(keys %skipped));
  printf("# errors: %d\n", $nerr);
}

exit($nerr ? 1 : 0);
