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

  acc2omp.pl [-o out.F90] in.F90
  acc2omp.pl --check [--pairs] [--summary] file.F90 ...

=head1 DESCRIPTION

Translation mode (default): write the OpenMP Target version of a free-form
Fortran file, on standard output or in the file given with -o. Each OpenACC
directive of the subset (below) is replaced by its OpenMP translation, on
the same lines: the output has exactly the same number of lines as the
input, so that line numbers (and the kernel names, which contain them) do
not change. If the OpenMP text needs fewer lines, the remaining lines are
empty; if it is longer, its lines are longer. A directive that has no
OpenMP equivalent (WAIT, LOOP SEQ) becomes empty lines. Directives inside
hand-written ACCGPU/OMPGPU pairs (see below) and all other lines are copied
unchanged. Any directive outside the subset is an error (as in check
mode), and then no output is written.

Translation table (synchronous: ASYNC and WAIT are dropped):

  DATA                     TARGET DATA
  END DATA                 END TARGET DATA
  ENTER DATA               TARGET ENTER DATA
  EXIT DATA                TARGET EXIT DATA
  UPDATE                   TARGET UPDATE
  HOST_DATA                TARGET DATA
  END HOST_DATA            END TARGET DATA
  PARALLEL LOOP            TARGET TEAMS DISTRIBUTE PARALLEL DO
  END PARALLEL LOOP        END TARGET TEAMS DISTRIBUTE PARALLEL DO
  LOOP SEQ, WAIT           (empty)

  PRESENT(x)               MAP(PRESENT,ALLOC:x)
  COPYIN(x)                MAP(TO:x)
  COPYOUT(x)               MAP(FROM:x)
  COPY(x)                  MAP(TOFROM:x)
  CREATE(x)                MAP(ALLOC:x)
  DELETE(x)                MAP(DELETE:x)
  HOST(x) / DEVICE(x)      FROM(x) / TO(x)
  USE_DEVICE(x)            USE_DEVICE_ADDR(x)
  TILE(a,b,...)            COLLAPSE(n), n = number of TILE arguments
  COLLAPSE, PRIVATE,
  FIRSTPRIVATE, IF         unchanged
  DEFAULT(NONE), ASYNC     (dropped)

A name in PRESENT(...) on a PARALLEL LOOP is dropped when an enclosing DATA
construct has an IF clause containing PRESENT(name) (an OPTIONAL argument
that may be absent): OpenACC ignores a data clause on an absent optional
argument, but MAP(PRESENT) would fail at run time. The enclosing DATA
constructs include those of the hand-written pairs (their ACCGPU half).

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
      $amp = '&' if ($body =~ s/^\s*&//o);   # "!$ACC &" continuation
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

# Split a clause argument list at the commas of the outer level; return the trimmed items
sub splitArgs
{
  my $args = shift;
  my (@items, $cur, $depth);
  ($cur, $depth) = ('', 0);
  for my $c (split(//, $args))
  {
    $depth++ if ($c eq '(');
    $depth-- if ($c eq ')');
    if ($c eq ',' && $depth == 0)
    {
      push @items, $cur; $cur = '';
      next;
    }
    $cur .= $c;
  }
  push @items, $cur;
  s/^\s+|\s+$//go for @items;
  return grep { length($_) } @items;
}

# Names of the data clauses of a DATA directive whose IF clause contains PRESENT(name)
sub conditionalNames
{
  my $clauses = shift;
  my %cond;
  for (@$clauses)
  {
    my ($c, $args) = @$_;
    next unless ($c eq 'IF');
    $cond{uc($1)} = 1 while ($args =~ m/\bPRESENT\s*\(\s*(\w+)\s*\)/gio);
  }
  return \%cond;
}

my %MAP = (PRESENT => 'PRESENT,ALLOC', COPYIN => 'TO', COPYOUT => 'FROM', COPY => 'TOFROM',
           CREATE => 'ALLOC', DELETE => 'DELETE');

my %OMPNAME = (
  'DATA' => 'TARGET DATA', 'END DATA' => 'END TARGET DATA',
  'ENTER DATA' => 'TARGET ENTER DATA', 'EXIT DATA' => 'TARGET EXIT DATA',
  'UPDATE' => 'TARGET UPDATE',
  'HOST_DATA' => 'TARGET DATA', 'END HOST_DATA' => 'END TARGET DATA',
  'PARALLEL LOOP' => 'TARGET TEAMS DISTRIBUTE PARALLEL DO',
  'END PARALLEL LOOP' => 'END TARGET TEAMS DISTRIBUTE PARALLEL DO',
  'LOOP' => undef, 'WAIT' => undef,
);

# Translate a parsed directive (of the subset); return the list of OpenMP tokens
# (directive name, then clauses), or an empty list if it has no translation.
# $cond: names that are conditionally present (enclosing DATA ... IF(PRESENT(name)))
sub translate
{
  my ($name, $clauses, $cond) = @_;
  my $omp = $OMPNAME{$name};
  return () unless defined($omp);
  my @tok = ($omp);
  for (@$clauses)
  {
    my ($c, $args) = @$_;
    if ($c eq 'ASYNC' || $c eq 'DEFAULT' || $c eq '(wait-arg)')
    {
      next;
    }
    elsif (exists $MAP{$c})
    {
      my @x = &splitArgs($args);
      @x = grep { ! $cond->{uc($_)} } @x if ($c eq 'PRESENT' && $name eq 'PARALLEL LOOP');
      push @tok, "MAP($MAP{$c}:" . join(',', @x) . ")" if @x;
    }
    elsif ($c eq 'HOST' || $c eq 'DEVICE')
    {
      push @tok, ($c eq 'HOST' ? 'FROM' : 'TO') . "(" . join(',', &splitArgs($args)) . ")";
    }
    elsif ($c eq 'USE_DEVICE')
    {
      push @tok, "USE_DEVICE_ADDR(" . join(',', &splitArgs($args)) . ")";
    }
    elsif ($c eq 'TILE')
    {
      push @tok, "COLLAPSE(" . scalar(&splitArgs($args)) . ")";
    }
    else   # COLLAPSE PRIVATE FIRSTPRIVATE IF
    {
      push @tok, "$c(" . join(',', &splitArgs($args)) . ")";
    }
  }
  return @tok;
}

# Lay out OpenMP tokens on exactly $n lines, with the given indentation; empty lines at the end
sub layout
{
  my ($tok, $n, $indent) = @_;
  return ('') x $n unless @$tok;
  my @groups = map { [] } 1 .. $n;
  my $total = 0; $total += length($_) + 1 for @$tok;
  my ($g, $len) = (0, 0);
  for my $i (0 .. $#$tok)
  {
    my $t = $tok->[$i];
    # start a new line when this one is full enough, if lines remain
    if (@{ $groups[$g] } && $g < $n - 1 && $len + length($t) > $total / $n)
    {
      $g++; $len = 0;
    }
    push @{ $groups[$g] }, $t;
    $len += length($t) + 1;
  }
  my @used = grep { @$_ } @groups;
  my @out;
  for my $i (0 .. $#used)
  {
    my $l = $indent . ($i ? '!$OMP& ' : '!$OMP ') . join(' ', @{ $used[$i] });
    $l .= ' &' if ($i < $#used);
    push @out, $l;
  }
  push @out, '' while (@out < $n);
  return @out;
}

# Translate a file; return (output lines, error messages [line, message, text])
sub translateFile
{
  my $file = shift;
  open(my $fh, '<', $file) or die "$file: $!\n";
  my @lines = <$fh>;
  close($fh);
  my $eol = (@lines && $lines[-1] !~ m/\n$/o) ? 0 : 1;   # final newline, kept as in the input
  chomp(@lines);

  my ($dir, $err) = &readDirectives($file);
  my @msg = @$err;
  my @stack;   # enclosing DATA constructs: hash of conditionally present names, or undef
  for my $d (@$dir)
  {
    my ($name, $clauses);
    eval { ($name, $clauses) = &parseDirective($d->{text}) };
    if ($@)
    {
      next if $d->{skip};
      chomp(my $m = $@);
      push @msg, [$d->{line}, $m, $d->{text}];
      next;
    }
    # enclosing DATA constructs, from all OpenACC directives (pairs included)
    push(@stack, &conditionalNames($clauses)) if ($name eq 'DATA');
    pop(@stack) if ($name eq 'END DATA');
    next if $d->{skip};

    my @m = &checkDirective($name, $clauses);
    if (@m)
    {
      push @msg, [$d->{line}, $_, $d->{text}] for @m;
      next;
    }
    my %cond = map { %$_ } @stack;
    my @tok = &translate($name, $clauses, \%cond);
    my ($indent) = ($lines[$d->{line} - 1] =~ m/^(\s*)/o);
    my @new = &layout(\@tok, $d->{last} - $d->{line} + 1, $indent);
    splice(@lines, $d->{line} - 1, scalar(@new), @new);
  }
  return (\@lines, \@msg, $eol);
}

my ($check, $pairs, $summary, $output);
my $usage = "usage: acc2omp.pl [-o out.F90] in.F90\n       acc2omp.pl --check [--pairs] [--summary] file.F90 ...\n";
GetOptions('check' => \$check, 'pairs' => \$pairs, 'summary' => \$summary, 'o=s' => \$output)
  or die $usage;

if (! $check)
{
  die $usage unless (@ARGV == 1 && ! $pairs && ! $summary);
  my $file = $ARGV[0];
  my ($lines, $msg, $eol) = &translateFile($file);
  if (@$msg)
  {
    for (sort { $a->[0] <=> $b->[0] } @$msg)
    {
      my ($ln, $m, $text) = @$_;
      $text =~ s/\s+/ /go; $text =~ s/^ //o;
      print STDERR "$file:$ln: error: $m: !\$ACC $text\n";
    }
    exit(1);
  }
  my $out = \*STDOUT;
  if ($output)
  {
    open($out, '>', "$output.tmp") or die "$output.tmp: $!\n";
  }
  print $out join("\n", @$lines), ($eol && @$lines ? "\n" : '');
  if ($output)
  {
    close($out) or die "$output.tmp: $!\n";
    rename("$output.tmp", $output) or die "$output: $!\n";
  }
  exit(0);
}

die $usage unless @ARGV;

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
