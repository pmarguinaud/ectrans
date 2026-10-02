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
# Add an explicit PRESENT(...) clause to every !$ACC PARALLEL LOOP that has
# no data clause other than PRESENT (an existing PRESENT list is
# completed), using the variables of the enclosing !$ACC DATA constructs
# that are referenced in the loop nest. Running it twice changes nothing.
#
#   add_present.pl [--write] [--conditional] file.F90 ...
#
# Without --write, only the report is printed. With --write, the files are
# modified in place. The clause is appended to the LAST line of the
# directive, so no line is added (line numbers and kernel names are kept).
#
# Report, one line per PARALLEL LOOP without data clause:
#   file:line: added PRESENT(...)            [notes]
# Notes:
#   kept implicit  the variable comes from a DATA construct with an IF clause:
#                it is NOT added (it may be absent when the loop runs),
#                unless --conditional is given (then noted "conditional")
#   nodefault    the loop has no DEFAULT(NONE): other variables may be
#                implicitly present/copied; listed as "unresolved" if they
#                are referenced with an index and are not in the list
#   nochange     no variable of the enclosing DATA is used: nothing added
#
use strict;
use warnings;
use Getopt::Long;

my ($write, $conditional);
GetOptions('write' => \$write, 'conditional' => \$conditional) or die "usage: add_present.pl [--write] [--conditional] file.F90 ...\n";

# Split a clause argument list at top-level commas
sub splitArgs
{
  my $s = shift;
  my (@a, $cur, $depth);
  ($cur, $depth) = ('', 0);
  for my $c (split(//, $s))
  {
    $depth++ if ($c eq '(');
    $depth-- if ($c eq ')');
    if ($c eq ',' && $depth == 0)
    {
      push @a, $cur; $cur = '';
      next;
    }
    $cur .= $c;
  }
  push @a, $cur if ($cur =~ m/\S/);
  s/^\s+|\s+$//g for @a;
  return @a;
}

# Clauses of a directive text: list of [NAME, args]
sub clauses
{
  my $text = shift;
  my @c;
  while ($text =~ m/\b(\w+)\s*(\((?:[^()]++|(?2))*\))?/g)
  {
    my ($n, $a) = (uc($1), $2);
    $a =~ s/^\(|\)$//g if defined($a);
    push @c, [$n, $a];
  }
  return @c;
}

my $nloop = 0;
for my $file (@ARGV)
{
  open(my $fh, '<', $file) or die "$file: $!\n";
  my @lines = <$fh>;
  close($fh);

  my @cpp;       # 1 if inside an OMPGPU block
  my @data;      # stack of [ [var names], if-clause or undef, line ]
  my $changed = 0;

  my $i = 0;
  while ($i <= $#lines)
  {
    my $s = $lines[$i];
    if ($s =~ m/^\s*#\s*(\w+)(.*)/)
    {
      my ($kw, $rest) = (lc($1), $2);
      push @cpp, ($rest =~ m/\bOMPGPU\b/ ? 1 : 0) if ($kw =~ m/^if/);
      pop @cpp if ($kw eq 'endif');
      $i++;
      next;
    }
    if ($s !~ m/^\s*!\$acc\b/i || grep { $_ } @cpp)
    {
      $i++;
      next;
    }

    # Join the directive
    my ($first, $last, $text) = ($i, $i, '');
    while (1)
    {
      (my $b = $lines[$last]) =~ s/^\s*!\$acc&?//i;
      $b =~ s/!.*$//; $b =~ s/\s+$//;
      my $cont = ($b =~ s/&$//);
      $text .= ' ' . $b;
      last unless $cont;
      $last++;
    }
    $text =~ s/\s+/ /g; $text =~ s/^ //;
    $i = $last + 1;

    if ($text =~ m/^DATA\b/i)
    {
      my (@v, $if);
      for my $c (&clauses($text))
      {
        my ($n, $a) = @$c;
        if ($n =~ m/^(PRESENT|COPY|COPYIN|COPYOUT|CREATE)$/)
        {
          push @v, map { (my $x = $_) =~ s/\s*\(.*//; $x } &splitArgs($a);
        }
        $if = $a if ($n eq 'IF');
      }
      push @data, [\@v, $if, $first + 1];
      next;
    }
    if ($text =~ m/^END\s+DATA\b/i)
    {
      die "$file:" . ($first + 1) . ": END DATA without DATA\n" unless @data;
      pop @data;
      next;
    }
    next unless ($text =~ m/^PARALLEL\s+LOOP\b/i);

    my @c = &clauses($text);
    # Loops with other data clauses than PRESENT are left alone
    next if (grep { $_->[0] =~ m/^(COPY|COPYIN|COPYOUT|CREATE|DEVICEPTR)$/ } @c);
    my %have = map { (uc($_) => 1) } map { &splitArgs($_->[1]) } grep { $_->[0] eq 'PRESENT' } @c;
    $nloop++;

    my %private = map { (uc($_) => 1) }
                  map { &splitArgs($_->[1]) } grep { $_->[0] =~ m/^(PRIVATE|FIRSTPRIVATE)$/ } @c;
    my $default = grep { $_->[0] eq 'DEFAULT' } @c;

    # Loop nest body: from the first DO after the directive to its END DO
    my ($j, $depth, $body) = ($last + 1, 0, '');
    for (; $j <= $#lines; $j++)
    {
      my $l = $lines[$j];
      next if ($l =~ m/^\s*#/ || $l =~ m/^\s*!/);
      (my $code = $l) =~ s/!.*$//;
      $depth++ if ($code =~ m/^\s*(\w+\s*:\s*)?DO\b/i);
      $depth-- if ($code =~ m/^\s*END\s*DO\b/i);
      $body .= $code;
      last if ($depth == 0);
    }
    die "$file:" . ($first + 1) . ": no DO loop after PARALLEL LOOP\n" unless $body =~ m/\bDO\b/i;

    # Identifiers of the body, not preceded by % (components)
    my %used;
    while ($body =~ m/(%\s*)?\b([A-Za-z]\w*)\b/g)
    {
      $used{uc($2)} = 1 unless $1;
    }

    my (@add, %seen, @notes);
    for my $d (@data)
    {
      my ($vars, $if, $dline) = @$d;
      for my $v (@$vars)
      {
        next if ($v =~ m/%/);
        next if ($seen{uc($v)}++);
        next unless ($used{uc($v)});
        next if ($private{uc($v)} || $have{uc($v)});
        if (defined($if))
        {
          if (! $conditional)
          {
            push @notes, "kept implicit: $v (DATA line $dline IF($if))";
            next;
          }
          push @notes, "conditional: $v (DATA line $dline IF($if))";
        }
        push @add, $v;
      }
    }

    if (! $default)
    {
      my %in = map { (uc($_) => 1) } @add;
      my %idx;
      while ($body =~ m/(%\s*)?\b([A-Za-z]\w*)\s*\(/g)
      {
        $idx{uc($2)} = 1 unless ($1 || $in{uc($2)} || $private{uc($2)});
      }
      my @unres = grep { !/^(MAX|MIN|ABS|SQRT|REAL|INT|NINT|MOD|SIZE|LBOUND|UBOUND|CMPLX|AIMAG|DBLE|EXP|LOG|SIGN|MERGE|HUGE|TINY|IF|PRESENT|TRANSFER|ATAN2|COS|SIN|FLOOR|CEILING)$/ } sort keys %idx;
      push @notes, "nodefault" . (@unres ? ", unresolved: @unres" : '');
    }

    my $where = "$file:" . ($first + 1);
    if (! @add)
    {
      print "$where: nochange @notes\n" if @notes;
      next;
    }
    my $list = join(',', @add);
    print "$where: " . (%have ? "completed PRESENT(...,$list)" : "added PRESENT($list)") . (@notes ? "   [@notes]" : '') . "\n";

    # Extend a PRESENT clause ending the last line of the directive, or
    # append a new one, before a trailing comment
    if ($lines[$last] =~ m/^(\s*!\$acc&?)(.*?)(\s*!.*)?(\r?\n)$/i)
    {
      my ($sent, $rest, $comment, $nl) = ($1, $2, $3 // '', $4);
      $rest =~ s/\s+$//;
      $rest .= " PRESENT($list)" unless ($rest =~ s/(\bPRESENT\s*\([^()]*)\)$/$1,$list)/i);
      $lines[$last] = "$sent$rest$comment$nl";
      $changed++;
    }
    else
    {
      die "$where: cannot edit line " . ($last + 1) . "\n";
    }
  }
  warn "$file: DATA not closed\n" if @data;

  if ($write && $changed)
  {
    open(my $out, '>', $file) or die "$file: $!\n";
    print $out @lines;
    close($out);
  }
}
print "# PARALLEL LOOPs examined (no data clause other than PRESENT): $nloop\n";
