#! /usr/bin/perl -w

use strict;

my $thorn = $ARGV[0];
print "Post-processing thorn $thorn\n";

my $filename = "$thorn/src/Boundaries.c";

open (FILE, "< $filename") or die "Cannot read file \"$filename\"";
my @lines = <FILE>;
close FILE or die;

#my $n=0;
#foreach my $line (@lines) {
#    print "$n   $line";
#    ++$n;
#}

$lines[33] = <<EOF;
  ierr = Boundary_SelectGroupForBC (cctkGH, CCTK_ALL_FACES, boundary_width, -1, "${thorn}::ML_Jadm", "scalar");
  if (ierr<0) CCTK_WARN (CCTK_WARN_ABORT, "Failed to select boundary condition for ${thorn}::ML_Jadm");

  ierr = Boundary_SelectGroupForBC (cctkGH, CCTK_ALL_FACES, boundary_width, -1, "${thorn}::ML_Madm", "scalar");
  if (ierr<0) CCTK_WARN (CCTK_WARN_ABORT, "Failed to select boundary condition for ${thorn}::ML_Madm");
EOF

#my $n=0;
#foreach my $line (@lines) {
#    print "$n   $line";
#    ++$n;
#}

open (FILE, "> $filename") or die "Cannot write file \"$filename\"";
print FILE join '', @lines;
close FILE or die;
