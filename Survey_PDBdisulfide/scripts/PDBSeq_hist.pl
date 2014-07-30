#!/usr/bin/env perl
# Copyright (C) 2014  Demian Riccardi <demianriccardi@gmail.com>
#
use Modern::Perl;
use HackaMol;
use YAML::XS qw(Dump DumpFile);
use MCE::Map;
use Time::HiRes qw(time);

my $t1 = time;

my $type = shift or die "pass type: xtals nmrs";
my $work = HackaMol->new( data => $type, scratch => "results/$type" );
my @pdbs = $work->data->children(qr/\.pdb/);
$work->scratch->mkpath unless $work->scratch->exists;

my @aas = mce_map { split }
          mce_map {
                map  { unpack "x17A53" }
                grep { /^SEQRES/ } $_->lines;
          } @pdbs;

my %aa;
$aa{$_}++ foreach @aas;

my %aa1 = map {(HackaMol::Atom->new(Z=>1,resname=>$_)->aa321, $aa{$_})} keys %aa;

my @a = qw/
L V I Y F W A P G M 
S T Q N H K R D E C
/;

my $t = 0;
$t += $aa1{$_} foreach @a;
printf( "%8.3f\n", $aa1{$_}/$t) foreach @a;
#printf( "%5s %8.6f %5i\n", $_, $aa1{$_}/$t,  $aa1{$_}) foreach @a;

