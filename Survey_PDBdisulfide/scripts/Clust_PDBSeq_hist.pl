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
my $work = HackaMol->new( data    => "$type\_clusters", 
                          scratch => "results/$type" );

my @pdbs = $work->data->children(qr/\.pdb/);

$work->scratch->mkpath unless $work->scratch->exists;

my @aas = mce_map{
    my $fpdb = $_;
    my $hack = HackaMol->new( name => "hackitup" );
    my @aas  = map  {$_->aa321}
               grep {$_->name eq 'CA' } $hack->read_file_atoms($fpdb);
}@pdbs;

my %aa ;
$aa{$_}++ foreach @aas;
$aa{C} -= 2*scalar(@pdbs); # do not include the disulfide itself

my @a = qw/
L V I Y F W A P G M 
S T Q N H K R D E C
/;

my $t = 0;

$t += $aa{$_} foreach @a;


printf( "%8.3f\n", $aa{$_}/$t) foreach @a;
#printf( "%5s %8.6f %5i\n", $_, $aa1{$_}/$t,  $aa1{$_}) foreach @a;

