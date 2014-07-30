#!/usr/bin/env perl
# Copyright (C) 2014  Demian Riccardi <demianriccardi@gmail.com>
#
use Modern::Perl;
use HackaMol;
use YAML::XS qw(Dump DumpFile);
use MCE::Map;
use Time::HiRes qw(time);

my $type = shift or die "pass type: xtals nmrs";
my $work = HackaMol->new(data=>"$type\_clusters", scratch=>"results/$type");
my @pdbs = $work->data->children(qr/\.pdb/);
$work->scratch->mkpath unless $work->scratch->exists;

my $t1 = time;

my @tchain = mce_map {
    my $fpdb = $_;
    my $hack = HackaMol->new( name => "hackitup" );
    my @atoms   = $hack->read_file_atoms($fpdb);
    my ($hg)    = grep {$_->symbol eq 'Hg' } @atoms; 
    my ($ss)    = grep { abs($_->COM-$hg->xyz) <= 0.01 } #only the one we want
                  $hack->find_disulfide_bonds(@atoms);
    my ($s1,$s2)= $ss->all_atoms;
    my $chain1  = $s1->chain;
    my $chain2  = $s2->chain;
    my $dres    = abs($s1->resid-$s2->resid);

    my $key = 'intra';
    my $val = [$fpdb->basename('.pdb'), $dres];

    if ($chain1 ne $chain2){
      $key = 'inter';
    }

    [$key,$val];

} @pdbs;

#convert array to hash for yaml dump
my %tchain;
push @{$tchain{$_->[0]}},$_->[1] foreach @tchain;

my $fyaml = $work->scratch->child("$type\_intra_inter.yaml");
DumpFile($fyaml,\%tchain);

my $t2 = time;

printf( "calculation time: %10.1f\n", $t2 - $t1 );
