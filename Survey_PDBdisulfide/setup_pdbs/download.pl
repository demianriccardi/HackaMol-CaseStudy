#!/usr/bin/env perl
# Copyright (C) 2014  Demian Riccardi <demianriccardi@gmail.com>
# This script pulls a json file of pdbids and downloads those not
# found in the corresponding directory.  See script in setup_pdbs
# for converting list of pdbs (from pdb.org) to json file. 
#
use Modern::Perl;
use Path::Tiny;
use JSON::XS;
use YAML::XS;
use File::Slurp;
use List::Compare;
use File::chdir;

my $type = shift || 'xtals';
my $text    = read_file( "setup_pdbs/$type.json", { binmode => ':raw' } );

my $json    = new JSON::XS;
$json->incr_parse($text);
my $stor = $json->incr_parse;

my @dir_pdbs = map{$_->basename(qr/\.pdb/)} path($type)->children(qr/\.pdb/);

my $lc = List::Compare->new(\@dir_pdbs,$stor->{pdbids});

my @sync_list = $lc->get_complement;

$CWD = $type;

foreach my $pdbid (@sync_list){
  say ("wget http://pdb.org/pdb/files/$pdbid.pdb"); 
  system ("wget http://pdb.org/pdb/files/$pdbid.pdb"); 
}

