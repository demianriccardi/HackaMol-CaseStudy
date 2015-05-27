#!/usr/bin/env perl
# Copyright (C) 2014  Demian Riccardi <demianriccardi@gmail.com>
use Modern::Perl;
use HackaMol;
use MCE::Loop max_workers => 4, chunk_size => 1;
use MCE::Subs qw( :worker );
use Time::HiRes qw(time);

my $type = shift or die 'pass nmrs xtals';

my $work = HackaMol->new(data=>$type, scratch=>"$type\_clusters");
my @pdbs = $work->data->children(qr/\.pdb/);
$work->scratch->mkpath unless $work->scratch->exists;

my $t1 = time;

mce_loop_s {

    my $fpdb  = $pdbs[$_];
    my $bld   = HackaMol->new( name => "local_builder" );

    my @atoms = grep { $_->Z != 1 } #no hydrogens
                $bld->read_file_atoms($fpdb);
    
    my @ss = $bld->find_disulfide_bonds(@atoms);

    unless (@ss){
      print STDERR $fpdb->basename(qr/\.pdb/) . " has no disulfide bonds\n";
    } 
    #loop over all disulfides
    my $i = 0;
    foreach my $ss (@ss) {

        my @cys_s = $ss->all_atoms;
        my @cut5  = grep {
                 $cys_s[0]->distance($_) <= 5.0
              or $cys_s[1]->distance($_) <= 5.0
        } @atoms;

        # throw out if disorder nearby
        next if ( grep { $_->occ != 1.0 } @cut5 );

        #create filter to pull all atoms in cut5 residues
        my %resid;
        foreach my $at (@cut5) {
            $resid{ $at->resid }{ $at->resname }{ $at->chain }++;
        }
        #apply filter
        my @bigcut =
          grep { exists( $resid{ $_->resid }{ $_->resname }{ $_->chain } ) }
          @atoms;

        #create molecule with bigcut
        my $mol = HackaMol::Molecule->new(
            atoms => [@bigcut],
        );

        #add a mercury atom between the sulfurs for convenience
        my $hg = HackaMol::Atom->new(
            name        => 'HG2',
            Z           => 80,
            resname     => "HG2",
            record_name => "HETATM",
            serial      => $bigcut[-1]->serial + 1,
        );

        #t loops are useful for NMR structures with multiple models
        foreach ( 0 .. $mol->tmax ) {

            # setting t from mol sets t for all atoms, including ss
            $mol->t($_);

            # fill Hg with coordinates from SS configurations
            $hg->push_coords( $ss->COM );
        }
        $mol->push_atoms($hg);    # Hg is now added to mol

        $mol->t(0);               #return to first t

        my $fname = $work->scratch->child($fpdb->basename(qr/\.pdb/))->absolute;
        #first print returns the filehandle for future writing
        my $fh = $mol->print_pdb_ts([0 .. $mol->tmax], $fname . "_$i.pdb" );
        $i++;
    }
} 0 , $#pdbs;

my $t2 = time;

printf( "calculation time: %10.1f\n", $t2 - $t1 );
