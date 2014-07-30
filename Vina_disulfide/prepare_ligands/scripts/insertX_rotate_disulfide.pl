#!/usr/bin/env perl
# Demian Riccardi, June 4, 2014
#
# This script will optionally insert an atom into a disulfide 
# and optionally set it's dihedral angle.  Probably will not work
# for molecules with > 2 sulfur atoms
#
# add Hg and rotate to 90:
# perl insertX_rotate_disulfide.pl GSSG.pdb Hg 90 
#
# rotate to 90
# perl insertX_rotate_disulfide.pl GSSG.pdb XX 90 
#
# add Hg
# perl insertX_rotate_disulfide.pl GSSG.pdb Hg
#
use Modern::Perl;
use HackaMol;
use Path::Tiny;

my $bldr = new HackaMol;
my $file = path(shift);
my $sym  = shift || 'Hg';
my $rot  = shift; # optional argument to rotate the disulfide
my $name = $file->basename(qr/\.\w+/);
my $mol  = $bldr->read_file_mol($file);
my @atoms = $mol->all_atoms;

my @ss = $bldr->find_disulfide_bonds( @atoms );

# set up the dihedral
my @Ss = grep { $_->symbol eq 'S' } @atoms;
my @Cs = grep { $_->symbol eq 'C' } @atoms;

#find S-C bonds
my @SCs = $bldr->find_bonds_brute(
    bond_atoms => [@Ss],
    candidates => [@Cs],
    fudge      => 0.45,
);

#bond_atoms (S) are first in the group! wanted: C-S -- S-C
my ($dihe) =
  $bldr->build_dihedrals( reverse( $SCs[0]->all_atoms ), $SCs[1]->all_atoms );

##############################################################################
#          Find atoms to rotate about dihedral for scan: qrotatable          #
##############################################################################

my $init = {
    $dihe->get_atoms(1)->iatom => 1,
    $dihe->get_atoms(2)->iatom => 1,
};
my $atoms_rotate = qrotatable( $mol->atoms, $dihe->get_atoms(2)->iatom, $init );
delete $atoms_rotate->{$dihe->get_atoms(1)->iatom};

my $group_rotate =
  HackaMol::AtomGroup->new( atoms => [ @atoms[ keys %{$atoms_rotate} ] ] );

#$mol->print_xyz;
if (defined($rot)){
  $mol->dihedral_rotate_groups( $dihe, $dihe->dihe_deg - $rot, $group_rotate );
}
#$mol->print_xyz;
if ($sym ne 'XX'){
  my $atom = HackaMol::Atom->new(
                                 name=>$sym ,
                                 resname=>uc($sym).2, 
                                 symbol=>$sym,
                                 chain=>$Ss[0]->chain
                                );
  my $targ = 2*$Ss[0]->covalent_radius + $atom->covalent_radius + 1.3; # 1.3 is fudge to make S-Hg-S 4.7...
  $group_rotate->translate($ss[0]->bond_vector->versor*($targ-$ss[0]->bond_length));
  $atom->push_coords($ss[0]->COM);
  $mol->push_atoms($atom);
  $mol->fix_serial(1);
}
$mol->print_pdb;

sub qrotatable {
    my $atoms   = shift;
    my $iroot   = shift;
    my $visited = shift;

    $visited->{$iroot}++;

    my @cands;
    foreach my $at (@$atoms) {
        push @cands, $at unless ( grep { $at->iatom == $_ } keys %{$visited} );
    }

    #find S-C bonds
    my @bonds = $bldr->find_bonds_brute(
        bond_atoms => [ $atoms->[$iroot] ],
        candidates => [@cands],
        fudge      => 0.45,
    );

    foreach my $cand ( map { $_->get_atoms(1) } @bonds ) {
        next if $visited->{ $cand->iatom };
        my $visited = qrotatable( $atoms, $cand->iatom, $visited );
    }
    return ($visited);
}

