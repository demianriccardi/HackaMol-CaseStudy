#!/usr/bin/env perl
# Demian Riccardi, June 3, 2014
use Modern::Perl;
use HackaMol;
use HackaMol::X::Calculator;
use lib 'basis_sets';
use YAML::XS;
use Path::Tiny;
use AugCcPv_t_z_pp;
use AugCcPv_td_z;

##############################################################
# load in the molecule and initialize charge and multiplicity#
##############################################################
my $bldr = new HackaMol;
my $file = path(shift); # or die "pass xyz/pdb file";
my $dang = shift or die "pass target dihedral";
my $name = $file->basename(qr/\..*/); 
my $mol  = $bldr->read_file_mol($file);
$mol->multiplicity(1);
$mol->push_charges(0);
$mol->fix_serial(1);
my @atoms = $mol->all_atoms;

##############################################################################
#                 Set up CS -- SC dihedral                                   #
# Designed for molecules with only two S atoms in disulfide or separated by  #
# atom in between:                                                           #
#    R-C-S.[x].S-C-R where [x] is optional atom                              #
##############################################################################

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
delete( $atoms_rotate->{ $Ss[0]->iatom } );
delete( $atoms_rotate->{ $Ss[1]->iatom } );

my $group_rotate =
  HackaMol::AtomGroup->new( atoms => [ @atoms[ keys %{$atoms_rotate} ] ] );


##############################################################################
#          Set up basis sets and input parameters for G09 input              #
##############################################################################
my $opt = 'opt=modredun freq'; # set to space if single point
my $thr = 'b3pw91';

my @basis =
  map { HackaMol::Atom->new( symbol => $_ ) } keys %{ $mol->bin_atoms };


foreach my $at (@basis) {
    if ( $at->symbol eq 'S' ) {
        $at->basis( $AugCcPv_td_z::basis{ $at->symbol } )
          if $AugCcPv_td_z::basis{ $at->symbol };

        $at->ecp( $AugCcPv_td_z::ecp{ $at->symbol } )
          if $AugCcPv_td_z::ecp{ $at->symbol };
    }
    else {
        $at->basis( $AugCcPv_t_z_pp::basis{ $at->symbol } )
          if $AugCcPv_t_z_pp::basis{ $at->symbol };

        $at->ecp( $AugCcPv_t_z_pp::ecp{ $at->symbol } )
          if $AugCcPv_t_z_pp::ecp{ $at->symbol };
    }
}

my $nbasis = 'aTZ';

#throw out the big basis above for optimization...
if ($opt =~ m/opt/){
  $_->basis('6-31+G**') foreach @basis;
  do{
     $_->basis('SDD'); 
     $_->ecp('SDD');
    } foreach grep {$_->has_ecp} @basis;
  $nbasis = '631+gss_opt';
}

# set the modred entry for gaussian
my @modred =
  ( sprintf( "D %i %i %i %i F", map { $_->serial } $dihe->all_atoms ) );


my %g09_param = (
    ppn     => 1,
    nodes   => 1,
    memory  => '250mw',
    job     => "#$thr/gen nosymm $opt test",
    Chk     => '0.chk',
    message => "sulfide dihedral set: $dang",
    basis   => \@basis,
    modred  => \@modred,
);

if ( grep { $_->has_ecp } @basis ) {
    $g09_param{job} = "#$thr/gen pseudo=read nosymm $opt test";
}

##############################################################################
#          Set up calculator for generating inputs                           #
##############################################################################

# instance of Calculator for setting up inputs
my $Calc = HackaMol::X::Calculator->new(
    mol     => $mol,
    scratch => "carver_SS.g09/gas/$thr",
    in_fn   => 'tmp.inp',
    map_in  => \&g09_input,    # pass anonymous subroutine
);


##############################################################################
#    print out input for a single  dihedral angle                            #
##############################################################################

    my $pname = sprintf("$name\_%03d\.$nbasis",$dang);
    $mol->dihedral_rotate_groups( $dihe, $dihe->dihe_deg - $dang,
        $group_rotate );
    $mol->print_xyz;
    $Calc->in_fn("$pname.inp");
    $g09_param{Chk}="$pname.chk";
    $Calc->map_input( \%g09_param );

exit;

sub g09_input {
    my $calc  = shift;
    my $param = shift;
    my $mol   = $calc->mol;

    my $fh = $calc->in_fn->openw_raw;
    print $fh "\%NprocShared = $param->{ppn}\n";
    print $fh "\%NprocLinda  = $param->{nodes}\n";
    print $fh "\%Mem = $param->{memory}\n";
    print $fh "\%Chk = $param->{Chk}\n";
    print $fh $param->{job} . "\n\n";
    print $fh $param->{message} . "\n\n";
    printf $fh ( "%-i %i \n", $mol->charge, $mol->multiplicity );
    printf $fh ( "%-3i %10.6f %10.6f %10.6f\n", $_->Z, @{ $_->xyz } )
      foreach $mol->all_atoms;
    print $fh "\n";

    if ( @{ $param->{modred} } ) {
        print $fh $_ . "\n" foreach @{ $param->{modred} };
    }

    print $fh "\n";

    foreach my $at ( @{ $param->{basis} } ) {
        print $fh $at->symbol . " 0 \n";
        print $fh $at->basis . "\n";
        print $fh "****\n";
    }
    print $fh "\n";
    foreach my $at ( grep { $_->has_ecp } @{ $param->{basis} } ) {
        print $fh $at->symbol . " 0 \n";
        print $fh $at->ecp . "\n";
    }
    print $fh "\n"

}

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
