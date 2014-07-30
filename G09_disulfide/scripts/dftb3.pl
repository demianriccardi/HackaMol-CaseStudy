#!/usr/bin/env perl
# DMR June 3, 2014
#
#   perl examples/dftd3.pl
#
# reads reads examples/xyzs
# and runs dftd3 xyz -func b3pw91 -bj
# and pulls the dispersion energy
#
# uses the DFT-D3 program
# by Stefan Grimme, Jens Antony, Stephan Ehrlich, and Helge Krieg
# J. Chem. Phys. 132, 154104 (2010); DOI:10.1063/1.3382344
#
# with BJ-damping,
# Stefan Grimme, Stephan Ehrlich and Lars Goerigk
# J. Comput. Chem. 32, 1456 (2011); DOI:10.1002/jcc.21759
#
# wget http://www.thch.uni-bonn.de/tc/downloads/DFT-D3/data/dftd3.tgz
#
# cd into directory and install with fortran compiler  (e.g. gfortran or intel's ifort)
# may have to edit make file to use whichever compiler you have

use Modern::Perl;
use HackaMol;
use HackaMol::X::Calculator;
use Path::Tiny;

my $eshft = shift || 0;
my $hack = HackaMol->new( data => "xyz_opts", );
my $i = 0;

my $scratch = path('tmp');

foreach my $xyz (grep {!/symbol/} $hack->data->children(qr/\.xyz$/) ) {
    my $mol = $hack->read_file_mol($xyz);

    my $Calc = HackaMol::X::Calculator->new(
        mol        => $mol,
        scratch    => $scratch,
        in_fn      => "bah$i.xyz",
        out_fn     => "calc-$i.out",
        map_in     => \&input_map,
        map_out    => \&output_map,
        exe        => '~/bin/dftd3',
        exe_endops => '-func b3pw91 -bj',

    );
    $Calc->map_input;
    $Calc->capture_sys_command;

    my $energy = $Calc->map_output(627.51);

    printf( "%-10.6f\n", $energy+$eshft );
    #printf( "%25s : %10.6f\n", $xyz->basename, $energy );

    $i++;

}

$scratch->remove_tree;

#  our functions to map molec info to input and from output
sub input_map {
    my $calc = shift;
    $calc->mol->print_xyz( $calc->in_fn );
}

sub output_map {
    my $calc = shift;
    my $conv = shift;
    my $re   = qr/-\d+.\d+/;
    my ($energy)  = $calc->out_fn->slurp =~ m /Edisp \/kcal,au:\s+${re}\s+(${re})/;
    return ( $energy * $conv );
}

