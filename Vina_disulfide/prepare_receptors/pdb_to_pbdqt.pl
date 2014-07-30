use Modern::Perl;
use MCE::Flow max_workers => 8, chunk_size => 1;
use MCE::Subs qw( :worker );
use File::chdir;
use File::Slurp;
use Capture::Tiny ':all';

my @pdbs = glob("xtals/*.pdb");
my $exe = "/Library/MGLTools/latest/bin/pythonsh prepare_receptor4.py -r";

mce_flow sub {

    my ($self, $chunk_ref, $chunk_id) = @_;
    my $pdb = $chunk_ref->[0];
    mce_say "$exe $pdb -A bonds_hydrogens";
    my ($stdout, $stderr, $exit) = capture { 
      system("$exe $pdb -A bonds_hydrogens");
    }

}, @pdbs;


