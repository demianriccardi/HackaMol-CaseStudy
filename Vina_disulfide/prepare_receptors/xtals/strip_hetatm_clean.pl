use Modern::Perl;
use HackaMol;
use Bio::PDB::Structure; # use this to write out a pdb that MGLTools likes
use List::Util qw(min);

foreach my $file (glob ("*.pdb")){

  my $mol = HackaMol::Molecule -> new (
      atoms=>[
              grep {$_->altloc ne 'B'}
              grep {
                   $_->record_name ne 'HETATM'
                   } HackaMol -> new -> read_file_atoms($file)
             ]
  );

  $mol->print_pdb($file); 

  my $mol1= Bio::PDB::Structure::Molecule->new;
  $mol1->read($file);
  unlink($file);
  $mol1->print($file);
}
