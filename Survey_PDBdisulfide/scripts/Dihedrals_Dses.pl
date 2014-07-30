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
my $work = HackaMol->new(data=>"$type\_clusters", scratch=>"results/$type");
my @pdbs = $work->data->children(qr/\.pdb/);
$work->scratch->mkpath unless $work->scratch->exists;

my %coords = mce_map {

  my $fpdb = $_;
  my $name = $fpdb->basename('.pdb');
  my $mol   = $work->read_file_mol($fpdb);
  my ($hg)  = grep {$_->symbol eq 'Hg' } $mol->all_atoms;
  my ($ss)  = grep {abs($_->COM-$hg->xyz) <= 0.001}
              $work->find_disulfide_bonds($mol->all_atoms);
  my @lcor;

  foreach my $t (0 .. $mol->tmax){
    $mol->t($t);
    my @SScor = CysCys_dse($ss, $mol->all_atoms);
    push @lcor, @SScor unless (grep {/MISSING/} @SScor);
  }
  ($name,\@lcor);

} @pdbs;

my $fyaml = $work->scratch->child("$type\_dihedrals_dses.yaml");
DumpFile($fyaml,\%coords);

my $t2 = time;

printf("Time %10.2f\n",$t2-$t1);

sub CysCys_dse{
#0    N CYS     
#1   CA CYS    
#2    C CYS     
#3    O CYS     
#4   CB CYS    
#5   SG CYS    
#
#
# chi3
#4a 5a 5b 4b
#
# chi2 
#1a 4a 5a 5b
# chi2p 
#1b 4b 5b 5a
#
# chi1
#2a 1a 4a 5a
# chi1p
#2b 1b 4b 5b

  my $ss    = shift;
  my $sa    = $ss->get_atoms(0);
  my $sb    = $ss->get_atoms(1);

  my @atoms = @_;

  my @cysa  = grep{ 
                    $_->resid eq $sa->resid 
                and $_->chain eq $sa->chain
                  } @atoms;

  my @cysb  = grep{ 
                    $_->resid eq $sb->resid
                and $_->chain eq $sb->chain
                  } @atoms;

  return ("MISSING_ATOMS") if (scalar(@cysa) < 6 or
                               scalar(@cysb) < 6);  
  my $hack = HackaMol->new(name=>"builder");

  my @dihes; 
  push @dihes,  $hack->build_dihedrals( @cysa[2,1,4,5] ); #chi1
  push @dihes,  $hack->build_dihedrals( @cysb[2,1,4,5] ); #chi1'
  push @dihes,  $hack->build_dihedrals( ( @cysa[1,4,5], $cysb[5]   ) ); #chi2
  push @dihes,  $hack->build_dihedrals( ( @cysb[1,4,5], $cysa[5]   ) ); #chi2'
  push @dihes,  $hack->build_dihedrals( ( @cysa[4,5],   @cysb[5,4] ) ); #chi3
#  push @dihes,  $hack->build_dihedrals( ( @cysa[4,5],   @cysb[5,4] ) ); #chi3

  #$dihes[0]->dihe_fc(8.37)  ;  $dihes[0]->dihe_mult(3) ;  
  #$dihes[1]->dihe_fc(8.37)  ;  $dihes[1]->dihe_mult(3) ;
  #$dihes[2]->dihe_fc(4.18)  ;  $dihes[2]->dihe_mult(3) ;
  #$dihes[3]->dihe_fc(4.18)  ;  $dihes[3]->dihe_mult(3) ;
  $dihes[0]->dihe_fc(2)  ;  $dihes[0]->dihe_mult(3) ;
  $dihes[1]->dihe_fc(2)  ;  $dihes[1]->dihe_mult(3) ;
  $dihes[2]->dihe_fc(1)  ;  $dihes[2]->dihe_mult(3) ;
  $dihes[3]->dihe_fc(1)  ;  $dihes[3]->dihe_mult(3) ;
  $dihes[4]->torsion_efunc(\&chi3_torsion_efunc); #defined below

#  $dihes[4]->dihe_fc(14.64) ;  $dihes[4]->dihe_mult(2) ;
#  $dihes[5]->dihe_fc(2.51)  ;  $dihes[5]->dihe_mult(3) ;
  
  my $dse1 = $dihes[0]->torsion_energy + $dihes[1]->torsion_energy; 
  my $dse2 = $dihes[2]->torsion_energy + $dihes[3]->torsion_energy;
  my $dse3 = $dihes[4]->torsion_energy;# + $dihes[5]->torsion_energy;
  my $tdse = $dse1+$dse2+$dse3;

  return (  { 
              dihedral => { 
                chi1  => $dihes[0]->dihe_deg, 
                chi1p => $dihes[1]->dihe_deg, 
                chi2  => $dihes[2]->dihe_deg, 
                chi2p => $dihes[3]->dihe_deg, 
                chi3  => $dihes[4]->dihe_deg, 
              },
              dse => {
                total => $tdse,
                chi1  => $dse1,
                chi2  => $dse2,
                chi3  => $dse3,
              },
              chain => dseq_inter_intra($sa,$sb),
            }
         );
}

sub chi3_torsion_efunc{
  my $dihe   = shift;
  my $ang    = $dihe->dihe_rad;
  my $energy = 3.5*(1+cos(2*$ang)) + 0.6*(1+cos(3*$ang)) ;
  return $energy;
}

sub dseq_inter_intra {
  my ($sa,$sb) = (shift,shift);
  my $dseq = abs($sa->resid - $sb->resid);
  my ($qinter,$qintra) = (0,1);
  ($qinter,$qintra)    = (1,0) unless ($sa->chain eq $sb->chain);
  return (
          {
           dseq   => $dseq,
           qinter => $qinter,
           qintra => $qintra,
          },
  );
}
