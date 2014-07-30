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
    my @SScor = CysCys_intcoords($ss, $mol->all_atoms);
    push @lcor, @SScor unless (grep {/MISSING/} @SScor);
    #push @lcor, [@SScor] unless (grep {/MISSING/} @SScor);
  }
  ($name,\@lcor);
} @pdbs;

my $fyaml = $work->scratch->child("$type\_distances_angles.yaml");
DumpFile($fyaml,\%coords);

my $t2 = time;

printf("Time %10.2f\n",$t2-$t1);

sub CysCys_intcoords{
# Depends on convention of PDB ordering 
#0    N CYS     
#1   CA CYS    
#2    C CYS     
#3    O CYS     
#4   CB CYS    
#5   SG CYS    
#
  my $ss    = shift;
  my $sa    = $ss->get_atoms(0);
  my $sb    = $ss->get_atoms(1);

  my $q_intra = 0;
  $q_intra = abs($sa->resid - $sb->resid)-1  if ($sa->chain eq $sb->chain);

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

  my $hack   = HackaMol->new(name=>"builder");

  my ($CACA)   = $hack->build_bonds(  $cysa[1], $cysb[1]  ) ;
  my ($CACBS1) = $hack->build_angles( @cysa[1,4,5]        ) ;
  my ($CACBS2) = $hack->build_angles( @cysb[1,4,5]        ) ;
  my ($CBS1S2) = $hack->build_angles( @cysa[4,5],$cysb[5] ) ;
  my ($CBS2S1) = $hack->build_angles( @cysb[4,5],$cysa[5] ) ;

  return (  
          { distance => {
              ss   => $ss->bond_length   ,
              caca => $CACA->bond_length , 
            },
            angle => {
              cacbs1 => $CACBS1->ang_deg   , 
              cacbs2 => $CACBS2->ang_deg   ,  
              cbs1s2 => $CBS1S2->ang_deg   , 
              cbs2s1 => $CBS2S1->ang_deg   , 
            },
            chain => dseq_inter_intra($sa,$sb),
          }
  );
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

