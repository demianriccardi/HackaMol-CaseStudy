use Modern::Perl;
use YAML::XS qw(LoadFile Dump);
use List::Util qw(max min sum);
use Math::SimpleHisto::XS;

my $type = shift or die "pass type: xtals nmrs";
my $yaml = LoadFile ("$type/$type\_dihedrals_dses.yaml");

my $nbin = 50;
my $min  = -180;
my $max  = 180;

my @chi3;
my @chi2;
my @chi1;
foreach my $cluster (keys %$yaml){
  foreach my $t (@{$yaml->{$cluster}}){
    push @chi3, $t->{dihedral}{chi3};
    push @chi2, $t->{dihedral}{chi2};
    push @chi2, $t->{dihedral}{chi2p};
    push @chi1, $t->{dihedral}{chi1};
    push @chi1, $t->{dihedral}{chi1p};
  }
}

  my $df = ($max-$min)/ $nbin;

my $cnt = 3;
foreach my $dihe (\@chi3,\@chi2,\@chi1){
  say "BEGIN chi$cnt";
  my $hist = Math::SimpleHisto::XS->new(
          bins => [ map{ $min + $df*$_ } (0 .. $nbin)],
  );

  $hist->fill($dihe);
  $hist->normalize($hist->total/scalar(@$dihe));
  my $data_bins   = $hist->all_bin_contents;
  my $bin_centers = $hist->bin_centers;

  my $sum = 0;
  foreach my $i (0 .. $#{$data_bins}){
    $sum += $data_bins->[$i];
    printf ("%10.4f\n", $data_bins->[$i]);
    #printf ("%10.4f %10.4f\n", $bin_centers->[$i], $data_bins->[$i]);
  }

  say "\n\nsum: $sum\n";
  say "END chi$cnt";
  $cnt--;
}
