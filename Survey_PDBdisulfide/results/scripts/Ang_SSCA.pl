use Modern::Perl;
use YAML::XS qw(LoadFile Dump);
use List::Util qw(max min sum);
use Math::SimpleHisto::XS;

my $type = shift or die "pass type: xtals nmrs";
my $yaml = LoadFile ("$type/$type\_distances_angles.yaml");

my $nbin = 50;
my $min  = 97;
my $max  = 112;

my @vals;
foreach my $cluster (keys %$yaml){
  foreach my $t (@{$yaml->{$cluster}}){
    push @vals, $t->{angle}{cbs1s2};
    push @vals, $t->{angle}{cbs2s1};
  }
}

my $df = ($max-$min)/ $nbin;
my $hist = Math::SimpleHisto::XS->new(
          bins => [ map{ $min + $df*$_ } (0 .. $nbin)],
);

$hist->fill(\@vals);
$hist->normalize($hist->total/scalar(@vals));
my $data_bins   = $hist->all_bin_contents;
my $bin_centers = $hist->bin_centers;

my $sum = 0;
foreach my $i (0 .. $#{$data_bins}){
  $sum += $data_bins->[$i];
  printf ("%10.4f %10.4f\n", $bin_centers->[$i], $data_bins->[$i]);
}

say "\n\nsum: $sum\n";
