use Modern::Perl;
use YAML::XS qw(DumpFile);
use File::Slurp;

my $file = shift;
my @pdbids = read_file($file);
chomp @pdbids;
my $yaml = $file =~ s/\.txt/\.yaml/r;
die "$yaml and $file are same" if $yaml eq $file;
my $hash = {pdbids => \@pdbids};
$hash->{date} = localtime();
DumpFile($yaml, $hash);
