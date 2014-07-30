use Modern::Perl;
use JSON::XS;
use File::Slurp;

my $file = shift;
my @pdbids = read_file($file);
chomp @pdbids;
my $json = $file =~ s/\.txt/\.json/r;
die "$json and $file are same" if $json eq $file;


my $hash = {pdbids => \@pdbids};
$hash->{date} = localtime();
write_file ($json, encode_json $hash);

