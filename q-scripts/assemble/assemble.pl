use strict;
use warnings;

my $html = slurpFile('html');
my $css  = slurpFile('css');
my $js   = slurpFile('js');
sub slurpFile {
    my ($ext) = @_;
    open my $fh, '<', "$ENV{SLAVES_DIR}/assemble/assemble.$ext" or die;
    $/ = undef;
    my $data = <$fh>;
    close $fh;
    return $data;
}
$html =~ s/__ASSEMBLE_CSS__/$css/;
$html =~ s/__ASSEMBLE_JS__/$js/;

my @li;
my @smp = split(" ", $ENV{SAMPLES});
my @als = split(" ", $ENV{ALIASES});
foreach my $i(sort {$als[$a] cmp $als[$b]} 0..$#smp){
    push @li, "<li data-sample=\"$smp[$i]\"  data-alias=\"$als[$i]\">$als[$i]</li>"  
}
my $li = join("", @li);
$html =~ s/__SAMPLE_LIST__/$li/;

print $html;
