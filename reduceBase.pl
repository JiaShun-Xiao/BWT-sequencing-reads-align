#!/usr/bin/perl
open FILE, '/2_disk/longyk/human/chr21.fasta';##換行的。
open FILE2, '>/2_disk/longyk/human/redChr21.fa';
open FILE3, '>/2_disk/longyk/human/chr21Location.txt';
=pod
    open FILE, 'E:\naive.fa';
open FILE2, '>E:\rednaive.fa';
open FILE3, '>E:\naivelocation.txt';
=cut
$_ = <FILE>;
chomp $_;
$chrom = substr($_,1);
print FILE2 "$chrom\n";
my $i = 0;
while (<FILE>) {
    chomp;
    my $length = length $_;
    my $location = $i*$length;
    if($_ !~ /N/ && $_ !~ /[a-z]+/ ){
print FILE2 "$_";
print FILE3 "$location\n";
}
    $i++;
}
print "$i";

