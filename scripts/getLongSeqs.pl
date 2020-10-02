#! /usr/bin/perl
#getLongSeqs.pl

my $usage = "getLongSeqs.pl, retrieve sequence longer than specified ".
            "length cutoff\n".
            "usage: remove_short_seqs.pl 1:input file(fasta) 2:minimum ".
            " length 3: output file\n";

my $in    = shift or die "$usage";
my $minlen= shift or die "$usage";
my $out   = shift or die "$usage";

open (IN, "$in") or die "Cannot open input file: $in"; 
open (OUT, ">$out"); 

my %hash=(); 
my $gene;
my $seqno=0;
my $pass=0;
 
# Loop through input fasta file and store
# sequences into dictionary
while (<IN>){
   if (/\>(\S+)/){
       ++$seqno;
       $gene=$1;
       $hash{$gene}->{'num'}=$seqno;
       $hash{$gene}->{len}=0;
       $hash{$gene}->{seq}=();
   } else {
       next if /^\s*$/;
       $seqbk=$_;
       chomp;
       s/\-|\s|\r//g;   
       my $seq=$_;
       my $len = length ($seq);
       $hash{$gene}->{seq} .= $seqbk;
       $hash{$gene}->{len} += $len;
   }
}


# foreach (keys %hash) {
#     push ( @len_sum, $hash{$_}->{len} );
# }

# Loop through %hash and export sequence longer than cutoff 
foreach (sort {$hash{$a}->{'num'} <=> $hash{$b}->{'num'}} keys %hash) {
    if ($_ !~ /plugin/i) {
       if ($hash{$_}->{len} >= $minlen) {
          ++$pass;
          print OUT '>'.$_."\n",
            $hash{$_}->{seq};
       }
    }
 }

print "$pass out of $seqno seqs in input data are longer".
      " than $minlen (not including gaps).\n";
