#!/usr/bin/env perl
# countCodons.pl
# A script to count codon for each sequence included in a fasta file
# (when codons are missing in a sequence, zero ("0") is assign.)

# Test command line argument
my $argm=join ("\*", @ARGV);
my $cmarg=$#ARGV+1;
my $usage = "Usage: countCodons.pl input_file(fasta)\n";
if ($argm =~ /\-help/i || $cmarg < 1) {
    die $usage;
}

# Get input file and create file handlers
$infile = shift or die $usage;
open (IN, "$infile") or die "Can't open file: $infile";
open (OUTg, ">$infile.af") or die "Can't write to output file: $infile.af";
# Write header to output file
print OUTg "\tttt\tttc\ttta\tttg\t",
             "ctt\tctc\tcta\tctg\t",
             "tct\ttcc\ttca\ttcg\t", 
             "agt\tagc\t",
             "tat\ttac\t",
             "tgt\ttgc\t",
             "cct\tccc\tcca\tccg\t",
             "cat\tcac\tcaa\tcag\t",
             "cgt\tcgc\tcga\tcgg\t", 
             "aga\tagg\t",
             "att\tatc\tata\t", 
             "act\tacc\taca\tacg\t", 
             "aat\taac\taaa\taag\t",
             "gtt\tgtc\tgta\tgtg\t",
             "gct\tgcc\tgca\tgcg\t",
             "gat\tgac\tgaa\tgag\t",
             "ggt\tggc\tgga\tggg\n";

# Loop through sequences and write out codon counts for each sequence
my $inseq  = "";
my $status = 0;
while (<IN>) {
    # Count codon and write to output file and reset relevant variables
    if (/>/) {
        $fseq="";
        if ($inseq) {
            $fseq = &Count_codon($inseq);
            print OUTg $def,"\t", $fseq, "\n";
        }
        $inseq="";
        s/\s|>//g;
        $def=$_;
    } else {
        s/\s//g;
        $inseq .= $_;
    }
}
# Count codon for the last sequence and write to output file
if ($inseq) {
    $fseq = &Count_codon ($inseq);
    print OUTg $def,"\t", $fseq, "\n";
}


sub Count_codon {
    # Count codons for a sequence
    # Argument: a string containing nucleotide sequence

    # Get input sequence and test if it is comprised of triplet codons
    my $seq=shift;
    my $len=length ($seq);
    my $codons=$len/3;
    die "Sequence length is not divisible by three.\n" if $len%3;
  
    # Initialize codon variables
    my ($ttt,$ttc,$tta,$ttg,$tct,$tcc,$tca,$tcg,$tat,$tac,$taa,$tag,$tgt,
        $tgc,$tga,$tgg,$ctt,$ctc,$cta,$ctg,$cct,$ccc,$cca,$ccg,$cat,$cac,
        $caa,$cag,$cgt,$cgc,$cga,$cgg,$att,$atc,$ata,$atg,$act,$acc,$aca,
        $acg,$aat,$aac,$aaa,$aag,$agt,$agc,$aga,$agg,$gtt,$gtc,$gta,$gtg,
        $gct,$gcc,$gca,$gcg,$gat,$gac,$gaa,$gag,$ggt,$ggc,$gga,$ggg
       ) = (0) x 64;
    my (
        $phe,$leu,$ser,$tyr,$ter,$cys,$trp,$pro,$his,$gln,$arg,$ile,$met,
        $thr,$asn,$lys,$arg,$val,$ala,$asp,$glu,$gly
       )=(0) x 22;
  
    # Loop through the length of sequence codon by codon and do the 
    # counting
    for my $i (1..$codons) {
        my $index=$i*3-3;
        my $cdn = substr($seq, $index, 3);
        
        if ($cdn=~/ttt/i)    {++$ttt;++$phe;}
        elsif ($cdn=~/ttc/i) {++$ttc;++$phe;}

        elsif ($cdn=~/tta/i) {++$tta;++$leu;} #leu
        elsif ($cdn=~/ttg/i) {++$ttg;++$leu;}
        elsif ($cdn=~/ctt/i) {++$ctt;++$leu;}
        elsif ($cdn=~/ctc/i) {++$ctc;++$leu;}
        elsif ($cdn=~/cta/i) {++$cta;++$leu;}
        elsif ($cdn=~/ctg/i) {++$ctg;++$leu;}

        elsif ($cdn=~/tct/i) {++$tct;++$ser;} #ser
        elsif ($cdn=~/tcc/i) {++$tcc;++$ser;}
        elsif ($cdn=~/tca/i) {++$tca;++$ser;}
        elsif ($cdn=~/tcg/i) {++$tcg;++$ser;}

        elsif ($cdn=~/cct/i) {++$cct;++$pro;} #pro
        elsif ($cdn=~/ccc/i) {++$ccc;++$pro;}
        elsif ($cdn=~/cca/i) {++$cca;++$pro;}
        elsif ($cdn=~/ccg/i) {++$ccg;++$pro;}

        elsif ($cdn=~/tat/i) {++$tat;++$tyr;} #tyr
        elsif ($cdn=~/tac/i) {++$tac;++$tyr;}
        elsif ($cdn=~/taa/i) {++$taa;++$ter;} #ter
        elsif ($cdn=~/tag/i) {++$tag;++$ter;}

        elsif ($cdn=~/tgt/i) {++$tgt;++$cys;} #cys
        elsif ($cdn=~/tgc/i) {++$tgc;++$cys;}
        elsif ($cdn=~/tga/i) {++$tga;++$ter;} #ter
        elsif ($cdn=~/tgg/i) {++$tgg;++$trp;} #trp

        elsif ($cdn=~/cat/i) {++$cat;++$his;} #his
        elsif ($cdn=~/cac/i) {++$cac;++$his;}
        elsif ($cdn=~/caa/i) {++$caa;++$gln;} #gln
        elsif ($cdn=~/cag/i) {++$cag;++$gln;}

        elsif ($cdn=~/cgt/i) {++$cgt;++$arg;} #arg
        elsif ($cdn=~/cgc/i) {++$cgc;++$arg;}
        elsif ($cdn=~/cga/i) {++$cga;++$arg;}
        elsif ($cdn=~/cgg/i) {++$cgg;++$arg;}

        elsif ($cdn=~/att/i) {++$att;++$ile;} #ile
        elsif ($cdn=~/atc/i) {++$atc;++$ile;}
        elsif ($cdn=~/ata/i) {++$ata;++$ile;}
        elsif ($cdn=~/atg/i) {++$atg;++$met;} #met

        elsif ($cdn=~/act/i) {++$act;++$thr;} #thr
        elsif ($cdn=~/acc/i) {++$acc;++$thr;}
        elsif ($cdn=~/aca/i) {++$aca;++$thr;}
        elsif ($cdn=~/acg/i) {++$acg;++$thr;}

        elsif ($cdn=~/aat/i) {++$aat;++$asn;} #asn
        elsif ($cdn=~/aac/i) {++$aac;++$asn;}
        elsif ($cdn=~/aaa/i) {++$aaa;++$lys;} #lys
        elsif ($cdn=~/aag/i) {++$aag;++$lys;}

        elsif ($cdn=~/agt/i) {++$agt;++$ser;} #ser
        elsif ($cdn=~/agc/i) {++$agc;++$ser;}
        elsif ($cdn=~/aga/i) {++$aga;++$arg;} #arg
        elsif ($cdn=~/agg/i) {++$agg;++$arg;}

        elsif ($cdn=~/gtt/i) {++$gtt;++$val;} #val
        elsif ($cdn=~/gtc/i) {++$gtc;++$val;}
        elsif ($cdn=~/gta/i) {++$gta;++$val;}
        elsif ($cdn=~/gtg/i) {++$gtg;++$val;}

        elsif ($cdn=~/gct/i) {++$gct;++$ala;} #ala
        elsif ($cdn=~/gcc/i) {++$gcc;++$ala;}
        elsif ($cdn=~/gca/i) {++$gca;++$ala;}
        elsif ($cdn=~/gcg/i) {++$gcg;++$ala;}

        elsif ($cdn=~/gat/i) {++$gat;++$asp;} #asp
        elsif ($cdn=~/gac/i) {++$gac;++$asp;}
        elsif ($cdn=~/gaa/i) {++$gaa;++$glu;} #glu
        elsif ($cdn=~/gag/i) {++$gag;++$glu;}

        elsif ($cdn=~/ggt/i) {++$ggt;++$gly;} #gly
        elsif ($cdn=~/ggc/i) {++$ggc;++$gly;}
        elsif ($cdn=~/gga/i) {++$gga;++$gly;}
        elsif ($cdn=~/ggg/i) {++$ggg;++$gly;}
    }
  
    my $vector="";        
    $vector .= ($ttt) . "\t";
    $vector .= ($ttc) . "\t";
    $vector .= ($tta) . "\t";
    $vector .= ($ttg) . "\t";
    $vector .= ($ctt) . "\t";
    $vector .= ($ctc) . "\t";
    $vector .= ($cta) . "\t";
    $vector .= ($ctg) . "\t"; 
    $vector .= ($tct) . "\t"; #ser
    $vector .= ($tcc) . "\t";
    $vector .= ($tca) . "\t";
    $vector .= ($tcg) . "\t";
    $vector .= ($agt) . "\t"; #ser
    $vector .= ($agc) . "\t";
    $vector .= ($tat) . "\t"; #
    $vector .= ($tac) . "\t";
    $vector .= ($tgt) . "\t"; #cys
    $vector .= ($tgc) . "\t";
    $vector .= ($cct) . "\t";
    $vector .= ($ccc) . "\t";
    $vector .= ($cca) . "\t";
    $vector .= ($ccg) . "\t";
    $vector .= ($cat) . "\t";
    $vector .= ($cac) . "\t";
    $vector .= ($caa) . "\t";
    $vector .= ($cag) . "\t";
    $vector .= ($cgt) . "\t"; # arg
    $vector .= ($cgc) . "\t";
    $vector .= ($cga) . "\t";
    $vector .= ($cgg) . "\t";
    $vector .= ($aga) . "\t";
    $vector .= ($agg) . "\t";
    $vector .= ($att) . "\t"; # ile
    $vector .= ($atc) . "\t";
    $vector .= ($ata) . "\t";
    $vector .= ($act) . "\t"; #thr
    $vector .= ($acc) . "\t";
    $vector .= ($aca) . "\t";
    $vector .= ($acg) . "\t";
    $vector .= ($aat) . "\t"; #
    $vector .= ($aac) . "\t";
    $vector .= ($aaa) . "\t";
    $vector .= ($aag) . "\t";
    $vector .= ($gtt) . "\t"; #val
    $vector .= ($gtc) . "\t";
    $vector .= ($gta) . "\t";
    $vector .= ($gtg) . "\t";
    $vector .= ($gct) . "\t"; #ala
    $vector .= ($gcc) . "\t";
    $vector .= ($gca) . "\t";
    $vector .= ($gcg) . "\t";
    $vector .= ($gat) . "\t"; #asp
    $vector .= ($gac) . "\t";
    $vector .= ($gaa) . "\t"; #glu
    $vector .= ($gag) . "\t";
    $vector .= ($ggt) . "\t"; #gly
    $vector .= ($ggc) . "\t";
    $vector .= ($gga) . "\t";
    $vector .= ($ggg) ;
  
    return $vector;
}
