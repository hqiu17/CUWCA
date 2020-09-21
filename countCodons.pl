#! /usr/bin/perl 
# missing value are given value zero "0"
# 
$argm=join ("\*", @ARGV);
$cmarg=$#ARGV+1;
if ($argm =~ /\-help/i||$cmarg<1) {
   die "usage: XXX.pl 1_infile(fasta)\n";
} 

$infile = shift or die "input is wrong \?";

open (IN, "$infile") or die "cof";

open (OUTg, ">$infile.af")or die "cof";

$inseq="";
$status=0;
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


while (<IN>) {
  if (/>/) {
     $fseq="";
     if ($inseq) {
         $fseq = &RUSC ($inseq);
         print OUTg $def,"\t", $fseq, "\n";
     }
     $inseq="";
     s/\s|>//g;
     $def=$_;
  } else {
     s/\s//g;
     $inseq .= $_;
  }
#last if $.==100;

}

if ($inseq) {
    $fseq = &RUSC ($inseq);
    print OUTg $def,"\t", $fseq, "\n";
}



sub RUSC {
  my $seq=shift;
  my $len=length ($seq);
  my $codons=$len/3;
  die "Sequence length is divisible by three.\n" if $len%3;
  
  my ($ttt,$ttc,$tta,$ttg,$tct,$tcc,$tca,$tcg,$tat,$tac,$taa,$tag,$tgt,
      $tgc,$tga,$tgg,$ctt,$ctc,$cta,$ctg,$cct,$ccc,$cca,$ccg,$cat,$cac,
      $caa,$cag,$cgt,$cgc,$cga,$cgg,$att,$atc,$ata,$atg,$act,$acc,$aca,
      $acg,$aat,$aac,$aaa,$aag,$agt,$agc,$aga,$agg,$gtt,$gtc,$gta,$gtg,
      $gct,$gcc,$gca,$gcg,$gat,$gac,$gaa,$gag,$ggt,$ggc,$gga,$ggg)=0;
  my ($phe,$leu,$ser,$tyr,$ter,$cys,$trp,$pro,$his,$gln,$arg,$ile,$met,
      $thr,$asn,$lys,$arg,$val,$ala,$asp,$glu,$gly)=0;
  
  for my $i (1..$codons) {
    my $index=$i*3-3;
    my $cdn=substr ($seq, $index, 3);
    if ($cdn=~/ttt/i) {++$ttt;++$phe;}
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
  my $index=($codons * 64);
  
  if ($phe) {                                   
    $vector .= ($ttt)/$index; $vector.="\t";
    $vector .= ($ttc)/$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}
  
  if ($leu) {
    $vector .= ($tta) /$index; $vector.="\t";
    $vector .= ($ttg) /$index; $vector.="\t";
    $vector .= ($ctt) /$index; $vector.="\t";
    $vector .= ($ctc) /$index; $vector.="\t";
    $vector .= ($cta) /$index; $vector.="\t";
    $vector .= ($ctg) /$index; $vector.="\t"; 
  } else {$vector .= "0\t0\t0\t0\t0\t0\t";}
  
  if ($ser) {
    $vector .= ($tct) /$index; $vector.="\t"; #ser
    $vector .= ($tcc) /$index; $vector.="\t";
    $vector .= ($tca) /$index; $vector.="\t";
    $vector .= ($tcg) /$index; $vector.="\t";
    $vector .= ($agt) /$index; $vector.="\t"; #ser
    $vector .= ($agc) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t0\t0\t0\t0\t";}
  
  if ($tyr) {
    $vector .= ($tat) /$index; $vector.="\t"; #
    $vector .= ($tac) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}

  if ($cys) {
    $vector .= ($tgt) /$index; $vector.="\t"; #cys
    $vector .= ($tgc) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}

  if ($pro) {
    $vector .= ($cct) /$index; $vector.="\t";
    $vector .= ($ccc) /$index; $vector.="\t";
    $vector .= ($cca) /$index; $vector.="\t";
    $vector .= ($ccg) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t0\t0\t";}

  if ($his) {
    $vector .= ($cat) /$index; $vector.="\t";
    $vector .= ($cac) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}
  
  if ($gln) {
    $vector .= ($caa) /$index; $vector.="\t";
    $vector .= ($cag) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}
  
  if ($arg) {
    $vector .= ($cgt) /$index; $vector.="\t"; # arg
    $vector .= ($cgc) /$index; $vector.="\t";
    $vector .= ($cga) /$index; $vector.="\t";
    $vector .= ($cgg) /$index; $vector.="\t";
    $vector .= ($aga) /$index; $vector.="\t";
    $vector .= ($agg) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t0\t0\t0\t0\t";}
  
  if ($ile) {
    $vector .= ($att) /$index; $vector.="\t"; # ile
    $vector .= ($atc) /$index; $vector.="\t";
    $vector .= ($ata) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t0\t";}

  if ($thr) {
    $vector .= ($act) /$index; $vector.="\t"; #thr
    $vector .= ($acc) /$index; $vector.="\t";
    $vector .= ($aca) /$index; $vector.="\t";
    $vector .= ($acg) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t0\t0\t";}
  
  if ($asn) {
    $vector .= ($aat) /$index; $vector.="\t"; #
    $vector .= ($aac) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}
  
  if ($lys) {
    $vector .= ($aaa) /$index; $vector.="\t";
    $vector .= ($aag) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}
  
  if ($val) {
    $vector .= ($gtt) /$index; $vector.="\t"; #val
    $vector .= ($gtc) /$index; $vector.="\t";
    $vector .= ($gta) /$index; $vector.="\t";
    $vector .= ($gtg) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t0\t0\t";}
  
  if ($ala) {
    $vector .= ($gct) /$index; $vector.="\t"; #ala
    $vector .= ($gcc) /$index; $vector.="\t";
    $vector .= ($gca) /$index; $vector.="\t";
    $vector .= ($gcg) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t0\t0\t";}
  
  if ($asp) {
    $vector .= ($gat) /$index; $vector.="\t"; #asp
    $vector .= ($gac) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}
  
  if ($glu) {
    $vector .= ($gaa) /$index; $vector.="\t"; #glu
    $vector .= ($gag) /$index; $vector.="\t";
  } else {$vector .= "0\t0\t";}
  
  if ($gly) {
    $vector .= ($ggt) /$index; $vector.="\t"; #
    $vector .= ($ggc) /$index; $vector.="\t";
    $vector .= ($gga) /$index; $vector.="\t";
    $vector .= ($ggg) /$index;#$vector.="\t";
  } else {$vector .= "0\t0\t0\t0";}
  
  return $vector;
}
