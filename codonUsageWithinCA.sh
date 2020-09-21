#set -x
#################################################################
#                    codonUsageWithinCA                         #
# Analyze codon usage using within-CA (correspondence analysis) #
#################################################################

# Check command arguments
if [[ ! $# == 2 ]]; then
    echo "Usage:"
    echo "codonUsage_withinCA.sh codingSeqs(fasta) length_cutoff(integer)"
    exit 1
fi

# set
file=$1
lenCut=$2

# Examine input file
if [[ ! -e $file ]]; then
    echo "Input file ($file) does not exist"
    exit 1
fi

# get input file basename
if [[ $file =~ (.+).txt ]]; then
	fileName=${BASH_REMATCH[1]}
elif [[ $file =~ (.+).fasta ]]; then
	fileName=${BASH_REMATCH[1]}
else
	fileName=$file
fi

# Remove short sequences
longSequences=${fileName}_len${lenCut}
./getLongSeqs.pl $file $lenCut $longSequences >/dev/null
echo "Done with sequence filtering"

# Count codon frequencies for each sequence
./countCodons.pl $longSequences $longSequences.af
mv $longSequences $longSequences.fasta
echo "Done with counting codons"

# Do within-correspondent analysis using codon count table
mv $longSequences.af $longSequences
./withinCA.R $longSequences >/dev/null
echo "Done with correspondence analysis and plotting"

# remove intermediate files
rm $longSequences
