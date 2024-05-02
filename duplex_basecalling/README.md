When using the duplex command, two types of DNA sequence results will be produced: 'simplex' and 'duplex'. Any specific position in the DNA which is in a duplex read is also seen in two simplex strands (the template and complement). So, each DNA position which is duplex sequenced will be covered by a minimum of three separate readings in the output.

The dx tag in the BAM record for each read can be used to distinguish between simplex and duplex reads:

    dx:i:1 for duplex reads.
    dx:i:0 for simplex reads which don't have duplex offsprings.
    dx:i:-1 for simplex reads which have duplex offsprings.

Dorado will report the duplex rate as the number of nucleotides in the duplex basecalls multiplied by two and divided by the total number of nucleotides in the simplex basecalls. This value is a close approximation for the proportion of nucleotides which participated in a duplex basecall.
