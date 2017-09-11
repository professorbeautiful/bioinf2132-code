###   sequence logos  ()

##T. D. Schneider and R. M. Stephens. Sequence logos: a new way to display consensus sequences. Nucleic Acids Research, Vol. 18, no. 20, pp 6097-6100, 1990. (Also check out Tom Schneiders page.)
library(seqLogo)
sequenceData = 
"CCAGAGGCCCAACUGGUAAACGGGC
CCG-AAGCUCAACGGGAUAAUGAGC
CCG-AAGCCGAACGGGAAAACCGGC
CC-CAAGCGC-AGGGGAGAA-GCGC
CCG-ACGCCA-ACGGGAGAA-UGGC
CCGUUUUCAG-UCGGGAAAAACUGA
CCGUUACUCC-UCGGGAUAAAGGAG
CCGUAAGAGG-ACGGGAUAAACCUC
CCG-UAGGAG-GCGGGAUAU-CUCC
CCG--UGCCG-GCGGGAUAU-CGGC
CCG-AACUCG-ACGGGAUAA-CGAG
CCG--ACUCG--CGGGAUAA-CGAG"
sequenceData = strsplit(sequenceData, "\n")[[1]]
library(plyr)
sequenceData = laply(sequenceData, strsplit, split="")
apply(Reduce(rbind, sequenceData), 1, paste0, collapse="") ###Check: ok? yes.
sequenceData = Reduce(rbind, sequenceData)

sequenceDataTableList = apply(sequenceData, 2, table)
symbols = unique( unlist(sapply(sequenceDataTableList, names)) ) 

sequenceDataTable = matrix(NA, nrow=length(symbols), ncol = length(sequenceDataTableList))
rownames(sequenceDataTable) = symbols

for(sym in symbols) sequenceDataTable[sym, ] = colSums(sequenceData==sym)
# This package can't handle skips "-".
sequenceDataTable = sequenceDataTable[c("C","U","A","G"), ]
normalize = function(v) v/sum(v)
sequenceDataProbs = apply(sequenceDataTable, 2, normalize)
pwm = makePWM(sequenceDataProbs)
plot(pwm)
getMethod("plot", "pwm")
seqLogo

## Information content.

