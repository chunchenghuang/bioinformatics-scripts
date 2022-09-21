{

if($3 == "gene"){gN+=1;gL+=($5-$4+1);gM=gL/gN}
else if($3 == "transcript"){tN+=1;tL+=($5-$4+1);tM=tL/tN}
else if($3 == "mRNA"){mN+=1;mL+=($5-$4+1);mM=mL/mN}
else if($3 == "start_codon"){stN+=1;stL+=($5-$4+1);stM=stL/stN}
else if($3 == "stop_codon"){spN+=1;spL+=($5-$4+1);spM=spL/spN}
else if($3 == "exon"){eN+=1;eL+=($5-$4+1);eM=eL/eN}
else if($3 == "CDS"){cN+=1;cL+=($5-$4+1);cM=cL/cN}
else if($3 == "intron"){iN+=1;iL+=($5-$4+1);iM=iL/iN}

}

END{

print "Genes: Total number = "gN"; Total length = "gL" bp; Mean length = "gM" bp"
print "Transcripts: Total number = "tN"; Total length = "tL"bp; Mean length = "tM" bp"
print "mRNAs: Total number = "mN", Total length = "mL" bp, Mean length = "mM" bp"
print "Start Codons: Total number = "stN", Total length = "stL" bp, Mean length = "stM" bp"
print "Stop Codons: Total number = "spN", Total length = "spL" bp, Mean length = "spM" bp"
print "Exons: Total number = "eN", Total length = "eL" bp, Mean length = "eM" bp"
print "CDSs: Total number = "cN", Total length = "cL" bp, Mean length = "cM" bp"
print "Introns: Total number = "iN", Total length = "iL" bp, Mean length = "iM" bp"

}