#!/bin/bash
dir="."
##Running RepeatModeler
##Build a database for your fasta file first
BuildDatabase -name Myrmecridium_db -engine ncbi contigs.fa
##Now run the modeler to create classified repeats output
RepeatModeler -database Myrmecridium_db -engine ncbi -pa 20 1> Myrmecridium_RepeatModeler.out 2> Myrmecridium_RepeatModeler.err
##Combining repeat libraries (optional) 
##If you are working with a species that has several repeat families already annotated in the RepBase library
##you may want to combine the new RepeatModeler generated sequences with the RepBase library
##This was the case for me as there were several well-annotated, hand-curated, cichlid-specific SINEs and LINEs already in the RepBase library
##and they are probably more accurate than the ones that RepeatModeler created.
##Even if this isn't the case for your project, you might want to combine these repeat libraries as
##RepeatModeler may not detect older ancestral repeats that may be present in your genome of interest.
cat /opt/RepeatMasker/Libraries/RepeatMasker.lib $dir/RM*/consensi.fa.classified > combined_repeat_libs.fasta

##then run the actually repeat masker
RepeatMasker -pa 5 -s -gff -xsmall -dir $dir -lib $dir/combined_repeat_libs.fasta contigs.fa 1> Myrmecridium_RepeatMasker.out 2> Myrmecridium_RepeatMasker.err

##each RMBlast process uses 4 threads. So if you have 20 total cores available, only specify -pa 5 as shown here.
##-s is a more sensitive approach compared to default
##default uses hard masking but return Xs as Ns
##Hard mask: Masked sequence is converted to "X" (-x)
##Soft mask: Masked sequence is converted to lower-case ATCG (-xsmall)
##-gff creates an additional Gene Feature Finding format output