# RNA-seq-scRipts
The repository gathers a set of (very) basic scripts in R to analyse or parse data related to RNA-seq. 

--------------
getcoverage.R

After transcript assembly, with "stringTie --merge" or cuffmerge, the final annotation file (.GTF) can be compared to original genome annotation using gffcompare utility (https://github.com/gpertea/gffcompare). This will classify transcripts as they relate to reference transcripts and collapse contained transfrags (intron-redundant) using -C option. To further evaluate the coverage of the assembled transcripts, the new annotation file is used as reference to generate coverage data [using Tablemaker (for Cufflinks annotation, https://github.com/leekgroup/tablemaker) or StringTie with "-e" and "-B" options (for StringTie annotation)]. 

getcoverage.R contains an function [getCoverage(...)] that generates histograms for transcript coverage, grouping transcripts according to gffcompare classes ("=", "j", "e", "i", "o", "p", "s", "u", "x", "c", see http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/index.html )

USAGE:
getCoverage(whole_tx_table, transcript_cov,
                 "/data/pedro/alternativeSplicing/gffcompare_run/cuffcmp-C.ref.transcriptID.txt",
                 "cuff-C.ref_repeat", complete = FALSE)

Default (complete = FALSE) run will generate histograms and cumulative density plots for most abundant classes of transcripts (conserved "="; novel isoforms "j"; unkown "u"; and generic exonic overlap "o")

To create histograms and cumulative density plots for all  transcript classes use "complete = FALSE"
