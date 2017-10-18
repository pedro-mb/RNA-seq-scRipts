# RNA-seq-scRipts
The repository gathers a set of (very) basic scripts in R to analyse or parse data related to RNA-seq. 

--------------
getcoverage.R

After transcript assembly, with "stringTie --merge" or cuffmerge, the final annotation file (.GTF) can be compared to original genome annotation using gffcompare utility (https://github.com/gpertea/gffcompare). This will classify transcripts as they relate to reference transcripts and collapse contained transfrags (intron-redundant) using -C option. To further evaluate the coverage of the assembled transcripts, the new annotation file is used as reference to generate coverage data [using Tablemaker (for Cufflinks annotation, https://github.com/leekgroup/tablemaker) or StringTie with "-e" and "-B" options (for StringTie annotation)]. 

getcoverage.R contains an function [getCoverage(...)] that generates histograms for transcript coverage, grouping transcripts according to gffcompare classes ("=", "j", "e", "i", "o", "p", "s", "u", "x", "c", see http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/index.html )

USAGE:

coverageTable = getCoverage(whole_tx_table, transcript_cov,
                 "/path/for/ref.transcriptID.txt", "outputID", complete = TRUE, xmax = 8000, stp = 20)
                 
path_transcriptsID: transcriptsID is a file generated after parsing the ".tracking" output file from gffcompare using the collowing bash command:
$ cat cuffcmp-C.tracking | cut -f1,2,4 > cuffcmp-C.ref.transcriptID.txt

whole_tx_table and transcript_cov are dataframes obtained using ballgown, as below (copy from Ballgown tutorial) :
# library(ballgown)
# data_directory <- "stringTie_ballgown/stringTie-cafj10.ref-T10-comb" # path to coverage data obtained using stringtie
# pattern <- "ballgown-BGI_ILU_" 
# bg = ballgown(dataDir=data_directory, samplePattern=pattern, meas='all')
# transcript_cov = texpr(bg, 'cov')
# whole_tx_table = texpr(bg_filt, 'all') 

### xmax and stp are plotting options, defining length of x axis (xmax) and setp (stp)           
Default (complete = FALSE) run will generate histograms and cumulative density plots for most abundant classes of transcripts (conserved "="; novel isoforms "j"; unkown "u"; and generic exonic overlap "o")

To create histograms and cumulative density plots for all  transcript classes use "complete = FALSE"
