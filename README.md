# RNA-seq-scRipts
This repository gathers a set of (very) basic scripts in R to analyse or parse data related to RNA-seq. 

### getcoverage.R 

After transcript assembly, with "stringTie --merge" or cuffmerge, the final annotation file (.GTF) can be compared to original genome annotation using gffcompare utility (https://github.com/gpertea/gffcompare). This will classify transcripts as they relate to reference transcripts and collapse contained transfrags (intron-redundant) using -C option. To further evaluate the coverage of the assembled transcripts, the new annotation file is used as reference to generate coverage data [using Tablemaker (for Cufflinks annotation, https://github.com/leekgroup/tablemaker) or StringTie with "-e" and "-B" options (for StringTie annotation)]. 

getcoverage.R contains an function [getCoverage(...)] that generates histograms for transcript coverage, grouping transcripts according to gffcompare classes ("=", "j", "e", "i", "o", "p", "s", "u", "x", "c", see http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/index.html )

USAGE:

coverageTable = getCoverage(whole_tx_table, transcript_cov,
                 PATH, outputID , complete = TRUE, xmax = int, stp = int)
                 
* **whole_tx_table** and **transcript_cov** are data frames obtained using ballgown, see comment section in the script-

* **PATH** is a string with full path to the transcriptsID file, which needs to be created by parsing the ".tracking" output file from gffcompare. use the collowing bash command:
$ cat file.tracking | cut -f1,2,4 > transcriptID.txt

* **outputID** is a string with a user ID to be added to plots filenames

* **xmax** and **stp** are integers providing plotting options. They define length of x axis (xmax) and setp (stp) [default xmax = 4000, stp = 10]         

* **complete** is a boolean option (default FALSE). Default run will generate histograms and cumulative density plots for most abundant classes of transcripts (conserved "="; novel isoforms "j"; unkown "u"; and generic exonic overlap "o"). To create histograms and cumulative density plots for all  transcript classes use "complete = TRUE"

