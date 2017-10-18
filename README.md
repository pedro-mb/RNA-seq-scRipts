# RNA-seq-scRipts
The repository gathers a set of (very) basic scripts in R to analyse or parse data related to RNA-seq. 


...............
getcoverage.R
...............
After transcript assembly, with "stringTie --merge" or cuffmerge, the final annotation file (.GTF) can be compared to original genome annotation using gffcompare utility (\url{https://github.com/gpertea/gffcompare}). This will classify transcripts as they relate to reference transcripts and collapse contained transfrags (intron-redundant) using -C option. To further evaluate the coverage of the assembled transcripts, the new annotation file is used as reference to generate coverage data [using Tablemaker (for Cufflinks annotation, \url{https://github.com/leekgroup/tablemaker}) or StringTie with "-e" and "-B" options (for StringTie annotation)]. Afterwards, histograms for transcript coverage were obtained using a custom built R script.
