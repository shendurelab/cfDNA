Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA
==================================================================================

The following sections provide details for the analyses performed in:

Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA Comprises 
an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin. Cell. 2016 Jan 
14;164(1-2):57-68. doi: 10.1016/j.cell.2015.11.050. PubMed PMID: [26771485](http://www.ncbi.nlm.nih.gov/pubmed/26771485)

*All scripts and binaries are provided as is, without any warrenty and for use at your own risk. This is not the release of a software package. We are only providing this information and code in addition to a description of methods for making it easier to reproduce are analyses. We are* _not_ *providing any support for these scripts.* 

Primary data processing of sequencing data
------------------------------------------

Barcoded paired-end (PE) sequencing data was split allowing up to one substitution in the barcode sequence. Fragments shorter than or equal to the read length were consensus-called and adapter-trimmed. Remaining consensus single-end reads (SR) and the individual PE reads were aligned to the human reference genome (GRCh37, 1000 Genomes phase 2 technical reference, ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/) using the ALN algorithm in BWA v0.7.10 (Li and Durbin, 2010). PE reads were further processed with BWA SAMPE to resolve ambiguous placement of read pairs or to rescue missing alignments by a more sensitive alignment step around the location of one placed read end. Aligned SR and PE data were stored in BAM format using the samtools API (Li et al., 2009). BAM files for each sample were merged across lanes and sequencing runs. BAM files for each sample were deposited in the NCBI Gene Expression Omnibus (GEO) with accession GSE71378 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378).

Simulated reads and nucleotide frequencies
------------------------------------------

Aligned sequencing data was simulated (SR if shorter than 45 bp, PE 45 bp otherwise) for all major chromosomes of the human reference (GRC37h). Dinucleotide frequencies were determined from real data on both fragment ends and both strand orientations, and for the reference genome on both strands. The insert size distribution of the real data was extracted for the 1-500 bp range. Reads were simulated procedurally: at each step (i.e., at least once at each genomic coordinate, depending on desired coverage), (1) the strand is randomly chosen, (2) the ratio of the dinucleotide frequency in the real data to that in the reference sequence is used to randomly decide whether the initiating dinucleotide is considered, (3) an length is sampled from the insert size distribution, and (4) the frequency ratio of the terminal dinucleotide is used to randomly decide whether the generated alignment is reported. The simulated coverage was matched to that of the original data after PCR duplicate removal.



Analysis of nucleotide composition of 167 bp fragments
------------------------------------------------------

Fragments with inferred lengths of exactly 167 bp were filtered within samples to remove duplicates. Dinucleotide frequencies were calculated in a strand-aware manner, using a sliding 2 bp window and reference alleles at each position, beginning 50 bp upstream of one fragment endpoint and ending 50 bp downstream of the other endpoint. Observed dinucleotide frequencies at each position were compared to expected dinucleotide frequencies determined from a set of simulated reads reflecting the same cleavage biases calculated in a library-specific manner. 

Coverage, fragment endpoints, and windowed protection scores
------------------------------------------------------------

Fragment endpoint coordinates were extracted from BAM files with the SAMtools API. Both outer alignment coordinates of PE data were extracted for properly paired reads. Both end coordinates of SR alignments were extracted when PE data was collapsed to SR data by adapter trimming. A fragment’s coverage is defined as all positions between the two (inferred) fragment ends, inclusive of endpoints. We define the Windowed Protection Score (WPS) of a window of size k as the number of molecules spanning the window minus those with an endpoint within the window. We assign the determined WPS to the center of the window. For 35-80 bp fragments (short fraction, S-WPS), k=16; for 120-180 bp fragments (long fraction, L-WPS), k=120.

Short fraction:
```bash
./samtools view -u -m 120 -M 180 BAMFILE.bam $EXTREGION | ./FilterUniqueBAM.py -p | ./extractReadStartsFromBAM2Wig.py -p -r $REGION -w 120 -c COVERAGE.wig -n WPS.wig -s STARTS.wig
```

Long fraction:
```bash
./samtools view -u -m 35 -M 80 BAMFILE.bam $EXTREGION | ./FilterUniqueBAM.py -p | ./extractReadStartsFromBAM2Wig.py -p -r $REGION -w 16 -c COVERAGE.wig -n WPS.wig -s STARTS.wig
```

Nucleosome peak calling
-----------------------

L-WPS is locally adjusted to a running median of zero in 1 kb windows and smoothed using a Savitzky-Golay filter (Savitzky and Golay, 1964) (window size 21, 2nd order polynomial). The L-WPS track is then segmented into above-zero regions (allowing up to 5 consecutive positions below zero). If the resulting region is 50-150 bp, we identify the median L-WPS value of that region and search for the maximum-sum contiguous window above the median. We report the start, end and center coordinates of this window as the “peak,” or local maximum of nucleosome protection. All calculations involving distances between peaks are based on these center coordinates. A score for each peak is determined as the distance between maximum value in the window and the average of the two adjacent L-WPS minima neighboring the region. If the identified region is 150-450 bp, we apply the same above median contiguous window approach, but only report those windows that are 50-150 bp. For score calculation of multiple windows derived from 150-450 bp regions, we set the neighboring minima to zero. We discard regions <50 bp or >450 bp.

Analysis of TFBSs and genomic features
--------------------------------------

We started with clustered FIMO (motif-based) intervals (Grant et al., 2011; Maurano et al., 2012) defining a set of computationally predicted TFBSs. For a subset of clustered TFs (AP-2-2, AP-2, CTCF_Core-2, E2F-2, EBF1, Ebox-CACCTG, Ebox, ESR1, ETS, MAFK, MEF2A-2, MEF2A, MYC-MAX, PAX5-2, RUNX2, RUNX-AML, STAF-2, TCF-LEF, YY1), we retained only predicted TFBSs that overlap with ENCODE ChIP-seq peaks (TfbsClusteredV3 set downloaded from http://hgdownload.cse.ucsc.edu/
goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/).

WPS values for CH01 and the corresponding simulation were extracted for each position in a 5 kb window around the start coordinate of each TFBS, and was aggregated within each TF cluster. The mean WPS of the first and last 500 bp (which is predominantly flat and represents a mean offset) of the 5 kb window was subtracted from the original WPS at each position. For L-WPS only, a sliding window mean is calculated using a 200 bp window and subtracted from the original signal. Finally, the corrected WPS profile for the simulation is subtracted from the corrected WPS profile for CH01 to correct for signal that is a product of fragment length and ligation bias. This final profile is plotted and termed the Adjusted WPS. In figures, CTCF binding sites are shifted such that the zero coordinate on the x-axis is the center of its 52 bp binding footprint (Ong and Corces, 2014). Genomic coordinates of transcription start sites, transcription end sites, start codons, and splice donor and acceptor sites were obtained from Ensembl Build version 75. Adjusted WPS surrounding these features was calculated as for TFBSs.

Analysis of CTCF sites
----------------------

CTCF sites first included clustered FIMO binding site predictions (described above). This set was intersected with ENCODE ChIP-seq peaks (TfbsClusteredV3, described above), and then further intersected with a set of CTCF binding sites experimentally observed to be active across 19 tissues (Wang et al., 2012), to produce three increasingly stringent sets. For each CTCF site, distances between each of 20 flanking nucleosomes  (10 upstream and 10 downstream) were calculated. The mean S-WPS and L-WPS at each position relative to the center of the CTCF binding motif were also calculated within bins defined by spacing between -1 and +1 nucleosomes (>160 bp, 161-200 bp, 201-230 bp, 231-270 bp, 271-420 bp, 421-460 bp, and >460 bp). 

Analysis of DHS sites
---------------------

DHS peaks for 349 primary tissue and cell line samples were downloaded from http://www.uwencode.org/
proj/Science_Maurano_Humbert_et_al/data/all_fdr0.05_hot.tgz. Samples derived from fetal tissues (233 of 349) were removed from the analysis as they behaved inconsistently within tissue type, possibly because of unequal representation of cell types within each tissue sample. 116 DHS callsets from a variety of cell lineages were retained for analysis. For the midpoint of each DHS peak in a particular set, the nearest upstream and downstream calls in the CH01 callset were identified, and the distance between the centers of those two calls was calculated. The distribution of all such distances was visualized for each DHS peak callset using a smoothed density estimate calculated for distances between 0 and 500 bp.

Gene expression analysis
------------------------

###Human Protein Atlas

FPKM gene expression (GE) values measured for 20,344 Ensembl gene identifiers in 44 human cell lines and 32 primary tissues by the Human Protein Atlas (Uhlén et al., 2015) was downloaded from http://www.proteinatlas.org/download/rna.csv.zip. Genes with 3 or more non-zero expression values were retained (n=19,378 genes). The GE data set is provided with one decimal precision for the FPKM values. Thus, a zero GE value (0.0) indicates expression in the interval [0, 0.05) Unless otherwise noted, we set the minimum GE value to 0.04 FPKM before log2-transformation..

###Fast Fourier transformation (FFT) and smoothing of trajectories

We use parameters to smooth (3 bp Daniell smoother; moving average giving half weight to the end values) and de-trend the data (i.e. subtract the mean of the series and remove a linear trend). A recursive time series filter implemented in R was used to remove high frequency variation from trajectories. 24 filter frequencies (1/seq(5,100,4)) were used, and the first 24 values of the trajectory were taken as init values. The 24-value shift in the resulting trajectories was corrected by repeating the last 24 values of the trajectory.

###FFT intensity correlation with expression

L-WPS was used to calculate periodograms of genomic regions using Fast Fourier Transform (FFT, spec.pgram in R) with frequencies between 1/500 and 1/100 bases. Intensity values for the 120-280 bp frequency range were determined from smooth FFT periodograms. S-shaped Pearson correlation between GE values and FFT intensities was observed around the major inter-nucleosome distance peak, along with a pronounced negative correlation in the 193-199 bp frequency range. The mean intensity in this frequency range was correlated with the average intensity with log2-transformed GE values for downstream analysis. 

###Running the analysis

```bash

zcat transcriptAnno-GRCh37.75.tsv.gz | awk 'BEGIN{ FS="\t"; OFS="\t" }{ if ($5 == "+") { print $1,$2,$3-10000,$3,$5 } else { print $1,$2,$4,$4+10000,$5 } }' > transcriptAnno-GRCh37.75.upstream.tsv                                                   
zcat transcriptAnno-GRCh37.75.tsv.gz | awk 'BEGIN{ FS="\t"; OFS="\t" }{ if ($5 == "+") { print $1,$2,$4,$4+10000,$5 } else { print $1,$2,$3-10000,$3,$5 } }' > transcriptAnno-GRCh37.75.downstream.tsv                                                 
zcat transcriptAnno-GRCh37.75.tsv.gz | awk 'BEGIN{ FS="\t"; OFS="\t" }{ if ($5 == "+") { print $1,$2,$3-1,$3-1+10000,$5 } else { print $1,$2,$4-1-10000,$4-1,$5 } }' > transcriptAnno-GRCh37.75.body.tsv

# Extract counts for a sample
mkdir -p body/$SAMPLE
./extractReadStartsFromBAM_Region_WPS.py --minInsert=120 --maxInsert=180 -i transcriptAnno-GRCh37.75.body.tsv -o 'body/$SAMPLE/block_%s.tsv.gz' BAMFILE.bam

# Run FFT & convert to summary tbl
mkdir -p /tmp/body/$SAMPLE/fft
( cd body/$SAMPLE/counts; ls block_*.tsv.gz ) | xargs -n 500 Rscript fft_path.R body/$SAMPLE/counts /tmp/body/$SAMPLE/fft
mkdir -p body/fft_summaries/
./convert_files.py -a transcriptAnno-GRCh37.75.body.tsv -t /tmp/ -r  -p body -i $SAMPLE
rm -fR /tmp/body/$SAMPLE/fft

# Correlate the intensities with the expression data and generate PDFs in R
R --vanilla --quiet < plots.R
```
