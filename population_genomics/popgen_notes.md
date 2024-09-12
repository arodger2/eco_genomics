# Coding and data notes for the population genomics module

## Author: Aly R

### 09-10-2024- Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS data from 3 regions (EU, NE, PNW) starting today with variant call format files (VCFs)

-   BASH has commands ex. ls- list

-   gz means g zip

-   zcat lets you look at file.gz\|head

-   we can specify the number of lines by giving head a number ex \|head4

-   CCACA is the barcode

-   In the fastq file, the 1st line is the sequence itself, second line is Q scores

-   Q scores range from 0-40+, letters like I and G mean high confidence in that base pair (mistakes every 1,000-10,000 bp)
