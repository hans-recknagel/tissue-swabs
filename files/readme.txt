### R script
Proteus_swabs_script.R - Script for generating figures, tables and statistics for the analysis.

### DNA and sequencing Quality 
qubit_values_swabs_mod2.csv - DNA quantity from swab samples and Proteus length measurements.
Pi_tissue_vs_swabs_outlier.csv - Coverage and nucleotide diversity (Pi) estimates for sequenced swab and tissue samples.
coverage_vs_loci.csv - Coverage compared to number of loci and type of source material.

### Error rate estimation
replicates_66_p8_r1.plink.raw - Replicate samples genotype matrix generated from the unguided analysis.
swab_tissue_guided_p66_sorted.plink.raw - Replicate samples genotype matrix generated from the tissue-guided analysis.
proteus_refmap_p8_r1.raw - Replicate samples genotype matrix generated from the genome-guided analysis.

### Tissue-guided analysis vs genome-guided
swab_catalogue_loci_bt2_local_results.csv - Proportion of swab loci that mapped to the tissue catalogue or the refrence genome.

### Locus annotation
annotation_loci_blast.csv - Annotation of uniquely mapping RAD loci in the olm's genome and if they fall into gene or repetitive regions.
replicate_filtered_loci_annotation.csv - Annotation of filtered uniquely mapping RAD loci in the olm's genome and if they fall into gene or repetitive regions.

### exogenous DNA
superkingdoms.csv - Proportion of individual swab and tissue samples that mapped to different superkingdoms. Limited to those samples that did not map to the olm's genome.
bacteria_order.csv - Proportion of individual swab and tissue samples that mapped to different bacterial orders. Limited to those samples that did not map to the olm's genome.
chordata_class.csv - Proportion of individual swab and tissue samples that mapped to different chordate classes. Limited to those samples that did not map to the olm's genome.

### conspecific DNA
heterozygosity_swab_tissue_genome_ref_vcftools.csv - Heterozygosity of replicate swab and tissue samples 