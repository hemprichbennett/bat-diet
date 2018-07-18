# bat-diet
# bat-diet
This is the R-project for all of the code for chapter 4 of my PhD thesis, using DNA metabarcoding to look at the dietary ecology of bats in the rainforests of Sabah. All data will be published with the upcoming paper. Please direct any questions to me at hemprich.bennett@gmail.com

betalink analysis of beta-diversity uses the script 'scripts/r/betalink_analysis.R'

iNEXT analysis of the data uses the script 'scripts/r/all_bats_inext.R'

For the OTU analysis data was first generated using the script 'scripts/r/otu_conclusion_check_generation.R', then analysed with 'scripts/r/otu_conclusion_check_output_use.R'


The important rarefaction analysis using LOTUS is done using the script array_site_motu95_confidence_intervals.R, initiated on apocrita using in_site_otu_95.sh. The output is then visualised locally using ranges_analysis.R