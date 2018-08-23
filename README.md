# bat-diet

This is the R-project for all of the code for chapter 4 of my PhD thesis, using DNA metabarcoding to look at the dietary ecology of bats in the rainforests of Sabah. All data will be published with the upcoming paper. Please direct any questions to me at hemprich.bennett@gmail.com

r_network_gen.r is a function used to make all the networks, using my old function The.matrix.reloader.R

betalink analysis of beta-diversity uses the script 'scripts/r/betalink_analysis.R'

iNEXT analysis of the data uses the script 'scripts/r/all_bats_inext.R'

For the OTU analysis data was first generated using the script 'scripts/r/otu_conclusion_check_generation.R', then analysed with 'scripts/r/otu_conclusion_check_output_use.R'


The important rarefaction analysis using LOTUS is done using the script array_site_motu95_confidence_intervals.R, initiated on apocrita using in_site_otu_95.sh. The output is then visualised locally using ranges_analysis.R

Degree analysis took place using all_species_sitewise_individual_analysis.R, plotting the degree of each bat and using a fixed effects model to look at how each site/habitat type/species varied.

netreducing occurs using the script 'netreducing.R', initiated on apocrita with the script 'in_a_netreducing.sh'

family_diets.R makes the family-level diet plots for each bat species

lulu_bold_querying.R is an array script ran on apocrita, initiated with in_array_bold_querying.sh to query BOLD for each MOTU, with the output then being analysed in bold_output_analysis.R. Or is it bold_querying.R?

all_species_sitewise_individual_analysis.R is it this that does corrs or neg_corrs.R

array_lulu.R is an array script, using LULU on datasets for each clustering level used, in_array_lulu.sh starts it

family_diets.R is for plotting the families

in_a_netreducing.sh is used to start netreducing script, which rarifies the networks, a separate script was used to start modularity (in_netreducing_JUST_MODULARITY.sh) because this array iteration was much more demanding than the others. The output is plotted using post_netreducing_plotting.R

in_metric_comparisons.sh

n_motu.R was used to count the MOTU in the networks at 95% similarity, post-lulu

sample_size.R was used to calculate the number of samples per site and species

shit_rates.R calculates the defecation rate of each bat species

sitewise_otu_conclusion_check_generation.R is the standard LOTUS metcalcs calculations for the network-level metrics


species_removal.R calculates how the networks change when we remove each species