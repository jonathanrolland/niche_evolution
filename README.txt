The script called "script_alt_lat_shared.r" is the central script to run the analysis. It contains all the functions related to the inference of ancestral characters and the code to run these functions on the empirical datasets.

The tables "birds_data.txt", "amphibians_data.txt", "squamates_data.txt" and "mammals_data.txt" contain the names of the species, the mean, minimum and maximum values of altitude (m), latitude (°) and temperature. Temperatures are in °C and are currently multiplied by 10.

The phylogenies called "birds.tre", "amphibians.tre", "mammals.tre", "squamates.tre" are the phylogenies from Jetz et al. 2012, Pyron & Wiens 2011, Pyron & Burbrink 2014, Bininda-emonds et al. 2007 (modified by Fritz and Khun et al. 2011), respectively. We modified these phylogenies to run the analyses in the best conditions, e.g. we removed the names of species for which altitude or latitude data were not available. See Materials and Methods section for more details.

The files "fossil_amphibia.txt", "fossil_birds.txt", "fossil_squamata.txt" and "fossil.mammals.txt" correspond to the raw data we downloaded from the paleobiology database website for each groups.

The files "fossil_tbl_amphibians_save.txt", "fossil_tbl_birds_save.txt", "fossil_tbl_mammals_save.txt", "fossil_tbl_squamata_save.txt", correspond to intermediate tables obtained after a first run of the script "script_alt_lat_shared.r" (see line 494 with "###"). These table contains the prior information necessary to constrain some nodes of the phylogenetic tree during the latitude reconstruction. As it is quite CPU intensive, we provide this to the users.


