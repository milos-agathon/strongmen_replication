# Replication files for "Strongmen Cry Too: The Effect of Aerial Bombing on Voting for The Incumbent in Competitive Autocracies"
This repo provides data and code for replicating JPR article "Strongmen Cry Too: The Effect of Aerial Bombing on Voting for The Incumbent in Competitive Autocracies" that was accepted for publication at Journal of Peace Research.

The replication material includes the following folders:

- ``R`` includes code.r and robust_summary.r files. The code file includes the replication code, which I originally ran in R version 3.5.1 (2018-07-02) -- "Feather Spray". robust_summary.r is the function developed by [Isidore Beautrelet ](https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R), which computes clustered robust standard errors. You will need to run this function using a line of code in code.r to replicate my results.

- ``data`` encompasses three data frames in .csv format. First, did_elections_final.csv is the panel dataset for the replication of the figures and tables in the article and appendix. Second, gradovi.csv includes the names and geo coordinates in EPSG4326 (WGS84) that are used to produce Figure 1. Finally, polls.csv consists of aggregated polling results for Milosevic's party in the 1990s and come from the Institute of Social Science in Belgrade. The raw polling data and codebooks are available upon request.

- 
