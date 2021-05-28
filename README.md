# Replication files for "Strongmen Cry Too: The Effect of Aerial Bombing on Voting for The Incumbent in Competitive Autocracies"

The NATO bombing of Yugoslavia, which lasted from March 24, 1999 until June 10, 1999, was the largest air campaign in Europe since the bombing of Britain and Germany in World War II. The air raids lasted for 78 days and hit 108 out of 160 municipalities, excluding Kosovo and Montenegro. The bombing was spread out and largely aimed at military barracks, industrial facilities, transportation networks, and communication lines. This repo provides a novel dataset with information on over 1,000 targets in the Federal Republic of Yugoslavia, including the date, location, target type, and fatalities. Included is also R code for the replication of my article "Strongmen Cry Too: The Effect of Aerial Bombing on Voting for The Incumbent in Competitive Autocracies" that was accepted for publication at Journal of Peace Research.

# Files

In a nutshell, the replication material includes the following folders:

- ``R`` includes code.r and robust_summary.r files. The code file includes the replication code, which I originally ran in R version 3.5.1 (2018-07-02) -- "Feather Spray". robust_summary.r is the function developed by [Isidore Beautrelet ](https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R), which computes clustered robust standard errors. You will need to run this function using a line of code in code.r to replicate my results.

- ``data`` encompasses three data frames in .csv format. First, did_elections_final.csv is the panel dataset for the replication of the figures and tables in the article and appendix. Second, gradovi.csv includes the names and geo coordinates in EPSG4326 (WGS84) that are used to produce Figure 1. Finally, polls.csv consists of aggregated polling results for Milosevic's party in the 1990s and come from the Institute of Social Science in Belgrade. The raw polling data and codebooks are available upon request.

- ``raw-data`` includes nato_bombing.csv, which is the dataset that includes information on the location, coordinates, and date of NATO bombing. The column type designates the nature of the target (if applicable) while descr column specifies the target subcategory (if applicable). The information on deaths from the bombing (deaths) and their categorization (death_stat) originates from [The Humanitarian Law Center's (HLC) data](http://www.hlc-rdc.org/db/nato_en/index.html). Finally the source of the bombing entry is provided under column source. 

This dataset comes with information on over 1,000 targets in the Federal Republic of Yugoslavia, including the date, location, target type, and fatalities. The majority of targets were military objects and forces (63%) followed by the industry (13%), transport infrastructure (9%), civilian (7%), communications facilities (7%), and other targets (1%). It is, by far, one of the most comprehensive and precise datasets on the NATO bombing of Yugoslavia. Other essential datasets include [Human Rights Data Analysis Group's dataset on killings in Kosovo](https://hrdag.org/kosovo-data-killings/), and [HLC's database of NATO bombing victims](http://www.hlc-rdc.org/db/nato_en/index.html). The dataset was manually coded and includes information on the location of bombings as reported in the media from March 24 until June 10, 1999. In specific, the information on bombed municipalities mainly comes from then pro-opposition Serbian daily (Glas Javnosti) and two major Serbian weeklies (NIN and Vreme). Reports from the state-owned news agency Tanjug, the Human Rights Watch and the Database on casualties of the HLC in Belgrade were used for data triangulation as well as the identification of under-reported strikes against the army.

While the dataset does not necessarily include the universe of all strikes, the level of under-reporting is restricted to targets in Kosovo and especially attacks on military forces. The under-reporting in Kosovo originates from limited media presence in the province since the outset of the conflict. Because Kosovo is excluded from the analysis, this issue does not affect the results of this study. The under-reporting of attacks on military units is due to restricted access to the placement of Yugoslav security forces in the field as well as their mobility in attempts to avoid NATO strikes. Some of these attacks were not reported by the media while others were reported but lacked information on the exact location. In the latter's case, the media vaguely referred to strikes in a ''wider area'' of a region. Such cases were omitted from the analysis. Fortunately, such occurrences were rare, less than 10 or equivalently less than 1\% of all strikes, and I was able to pin down a few unreported locations using the HLC database of casualties. Any remaining bias from the omission of strikes on military forces may not necessarily affect the inference because the focus of this article is on the bombing events that directly impact on the lives of voters. 

Below is the codebook for the NATO bombing dataset.

| column     | description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| loc        | target location denoted by<br>object name and/or settlement<br>name                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| long       | longitude of the target in EPSG:4326 (WGS84) coordinate                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| lat        | latitude of the target in EPSG:4326 (WGS84) coordinate                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| region     | denoting a wider area where the target is based including<br>Belgrade, Kosovo, Montenegro, Serbia (central) and Vojvodina                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| date       | date of the strike(s) on a target as reported in a source or sources.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| type       | target type based on the nature of the targeted object (not resulting <br>fatalities):<br><br><em>civilian</em> (civilian object like hospital, house or private farm) <br><em>communications</em> (media station, broadcast relay station)<br><em>government</em> (government objects such as ministry, excluding military and<br>police) <br><em>industry</em> (depot, factory, state-owned farm, plant, refinery)<br><em>infrastructure</em> (bridge, power grid, railway)<br><em>military facility</em> (airport, barracks, police station)<br><em>military units</em> (strikes against the military objects or <br>instances with unidentified objects)                                                                                                                                                                                                             |
| descr      | further description of the target if applicable                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| deaths     | fatality count resulting from the aerial strikes on the target, according to <br>HLC's database of NATO bombing victims                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| death_stat | following HLC classification, denoting whether a victim or victims were<br>civilian, military, and/or police                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| source     | online media source used to determine a strike or strikes on a target on a specific date. The <em>Glas Javnosti</em> newspaper was denoted using the link to its archived [main page](http://arhiva.glas-javnosti.rs/arhiva/1999/03/25/gv-01-1999-03-25.html) dedicated to the start of the bombing because link to the pages for specific dates sometimes failed to load or led back to this page. <em>Bela Knjiga</em> denotes the White Book of the Yugoslav government enlisting the targets and dates of strikes. <em>FHP</em> stands for the Humanitarian Law Center's database above. <em>HRW</em> refers to Human rights watch and their online reports. Additional sources include NIN, Vreme, book "Agresija NATO" by Smiljanic (referenced in the paper) as well as local online sources. While I massively relied on the Glas Javnosti reports, I aimed to use multiple sources where possible.  |

- ``shp`` includes shapefiles neccessary for the replciation of the maps. SRB_adm2.shp are municipality polygons for Serbia proper in the 1990s and I created them from the scratch using the latest settlement shapefile from [Geosrbija](https://geosrbija.rs/) and the 1991 census. Following list of municipalities and their settlements from the 1991 census, I aggregated the settlement level shapefile to the 1991 municipality level shapefile. Other shapefiles in this folder are administrative level 0 polygons for Montenegro (MNE_adm0), Kosovo (XKO_adm0) and Serbia proper (SRB_adm0), which I downloaded from [OSM Geofabrik](https://download.geofabrik.de/europe.html).

- ``output`` contains plots that are part of the publication or online appendix.

- There is also ``article.pdf`` in the root directory, which is the last version of the paper that was conditionally accepted by Journal of Peace Research. At page 35, you'll find the online supplementary file that includes information on the datasets used such as the codebook for did_elections_final.csv as well as additional robustness tests.

# Recommended citation

If any feature from the original NATO bombing dataset (``nato_bombing.csv``) or the dataset used for the analysis (``did_elections_final.csv``) appears in your project, please use the citations as follows:

Milos Popovic. (2021). Replication files for: Strongmen Cry Too: The Effect of Aerial Bombing on Voting for The Incumbent in Competitive Autocracies [Data set]. Zenodo. http://doi.org/10.5281/zenodo.4814191

Milos Popovic. (2021). Strongmen Cry Too: The Effect of Aerial Bombing on Voting for The Incumbent in Competitive Autocracies. Journal of Peace Research, accepted & forthcooming


![alt text](https://github.com/milos-agathon/strongmen_replication/blob/main/bombing.png?raw=true)

