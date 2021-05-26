# Replication code for "Strongmen Cry Too: The Effect of Aerial Bombing on Voting for The Incumbent in Competitive Autocracies"
# Milos Popovic
# 2021/5/25

# load libraries
library(plyr, quietly=T)
library(tidyr, quietly=T)
library(gridExtra, quietly=T)
library(grid, quietly=T)
library(ggplot2, quietly=T)
library(RCurl, quietly=T)
library(ggcorrplot, quietly=T)
library(ggrepel, quietly=T)
library(scales, quietly=T)
library(marmap, quietly=T)
library(ggeffects, quietly=T)
library(raster, quietly=T)
library(gdalUtils, quietly=T)
library(rgdal, quietly=T)
library(reshape2, quietly=T)
library(ascii, quietly=T)
library(stargazer, quietly=T)

################################################
#                                              #
#                  ANALYSIS                    #
################################################

# load function for clustered robust SEs
# source: https://raw.githubusercontent.com/IsidoreBeautrelet/economictheoryblog/master/robust_summary.R
source("robust_summary.r")

# load dataset
a <- read.csv("did_elections_final.csv")
b <- a

#tranform vars
b$pop_per_doctor <-  b$pop_per_doctor/1000
b$employed_per_100000inh <-  b$employed_per_1000inh/100
b$vote_milosevic <- b$vote_milosevic/100
b$vote_srs <- b$vote_srs/100
b$vote_opposition <- b$vote_opposition/100
b$turnout <- b$turnout/100
b$okrug_id <-as.factor(b$okrug_id)
b$year92 <- NULL
b$year92[b$year==1992] <- 1
b$year92[b$year!=1992] <- 0
b$year96 <- NULL
b$year96[b$year==1996] <- 1
b$year96[b$year!=1996] <- 0

# Base model, DV = vote-share for Milosevic
# Table 1, column 1
m1a <- lm(vote_milosevic ~ okrug_id + year96 +
          bombed*year00,
     data=b)

summary(m1a, cluster = c("id")) # clustered robust SEs

# Full model: DV = vote-share for Milosevic
# Table 1, column 2
m1b <- lm(vote_milosevic ~ okrug_id + year96 +
          bombed*year00 + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91,
     data=b)
summary(m1b, cluster = c("id")) # clustered robust SEs

# Placebo test (take only 1992 and 1996): 
# Table 1, column 3
c <- subset(b, year!=2000) # keep only 1992 and 1996

m2 <- lm(vote_milosevic ~ okrug_id + 
          bombed*year96 + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91,
     data=c)
 summary(m2, cluster = c("id")) # clustered robust SEs

# Table 1, column 4: Did the bombing favor the radicals?
# DV = vote-share for the radicals
m3 <- lm(vote_srs ~ okrug_id + year96 + 
          year00*bombed + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91,
     data=b)
 summary(m3, cluster = c("id")) # clustered robust SEs

# Table 2, column 1 
# Alternative 1: Did the bombing affect voter turnout?
# DV = turnout
m4 <- lm(turnout ~ okrug_id + year96 + 
          year00*bombed + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91,
     data=b)
 summary(m4, cluster = c("id")) # clustered robust SEs

 # Table 2, column 2 
 # did the bombing influence turnout among Milosevic's voters?
 # DV = votes for Milosevic * turnout
 b$vote_milosevic_turn <- b$vote_milosevic * b$turnout

 m5 <- lm(vote_milosevic_turn ~ okrug_id + year96 + 
          year00*bombed + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91,
     data=b)
 summary(m5, cluster = c("id")) # clustered robust SEs

# Table 2, column 3 
# Alternative 2: Did the electorate change affect vote-share?
# added change in population between 98 and 2002 as a covariate
m6 <- lm(vote_milosevic ~ okrug_id + year96 + 
          year00*bombed + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91 + pop_chng_98_02,
     data=b)
 summary(m6, cluster = c("id")) # clustered robust SEs

 # Table 2, column 4 
 # And what about the influx of displaced persons from Kosovo
 # added share of displaced persons from Kosovo in 2000 as a covariate
m7 <- lm(vote_milosevic ~ okrug_id + year96 + 
          year00*bombed + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91 + irp,
     data=b)
 summary(m7, cluster = c("id")) # clustered robust SEs

###############################################
#              FIGURES                        #
#             MAIN TEXT                       #
###############################################

# Figure 1
# load the bombing dataset
n <- read_csv("nato_bombing.csv")

# load names of major towns
grad <- read_csv("gradovi.csv")

# load shapefiles
# kosovo
ks <- readOGR(getwd(), "XKO_adm0", 
 verbose = TRUE, stringsAsFactors = FALSE)
# serbia proper - municipalities
rs <- readOGR(getwd(), "SRB_adm2", 
 verbose = TRUE, stringsAsFactors = FALSE)
# serbia proper
rs0 <- readOGR(getwd(), "SRB_adm0", 
 verbose = TRUE, stringsAsFactors = FALSE)
# montenegro
me <- readOGR(getwd(), "MNE_adm0", 
 verbose = TRUE, stringsAsFactors = FALSE)

# merge shapefiels
sc <- rbind(rs0, ks, me)
ss <- fortify(sc, "id")

# download SRTM data with spatial resolution at 30 meter on the line of the equator
# source: https://srtm.csi.cgiar.org/srtmdata/
url1 <- "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_41_04.zip" # part1
url2 <- "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_40_04.zip" # part 2
url3 <- "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_41_03.zip" # part 3
url4 <- "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/srtm_40_03.zip" # part 4
download.file(url1, basename(url1), mode="wb")
download.file(url2, basename(url2), mode="wb")
download.file(url3, basename(url3), mode="wb")
download.file(url4, basename(url4), mode="wb")
unzip("srtm_40_04.zip") # unzip part1
unzip("srtm_41_04.zip") # unzip part2
unzip("srtm_40_03.zip") # unzip part3
unzip("srtm_41_03.zip") # unzip part4

# load into R and merge raster files
rasters <- c("srtm_40_04.tif", "srtm_41_04.tif",
"srtm_40_03.tif", "srtm_41_03.tif")
ras <- mosaic_rasters(gdalfile=rasters,
              dst_dataset="raster.tif",
              of="GTiff")
ras <- raster("raster.tif")
crs(ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# crop raster
j <- crop(x = ras, y = sc) # crop by FRY
jm <- mask(x=j, mask=sc)
jspdf <- as(jm, "SpatialPixelsDataFrame")
jj <- as.data.frame(jspdf)


# aggregate attacks by location
i <- ddply(n, "loc", 
  summarise, 
  attacks = length(loc), 
  lat = max(lat), 
  long = max(long))

#plot
g1 <- ggplot() +
geom_raster(data = jj, 
    aes_string(x = "x", 
               y = "y", 
               alpha = "raster")
               ) +
geom_path(data=rs, aes(x = long, 
              y = lat, 
              group = group),
              fill = NA, 
              color="#066598", 
              size=0.25,
              alpha=1) +
coord_fixed(ratio = 1.35) +
  geom_point(data=i, 
    aes(x=long, y=lat, size=attacks), 
    fill="#d11135", 
    col="#d11135", 
    alpha=0.65, 
    stroke=.5) +
  geom_point(data=grad, 
    aes(x=long, 
        y=lat), 
    fill="grey10", 
    col="grey10",
    pch=15,
    alpha=1,
    size=.75)+
  geom_text_repel(data=grad, 
    aes(x=long, 
      y=lat, 
      label = place),
    size = 3, 
    hjust = 0.5, 
    vjust=-0.25,
    col="grey10",
    alpha=1)+
  scale_size(name="Count",
    breaks=c(1,2,5,10,26), 
    range=c(1,8), 
    limits=c(min(i$attacks), max(i$attacks)))+
  scale_alpha(name = "", 
    range = c(.15, .85), 
    guide = F) +
guides(fill=F,
size= guide_legend(override.aes = list(alpha=1),
            direction = "vertical",
            title.position = 'top',
            title.hjust = 0.5,
            label.hjust = 0,
            nrow = 6,
            byrow = T,
            reverse = F,
            label.position = "right",
      family="georg")
  ) +
  theme_minimal() +
    theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title=element_text(size=11, face="bold"),
  plot.title = element_text(size=14, hjust=.5),
  plot.caption = element_text(size=11, hjust=.65, vjust=5),
    legend.position = c(.8, .8),
    legend.text = element_text(size=10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    panel.border = element_blank()
    ) +
labs(x=NULL, y=NULL,title="",
caption="")
# export
ggsave(filename="fig11.png", width=7, height=8, dpi = 600, device='png', g1)

  # Figure 2
  rs$bombed <- as.factor(rs$bombed)
  rs$bombed <- factor(rs$bombed, levels=c("1", "0"))

  f <- as.data.frame(rs)
  r <- fortify(rs, region = "id")
  d <- r %>% left_join(f, by = "id")

g2 <- ggplot() +
      geom_polygon(data=d, aes(x = long, 
              y = lat, fill=bombed, 
              group = group)) +
      geom_path(data=d, aes(x = long, 
              y = lat, 
              group = group), color="white", size=0.25)+
      geom_path(data=ks, aes(x = long, 
              y = lat, 
              group = group), color="grey60", size=0.25) +
      geom_path(data=me, aes(x = long, 
              y = lat, 
              group = group), color="grey60", size=0.25) +
      coord_map() +
      scale_fill_manual(values=c('#d11135', '#066598'),
                labels=c("bombed", "not bombed"),
                name="") +
      guides(fill=guide_legend(
            direction = "vertical",
            keyheight = unit(5, units = "mm"),
            keywidth = unit(5, units = "mm"),
            title.position = 'top',
            title.hjust = 0.5,
            label.hjust = 0,
            nrow = 2,
            byrow = T,
            reverse = F,
            label.position = "right"
          )
      ) +
      theme_minimal() +
      theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.title=element_text(size=7),
            plot.title = element_text(size=14, hjust=.25),
            legend.position = c(.8, .9),
            legend.text = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            legend.background = element_blank(),
            plot.margin=unit(c(0,0,-10,0),"mm"),
            panel.border = element_blank()
      ) +
labs(x=NULL, y=NULL,title="",
caption="")

ggsave("fig2.png", height=8, width=7, g2)

# Figure 3
did <- a
df <- did[,c("id", "okrug_id", "year", "vote_milosevic", "bombed")]

# parallel trends assumption
df1 <- df
df1$year <- as.factor(df1$year)
df1$year <- relevel(df1$year, ref = 1)

m0 <- lm(vote_milosevic ~ okrug_id + 
  bombed*year, 
data = df1)

summary(m0, cluster = c("id")) # clustered robust SEs

pr <- as.data.frame(ggpredict(m0, c("bombed", "year"),
ci.lvl = 0.95
))
pr <- pr[,-3]
names(pr) <- c("bombed", "predicted", "conf.low",  "conf.high", "year")       
pr$year <- as.numeric(levels(pr$year))[pr$year]
pr$bombed <- factor(pr$bombed, levels=c(1, 0))

arrow = arrow(length = unit(1.5, "mm"), angle=30, type = "closed")

p2 <- ggplot(data=pr, aes(x=year, 
  y=predicted, 
  col=bombed, 
  group=bombed)) +
geom_line(size=1.25) +
geom_pointrange(aes(ymin = conf.low,
                    ymax = conf.high), 
                    lwd = 0.55, 
                    position = position_dodge(width = 0)) +
geom_vline(xintercept = 1999, 
           colour = "grey30", 
           lty = 2) +
scale_color_manual(values=c('#d11135', '#066598'),
                labels=c("bombed", "not bombed"),
                name="Municipality") +
scale_y_continuous(limits=c(25, 55),
  breaks=c(25, 35, 45, 55), 
  labels=c("25%", "35%", "45%", "55%")) +
scale_x_continuous(breaks=c(1992, 1996, 1999, 2000)) +
annotate(geom="curve", 
  xend=1999, 
  x=1998.5, 
  y=52, 
  yend=50, 
  curvature=.2, 
  colour="grey30", 
  arrow = arrow)+
geom_text(
    x = 1997.5,
    y = 52.5,
    inherit.aes = FALSE,
    label = "NATO bombing",
    check_overlap = TRUE,
    hjust = 0.5,
    fontface = 'bold',
    colour="grey30",
    family = "georg",
    size = 4
  ) +
  guides(group=F,
  col = guide_legend(override.aes = list(size=1, linetype=0, alpha=1))) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major=element_blank(),
        panel.border=element_blank(),
        axis.title.x = element_text(size=11,colour="grey20", face="bold"),
        axis.title.y = element_text(size=11, colour="grey20", face="bold", angle=0, vjust=.5),
        axis.line   = element_line(colour=NA),
        axis.line.x = element_line(colour="grey20"),
        axis.line.y = element_line(colour="grey20"),
        axis.text.x = element_text(size=11,colour="grey20", hjust = 1),
        axis.text.y = element_text(size=11,colour="grey20"),
        plot.title = element_text(size=13,colour="grey20", hjust=0.5),          
        plot.caption = element_text(size=8,colour="grey20", hjust=0), 
        legend.title = element_text(size=12,colour="grey20", hjust=0.5, face="bold"), 
        legend.key=element_rect(fill="white", color="white"),
        legend.text = element_text(size=11,colour="grey20"),
        legend.position=c(.2, .8)) +
labs(title = "",
        x="Year", 
              y="Predicted\nShare of Votes\nfor SPS/JUL")
#print(p2)

ggsave(file="fig3.png", p2)

###############################################
#              TABLE 1                        #
#             MAIN TEXT                       #
###############################################

t <- b[,c(6:10, 12:16)]
stargazer(t, summary.stat = c("n", "mean", "sd", "min", "max"), 
             type="text")

###############################################
#              TABLES                         #
#             APPENDIX                        #
###############################################

# Table A1
ct <- a[,c(10, 12:16, 19:21)]

ct1 <- ddply(ct, "bombed", summarise,
  employed_per_1000inh=mean(employed_per_1000inh, na.rm=T),
  pop_per_doctor=mean(pop_per_doctor, na.rm=T),
  refugees_perc_pop=mean(refugees_perc_pop, na.rm=T),
  minorpc91=mean(minorpc91, na.rm=T),
  femalepc91=mean(femalepc91, na.rm=T),
  urbanpc91=mean(urbanpc91, na.rm=T),
  age20_34pc=mean(age20_34pc, na.rm=T),
  illiterate_pc=mean(illiterate_pc, na.rm=T)) 

ct2 <- melt(ct1, id.vars=c("employed_per_1000inh", 
  "pop_per_doctor", 
  "refugees_perc_pop", 
  "minorpc91", 
  "femalepc91", 
  "urbanpc91",
  "age20_34pc", 
  "illiterate_pc"))
ct3 <- as.data.frame(t(ct2))

ascii(ct3)

# Table A2, column 1
 r1a <- lm(vote_opposition ~ okrug_id + year96 + 
          year00*bombed + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91,
     data=b)
# clustered robust SEs
 summary(r1a, cluster = c("id"))

 # Table A2, column 2
b$vote_opposition_turn <- b$vote_opposition * b$turnout

r1b <- lm(vote_opposition_turn ~ okrug_id + year96 + 
          year00*bombed + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91,
    data=b)
# clustered robust SEs
 summary(r1b, cluster = c("id"))

# Table A3
 r2 <- lm(vote_milosevic ~ okrug_id + year96 + 
          year00*allen_vincent_2011 + employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          minorpc91 + femalepc91,
     data=b)
# clustered robust SEs
summary(r2, cluster = c("id"))

 # Table A4
c1 <- subset(b, year==2000)
r3 <- glm(bombed ~ employed_per_100000inh +  
          pop_per_doctor + refugees_perc_pop +
          urbanpc91 + illiterate_pc,
     family="binomial",
     data=c1)
# clustered robust SEs
 summary(r3, cluster = c("id"))

 # Table A5
 r4 <- lm(vote_milosevic ~ okrug_id + year96 +
          bombed*year00 + employed_per_100000inh +
          refugees_perc_pop +
          minorpc91 + femalepc91,
     data=b)
 # clustered robust SEs
 summary(r4, cluster = c("id"))

###############################################
#              FIGURES                        #
#             APPENDIX                        #
###############################################

# Figure A1
cr <- t[,c(5:10)]
cr1 <- na.omit(cr)

colnames(cr1) <- c("Bombed", "Employed per 100,000",
 "Population per doctor",       "Refugees (%)",
     "Minorities (%)", "Females (%)")

corr <- round(cor(cr1), 3)

p <- ggcorrplot(corr, hc.order = TRUE, type = "lower",
   lab = TRUE,  
outline.col = "white",
colors = c("#BB4444", "#FFFFFF", "#4477AA"),
ggtheme = ggplot2::theme_classic)

png("fig5.png", res=600, height=7, width=7, units="in")
print(p)
dev.off()

# Figure A2
polls <- read.csv("polls.csv")
polls$period <- as.factor(polls$period)
polls$period <- factor(polls$period, levels=c("1992 (I)",
"1992 (II)", "1992 (III)", "1993",          
"1995", "1996", "1997 March", "1997 May",
"1998", "1999 September", "1999 November",   
"2000 February", "2000 March", "2000 August",
"2000 September"))

arrow = arrow(length = unit(1.5, "mm"), angle=30, type = "closed")

g3 <- ggplot(data=polls, 
  aes(x=ord, 
    y=perc)) +
geom_line(size=1.25, col='#066598', alpha=.5) + 
geom_point(size=2, col='#066598', alpha=1) +
geom_text_repel(data=polls, 
    aes(x=ord, 
    y=perc, 
    label = paste0(perc, "%")),
    size = 3, 
    hjust = 0.5, 
    fontface = 'bold',
    col="#1D02B0")+
geom_vline(xintercept = 9.5, 
  colour = "#d11135", 
  lty = 2) +
scale_y_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60), 
  labels=c("0%", "10%", "20%", "30%","40%", "50%", "60%")) +
scale_x_continuous(breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
  11, 12, 13, 14, 15),
  labels=c("1992\n(I)",
"1992\n(II)", "1992\n(III)", "1993",          
"1995", "1996", "1997\nMarch", "1997\nMay",
"1998", "1999\nSep", "1999\nNov",   
"2000\nFeb", "2000\nMar", "2000\nAug",
"2000\nSep")) +
annotate(geom="curve", 
  xend=9.5, 
  x=10.5, 
  y=44.5, 
  yend=40, 
  curvature=-.6, 
  colour="grey30", arrow = arrow)+
geom_text(
    x = 9.75,
    y = 45,
    inherit.aes = FALSE,
    label = "NATO bombing",
    check_overlap = TRUE,
    hjust = 0,
    fontface = 'bold',
    colour="grey40",
    family = "georg",
    size = 4
  ) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major=element_blank(),
        panel.border=element_blank(),
        axis.title.x = element_text(size=11,colour="grey20", face="bold", vjust=-2),
        axis.title.y = element_text(size=11, colour="grey20", face="bold", angle=0, vjust=.5),
        axis.line   = element_line(colour=NA),
        axis.line.x = element_line(colour="grey20"),
        axis.line.y = element_line(colour="grey20"),
        axis.text.x = element_text(size=9,colour="grey20", hjust = .5),
        axis.text.y = element_text(size=9,colour="grey20"),
        plot.title = element_text(size=13,colour="grey20", hjust=0.5),          
        plot.caption = element_text(size=8,colour="grey20", hjust=0), 
        legend.title = element_text(size=12,colour="grey20", hjust=0.5, face="bold"), 
        legend.key=element_rect(fill="white", color="white"),
        legend.text = element_text(size=11,colour="grey20"),
        legend.position="none") +
labs(title = "",
        x="Poll", 
              y="Would vote for\nSPS/JUL (%)")

ggsave("fig4.png", height=7.5, width=7, g3)