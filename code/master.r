# ---
# title: "Supplementary Materials"
# author: Ben Marwick
# date: Saturday, March 22, 2014
# ---


```{r setup_environment, message=FALSE, echo=FALSE}
# This is necessary to direct knitr to find the 
# 'code', 'data', and other directories that contain
# files needed to execute this document
knitr::opts_knit$set(root.dir=normalizePath('../'))
```


```{r load_libraries_and_data}


# check that dependencies are locally available
# source("code/packrat.r")

suppressMessages(source("code/libraries.r"))

# get versions of everything
# sessionInfo()

## load data

source("code/load.r") 
```


### Site chronology: analysing the dates


```{r chronology}

# plot dates and label with phases

dates$cal_age_min <- as.numeric(gsub("\\D", "", dates$cal_age_min))
dates$cal_age_max <- as.numeric(gsub("\\D", "", dates$cal_age_max))
dates$mid <- with(dates, cal_age_min + ((cal_age_max - cal_age_min) /2) )
ggplot(dates, aes(mid, depth_bs)) +
  geom_point() +
  scale_y_reverse() + 
  geom_errorbarh(aes(xmin=cal_age_min,xmax=cal_age_max,  height = 0)) +
  theme_minimal() +
  #geom_smooth(span = 0.19, se = FALSE) +
  ylab("depth below surface (m)") +
  xlab("calibrated age (BP)") +
  geom_vline(xintercept = 42000, colour = 'grey') +
  annotate("text", x = 40000, y = 0.25, label = "Phase I", angle = 90) + 
  geom_vline(xintercept = 38000, colour = 'grey') +
  geom_vline(xintercept = 17000, colour = 'grey') +
  annotate("text", x = 13000, y = 0.25, label = "Phase II", angle = 90) + 
  geom_vline(xintercept = 9000, colour = 'grey') +
  geom_vline(xintercept = 6500, colour = 'grey') +
  annotate("text", x = 5900, y = 0.25, label = "Phase III", angle = 90) + 
  geom_vline(xintercept = 5500,colour = 'grey') +
  annotate("text", x = 1000, y = 0.25, label = "Phase IV", angle = 90) +
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_blank()  
  )

# From Science paper - four chonological divisions
# The chronological divisions are: 
# Phase I (spits 69-50) dated ~42 to 38 ka; 
# Phase II (spits 49-40) # dated between ~ 17 to 9 ka, 
# Phase III (spits 39-21) dated between 6.5 and 5.5 ka, 
# Phase IV (spits 20-3) dated ~5.5 - modern.

# My Q: why omit 17-38 ky BP? 
 
ggsave("figures/Jeremalai-dates.png")
```

                 
### Site chronology: analysing the lithic distribution over time

```{r lithis_over_time}
# omit rows with blanks or NAs
flakes <- flakes[!(flakes$Weight == "" | is.na(flakes$Weight)), ]

# put depths on lithic data
flakes$depth <- depths$Depth.bs..m[match(flakes$Spit,depths$Spit.no)]

# omit rows with blanks or NAs... again
flakes <- flakes[!(flakes$depth == "" | is.na(flakes$depth)), ]

# From Science paper - four chonological divisions
# The chronological divisions are: 
# Phase I (spits 69-50) dated ~42 to 38 ka; 
# Phase II (spits 49-40) # dated between ~ 17 to 9 ka, 
# Phase III (spits 39-21) dated between 6.5 and 5.5 ka, 
# Phase IV (spits 20-3) dated ~5.5 - modern.

flakes$phase <- ifelse(flakes$Spit > 0 & flakes$Spit <= 20, 4,
                       ifelse(flakes$Spit >= 21 & flakes$Spit <= 39, 3,
                              ifelse(flakes$Spit >= 40 & flakes$Spit <= 49, 2,
                                     ifelse(flakes$Spit >= 50 & flakes$Spit <= 69, 1, NA))))
# check if any NA
unique(flakes$phase)

# My version after Sue's comment
# The chronological divisions are: 
# Phase I (spits 69-56) dated ~42 to 38 ka; 
# Phase II (spits 55-49) # dated between ~ 38 to 17 ka, 
# Phase III (spits 48-40) # dated between ~ 17 to 9 ka, 
# Phase IV (spits 39-21) dated between 6.5 and 5.5 ka, 
# Phase V (spits 20-3) dated ~5.5 - modern.


# here's a function to assign phases based on spit numbers
  makephases <- function(x) {ifelse(x > 0 & x <= 20, 5,
                       ifelse(x >= 21 & x <= 39, 4,
                              ifelse(x >= 40 & x <= 48, 3,
                                     ifelse(x >= 49 & x <= 55, 2,
                                            ifelse(x >= 56 & x <= 69, 1, NA)))))}

flakes <- flakes[!(is.na(flakes$Spit)),]
flakes$phase  <- makephases(flakes$Spit)
# check if any NA, don't want any, should return all TRUE
is.na(unique(flakes$phase)) %in% FALSE
# get phase durations
phases <- data.frame(phase = 1:5,
                     start = c(42, 38, 17, 6.5, 5.5),
                     end =   c(38, 17, 9,  5.5, 0 ))
phases$duration <- with(phases, start - end)
```


### Lithics: analysing the raw material proportions and change over time


```{r lithics_raw_material_over_time}
# raw material
raw <- dcast(flakes, Material ~ phase) # change group to depth for hi res
# remove row with no raw material
raw <- raw[-10,]
rownames(raw) <- raw[,1]
# get rid of rows with no value
raw <- raw[rownames(raw) != "",]
# remove col of NA
raw <- raw[,colnames(raw) != "NA"]
raw <- raw[,-1]
# stat test on all
results <- chisq.test(raw, simulate.p.value=T, B=9999)
results$expected
prop.table(raw)
colSums(prop.table(raw))
rowSums(prop.table(raw))
# stat test on just dominant materials
dom <- raw[rowSums(raw) > 10,]
prop.table(as.matrix(raw[rowSums(raw) > 20,]), 2)
# a slight decrease in chert in Holocene
result <- chisq.test(dom, simulate.p.value=FALSE, B=999)
result
# so it's significant in the frequentist system, but what's the strength?
# Cramer's V: is a measure of association for nominal variables. Effectively it is the Pearson chi-square statistic rescaled to have values between 0 and 1,
assocstats(as.matrix(dom)) # very small V
# plot dominant raw materials
dom_m <- dom
dom_m$raw_type <- rownames(dom_m) 
dom_m <- melt(dom_m)
ggplot(dom_m, aes(variable, value, fill = raw_type)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("Depositional phase") +
  ylab("Artefact frequency") +
  scale_fill_discrete(name="Raw material")
# save plot
ggsave("figures/Jeremalai-raw-materials-group-counts.png")

# send out for bayesian contingency table test...
# reshape and give meaningful names...
colnames(dom) <- c("phase 1", "phase 2", "phase 3", "phase 4", "phase 5" )
data <- melt(as.matrix(dom),  varnames=c("raw_material", "phase"), value.name="Freq")
data <- data[data$Freq != 0 | !is.numeric(data$Freq),] # zeros break the model
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has a last
# column called 'Freq'
# source("code/BensPoissonExponentialJagsSTZ.R")

# all interactions include zero in the HDI, so no credible differences 
# in raw materials between phases. 

# table to html file to paste into word doc, include obsidian
dom_obs <- raw[rowSums(raw) > 10 | row.names(raw) == 'Obsidian',]
print(xtable(dom_obs), type = 'html', html.table.attributes = (border=0), file = "raw_tab.html")

raw_freqs <- data.frame(raw_types = rownames(raw), Freq = rowSums(raw))
ggplot(raw_freqs, aes(reorder(raw_types, Freq), Freq)) + 
  geom_bar(stat="identity")
# raw materials by site
# compute proportions per layer (col props)
all_tab <- data.frame()
for(i in seq(ncol(raw))){
  for(j in seq(nrow(raw))){
    all_tab[j,i] <- raw[j,i]/colSums(raw)[i]
  }
}
# check 
colSums(all_tab) # should == 1
colnames(all_tab) <- colnames(raw) 
rowSums(all_tab) # should be various
all_tab$raw_type <- rownames(raw) 
# get rid of raw materials that are not very abundant
all_tab <- all_tab[which(rowSums(all_tab[,1:ncol(all_tab)-1]) > 0.02) , ]

# get rid of NA column
all_tab <- all_tab[,names(all_tab) != 'NA']

# plot
all_tab_m <- melt(all_tab, id.var = 'raw_type')
ggplot(all_tab_m, aes(variable, value, fill = raw_type)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("Depth below surface (m)") +
  ylab("Proportion of depositional phase") +
  scale_fill_discrete(name="Raw material")
# save plot
ggsave("figures/Jeremalai-raw-materials-group.png")
```

### Lithics: analysing discard rates and change over time


```{r lithics_discard_over_time}
# discard rates
discard <- aggregate(Weight ~ depth + Spit, flakes, length)

# sediment volumes: put volumes on
discard$sedvol <- vols$Soil[match(discard$Spit, vols$X)]
# put spit thickesses on
discard$thick <- c(0.018, diff(discard$depth)) # add first value from depth_and_dates.xsl
# compute artefacts per kg of sediment
discard$kgsed <- with(discard, Weight / sedvol) # weight is count of artefacts that have a weight
# compute artefact per cubic meter (spit thickess)
discard$cubmet <- with(discard, Weight / thick)
# seems we have an unusually extreme value in spit 34
# omit - perhaps a data collection typo
discard <- discard[discard$Spit != 34, ]
# Plot
ggplot(discard, (aes(depth, cubmet))) +
  geom_point() +
  theme_minimal() +
  stat_smooth(span = 0.5, se = FALSE) +
  xlab("Depth below surface (m)") +
  ylab("Number of chert flakes per cubic meter of deposit") +
  coord_flip() +
  scale_x_reverse()
# save plot
ggsave("figures/Jeremalai-flake-discard.png")

# or put ages on the y-axis
x <- dates$depth_bs *100; y <-  dates$mid

# try loess... not ideal because of loops over some points...
span <- 0.3 # this has a big influence on the shape of the curve, experiment with it!
cal.date.lo <- loess(mid  ~ depth_bs, dates, span = span)
cal.date.pr <- predict(cal.date.lo, data.frame(depth_bs = seq(0, max(dates$depth_bs), 0.01)))
plot(cal.date.pr) # inspect to seee that it looks right
points(dates$depth_bs *100, dates$mid, pch = 19, col = "green") # check fit
cal.date.pr <- data.frame(age = unname(cal.date.pr), depth = as.numeric(names(cal.date.pr)))

# try polynomial regression, also not ideal because of looping over
polyfit <- function(i) x <- AIC(lm(y~poly(x,i)))
n <- as.integer(optimize(polyfit,interval = c(1,length(x)-1))$minimum)
fit <- lm(y ~ poly(x, 9))
summary(fit)
plot(x, predict(fit), type="l", col="red", lwd=2)
points(dates$depth_bs *100, dates$mid, col = 'blue') # check fit
xx <- seq(min(x),max(x), length=200)
plot(x, y, pch = 19, col = "green")
lines(xx, predict(fit, data.frame(x=xx)), col="red")
points(dates$depth_bs *100, dates$mid, col = 'blue', pch = 19) # check fit

# try approx, looks better, use this
plot(approx(x, y), col = 'red')
points(dates$depth_bs *100, dates$mid, col = 'blue', pch = 19) # check fit
# get output
discard$age <- approx(x/100, y, xout = discard$depth)$y
# omit spit 34 since it seems to be anomalously small
# perhaps a data collection typo
discard <- discard[discard$Spit != 34, ]
# assign spits 1 and 2 a nominal age of 100 BP
discard[discard$Spit == 1, ]$age <- 100
discard[discard$Spit == 2, ]$age <- 100

# plot discard density by age
ggplot(discard, (aes(age/1000, cubmet))) +
  geom_point() +
  theme_minimal() +
 #  stat_smooth(span = 0.7, se = FALSE) +
  xlab("Age ka BP") +
  ylab("Number of chert flakes per cubic meter of deposit") +
  coord_flip() +
  scale_x_reverse()
# save plot
ggsave("figures/Jeremalai-flake-discard-age.png")

# plot artefacts/cubic meter by phase, get phase number for each spit
discard$phase <- makephases(discard$Spit)
discard_agg <- aggregate(cubmet ~ phase, discard, mean)
ggplot(discard_agg, (aes(phase, cubmet))) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("Depositional phase") +
  ylab("Mean number of chert flakes per cubic meter of deposit") 
# save plot
ggsave("figures/Jeremalai-flake-discard-phase.png")

# plot artefacts/cubic meter/1000 years by phase, get phase number for each spit
# this is the most sensible option
discard$phase <- makephases(discard$Spit)
discard_agg <- aggregate(cubmet ~ phase, discard, mean)
discard_agg$cubmetperkyr <- discard_agg$cubmet / phases$duration
ggplot(discard_agg, (aes(phase, cubmetperkyr))) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("Depositional phase") +
  ylab("Mean number of chert flakes per \ncubic meter of deposit per thousand years") 
# save plot
ggsave("figures/Jeremalai-flake-discard-phase-m3.png")

# plot artefacts/kg sed by phase, get phase number for each spit
discard$phase <- makephases(discard$Spit)
discard_agg <- aggregate(sedvol ~ phase, discard, mean)
# per thousand years
discard_agg$sedvolky <- discard_agg$sedvol /  phases$duration
ggplot(discard_agg, (aes(phase, sedvolky))) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("Depositional phase") +
  ylab("Mean number of chert flakes \n per kg of deposit per 1000 yrs") 
ggsave("figures/Jeremalai-flake-discard-phase-kg.png")
```


### Lithics: artefact taphonomy


```{r lithic_taphonomy}
allchert <- all[all$Material == 'Chert', ]
allchert$phase <- makephases(allchert$Spit)
# make Artclass that is long and transv breaks
allchert$Artclas <- ifelse(allchert$Breaks == "", 
       as.character(allchert$Artclas), 
       paste(allchert$Artclas, allchert$Breaks, sep = "-"))

taph <- data.frame(table(allchert$Artclas))

# use regex to get broken flakes -b- 
broken <- allchert[grep("-b", allchert$Artclas), ]


# get counts of broken to complete per phase
# flake to -b-

breaks <- dcast(allchert, Artclas ~ phase)[-1,]
breaks <- breaks[breaks$Artclas =="Flake" | grepl("-b-", breaks$Artclas), ]
# total count of broken artefacts per phase
brk_tot <- colSums(filter(breaks, grepl("-b-", Artclas))[,-1])
flk_tot <- breaks[1,-1]
# table of broken to complete flakes
brk_flk <- data.frame(rbind(brk_tot = brk_tot, 
                      flk_tot = flk_tot))
colnames(brk_flk) <- c("1", "2", "3", "4", "5" )

# stat test
assocstats(as.matrix(brk_flk))

# plot
brk_flk_m <- melt(as.matrix(brk_flk),  varnames=c("breakt", "phase"), value.name="Freq")
ggplot(brk_flk_m, aes(phase, Freq, fill = breakt)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_discrete(name="flake type",
                      labels=c("broken flakes", "complete flakes"))
# save plot
ggsave("figures/Jeremalai-flake-broken-phase.png")

# send out for bayesian contingency table test...
# reshape and give meaningful names...
colnames(brk_flk) <- c("phase1", "phase2", "phase3", "phase4", "phase5" )
data <- melt(as.matrix(brk_flk),  varnames=c("breakt", "phase"), value.name="Freq")
data <- data[data$Freq != 0,] # zeros break the model
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has a last
# column called 'Freq'
# source("code/BensPoissonExponentialJagsSTZ.R")

# credible interactions in phase 4 and 5


# and transverse and longitude
allchert$Artclas <- tolower(allchert$Artclas)
allchert$breakt <- "" # create variable to fill

allchert$breakt[grep("trans", allchert$Artclas)] <- "trans"
allchert$breakt[grep("long", allchert$Artclas)] <- "long"

# per depositional phase

breakt <- dcast(allchert, breakt ~ phase)[-1,]
# add complete flake counts
breakt <- rbind( breakt , setNames( breaks[1, ] , names( breakt ) ) )
# shift rownames out and delete them
rownames(breakt) <- breakt[,1]
breakt <- breakt[,-1]
# stat test
assocstats(as.matrix(breakt))

# do bayesian contingency table test
colnames(breakt) <- c("phase 1", "phase 2", "phase 3", "phase 4", "phase 5" )
data <- melt(as.matrix(breakt),  varnames=c("breakt", "phase"), value.name="Freq")
data <- data[data$Freq != 0,] # zeros break the model
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has a last
# column called 'Freq'
# source("code/BensPoissonExponentialJagsSTZ.R")

# plot 
data$phase <- gsub("[[:alpha:]]*", "", data$phase)
ggplot(data, aes(phase, Freq, fill = breakt)) +
  ylab("Frequency") +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_discrete(name="flake type",
                      labels=c("complete flakes", 
                               "trans. broken flakes",
                               "long. broken flakes"))
ggsave("figures/Jeremalai-flake-broken-phase.png")
```

### Lithics: heat treatment

```{r lithics_heat_treatment}
sum(flakes$Heat, na.rm = TRUE) / nrow(flakes)
heat <- aggregate(Heat ~ phase, flakes, length)
total <- aggregate(Spit ~ phase, flakes, length)
heat$notheat <- total$Spit - heat$Heat
# show proportions that are heat-treated
heat$Heat / total$Spit
heat <- t(heat)

# frequentist stat test
assocstats(as.matrix(heat))

# do bayesian contingency table test
colnames(heat) <- c("phase 1", "phase 2", "phase 3", "phase 4", "phase 5" )
heat <- heat[-1,]
data <- melt(as.matrix(heat),  varnames=c("heat", "phase"), value.name="Freq")
data <- data[data$Freq != 0,] # zeros break the model
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has a last
# column called 'Freq'
# source("code/BensPoissonExponentialJagsSTZ.R")

# plot
data$phase <- gsub("[[:alpha:]]*", "", data$phase)
ggplot(data, aes(phase, Freq, fill = heat)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_discrete(name="heat treatment",
                      labels=c("heat treated", 
                               "unheated"))
ggsave("figures/Jeremalai-flake-heat-phase.png")
```

### Lithics: chert flake metrics 

```{r lithics_chert_metrics}
metrics <- flakes %.% group_by(phase) %.% summarise(mean(Weight), mean(Length), mean(Width), mean(Thick))

ggplot(flakes, aes(as.factor(Spit), Weight) ) +
  geom_point() +
  scale_y_log10() 
#ggsave("figures/Jeremalai-mass-spit.png")

ggplot(flakes, aes(as.factor(phase), Weight) ) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  xlab("depositional phase") +
  ylab("mass (g)")
ggsave("figures/Jeremalai-flake-mass-phase.png")

# ANOVA with Tukey's HSD
fit <- aov(Weight ~ as.factor(phase), flakes)
summary(fit)
tuk <- TukeyHSD(fit)
plot(tuk, las = 2, cex = 0.1)

# do bayesian ANOVA
data <- data.frame(phase = flakes$phase, mass = flakes$Weight)
data <- data[data$mass != 0,] # zeros break the model
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has col1
# as numbers indicating groups and col2 as the measurement
# source("code/BensANOVAonewayJagsSTZ.R")

# core mass
cores_mass <- allchert[allchert$Artclas == "core", ]
core_metrics <- cores_mass %.% group_by(phase) %.% summarise(mean(Weight), mean(Length), mean(Width), mean(Thick))

ggplot(cores_mass, aes(as.factor(Spit), Weight) ) +
  geom_point() +
  scale_y_log10() 
  #ggsave("figures/Jeremalai-mass-spit.png")
  
ggplot(cores_mass, aes(as.factor(phase), Weight) ) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  xlab("depositional phase") +
  ylab("mass (g)")
ggsave("figures/Jeremalai-core-mass-phase.png")

# ANOVA with Tukey's HSD
fit <- aov(Weight ~ as.factor(phase), cores_mass)
summary(fit)
tuk <- TukeyHSD(fit)
plot(tuk, las = 2, cex = 0.1)

# do bayesian ANOVA
data <- data.frame(phase = cores_mass$phase, mass = cores_mass$Weight)
data <- data[data$mass != 0,] # zeros break the model)
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has col1
# as numbers indicating groups and col2 as the measurement
# source("code/BensANOVAonewayJagsSTZ.R")

# plot flake and core mass by phase in one box plot
flakes_cores_weight <- allchert[(allchert$Artclas == "flake" | allchert$Artclas == "core"), c('Artclas', 'Weight', 'phase') ]

ggplot(flakes_cores_weight, aes( fill = Artclas, as.factor(phase), Weight)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  xlab('depostional phase') +
  ylab("mass (g)") +
  scale_fill_manual(name="artefact type",
                      labels=c("core", "complete flake"), 
                      values = c("black", "white"))
# save plot
ggsave("figures/Jeremalai-flake-core-mass-phase.png")
```

### Lithics: chert flake platform metrics 

```{r chert_flake_platform}
# flake platform
plat <- dcast(flakes, Plat ~ depth) 
rownames(plat) <- plat[,1]
# get rid of rows with no plat
plat <- plat[rownames(plat) != "",]
plat <- plat[,-1]
# stat test on all
# frequentist stat test
assocstats(as.matrix(plat))
chisq.test(plat, simulate.p.value=T, B=999)
# stat test on just dominant plats
chisq.test(plat[rowSums(plat) > 70,], simulate.p.value=T, B=999)
assocstats(as.matrix(plat[rowSums(plat) > 70,]))
plat_freqs <- data.frame(plat_types = rownames(plat), Freq = rowSums(plat))
plat_freqs <- plat_freqs[plat_freqs$Freq > 90,]
plat_freqs$plat_types <- c('2-scars', '3-scars', 'cortical', 'focalised', 'single')
# plot freq of platform types
main <- ggplot(plat_freqs, aes(reorder(plat_types, Freq), Freq)) + 
  geom_bar(stat="identity") +
  # theme_minimal() +
  xlab("platform type") +
  ylab("frequency") +
  # remove grid lines for subplot
  theme_update(panel.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())

# raw materials by spit
# compute proportions per layer (col props)
all_tab <- data.frame()
for(i in seq(ncol(plat))){
  for(j in seq(nrow(plat))){
    all_tab[j,i] <- plat[j,i]/colSums(plat)[i]
  }
}
# check 
colSums(all_tab) # should == 1
colnames(all_tab) <- colnames(plat) 
rowSums(all_tab) # should be various
all_tab$plat_type <- rownames(plat) 
# get rid of raw materials that are not very abundant
all_tab <- all_tab[which(rowSums(all_tab[,1:ncol(all_tab)-1]) > 1.5) , ]

# get rid of NA column
all_tab <- all_tab[,names(all_tab) != 'NA']
# stat test
chisq.test(all_tab[,-ncol(all_tab)])
assocstats(as.matrix(all_tab[,-ncol(all_tab)]))
# plot
all_tab_m <- melt(all_tab, id.var = 'plat_type')
ggplot(all_tab_m, aes(variable, value, fill = plat_type)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("Depth below surface (m)") +
  ylab("Proportion of spit") +
  scale_fill_discrete(name="Platform condition")
#ggsave("F:/My Documents/My Papers/conferences/IPPA14/Jeremalai/figures/Jeremalai-platform-depth.png")

# by phase
plat <- dcast(flakes, Plat ~ phase)
rownames(plat) <- plat[,1]
# get rid of rows with no plat
plat <- plat[rownames(plat) != "",]
plat <- plat[,-1]
# stat test on all
# frequentist stat test
assocstats(as.matrix(plat))
chisq.test(plat, simulate.p.value=T, B=999)
# stat test on just dominant plats
chisq.test(plat[rowSums(plat) > 20,], simulate.p.value=T, B=999)
assocstats(as.matrix(plat[rowSums(plat) > 20,]))

# raw materials by site
# compute proportions per phase (col props)
all_tab <- data.frame()
for(i in seq(ncol(plat))){
  for(j in seq(nrow(plat))){
    all_tab[j,i] <- plat[j,i]/colSums(plat)[i]
  }
}
# check 
colSums(all_tab) # should == 1
colnames(all_tab) <- colnames(plat) 
rowSums(all_tab) # should be various
all_tab$plat_type <- rownames(plat) 
# get rid of raw materials that are not very abundant
all_tab <- all_tab[which(rowSums(all_tab[,1:ncol(all_tab)-1]) > 0.15) , ]

# get rid of NA column
all_tab <- all_tab[,names(all_tab) != 'NA']
# stat test
chisq.test(all_tab[,-ncol(all_tab)])
# plot
all_tab_m <- melt(all_tab, id.var = 'plat_type')
ggplot(all_tab_m, aes(variable, value, fill = plat_type)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("depositional phase") +
  ylab("Proportion of phase") +
  scale_fill_discrete(name="Platform condition")
ggsave("figures/Jeremalai-platform-phase.png")


# plot distibution of platform sizes for each type
flakes$Platarea <- with(flakes, (Platthic * Platwid))
plat_area_type <- flakes[flakes$Plat %in% c("Single", "Focal", "2-scars", "Cort", '3-scars'),]

# make names a bit more readable
plat_area_type$Plat <- ifelse(plat_area_type$Plat == 'Cort', 'cortical',
                      ifelse(plat_area_type$Plat == 'Focal',  'focalised', 
                    ifelse(plat_area_type$Plat == 'Single',  'single', as.character(plat_area_type$Plat))))  

# put types in same order as frequency plot
plat_area_type$Plat <- factor(plat_area_type$Plat, 
                              levels = c('cortical', '3-scars','2-scars', 'focalised', 'single'), ordered = TRUE)
                               
sub <- ggplot(plat_area_type, aes(Plat, Platarea)) +
  geom_boxplot() +
  scale_y_log10() +
  ylab(as.expression(bquote('platform area (' * mm^{2} * ")" ))) +
  xlab("") +
  theme_minimal()
ggsave("figures/Jeremalai-platform-area-by-phase.png")

# plot freq of plat type and platform area together in one plot
  
vp <- viewport(width = 0.4, height = 0.5, 
               x = 0.5, y = 0.4, 
               just = c("right", "bottom"))

# combine plots, print and save (wont show in console)
png("figures/Jeremalai-platform-area-by-plat-type.png", w = 650, h = 350)
   print(main)
   print(sub + theme_bw(base_size = 10), vp = vp)
dev.off()

# test for diff in area between types
plat_type_area_test <- na.omit(flakes[flakes$Plat %in% c("Single", "Focal", "2-scars"),c('Plat', 'Platarea')])

# ANOVA with Tukey's HSD
fit <- aov(Platarea ~ as.factor(Plat), plat_type_area_test)
summary(fit)
tuk <- TukeyHSD(fit)
par(mar=c(5, 7, 5, 5))
plot(tuk, las = 2, cex = 0.1)
```

### Lithics: chert flake cortex  

```{r lithics_chert_cortex}

# dorsal cortex on flakes
ggplot(flakes, aes(as.factor(phase), Cortex)) + 
  geom_boxplot() +
  theme_minimal() +
  scale_y_log10() +
  xlab("depositional phase") +
  ylab("Percent of dorsal cortex")
ggsave("figures/Jeremalai-flake-dorsal-cortex-phase.png")

# do bayesian ANOVA
data <- data.frame(phase = flakes$phase, cortex = flakes$Cortex)
data <- data[data$cortex != 0,] # zeros break the model)
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has col1
# as numbers indicating groups and col2 as the measurement
# source("code/BensANOVAonewayJagsSTZ.R")

# Frequentist test also output from that script.
# no credible interactions under expectation of not independant 
```


### Lithics: chert flake and core cortex  

```{r lithics_core_cortex}

# dorsal cortex on cores
ggplot(cores_mass, aes(as.factor(phase), Cortex)) + 
  geom_boxplot() +
  theme_minimal() +
  scale_y_log10() +
  xlab("depositional phase") +
  ylab("Percent of dorsal cortex")
ggsave("figures/Jeremalai-core-cortex-phase.png")

# do bayesian ANOVA
data <- data.frame(phase = cores_mass$phase, cortex = cores_mass$Cortex)
data <- na.omit(data[data$cortex != 0 ,]) # zeros break the model)
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has col1
# as numbers indicating groups and col2 as the measurement
# source("code/BensANOVAonewayJagsSTZ.R")

# Frequentist test also output from that script.
# no credible interactions under expectation of not independant 

# plot flake and core cortex by phase in one box plot
flakes_cores_cortex <- allchert[(allchert$Artclas == "flake" | allchert$Artclas == "core"), c('Artclas', 'Cortex', 'phase') ]

ggplot(flakes_cores_cortex, aes( fill = Artclas, as.factor(phase), Cortex)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  xlab('depostional phase') +
  ylab("percent cortex") +
  scale_fill_manual(name="artefact type",
                    labels=c("core", "complete flake"), 
                    values = c("black", "white"))
# save plot
ggsave("figures/Jeremalai-flake-core-cortex-phase.png")
```

### Lithics: chert flake scars

```{r lithics_flake_scars}

# number of flake scars per flake
# remove outliers
flakes <- flakes[flakes$NoDS < 20 & !(is.na(flakes$phase)), ]
ggplot(flakes, aes(as.factor(phase), NoDS)) + 
  geom_boxplot() +
  theme_minimal() +
  scale_y_log10() +
  xlab("Depth below surface (m)") +
  ylab("Number of dorsal flake scars")
# save
ggsave("figures/Jeremalai-dorsal-cortex-group.png")

# tables of core and flake metrics, etc
# flake cortex, dorsal flake scars and OHR together
ag_fl <- aggregate(NoDS ~ phase, flakes, function(x) c(fl_mean = mean(x), fl_sd = sd(x)))
ag_co <- aggregate(NoDS ~ phase, cores_mass, function(x) c(co_mean = mean(x), co_sd = sd(x)))
# get proportion of all flakes with OHR in each phase
flakes$ohr <- ifelse(flakes$Overhang == 'Yes', 1, 0)
ohr <- aggregate(ohr ~ phase, flakes, function(x) c(ohr_mean = sum(x)/length(x)))
# output to format for doc
write.csv(t(cbind(ag_fl, ag_co[,-1], ohr[,-1])), "table.csv")


### come back to this ###
```

### Lithics: core mass by type

```{r lithic_core_mass}
core_types$Type_long <- with(core_types, ifelse(Type == "SPC", "Single plaform",
                                     ifelse(Type == "RC", "Radial",
                                      ifelse(Type == "BDC", "Bidirectional",
                                       ifelse(Type == "BiC", "Bipolar",
                                        ifelse(Type == "MPC", "Multi-platform",
                                         ifelse(Type == "LLC", "Levallois-like",
                                          ifelse(Type == "FFC", "Faceted flake", NA))))))) )

# plot
ggplot(core_types, aes(reorder(Type_long, -Mass,  FUN=median), Mass)) +
  geom_boxplot() +
  ylim(0,30) +
  theme_minimal() +
  xlab("Core type") +
  ylab("Core mass (g)")
# save
ggsave("figures/Jeremalai-core-by-type.png")

aggregate(Mass ~ Type_long, data = core_types, median)

```

### Lithics: retouch

```{r lithics_retouch}

# frequency of flakes with retouch per phase
rt <- flakes[flakes$Rtch == "Yes", ]
rt_tab <- t(data.frame(rt = aggregate(Retype ~ phase, rt, function(x) rt_count = length(x)), flakes =  aggregate(Retype ~ phase, flakes, function(x) fk_count = length(x))[,2] ))

# send out for bayesian contingency table test...
# reshape and give meaningful names...
colnames(rt_tab) <- c("phase 1", "phase 2", "phase 3", "phase 4", "phase 5" )
rt_tab <- rt_tab[-1,]
data <- melt(as.matrix(rt_tab),  varnames=c("retouch", "phase"), value.name="Freq")
data <- data[data$Freq != 0 | !is.numeric(data$Freq),] # zeros break the model
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has a last
# column called 'Freq'
# source("code/BensPoissonExponentialJagsSTZ.R")

# nothing from bayes, nothing from freq.
assocstats(as.matrix(rt_tab))

# get proportion
rt_tab <- data.frame(t(rt_tab))
rt_tab$prop <- with(rt_tab, (rt.Retype/flakes))
  
# amount of the flake that is retouched as measured by retouch length 

rt_len <- data.frame(phase = rt$phase, rt_len = rt$Retlen)
# omit NAs
rt_len <- rt_len[!(is.na(rt_len$rt_len)), ]
ggplot(rt_len, aes(as.factor(phase), rt_len)) +
  scale_y_log10() + 
  geom_boxplot()
# get mean lengths
aggregate(rt_len ~ phase, rt_len, mean)

# do bayesian ANOVA
data <- rt_len
data <- na.omit(data[data$rt_len != 0 ,]) # zeros break the model)
# This script is in the same folder as the current script
# and generates plots and a lot of data objects. The only 
# input is a long data.frame called 'data' that has col1
# as numbers indicating groups and col2 as the measurement
# source("code/BensANOVAonewayJagsSTZ.R")

# no credible differences

# plot retouch flake size and complete unretouch flake size
#  in one box plot
flakes_retouch_size <- allchert[(allchert$Artclas == "flake" | allchert$Artclas == "retf"), c('Artclas', 'Weight', 'Length', 'phase') ]

ggplot(flakes_retouch_size, aes( fill = Artclas, as.factor(phase), Weight)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  xlab('depostional phase') +
  ylab("mass (g)") +
  scale_fill_manual(name="artefact type",
                    labels=c("unretouched flake", "retouched flake"), 
                    values = c("black", "white"))
# save plot
ggsave("figures/Jeremalai-flake-retouchedflake-mass-phase.png")

# any difference between mass of retouched and unretouched, in total?
# Bayesian t-test of the mass of all retouched vs all unretouched flakes
# flakes_retouch_size_t_test <- BESTmcmc(flakes_retouch_size[flakes_retouch_size$Artclas == "flake", ]$Weight, flakes_retouch_size[flakes_retouch_size$Artclas == "retf", ]$Weight)
```

### Lithics: retouch location

```{r lithics_retouch_location}

# retouch locations
rt_loc <- data.frame(phase = rt$phase, rt_loc = rt$Retloc)
rt_loc <- rt_loc[!(is.na(rt_loc$phase)), ]
rt_sum <- as.data.frame.matrix(table(rt_loc))[,-1]
rt_sum$phase <- paste0("phase ", row.names(rt_sum))
colnames(rt_sum) <- c("Both margins", "Distal", "Perimeter", "Left lateral", "Medial", "Proximal", "Right lateral", "phase")
rt_sum_m <- melt(rt_sum)
# plot
ggplot(rt_sum_m, aes(variable, value)) +
  geom_bar(stat="identity") +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12)) +
  ylab("frequency") +
  xlab("retouch location") +
  facet_wrap(~phase, ncol = 3)
# save plot
ggsave("figures/Jeremalai-flake-retouched-flake-location-phase.png")
```

### Lithics: retouch indices


```{r lithics_retouch_indices}
retouch_indices[is.na(retouch_indices)] <- 0
retouch_indices$GIUR <- with(retouch_indices, t1/T1 + t2/T2 + t3/T3)/3
retouch_indices$perimeter_perc <- with(retouch_indices, length/perimeter * 100)
retouch_indices$II <- with(retouch_indices, ((X0.5 * 0.5) + (X1 * 1))/16)
# get mean and standard deviation for each index
retouch_indices_subset <- retouch_indices %>% select(GIUR, perimeter_perc, II) 
# sweep over the columns to compute mean and standard deviation
retouch_indices_means <- data.frame(t(round(apply(retouch_indices_subset, 2, mean, na.rm = TRUE),2)))
retouch_indices_sds <- data.frame(t(round(apply(retouch_indices_subset, 2, sd, na.rm = TRUE),2)))


GIUR = 0.41 ± 0.26, % perimeter retouched = 29 ± 14, II = 0.16 ± 0.13).

```

The retouch intensity can be summarised with the following metrics:
GIUR = `r retouch_indices_means$GIUR` +/- `r retouch_indices_sds$GIUR`
perimeter = `r retouch_indices_means$perimeter_perc` +/- `r retouch_indices_sds$ perimeter_perc`%
II = `r retouch_indices_means$II` +/- `r retouch_indices_sds$II`



### Lithics: technological types

```{r lithics_technological_types}

# all 

# combine

techno <- data.frame(cores, types, retouch, features, ground,  stringsAsFactors = FALSE)
# remove extra Spit cols (is those called 'Spit.1' etc)
techno <- techno[,-grep("Spit\\.", colnames(techno)) ]
# put depths on 
techno$depth <- depths$Depth.bs..m[match(techno$Spit,depths$Spit.no)]
# put phases on 
techno$phase <- flakes$phase[match(techno$Spit,flakes$Spit)]

# Chi Square based on % of spits in each temporal block that have a technological type present

techno[] <- lapply(techno, as.character) # change factor to character
techno[techno == 'x'] <-  1 # replace x with 1
techno[] <- lapply(techno, as.numeric) # change char to num
# get row sums
rs <- rowSums(techno[,2:32], na.rm = TRUE )
# get row sums by group
rs_phase <- data.frame(rs, phase = techno$phase)
# % that have a type present
# phases with no types at all
rs_phase[rs_phase$rs == 0,]
# set 1 or above to 1
rs_phase$rs <- ifelse(rs_phase$rs == 0, 0, 1)
aggregate(rs ~ phase, rs_phase[rs_phase$rs == 0,], length) # counts
t(apply(table(rs_phase), 2, function(x) x/sum(x))) # props
# just three spits with zero types... let's do it by class

## cores ##
# put depths on 
cores$depth <- depths$Depth.bs..m[match(cores$Spit,depths$Spit.no)]
# put groups on 
cores$phase <- flakes$phase[match(cores$Spit,flakes$Spit)]
cores[] <- lapply(cores, as.character) # change factor to character
cores[cores == 'x'] <-  1 # replace x with 1
cores[] <- lapply(cores, as.numeric) # change char to num
# get row sums
rs <- rowSums(cores[,2:(ncol(cores)-2)], na.rm = TRUE )
# get row sums by phase
rs_phase <- data.frame(rs, group = cores$phase)
# % that have a type present
# groups with no types at all
rs_phase[rs_phase$rs == 0,]
# set any non-zero to 1
rs_phase$rs <- ifelse(rs_phase$rs != 0, 1, rs_phase$rs)
# yes, more with zero here...
dc_core <- data.frame(t(apply(table(rs_phase), 2, function(x) x/sum(x))))
dc_core <- prop.table(as.matrix(dc_core), 1)
# chi-square test
chisq.test(dc_core, simulate.p.value=T, B=9999)


## types ##
# put depths on 
types$depth <- depths$Depth.bs..m[match(types$Spit,depths$Spit.no)]
# put groups on 
types$phase <- flakes$phase[match(types$Spit,flakes$Spit)]
types[] <- lapply(types, as.character) # change factor to character
types[types == 'x'] <-  1 # replace x with 1
types[] <- lapply(types, as.numeric) # change char to num
# get row sums
rs <- rowSums(types[,2:(ncol(types)-2)], na.rm = TRUE )
# get row sums by group
rs_phase <- data.frame(rs, phase = types$phase)
# % that have a type present
# groups with no types at all
rs_phase[rs_phase$rs == 0,]
# set any non-zero to 1
rs_phase$rs <- ifelse(rs_phase$rs != 0, 1, rs_phase$rs)
# yes, more with zero here...
dc_types <- data.frame(t(apply(table(rs_phase), 2, function(x) x/sum(x))))
dc_types <- prop.table(as.matrix(dc_types), 1)
# chi-square test
chisq.test(dc_types, simulate.p.value=T, B=9999)


## retouch ##
# put depths on 
retouch$depth <- depths$Depth.bs..m[match(retouch$Spit,depths$Spit.no)]
# put groups on 
retouch$phase <- flakes$phase[match(retouch$Spit,flakes$Spit)]
retouch[] <- lapply(retouch, as.character) # change factor to character
retouch[retouch == 'x'] <-  1 # replace x with 1
retouch[] <- lapply(retouch, as.numeric) # change char to num
# get row sums
rs <- rowSums(retouch[,2:(ncol(retouch)-2)], na.rm = TRUE )
# get row sums by group
rs_phase <- data.frame(rs, phase = retouch$phase)
# % that have a type present
# groups with no types at all
rs_phase[rs_phase$rs == 0,]
# set any non-zero to 1
rs_phase$rs <- ifelse(rs_phase$rs != 0, 1, rs_phase$rs)
# yes, more with zero here...
dc_retouch <- data.frame(t(apply(table(rs_phase), 2, function(x) x/sum(x))))
dc_retouch <- prop.table(as.matrix(dc_retouch), 1)
# chi-square test
chisq.test(dc_retouch, simulate.p.value=T, B=9999)


## features ##
# put depths on 
features$depth <- depths$Depth.bs..m[match(features$Spit,depths$Spit.no)]
# put groups on 
features$phase <- flakes$phase[match(features$Spit,flakes$Spit)]
features[] <- lapply(features, as.character) # change factor to character
features[features == 'x'] <-  1 # replace x with 1
features[] <- lapply(features, as.numeric) # change char to num
# get row sums
rs <- rowSums(features[,2:(ncol(features)-2)], na.rm = TRUE )
# get row sums by group
rs_phase <- data.frame(rs, phase = features$phase)
# % that have a type present
# groups with no types at all
rs_phase[rs_phase$rs == 0,]
# set any non-zero to 1
rs_phase$rs <- ifelse(rs_phase$rs != 0, 1, rs_phase$rs)
# yes, more with zero here...
dc_feat <- data.frame(t(apply(table(rs_phase), 2, function(x) x/sum(x))))
dc_feat <- prop.table(as.matrix(dc_feat), 1)
# chi-square test
chisq.test(dc_feat, simulate.p.value=T, B=9999)

## ground ##
# put depths on 
ground$depth <- depths$Depth.bs..m[match(ground$Spit,depths$Spit.no)]
# put phases on 
ground$phase <- flakes$phase[match(ground$Spit,flakes$Spit)]
ground[] <- lapply(ground, as.character) # change factor to character
ground[ground == 'x'] <-  1 # replace x with 1
ground[] <- lapply(ground, as.numeric) # change char to num
# get row sums
rs <- rowSums(ground[,2:(ncol(ground)-2)], na.rm = TRUE )
# get row sums by phase
rs_phase <- data.frame(rs, phase = ground$phase)
# % that have a type present
# phases with no types at all
rs_phase[rs_phase$rs == 0,]
# set any non-zero to 1
rs_phase$rs <- ifelse(rs_phase$rs != 0, 1, rs_phase$rs)
# yes, more with zero here...
dc_gr <- data.frame(t(apply(table(rs_phase), 2, function(x) x/sum(x))))
dc_gr <- prop.table(as.matrix(dc_gr), 1)
# chi-square test
chisq.test(dc_gr, simulate.p.value=T, B=9999)

# put them together
lst <- list(cores = dc_core,  retouch = dc_retouch, types = dc_types, features = dc_feat,  ground = dc_gr)
df <- ldply(lst, data.frame)
df$phase <- rep(seq_along(unique(na.omit(flakes$phase))), length(lst))

# plot proportion of spits having a techno-type present
df_m <- melt(df, id.var = c('phase', '.id'))
df_m$variable <- ifelse(df_m$variable == 'X0', 'Absent', 'Present')
ggplot(df_m, aes(as.factor(phase), value, fill = variable)) +
  geom_bar(stat="identity") +
  facet_grid(. ~ .id) +
  theme_minimal() +
  xlab("Depositional phase") +
  ylab("Proportion of spits") +
  scale_fill_discrete(name="")
ggsave("figures/Jeremalai-techno-types.png")

# now we can see, let's explore some of the minor patterns...
# what are the counts of each class in each phase?

# more about features
l <- lapply(features[,2:(ncol(features)-2)], function(i)  aggregate( i ~ phase, features, sum, na.rm = TRUE))
df <- do.call(rbind.data.frame, l)
df$name <- unlist(lapply(1:length(l), function(i) rep(names(l)[i], nrow(l[[i]]))))
dcast(phase ~ name, value.var = 'i', data = df)

# more about ground
l <- lapply(ground[,2:(ncol(ground)-2)], function(i)  aggregate( i ~ phase, ground, sum, na.rm = TRUE))
df <- do.call(rbind.data.frame, l)
df$name <- unlist(lapply(1:length(l), function(i) rep(names(l)[i], nrow(l[[i]]))))
dcast(phase ~ name, value.var = 'i', data = df)

# more about retouch
l <- lapply(retouch[,2:(ncol(retouch)-2)], function(i)  aggregate( i ~ phase, retouch, sum, na.rm = TRUE))
df <- do.call(rbind.data.frame, l)
df$name <- unlist(lapply(1:length(l), function(i) rep(names(l)[i], nrow(l[[i]]))))
dcast(phase ~ name, value.var = 'i', data = df)

# more about types
l <- lapply(types[,2:(ncol(types)-2)], function(i)  aggregate( i ~ phase, types, sum, na.rm = TRUE))
df <- do.call(rbind.data.frame, l)
df$name <- unlist(lapply(1:length(l), function(i) rep(names(l)[i], nrow(l[[i]]))))
dcast(phase ~ name, value.var = 'i', data = df)

# more about cores
l <- lapply(cores[,2:(ncol(cores)-2)], function(i)  aggregate( i ~ phase, cores, sum, na.rm = TRUE))
df <- do.call(rbind.data.frame, l)
df$name <- unlist(lapply(1:length(l), function(i) rep(names(l)[i], nrow(l[[i]]))))
dcast(phase ~ name, value.var = 'i', data = df)

# obsidian?
all$phase <- makephases(all$Spit)
table(data.frame(phase = all$phase, rm = all$Material))

# summary table of all techno-types. This will give a count of spits
# in each phase that contain at least one artefact in the category.
summaryt <- (t(cbind(
  ddply(retouch[,2:(ncol(retouch))], "phase", numcolwise(sum, na.rm = TRUE)),
  ddply(features[,2:(ncol(features))], "phase", numcolwise(sum, na.rm = TRUE)),
  ddply(ground[,2:(ncol(ground))], "phase", numcolwise(sum, na.rm = TRUE)),
  ddply(types[,2:(ncol(types))], "phase", numcolwise(sum, na.rm = TRUE)),
  ddply(cores[,2:(ncol(cores))], "phase", numcolwise(sum, na.rm = TRUE))
)))
summaryt <- summaryt[!(row.names(summaryt) %in% c("depth", "phase")),]
# compute proportions so we have the proportion of spits in each phase 
# that contains at least one of each class. 
spits_per_phase <- aggregate(Spit ~ phase, data = cores, length)
summaryt_props <- round(t(t(summaryt) / spits_per_phase$Spit),2)

# write table to csv to put into word doc
write.csv(summaryt_props, "techno_types_table.csv")
# this is table 5

```

#########################################################################
### summary tables by phase

# complete chert flakes
Mean±Standard deviation
Length (mm)
Width (mm)
Thickness (mm)
Mass (g)
dorsal flake scars
dorsal cortex
retouched (%)
retouch length


DT <- data.table(flakes)
fl_means <- DT[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Length', 'Width', 'Thick', 'Weight')] 
fl_sds <- DT[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Length', 'Width', 'Thick', 'Weight') ]
dfs_means <- DT[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('NoDS')]
dfs_sds <- DT[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('NoDS')]
ctx_means <- DT[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Cortex')]
ctx_sds <- DT[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Cortex')]
rt_means <- DT[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Retlen')]
rt_sds <- DT[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Retlen')]
rt_count <- aggregate(Retlen ~ phase, flakes, length)
rt_prop <- matrix(rt_tab['prop',], ncol = ncol(rt_tab))

# get proportion of all flakes with OHR in each phase
flakes$ohr <- ifelse(flakes$Overhang == 'Yes', 1, 0)
oh_count <- aggregate(ohr ~ phase, flakes, function(x) sum(x))
oh_prop <- aggregate(ohr ~ phase, flakes, function(x) c(ohr_mean = sum(x)/length(x)))

# platform dimensions
pw_means <- DT[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Platwid')]
pw_sds <- DT[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Platwid')]
pth_means <- DT[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Platthic')]
pth_sds <- DT[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Platthic')]

summary_flakes <- list(
t(fl_means),
t(fl_sds),
t(pw_means),
t(pw_sds),
t(pth_means),
t(pth_sds),
t(dfs_means),
t(dfs_sds),
t(ctx_means),
t(ctx_sds),
t(rt_means),
t(rt_sds),
t(rt_count)  ,
rt_prop, 
t(oh_count),
t(oh_prop)
)

# total n of flakes per phase
aggregate(Spit ~ phase, flakes, length)


# make one big table
flakes_table <- data.table(ldply(summary_flakes, data.frame))
flakes_table <- data.frame(sapply(flakes_table, as.character))
# write csv
write.csv(flakes_table, "flakes_table.csv")

# cores
all$phase <- makephases(all$Spit)
DTc <- data.table(all[all$Artclas == 'Core', ])
co_means <- DTc[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Length', 'Width', 'Thick', 'Weight')] 
co_sds <- DTc[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Length', 'Width', 'Thick', 'Weight') ]
coscrs_means <- DTc[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('NoDS')]
coscrs_sds <- DTc[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('NoDS')]
coctx_means <- DTc[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Cortex')]
coctx_sds <- DTc[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Cortex')]
corot_means <- DTc[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Corerot')]
corot_sds <- DTc[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Corerot')]
# platform dimensions
copw_means <- DTc[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Platwid')]
copw_sds <- DTc[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Platwid')]
copth_means <- DTc[order(phase), lapply(.SD, mean, na.rm=TRUE), by=phase, .SDcols=c('Platthic')]
copth_sds <- DTc[order(phase), lapply(.SD, sd, na.rm=TRUE), by=phase, .SDcols=c('Platthic')]

summary_cores <- list(
  t(co_means),
  t(co_sds),
  t(coscrs_means),
  t(coscrs_sds),
  t(coctx_means),
  t(coctx_sds),
  t(corot_means),
  t(corot_sds),
  t(copw_means),
  t(copw_sds),
  t(copth_means),
  t(copth_sds)
)
# make one big table
cores_table <- ldply(summary_cores, data.frame)
write.csv(cores_table, "core_table.csv")


# total n of flakes per phase
aggregate(Spit ~ phase, DTc, length)









