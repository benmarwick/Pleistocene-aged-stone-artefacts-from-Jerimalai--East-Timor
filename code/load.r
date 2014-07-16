## ---- load-data

# load data
 
# dates
dates <- read.csv("data/Jeremalai_dates.csv", as.is = TRUE)

# complete flake data
flakes <- read.csv("data/JB_Chert_Flakes_and_Retouch.csv")

# spit depths
depths <- read.csv("data/Jeremalai_spit_depths.csv")

# sediment volumes
vols <- read.csv("data/Artefact densities with soil volumes Sq B.csv", skip = 1)

# all artefacts
all <- read.csv("data/Jerimalai_All_Artefacts.csv")

# techno types
cores <- read.csv("data/Jerimalai_tech_table_cores.csv")
types <-  na.omit(read.csv("data/Jerimalai_tech_table_types.csv"))
retouch <- read.csv("data/Jerimalai_tech_table_retouch.csv")
features <- read.csv("data/Jerimalai_tech_table_features.csv")
ground <- read.csv("data/Jerimalai_tech_table_ground.csv")

# retouch indices
retouch_indices <- read.csv("data/Jerimalai_retouch_indices.csv")


# Below here refer to the original locations of these files on my computer,
# which might be good to know if I have conflicting versions

# # dates
# dates <- read.csv("F:/My Documents/My Papers/conferences/IPPA14/Jeremalai/Jeremalai_Science_dates.csv", as.is = TRUE)
# 
# # complete flake data
# flakes <- read.xls("F:/My Documents/My Papers/conferences/EurASEAA2010/Jeremalai/JB_Chert_Flakes_and_Retouch.xls")
# 
# # spit depths
# depths <- read.csv("F:/My Documents/My Papers/conferences/IPPA14/Jeremalai/Jeremalai_spit_depths.csv")
# 
# # sediment volumes
# vols <- read.xls("F:\\My Documents\\My Papers\\conferences\\Past conferences\\AAA conf\\AAA2007\\Sophie\\Artefact densities with soil volumes Sq B.xls", sheet = "Sheet1", skip = 1)
# 
# # all artefacts
# all <- read.xls("F:\\My Documents\\My Papers\\conferences\\Past conferences\\AAA conf\\AAA2007\\Sophie\\Jerimalai_All_Artefacts.xlsx", sheet = 'Square B')
# 
# # techno types
# 
# cores <- read.csv("F:/My Documents/My Papers/conferences/IPPA14/Jeremalai/Jerimalai_tech_table_cores.csv")
# types <-  na.omit(read.csv("F:/My Documents/My Papers/conferences/IPPA14/Jeremalai/Jerimalai_tech_table_types.csv"))
# retouch <- read.csv("F:/My Documents/My Papers/conferences/IPPA14/Jeremalai/Jerimalai_tech_table_retouch.csv")
# features <- read.csv("F:/My Documents/My Papers/conferences/IPPA14/Jeremalai/Jerimalai_tech_table_features.csv")
# ground <- read.csv("F:/My Documents/My Papers/conferences/IPPA14/Jeremalai/Jerimalai_tech_table_ground.csv")