### R code used to generate Figure 1 for Bacterial DNA invasion triggers transposable element bursts and genome expansion


setwd("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/")
#file.exists("/home/zach/Sitophiles_BUSCO_phylogeny/Synteny/SzeSor1004_tst1022.karyo.svg")

install.packages("rsvg")
library(rsvg)

### Ccar & Cnanu
rsvg_png("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/Ccar_CnanuRT_42426.karyo.svg",
         "/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/Ccar_CnanuRT_42426.karyo.png")
karyotype <- read.table("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/Ccar_CnanuRT_karyo_R.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
synteny <- read.table("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/CcCnRT_synteny.rtsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
head(karyotype)
head(synteny)
#colors for each fill:
unique_species <- unique(synteny$Species_1)
colors <- custom_palette(length(unique_species))
color_mapping <- setNames(colors, unique_species)
color_mapping[c("5", "2", "9", "3", "12", "8", "1", "4", "7", "6", "10", "13", "11")] <- c(
  "#FBB4AE",
  "#B3CDE3",
  "#CCEBC5",
  "#DECBE4",
  "#FED9A6",
  "#FFFFCC",
  "#E5D8BD",
  "#8DD3C7",
  "#FDDAEC",
  "#BEBADA",
  "#F2F2F2",
  "#FFFFB3",
  "#FB8072"
)

### Cglan & Cnanu

karyotype <- read.table("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/CnanuRT_Clan_R.rtsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
synteny <- read.table("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/CNanuRT_Cglan_SCO_synteny2.gff", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
View(karyotype)
View(synteny)
#colors for each fill:
unique_species <- unique(synteny$Species_1)
colors <- custom_palette(length(unique_species))
color_mapping <- setNames(colors, unique_species)
synteny$fill <- color_mapping[synteny$Species_1]
synteny$fill <- gsub("#", "", synteny$fill)
?ideogram
setwd("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/")
ideogram(karyotype = karyotype, synteny = synteny, output = "CnanuRT_Cglan.karyo.svg")
#file.exists("/home/zach/Sitophiles_BUSCO_phylogeny/Synteny/SzeSor1004_tst1022.karyo.svg")

install.packages("rsvg")
library(rsvg)
rsvg_png("/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/CnanuRT_Cglan.karyo.svg",
         "/Users/zcohen5/Desktop/FAU_ZPC/Curculio_genomes/CnanuRT_Cglan.karyo.png")
