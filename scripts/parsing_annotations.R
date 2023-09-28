### annotations
annotTab <- read_csv("/Users/mo/mmorando@ucsc.edu - Google Drive/My Drive/data/amplicon_review/all_studies/Thal_RwrkSpc/annotTab.csv")

annotationsRA <- annotTab

### phylo file
# phyloAll = read_csv('all_studies/master_annotation/phylogenies/nifH_Global_ALL_classification_Oct2022.txt')%>%
#   rename(AUID = ASV,
#          ConsClassKTK = `Consensus Classification`,
#          gamma_cluster = `gamma cluster`)


phyloNCD <- read_csv("all_studies/master_annotation/phylogenies/nifH_Global_NCD_classification_Oct2022.txt") %>%
    rename(
        AUID = ASV,
        ConsClassKTK = `Consensus Classification`,
        gamma_cluster = `gamma cluster`
    )


# phyloCyano = read_csv('all_studies/master_annotation/phylogenies/nifH_Global_1B_classification_Jan2023.txt') %>%
#   rename(AUID = ASV,
#          ConsClassKTK = `Consensus Classification`)

phyloCyano <- read_csv("all_studies/master_annotation/phylogenies/nifH_Global_1B_classification_25Jan2023.txt") %>%
    rename(
        AUID = ASV,
        ConsClassKTK = `Consensus Classification`
    )

# phyloNCD %>%
#   filter(ConsClassKTK=='AUID.17586') %>%
#   view()

### fix the file bc some names are incorrect
phyloNCD <- phyloNCD %>%
    mutate(
        ConsClassKTK = if_else(ConsClassKTK == "oligo1_A1_ATCTCGCTTCTTT", "UCYN-A1",
            if_else(ConsClassKTK == "oligo3_A2_ATTCTATTTTCTT", "UCYN-A2",
                if_else(ConsClassKTK == "oligo2_A3_GCTCTATCCTTCA", "UCYN-A3",
                    if_else(ConsClassKTK == "oligo4_A4_ACTCTATTTCCCT", "UCYN-A4", ConsClassKTK)
                )
            )
        ),
        ConsClassKTK = if_else(ConsClassKTK == "Gamma", "GammaA_Lang", ConsClassKTK),
        ConsClassKTK = if_else(ConsClassKTK == "AUID.17586", "unknown1G", ConsClassKTK)
    ) # %>%
# filter(ConsClassKTK=='AUID.17586') %>%
# view()


phyloCyano <- phyloCyano %>%
    mutate(
        ConsClassKTK = if_else(ConsClassKTK == "oligo1_A1_ATCTCGCTTCTTT", "UCYN-A1",
            if_else(ConsClassKTK == "oligo3_A2_ATTCTATTTTCTT", "UCYN-A2",
                if_else(ConsClassKTK == "oligo2_A3_GCTCTATCCTTCA", "UCYN-A3",
                    if_else(ConsClassKTK == "oligo4_A4_ACTCTATTTCCCT", "UCYN-A4", ConsClassKTK)
                )
            )
        ),
        ConsClassKTK = if_else(ConsClassKTK == "Gamma", "GammaA_Lang", ConsClassKTK),
        ConsClassKTK = if_else(ConsClassKTK == "AUID.17586", "unknown1G", ConsClassKTK)
    ) # %>%
# filter(ConsClassKTK=='AUID.17586') %>%
# view()


phyloAll <- phyloCyano %>%
    bind_rows(phyloNCD)
