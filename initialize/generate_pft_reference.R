# Assign a PFT to each species

library(dplyr)

# load data

soap <- readr::read_csv("data/raw/inventory/SOAP/veg_structure.csv") %>%
    mutate(growthForm_neon = ifelse(grepl("tree", growthForm), "tree",
                                   ifelse(grepl("sapling", growthForm), "tree",
                                          ifelse(grepl("shrub", growthForm), "shrub", growthForm)))) %>%
    dplyr::select(scientificName, taxonID, growthForm_neon, SOAP=siteID)

sjer <- readr::read_csv("data/raw/inventory/SJER/veg_structure.csv") %>%
    mutate(growthForm_neon = ifelse(grepl("tree", growthForm), "tree",
                                   ifelse(grepl("sapling", growthForm), "tree",
                                          ifelse(grepl("shrub", growthForm), "shrub", growthForm)))) %>%
    dplyr::select(scientificName, taxonID, growthForm_neon, SJER=siteID) 

teak <- readr::read_csv("data/raw/inventory/TEAK/veg_structure.csv") %>%
    mutate(growthForm_neon = ifelse(grepl("tree", growthForm), "tree",
                                   ifelse(grepl("sapling", growthForm), "tree",
                                          ifelse(grepl("shrub", growthForm), "shrub", growthForm)))) %>%
    dplyr::select(scientificName, taxonID, growthForm_neon, TEAK=siteID) 

trait_table <- readr::read_csv("data/raw/inventory/NEON_trait_table.csv") %>%
    dplyr::select(scientific, leafPhen=leaf.phen, growthForm_traittable=growth.form)

# carissa <- readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/carissa/neon_2023_field_work_final_ais_grouped_species.csv") %>%
#     mutate(scientificName=ifelse(Species=="OTHER",`If other species: describe here`,Species)) %>%
#     filter(!is.na(scientificName)) %>%
#     #rename(taxonID=) %>%
#     rename(siteID = `NEON Site`) %>%
#     #rename(growthForm=) %>%
#     dplyr::select(scientificName,  siteID) %>% #taxonID,growthForm
#     distinct() %>%
#     # Use NEON taxonID table taht I downloaded to raw/inventory
#     # https://data.neonscience.org/taxonomic-lists?taxonTypeCode=PLANT
#     mutate(taxonID = case_when(
#         scientificName == "Broad-leaved lupine (Lupinus latifolia var. columbianus)" | 
#             scientificName == "Lupinus sp." ~ "LUPINSPP",
#         scientificName == "Abies concolor (Gord. & Glend.) Lindl. ex Hildebr. var. lowiana (Gord. & Glend.) Lemmon (white fir)" |
#             scientificName == "Abies magnifica A. Murray bis (CA red fir)" ~ "ABIES",
#         scientificName == "White leaf manzanita (arctostaphulos viscida)" |
#             scientificName == "Arctostaphylos patula Greene (greenleaf manzanita)" |
#             scientificName == "White leaf manzanita (arctostaphulos viscida)"  |
#             scientificName == "Pinemat manzanita (Arctostaphylos nevadensis)"  |
#             scientificName == "Arctostaphylos nevadensis\n\n\n" |
#             scientificName == "Arctostaphylos patula Greene (greenleaf manzanita)" ~ "ARCTO3SPP",
#         scientificName == "Chrysolepis sempervirens (Kellogg) Hjelmqvist (Sierran Chinkapin)" |
#             scientificName == "Bush chinquapin (Chrysolepis (Castonopsis) sempervirens)" ~ "CHSE11",
#         scientificName == "Arroyo or Lemmon\xd5s willow" |
#             scientificName == "Lemmon\xd5s willow (Salix Lemmonii - see notes)" |
#             scientificName == "Pacific Willow or Arroyo Willow" ~ "SALIX",
#         scientificName == "Vaccinium sp." ~ "VACCI",
#         scientificName == "Sambucus nigra L. (elderberry)" ~ "SANI4",
#         scientificName == "Artemesia californica" ~ "ARCA11",
#         scientificName == "coffee berry" ~ "PSYCH",
#         scientificName == "Quercus douglasii Hook. & Arn. (blue oak)" ~ "QUDO",
#         scientificName == "Pinus sabiniana Douglas ex Douglas (CA foothill pine)" ~ "PISA2",
#         scientificName == "Ceanothus cuneatus (Hook.) Nutt. (buckbrush)" |
#             scientificName == "Ceanothus? thyrsiflorus - blue blossom" ~ "CECU",
#         scientificName == "Quercus wislizeni A. DC. (interior live oak)" ~ "QUWI2",
#         scientificName == "Ceanothus leucodermis Greene (chaparral whitethorn)" ~ "CELE2",
#         scientificName == "Aesculus californica (Spach) Nutt. (CA buckeye)" ~ "AECA",
#         scientificName == "Rhamnus ilicifolia Kellogg (hollyleaf redberry)" ~ "RHIL",
#         scientificName == "90% thimbleberry, 10% some kind of pea" ~ "RUBUS",
#         scientificName == "Scouler\xd5s willow \nSalix scouleriana" ~ "SASC",
#         scientificName == "Black oak \n\nQuecus kolleggii" |
#             scientificName == "Quercus kelloggii Newberry (CA black oak)" ~ "QUKE",
#         scientificName == "Incense cedar" |
#             scientificName == "Calocedrus decurrens (Torr.) Florin (CA incense cedar)" ~ "CADE27",
#         scientificName == "Quaking aspen (populous tremuloides)" ~ "POTR5",
#         scientificName == "Broad-leaved lotus (Lotus crassifolius)" ~ "LOCR",
#         scientificName == "Western juniper (Juniperous occidentalis)" ~ "JUOC",
#         scientificName == "Black Cottonwood (populus balsamifera)" ~ "POBA2",
#         scientificName == "Sequoiadendron giganteum" ~ "SEGI2",
#         scientificName == "Chamaebatia foliolosa Benth. (mountain misery)" ~ "CHFO",
#         scientificName == "Mountain alder (Alnus incana spp. tenuifolia)" ~ "ALIN2",
#         scientificName == "Wax currant (Ribes cereum)" ~ "RICE",
#         scientificName == "Pinus lambertiana Douglas (sugar pine)" ~ "PILA",
#         scientificName == "Prunus emarginata (Douglas ex Hook.) D. Dietr. (bitter cherry)" ~ "PREM",
#         scientificName == "Pinus jeffreyi Balf. (Jeffrey pine)" ~ "PIJE",
#         scientificName == "Pinus contorta Douglas ex Loudon (lodgepole pine)" ~ "PICO",
#         scientificName == "Fern" ~ "ASPLE",
#         .default="2PLANT")) %>%
#     mutate(TEAK=ifelse(siteID=="TEAK","TEAK",NA),
#            SJER=ifelse(siteID=="SJER","SJER",NA)) %>%
#     dplyr::select(-siteID)

# Join data
neon_sites_raw <- soap %>%
    full_join(sjer, by=join_by(scientificName, taxonID, growthForm_neon)) %>%
    full_join(teak, by=join_by(scientificName, taxonID, growthForm_neon)) %>%
    arrange(taxonID) %>%
    distinct() 

neon_sites <- neon_sites_raw %>%
    filter(!is.na(scientificName)) %>%
    mutate(scientific = word(scientificName, 1,2, sep=" ")) %>%
    dplyr::select(scientific, taxonID, SOAP, SJER, TEAK, growthForm_neon) %>%
    mutate(growthForm_neon = ifelse(taxonID=="CONU4","tree",
                                   ifelse(taxonID=="ARNE","shrub",growthForm_neon))) %>%
    # Otherwise get rid of rows with NA in growthForm_neon (these are duplicates)
    filter(!is.na(growthForm_neon)) %>%
    left_join(trait_table) %>%
    # By default use trait table growth form
    mutate(growthForm = ifelse(is.na(growthForm_traittable),
                               growthForm_neon,growthForm_traittable)) %>%
    # Filter out remaining unlikely growth forms - AIS search on Calscape
    mutate(growthForm = case_when(
        taxonID == "ABIES" |
            taxonID == "ABLO" |
            taxonID == "CADE27" ~ "tree",
        taxonID == "CECO" |
            taxonID == "FRCA6" |
            taxonID == "LOIN4" |
            taxonID == "RHIL" |
            taxonID == "RIBES" |
            taxonID == "RIRO" |
            taxonID == "TODI" ~ "shrub",
        taxonID == "AECA" |
            taxonID == "SALIX" ~ "shrub/tree",
        taxonID == "2PLANT" |
            taxonID == "2PLANT-S" |
            taxonID == "LUAL4" ~ NA,
        .default = growthForm  )) %>%
    dplyr::select(scientific,taxonID, SOAP, SJER, TEAK, growthForm) %>%
    distinct()
    
# Add in Carissa's dataset
neon_carissa <- neon_sites %>%
    full_join(carissa) %>%
    # Use NEON species name if it exists, otherwise Carissa's
    mutate(scientific = ifelse(is.na(scientific),scientificName,scientific)) %>%
    dplyr::select( scientific,   taxonID, SOAP, SJER, TEAK, growthForm) %>%
    distinct() %>% 
    tidyr::unite("siteID", c(SOAP,  SJER,  TEAK), na.rm = TRUE) %>%
    # Manually remove remaining duplicates
    mutate(siteID=ifelse(taxonID=="AECA", "SOAP_SJER", siteID),
           siteID=ifelse(taxonID=="QUWIF", "SJER_TEAK", siteID),
           
           siteID=ifelse(taxonID=="ARCTO3SPP", "SJER_TEAK", siteID),
           scientific=ifelse(taxonID=="ARCTO3SPP","Arctostaphylos sp.",scientific),
           
           siteID=ifelse(taxonID=="CECU", "SOAP_SJER_TEAK", siteID),
           scientific=ifelse(taxonID=="CECU","Ceanothus sp.",scientific),
           
           siteID=ifelse(taxonID=="CHFO", "SOAP_TEAK", siteID),
           scientific=ifelse(taxonID=="CHFO","Chamaebatia foliolosa",scientific),
           
           siteID=ifelse(taxonID=="LUPINSPP", "SJER_TEAK", siteID),
           scientific=ifelse(taxonID=="LUPINSPP","Lupinus sp.",scientific),
           
           siteID=ifelse(taxonID=="QUKE", "SJER_TEAK", siteID),
           scientific=ifelse(taxonID=="QUKE","Quercus kelloggii",scientific),
           
           siteID=ifelse(taxonID=="QUWI2", "SOAP_TEAK", siteID),
           scientific=ifelse(taxonID=="QUWI2", "Quercus wislizeni", scientific),
           
           scientific=ifelse(taxonID=="SANI4","Sambucus nigra",scientific),
           siteID=ifelse(scientific=="Sambucus nigra", "SOAP_SJER", siteID),
           taxonID=ifelse(scientific=="Sambucus nigra", "SANI4", taxonID),
           
           scientific=ifelse(taxonID=="SALE","Salix sp.",scientific)) %>%
    filter(!(taxonID == "CADE27" & siteID=="SOAP")) %>%
    filter(!(taxonID == "2PLANT" & siteID=="SJER")) %>%
    filter(!(taxonID == "RIRO" & siteID=="SOAP")) %>%
    filter(!(taxonID == "TODI" & siteID=="SJER")) %>%
    filter(!(taxonID == "CECU" & is.na(growthForm))) %>%
    filter(!(taxonID == "CHFO" & is.na(growthForm))) %>%
    filter(!(taxonID == "QUKE" & is.na(growthForm))) %>%
    filter(!(taxonID == "QUWI2" & is.na(growthForm))) %>%
    filter(!(taxonID == "SANI4" & is.na(growthForm))) %>%
    distinct()
    
    
pft_assignment_df <- neon_carissa %>%
    mutate(pft = case_when(
        taxonID == "ABIES" | taxonID == "ABCO" |
            taxonID == "ABLO" |
            taxonID == "ABMA" ~ "fir", #white or red fir
        taxonID == "CADE27" | taxonID == "SEGI2" ~ "cedar", #incense cedar
        taxonID == "PICO" | taxonID == "PIJE" |
            taxonID == "PILA" | taxonID == "PINACE" |
            taxonID == "PINACE" | taxonID == "PIPO" |
            taxonID == "PISA2" |
            taxonID == "PINUS" ~ "pine", #Ponderosa pine
        taxonID == "QUCH2" | taxonID == "QUKE" |
            taxonID == "QUDO" | taxonID == "QUERC" |
            taxonID == "QUWI2" ~ "oak", #target:canyon live oak (evergreen), also Black oak(deciduous), 
        taxonID == "ARNE" | taxonID == "ARPA" |
            taxonID == "ARPA6" | taxonID == "ARVIM" |
            taxonID == "CEANO" | taxonID == "CECA" |
            taxonID == "CECO" | taxonID == "CECU" |
            taxonID == "CEIN3" | taxonID == "CELE2" |
            taxonID == "CHFO" | taxonID == "CHSE11" |
            taxonID == "CRSE11" | taxonID == "DAWR2" |
            taxonID == "FRCA6" | taxonID == "FRCAC7" |
            taxonID == "LOIN4" | taxonID == "RHAMNA" |
            taxonID == "RHIL" | taxonID == "PREM" |
            taxonID == "RICE" | taxonID == "RIBES" |
            taxonID == "RIRO" | taxonID == "RIVI3" |
            taxonID == "SALIX" | taxonID == "SANI4" |
            taxonID == "CEMOG" | taxonID == "ARCTO3SPP" |
            taxonID == "SEFL3" | taxonID == "TODI" |
            taxonID == "SASC" | taxonID == "AECA" |
            taxonID == "CONU4" | taxonID == "LUPINSPP" ~ "montane_shrub", #Ceanothus sp.
        .default = "other"
    ))

readr::write_csv(pft_assignment_df, "/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/pft_reference.csv")	
