  
    # assign PFTs to Carissa's data
    site<-"carissa"
    year<-"2023"
    data_out_path<- "/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/intermediate"
    px_thresh <- 2
    pft_reference<-readr::read_csv("/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/pft_reference.csv")
    
# read in shp
shp_path <- "/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/raw/inventory/carissa/2023_NEON_species/carissa_shapes_raw_filteredbyAIS.shp"
shp_sf <- sf::read_sf(shp_path) %>%   
    filter(!is.na(species)) %>%
    # Use NEON taxonID table taht I downloaded to raw/inventory
    # https://data.neonscience.org/taxonomic-lists?taxonTypeCode=PLANT
    mutate(taxonID = case_when(
        species == "Broad-leaved lupine (Lupinus latifolia var. columbianus)" | 
            species == "Lupinus sp." ~ "LUPINSPP",
        species == "Abies concolor (Gord. & Glend.) Lindl. ex Hildebr. var. lowiana (Gord. & Glend.) Lemmon (white fir)" |
            species == "Abies magnifica A. Murray bis (CA red fir)" ~ "ABIES",
        species == "White leaf manzanita (arctostaphulos viscida)" |
            species == "Arctostaphylos patula Greene (greenleaf manzanita)" |
            species == "White leaf manzanita (arctostaphulos viscida)"  |
            species == "Pinemat manzanita (Arctostaphylos nevadensis)"  |
            species == "Arctostaphylos nevadensis\n\n\n" |
            species == "Arctostaphylos patula Greene (greenleaf manzanita)" ~ "ARCTO3SPP",
        species == "Chrysolepis sempervirens (Kellogg) Hjelmqvist (Sierran Chinkapin)" |
            species == "Bush chinquapin (Chrysolepis (Castonopsis) sempervirens)" ~ "CHSE11",
        species == "Arroyo or Lemmon\xd5s willow" |
            species == "Lemmon\xd5s willow (Salix Lemmonii - see notes)" |
            species == "Pacific Willow or Arroyo Willow" ~ "SALIX",
        species == "Vaccinium sp." ~ "VACCI",
        species == "Sambucus nigra L. (elderberry)" ~ "SANI4",
        species == "Artemesia californica" ~ "ARCA11",
        species == "coffee berry" ~ "PSYCH",
        species == "Quercus douglasii Hook. & Arn. (blue oak)" ~ "QUDO",
        species == "Pinus sabiniana Douglas ex Douglas (CA foothill pine)" ~ "PISA2",
        species == "Ceanothus cuneatus (Hook.) Nutt. (buckbrush)" |
            species == "Ceanothus? thyrsiflorus - blue blossom" ~ "CECU",
        species == "Quercus wislizeni A. DC. (interior live oak)" ~ "QUWI2",
        species == "Ceanothus leucodermis Greene (chaparral whitethorn)" ~ "CELE2",
        species == "Aesculus californica (Spach) Nutt. (CA buckeye)" ~ "AECA",
        species == "Rhamnus ilicifolia Kellogg (hollyleaf redberry)" ~ "RHIL",
        species == "90% thimbleberry, 10% some kind of pea" ~ "RUBUS",
        species == "Scouler\xd5s willow \nSalix scouleriana" ~ "SASC",
        species == "Black oak \n\nQuecus kolleggii" |
            species == "Quercus kelloggii Newberry (CA black oak)" ~ "QUKE",
        species == "Incense cedar" |
            species == "Calocedrus decurrens (Torr.) Florin (CA incense cedar)" ~ "CADE27",
        species == "Quaking aspen (populous tremuloides)" ~ "POTR5",
        species == "Broad-leaved lotus (Lotus crassifolius)" ~ "LOCR",
        species == "Western juniper (Juniperous occidentalis)" ~ "JUOC",
        species == "Black Cottonwood (populus balsamifera)" ~ "POBA2",
        species == "Sequoiadendron giganteum" ~ "SEGI2",
        species == "Chamaebatia foliolosa Benth. (mountain misery)" ~ "CHFO",
        species == "Mountain alder (Alnus incana spp. tenuifolia)" ~ "ALIN2",
        species == "Wax currant (Ribes cereum)" ~ "RICE",
        species == "Pinus lambertiana Douglas (sugar pine)" ~ "PILA",
        species == "Prunus emarginata (Douglas ex Hook.) D. Dietr. (bitter cherry)" ~ "PREM",
        species == "Pinus jeffreyi Balf. (Jeffrey pine)" ~ "PIJE",
        species == "Pinus contorta Douglas ex Loudon (lodgepole pine)" ~ "PICO",
        species == "Fern" ~ "ASPLE",
        .default="2PLANT")) %>%
    # join pft_reference
    left_join(pft_reference %>% dplyr::select(c(taxonID, pft)), by=join_by(taxonID)) #%>%
    
    #save shp again
sf::write_sf(shp_sf,"/Users/AISpiers/dev/RS-PRISMATIC/preprocessing/data/intermediate/carissa/2023/training/tree_crowns_training.shp")
