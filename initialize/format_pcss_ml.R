#==========================================================================================#
#==========================================================================================#
#     Script for formatting initial conditions in the FATES format.
#------------------------------------------------------------------------------------------#

#---~---
#   Fire pre-work.
#---~---
cat0("   - Read in the patch and cohort files.")
i_pss = read_csv( "/Users/AISpiers/Downloads/rs_random_plots/Inputs/patch_ic_rs_over_random.pss", show_col_types = FALSE )
i_css = read_csv( "/Users/AISpiers/Downloads/rs_random_plots/Inputs/cohort_ic_rs_over_random.css", show_col_types = FALSE )
#---~---

#---~---
#   AIS read in fire data
#---~---
cat0("   - Read in fire data.")
library(sf)
plots <- read_sf("/Users/AISpiers/Downloads/rs_random_plots/Inputs/plots_overlapping_fire.shp") %>%
    mutate(burned="y") %>% tibble() %>% rename(patch=PLOTID)
i_pss_burned <- i_pss %>% left_join(plots) %>% filter(burned=="y") %>% dplyr::select(-burned)
i_pss_not_burned <- i_pss %>% left_join(plots) %>% filter(is.na(burned)) %>% dplyr::select(-burned)
i_css_burned <- i_css %>% left_join(plots) %>% filter(burned=="y") %>% dplyr::select(-burned)
i_css_not_burned <- i_css %>% left_join(plots) %>% filter(is.na(burned)) %>% dplyr::select(-burned)

dummy = write_delim(x=i_pss_burned,file="/Users/AISpiers/Downloads/rs_random_plots/Inputs/ic_rs_over_random_burned.pss",delim=",",append=FALSE,quote="none")
dummy = write_delim(x=i_css_burned,file="/Users/AISpiers/Downloads/rs_random_plots/Inputs/ic_rs_over_random_burned.css",delim=",",append=FALSE,quote="none")

dummy = write_delim(x=i_pss_not_burned,file="/Users/AISpiers/Downloads/rs_random_plots/Inputs/ic_rs_over_random_unburned.pss",delim=",",append=FALSE,quote="none")
dummy = write_delim(x=i_css_not_burned,file="/Users/AISpiers/Downloads/rs_random_plots/Inputs/ic_rs_over_random_unburned.css",delim=",",append=FALSE,quote="none")

#---~---




#---~---
#   Reset session
#---~---
rm(list=ls())
graphics.off()
options(warn=0)
#---~---



#---~---
#   Main path settings.
#---~---
main_path   = "/Users/AISpiers/Downloads/rs_random_plots"                        # Main working path
input_path  = file.path(main_path,"Inputs")         # Path with input files.
util_path   = "/Users/AISpiers/Downloads/Util_xgao/RUtils"  # Path with utility functions.
output_path = file.path(main_path,"Outputs")        # Path with output files.
figure_path = file.path(main_path,"Figures")        # Path with output plots.
#---~---


#---~---
#   Prefix for output inventories.
#---~---
site_prefix  = "1x1pt-soaprootCAUS"
#---~---



#---~---
#   List of settings.
#---~---
initial_list = 
    tibble::tribble( ~desc               , ~in_suffix         , ~out_suffix       , ~year, ~colour  , ~f_net_def
                     #                  , "Inventory (Anna)"  , "ic_rs_over_inv"   , "anna_neon"       , 2019L, "#CC8829", 1.0
                     #                  , "Inventory (Marcos)", "ic_rsml"          , "marcos_neon"     , 2019L, "#47B2AA", 1.0
                     #, "Inventory (Field)" , "ic_inv"           , "inventory_neon"  , 2019L, "#CC8829", 1.0
                     #, "Inventory (RS)"    , "ic_rs_over_inv"   , "rsinventory_neon", 2019L, "#47B2AA", 1.0
                     , "Entire site (RS)"  , "ic_rs_over_random_burned", "prismatic_neon"  , 2019L, "#AC5CE5", 1.0 # 2.018285
                     , "Entire site (RS)"  , "ic_rs_over_random_unburned", "prismatic_neon"  , 2019L, "#AC5CE5", 1.0 # 2.018285
    )#end tibble::tribble
#---~---


#---~---
#   PFT translation.
#---~---
pft_lookup =
    tibble::tribble( ~as_name   , ~fates_id, ~a1     , ~a2  ,     ~l1,    ~l2,   ~SLA, ~dbh_max, ~colour  , ~desc
                     , "pine" , 1L        , 0.03237, 2.482, 0.01017, 2.0112,  9.351,      90., "#66CCAE", "Pine"
                     , "cedar", 2L        , 0.03237, 2.482, 0.01017, 2.0112,  9.639,      90., "#1B9E77", "Cedar"
                     , "fir"  , 3L        , 0.03237, 2.482, 0.13748, 1.4371,  9.549,      90., "#095941", "Fir"
                     , "montane_shrub", 4L        , 0.01289, 2.586, 0.01023, 1.7698, 18.500,     12.6, "#D95F02", "Shrub"
                     , "oak"  , 5L        , 0.04267, 2.488, 0.02879, 1.8926, 13.000,      65., "#7570B3", "Oak"
    )#end tibble:tribble
#---~---


#---~---
#   List of variables for generating comparisons between inventory and remote sensing.
#---~---
var_comp = 
    tibble::tribble( ~var    , ~desc                       , ~unit    , ~add0, ~mult
                     , "nplant", "Stem number density"       , "oneoha" , "0." ,        "ha.2.m2"
                     , "ba"    , "Basal area"                , "m2oha"  , "0." ,             "1."
                     , "agb"   , "Aboveground carbon density", "Mgcoha" , "0." , "kgom2.2.tonoha"
                     , "lai"   , "Leaf area index"           , "m2lom2" , "0." ,             "1."
    )#end tibble:tribble
#---~---



#---~---
#   DBH and height classes.
#---~---
dbh_cut  = c( 5.,10.,20.,35.,55.,80.)
hgt_cut  = c(2.5,5.0,10.,15.,20.,30.)
dbh0_std = 5
#---~---



#---~---
#   Plot settings.
#---~---
gg_device      = c("png") # Output devices to use (Check ggsave for acceptable formats)
gg_depth       = 300      # Plot resolution (dpi)
gg_ptsz        = 24       # Font size
gg_width       = 11.0     # Plot width (units below)
gg_height      = 8.5      # Plot height (units below)
gg_units       = "in"     # Units for plot size
gg_fleg        = 1./6.    # Fraction of plotting area dedicated for legend
#---~---




#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#            CHANGES BEYOND THIS POINT ARE FOR ADJUSTING THE INPUT FILE ONLY.              
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------




#---~---
#   Load functions and packages
#---~---
source(file.path(util_path,"load.everything.r"),chdir=TRUE)
#---~---



#---~---
#   Make output path
#---~---
dummy = dir.create(path=output_path,recursive=TRUE,showWarnings=FALSE)
dummy = dir.create(path=figure_path,recursive=TRUE,showWarnings=FALSE)
#---~---



#---~---
#   Find the number of file formats to process
#---~---
n_initial = nrow  (initial_list)
n_pft     = nrow  (pft_lookup  )
n_var     = nrow  (var_comp    )
n_device  = length(gg_device   )
#---~---


#---~---
#   DBH and height classes.
#---~---
dbh_cut    = c(0.,dbh_cut,+Inf)
hgt_cut    = c(0.,hgt_cut,+Inf)
n_dbh_cut  = length(dbh_cut)
n_hgt_cut  = length(hgt_cut)
dbh_at     = sequence(n_dbh_cut)-0.5
dbh_labels = c(dbh_cut[-n_dbh_cut],"infinity")
hgt_at     = sequence(n_hgt_cut)-0.5
hgt_labels = c(hgt_cut[-n_hgt_cut],"infinity")
#---~---


#---~---
#   Initialise output tibbles
#---~---
bypatch = NULL
bydbh   = NULL
byhgt   = NULL
#---~---


#---~---
#   Loop through output types, and fix files.
#---~---
for (i in sequence(n_initial)){
    #---~---
    #   Load information.
    #---~---
    io_desc   = initial_list$desc      [i]
    i_suffix  = initial_list$in_suffix [i]
    o_suffix  = initial_list$out_suffix[i]
    io_year   = initial_list$year      [i]
    f_net_def = initial_list$f_net_def [i]
    cat0(" + Process data for ",tolower(io_desc)," initialisation.")
    #---~---
    
    
    #---~---
    #   Set patch and cohort files.
    #---~---
    i_pss_base = paste0(i_suffix,".pss")#paste0("patch_" ,i_suffix,".pss")
    i_css_base = paste0(i_suffix,".css")#paste0("cohort_",i_suffix,".css")
    o_pss_base = paste0(site_prefix,"_",o_suffix,"_",io_year,".pss")
    o_css_base = paste0(site_prefix,"_",o_suffix,"_",io_year,".css")
    i_pss_file = file.path(input_path ,i_pss_base)
    i_css_file = file.path(input_path ,i_css_base)
    o_pss_file = file.path(output_path,o_pss_base)
    o_css_file = file.path(output_path,o_css_base)
    #---~---
    
    
    #---~---
    #   Read in the patch and cohort files.
    #---~---
    cat0("   - Read in the patch and cohort files.")
    i_pss = read_csv( i_pss_file, show_col_types = FALSE,col_names = T )
    i_css = read_csv( i_css_file, show_col_types = FALSE,col_names=T )
    #---~---

    
    #---~---
    #   Validate files
    #---~---
    cat0("   - Validate files.")
    is_pcss_match = all(i_css$patch %in% i_pss$patch)
    if (! is_pcss_match){
        invalid_patch = i_css$patch[! (i_css$patch %in% i_pss$patch)]
        
        cat0("---~---"                                                                 )
        cat0(" FATAL ERROR!"                                                           )
        cat0("---~---"                                                                 )
        cat0("  The following patches exist in cohort files but not in the patch file:")
        for (p in seq_along(invalid_patch)){
            cat0(" + ",invalid_patch[p],".")
        }#end for (p in seq_along(invalid_patch))
        cat0("---~---"                                                                 )
        stop(" All patches in cohort file must appear in the patch file.")
    }#end if (! is_pcss_match)
    #---~---
    
    
    #---~---
    #   Validate and standardise PFTs.
    #---~---
    cat0("   - Validate PFTs.")
    if (is.character(i_css$pft)){
        is_pft_fine = all(i_css$pft %in% pft_lookup$as_name)
        
        
        #---~---
        #   Report invalid PFTs.
        #---~---
        if (! is_pft_fine){
            invalid_pft = sort(unique(i_css$pft[! (i_css$pft %in% pft_lookup$as_name)]))
            
            cat0("---~---"                                                           )
            cat0(" FATAL ERROR!"                                                     )
            cat0("---~---"                                                           )
            cat0("  The following PFTs exist in cohort files but are not recognised:")
            for (p in seq_along(invalid_pft)){
                cat0(" + ",invalid_pft[p],".")
            }#end for (p in seq_along(invalid_pft))
            cat0("---~---"                                                           )
            stop(" All PFTs in cohort file must match a value in \"pft_lookup\"."    )
        }#end if (! is_pft_fine)
        #---~---
        
        #---~---
        #   Substitute names with the FATES id
        #---~---
        i_css$pft = pft_lookup$fates_id[match(i_css$pft,pft_lookup$as_name)]
        #---~---
        
    }else{
        is_pft_fine = all(i_css$pft %in% pft_lookup$fates_id)
        if (! is_pft_fine){
            invalid_pft = sort(unique(i_css$pft[! (i_css$pft %in% pft_lookup$fates_id)]))
            
            cat0("---~---"                                                           )
            cat0(" FATAL ERROR!"                                                     )
            cat0("---~---"                                                           )
            cat0("  The following PFTs exist in cohort files but are not recognised:")
            for (p in seq_along(invalid_pft)){
                cat0(" + ",invalid_pft[p],".")
            }#end for (p in seq_along(invalid_pft))
            cat0("---~---"                                                           )
            stop(" All PFTs in cohort file must match a value in \"pft_lookup\"."    )
        }#end if (! is_pft_fine)
    }#end if (is.character(icss$pft))
    #---~---
    
    
    #---~---
    #   Validate and standardise PFTs.
    #---~---
    cat0("   - Validate cohorts.")
    if (any(i_css$n == 0.)){
        #---~---
        #   Find invalid.
        #---~---
        sel_invalid = which(i_css$n == 0.)
        i_invalid   = i_css[sel_invalid,,drop=""]
        #---~---
        
        
        #---~---
        #   Report invalid PFTs.
        #---~---
        cat0 ("---~---"                                                            )
        cat0 (" FATAL ERROR!"                                                      )
        cat0 ("---~---"                                                            )
        cat0 ("   There are cohorts with zero stem density, and this is not valid.")
        cat0 (" The first few lines are printed here. Check variable \"i_invalid\"")
        cat0 (" for the full list."                                                )
        cat0 ("---~---"                                                            )
        print(i_invalid)
        cat0 ("---~---"                                                            )
        stop(" Cohort cannot have zero density."                                   )
        #---~---
    }#end if (any(i_css$n == 0.))
    #---~---
    
    
    #---~---
    #   Remove empty patches (they are not accepted in FATES).
    #---~---
    #sscat0("   - Remove empty patches.")
    #i_pss = i_pss %>%
    #   filter( patch %in% i_css$patch )
    #---~---
    
    
    #---~---
    #   Make sure that the patch area adds to 1
    #---~---
    cat0("   - Standardise patch area.")
    i_pss$area = i_pss$area / sum(i_pss$area)
    #---~---
    
    
    #---~---
    #   Standardise names...
    #---~---
    cat0("   - Standardise names.")
    if (! "nplant" %in% names(i_css)) i_css = i_css %>% rename( nplant = n     )
    if ("index"    %in% names(i_css)) i_css = i_css %>% rename( cohort = index )
    #---~---
    
    
    #---~---
    #   For patches and cohorts, we use either the original name or the hexadecimal code.
    #---~---
    if (is.numeric(i_pss$patch )) i_pss$patch  = sprintf("0x%3.3X",i_pss$patch )
    if (is.numeric(i_css$patch )) i_css$patch  = sprintf("0x%3.3X",i_css$patch )
    if (is.numeric(i_css$cohort)) i_css$cohort = sprintf("0x%3.3X",i_css$cohort)
    #---~---
    
    
    #---~---
    #   Apply correction factor if needed.
    #---~---
    i_css$nplant = i_css$nplant * f_net_def
    #---~---
    
    
    
    #---~---
    #   Discard patches with tiny populations.
    #---~---
    no_css = i_css %>% filter(i_css$nplant <  1.e-8)
    i_css  = i_css %>% filter(i_css$nplant >= 1.e-8)
    #---~---
    
    #---~---
    #   Sort patches and cohorts
    #---~---
    i_pss = i_pss %>% arrange(patch)
    i_css = i_css %>% arrange(patch,-dbh)
    #---~---
    
    
    #---~---
    #   Create output patch and cohort structures.
    #---~---
    cat0("   - Create output structures.")
    o_pss = tibble::tibble( time  = sprintf("%4.4i"  , io_year    )
                            , patch = sprintf("%s"     , i_pss$patch)
                            , trk   = sprintf("%5i"    , 2L         )
                            , age   = sprintf("%6.1f"  , 0.         )
                            , area  = sprintf("%16.14f", i_pss$area )
                            , water = sprintf("%5i    ", 0L         )
                            , fsc   = sprintf("%10.5f" , 0.         )
                            , stsc  = sprintf("%10.5f" , 0.         )
                            , stsl  = sprintf("%10.5f" , 0.         )
                            , ssc   = sprintf("%10.5f" , 0.         )
                            , psc   = sprintf("%10.5f" , 0.         )
                            , msn   = sprintf("%10.5f" , 0.         )
                            , fsn   = sprintf("%10.5f" , 0.         )
    )#end tibble::tibble
    #---~---
    
    
    #---~---
    #   Create output cohort structure.
    #---~---
    o_css = tibble::tibble( time   = sprintf("%4.4i"  , io_year     )
                            , patch  = sprintf("%s"     , i_css$patch )
                            , cohort = sprintf("%s"     , i_css$cohort)
                            , dbh    = sprintf("%9.3f"  , i_css$dbh   )
                            , height = sprintf("%9.3f"  , i_css$height)
                            , pft    = sprintf("%5i"    , i_css$pft   )
                            , nplant = sprintf("%16.10f", i_css$nplant)
                            , bdead  = sprintf("%9.3f"  , 0.          )
                            , balive = sprintf("%9.3f"  , 0.          )
                            , avgRg  = sprintf("%9.3f"  , 0.          )
    )#end tibble:tibble
    #---~---
    
    
    
    #---~---
    #   Write the output files
    #---~---
    cat0("   - Write formatted files.")
    dummy = write_delim(x=o_pss,file=o_pss_file,delim=" ",append=FALSE,quote="none")
    dummy = write_delim(x=o_css,file=o_css_file,delim=" ",append=FALSE,quote="none")
    #---~---
    
    
    
    #---~---
    #   Find derived quantities.
    #---~---
    cat0("   - Find derived quantities.")
    i_css$ba  = 0.25 * pi * i_css$dbh^2
    i_pft     = match(i_css$pft,pft_lookup$fates_id)
    i_css$agb = pft_lookup$a1[i_pft] * i_css$dbh ^ pft_lookup$a2[i_pft]
    i_css$lai = 
        ( i_css$nplant * pft_lookup$SLA[i_pft] * pft_lookup$l1[i_pft]
          * pmin(i_css$dbh,pft_lookup$dbh_max[i_pft])^pft_lookup$l2[i_pft] )
    #---~---
    
    
    
    #---~---
    #   Aggregate data by patch.
    #---~---
    cat0("   - Aggregate data by patch.")
    i_agg = i_css                                                                      %>%
        filter( dbh %ge% dbh0_std )                                                     %>%
        mutate( ba  = nplant * ba , agb = nplant * agb )                                %>%
        group_by(patch)                                                                 %>%
        summarise( nplant = sum(nplant), ba = sum(ba), agb = sum(agb), lai = sum(lai) ) %>%
        ungroup()
    i_match = match(i_agg$patch,i_pss$patch)
    #---~---
    
    
    #---~---
    #   Populate the patch properties.
    #---~---
    cat0("   - Assign aggregated data to patches.")
    i_pss = i_pss %>%
        mutate( nplant = 0., ba = 0., agb = 0., lai = 0.)
    i_pss$nplant[i_match] = i_agg$nplant
    i_pss$ba    [i_match] = i_agg$ba
    i_pss$agb   [i_match] = i_agg$agb
    i_pss$lai   [i_match] = i_agg$lai
    #---~---
    
    
    
    #---~---
    #   Create data for binning by DBH and height class.
    #---~---
    cat0("   - Assign DBH and height categories.")
    i_work         = i_css
    i_match        = match(i_work$patch,i_pss$patch)
    i_work$area    = i_pss$area[i_match]
    i_work$dbh_bin = as.integer(cut(i_work$dbh,breaks=dbh_cut,right=FALSE))
    i_work$hgt_bin = as.integer(cut(i_work$height,breaks=hgt_cut,right=FALSE))
    #---~---
    
    
    
    
    
    #---~---
    #   Aggregate data by patch, dbh class and height class
    #---~---
    cat0("   - Define aggregated tibbles.")
    i_bypatch = i_pss                                                      %>%
        mutate( init = i )                                                  %>%
        dplyr::select(c(init,patch,nplant,ba,agb,lai))
    i_bydbh   = i_work                                                     %>%
        mutate( nplant = nplant * area, lai = lai * area )                  %>%
        mutate( ba = nplant * ba, agb = nplant * agb )                      %>%
        group_by(pft,dbh_bin)                                               %>%
        summarise(nplant=sum(nplant),ba=sum(ba),agb=sum(agb),lai=sum(lai))  %>%
        ungroup()                                                           %>%
        mutate( init = i )                                                  %>%
        dplyr::select(c(init,pft,dbh_bin,nplant,ba,agb,lai))
    i_byhgt   = i_work                                                     %>%
        mutate( nplant = nplant * area, lai = lai * area )                  %>%
        mutate( ba = nplant * ba, agb = nplant * agb )                      %>%
        group_by(pft,hgt_bin)                                               %>%
        summarise(nplant=sum(nplant),ba=sum(ba),agb=sum(agb),lai=sum(lai))  %>%
        ungroup()                                                           %>%
        mutate( init = i )                                                  %>%
        dplyr::select(c(init,pft,hgt_bin,nplant,ba,agb,lai))
    #---~---
    
    
    #---~---
    #   Append data.
    #---~---
    cat0("   - Append data.")
    if (is.null(bypatch)){bypatch = i_bypatch}else{bypatch = rbind(bypatch,i_bypatch)}
    if (is.null(bydbh  )){bydbh   = i_bydbh  }else{bydbh   = rbind(bydbh  ,i_bydbh  )}
    if (is.null(byhgt  )){byhgt   = i_byhgt  }else{byhgt   = rbind(byhgt  ,i_byhgt  )}
    #---~---
}#end for (i in sequence(n_initial))
#---~---




#---~---
#   Turn categories into factors and standardise units.
#---~---
bypatch = bypatch %>%
    mutate( init   = factor(x=init,levels=sequence(n_initial),labels=initial_list$desc))
bydbh = bydbh %>%
    mutate( pft    = factor(x=pft ,levels=pft_lookup$fates_id,labels=pft_lookup$desc  )
            , init   = factor(x=init,levels=sequence(n_initial),labels=initial_list$desc) )
byhgt = byhgt %>%
    mutate( pft    = factor(x=pft ,levels=pft_lookup$fates_id,labels=pft_lookup$desc  )
            , init   = factor(x=init,levels=sequence(n_initial),labels=initial_list$desc) )
for (v in sequence(n_var)){
    #---~---
    #   Load variable transformation settings.
    #---~---
    v_var  = var_comp$var[v]
    v_add0 = eval(parse(text=var_comp$add0[v]))
    v_mult = eval(parse(text=var_comp$mult[v]))
    #---~---
    
    #---~---
    #   Apply transformations.
    #---~---
    bypatch[[v_var]] = v_add0 + v_mult * bypatch[[v_var]]
    bydbh  [[v_var]] = v_add0 + v_mult * bydbh  [[v_var]]
    byhgt  [[v_var]] = v_add0 + v_mult * byhgt  [[v_var]]
    #---~---
}#end for (v in n_var)
#---~---





#---~---
#   Define some plot settings for comparison.s
#---~---
init_colours = initial_list$colour; names(init_colours) = initial_list$desc
pft_colours  = pft_lookup$colour  ; names(pft_colours ) = pft_lookup$desc
#---~---



#---~---
#   Plot comparisons between lidar and inventory by patch.
#---~---
cat0(" + Plot patch distribution across initial condition types.")
bypatch_path = file.path(figure_path,"byPatch")
for (v in sequence(n_var)){
    #---~---
    #   Load variable transformation settings.
    #---~---
    v_var  = var_comp$var [v]
    v_desc = var_comp$desc[v]
    v_unit = var_comp$unit[v]
    #---~---
    
    
    #---~---
    #   Create output patch if needed
    #---~---
    dummy = dir.create(path=bypatch_path,recursive=TRUE,showWarnings=FALSE)
    #---~---
    
    
    #---~---
    #   Develop the histogram.
    #---~---
    gg_patch = ggplot(data=bypatch,mapping=aes(x=v_var,colour="init",fill="init"))
    gg_patch = gg_patch + geom_histogram( mapping  = aes(y=after_stat(density))
                                          #, position = "identity"
                                           , stat="count"
                                          , bins     = 10
                                          , alpha    = 1./3.
    )#end geom_historgram
    gg_patch = gg_patch + scale_colour_manual(name="Initial Conditions",values=init_colours)
    gg_patch = gg_patch + scale_fill_manual  (name="Initial Conditions",values=init_colours)
    gg_patch = gg_patch + xlab(desc.unit(desc=v_desc   ,unit=untab  [[v_unit]]))
    gg_patch = gg_patch + ylab("Density")
    gg_patch = gg_patch + theme_grey( base_size      = gg_ptsz
                                      , base_family    = "Helvetica"
                                      , base_line_size = 0.5
                                      , base_rect_size = 0.5
    )#end theme_grey
    gg_patch = gg_patch + theme( axis.text.x       = 
                                     element_text( size= gg_ptsz
                                                   , margin = unit(rep(0.35,times=4),"cm")
                                     )#end element_text
                                 , axis.text.y       = 
                                     element_text( size   = gg_ptsz
                                                   , margin = unit(rep(0.35,times=4),"cm")
                                     )#end element_text
                                 , axis.ticks.length = unit(-0.25,"cm")
                                 , legend.position   = "right"
                                 , legend.direction  = "vertical"
    )#end theme
    #---~---
    
    
    #---~---
    #   Save plots.
    #---~---
    for (d in sequence(n_device)){
        h_output = paste0("byPatch-",v_var,".",gg_device[d])
        dummy    = ggsave( filename = h_output
                           , plot     = gg_patch
                           , device   = gg_device[d]
                           , path     = bypatch_path
                           , width    = gg_width
                           , height   = gg_height
                           , units    = gg_units
                           , dpi      = gg_depth
        )#end ggsave
        
    }#end for (d in sequence(ndevice))
    #---~---
}#end for (v in n_var)
#------------------------------------------------------------------------------------------#





#---~---
#   Plot comparisons between lidar and inventory by patch.
#---~---
cat0(" + Plot size and PFT distribution across initial condition types.")
bydbh_path = file.path(figure_path,"byDBH")
for (v in sequence(n_var)){
    #---~---
    #   Load variable transformation settings.
    #---~---
    v_var  = var_comp$var [v]
    v_desc = var_comp$desc[v]
    v_unit = var_comp$unit[v]
    #---~---
    
    
    #---~---
    #   Create output patch if needed
    #---~---
    dummy = dir.create(path=bydbh_path,recursive=TRUE,showWarnings=FALSE)
    #---~---
    
    
    #---~---
    #   Develop the histogram.
    #---~---
    gg_dbh = ggplot( data    = bydbh
                     , mapping = aes( x      = "dbh_bin"
                                             , y      = v_var
                                             , colour = "pft"
                                             , fill   = "pft"
                                             , group  = "pft"
                     )#end aes_string
    )#end ggplot
    gg_dbh = gg_dbh + facet_wrap(~ init, ncol=n_initial)
    gg_dbh = gg_dbh + geom_bar(stat="identity")
    gg_dbh = gg_dbh + scale_colour_manual(values=pft_colours)
    gg_dbh = gg_dbh + scale_fill_manual  (values=pft_colours)
    gg_dbh = gg_dbh + scale_x_continuous ( breaks   = dbh_at
                                           , labels   = parse(text=dbh_labels)
    )#end scale_x_discrete
    gg_dbh = gg_dbh + labs( title = "Size and PFT distribution")
    gg_dbh = gg_dbh + xlab(desc.unit(desc="DBH classes",unit=untab$cm       ))
    gg_dbh = gg_dbh + ylab(desc.unit(desc=v_desc       ,unit=untab[[v_unit]]))
    gg_dbh = gg_dbh + theme_grey( base_size      = gg_ptsz
                                  , base_family    = "Helvetica"
                                  , base_line_size = 0.5
                                  , base_rect_size = 0.5
    )#end theme_grey
    gg_dbh = gg_dbh + theme( axis.text.x       = 
                                 element_text( size= gg_ptsz
                                               , margin = unit(rep(0.35,times=4),"cm")
                                 )#end element_text
                             , axis.text.y       = 
                                 element_text( size   = gg_ptsz
                                               , margin = unit(rep(0.35,times=4),"cm")
                                 )#end element_text
                             , axis.ticks.length = unit(-0.25,"cm")
                             , legend.position   = "right"
                             , legend.direction  = "vertical"
    )#end theme
    #---~---
    
    
    #---~---
    #   Save plots.
    #---~---
    for (d in sequence(n_device)){
        h_output = paste0("byDBH-",v_var,".",gg_device[d])
        dummy    = ggsave( filename = h_output
                           , plot     = gg_dbh
                           , device   = gg_device[d]
                           , path     = bydbh_path
                           , width    = gg_width *14/11
                           , height   = gg_height
                           , units    = gg_units
                           , dpi      = gg_depth
        )#end ggsave
        
    }#end for (d in sequence(ndevice))
    #---~---
}#end for (v in n_var)
#------------------------------------------------------------------------------------------#





#---~---
#   Plot comparisons between lidar and inventory by patch.
#---~---
cat0(" + Plot height and PFT distribution across initial condition types.")
byhgt_path = file.path(figure_path,"byHeight")
for (v in sequence(n_var)){
    #---~---
    #   Load variable transformation settings.
    #---~---
    v_var  = var_comp$var [v]
    v_desc = var_comp$desc[v]
    v_unit = var_comp$unit[v]
    #---~---
    
    
    #---~---
    #   Create output patch if needed
    #---~---
    dummy = dir.create(path=byhgt_path,recursive=TRUE,showWarnings=FALSE)
    #---~---
    
    
    #---~---
    #   Develop the histogram.
    #---~---
    gg_hgt = ggplot( data    = byhgt
                     , mapping = aes( x      = v_var
                                             , y      = "hgt_bin"
                                             , colour = "pft"
                                             , fill   = "pft"
                                             , group  = "pft"
                     )#end aes_string
    )#end ggplot
    gg_hgt = gg_hgt + facet_wrap(~ init, ncol=n_initial)
    gg_hgt = gg_hgt + geom_bar(stat="identity",orientation="y")
    gg_hgt = gg_hgt + scale_colour_manual(values=pft_colours)
    gg_hgt = gg_hgt + scale_fill_manual  (values=pft_colours)
    gg_hgt = gg_hgt + scale_y_continuous ( breaks   = hgt_at
                                           , labels   = parse(text=hgt_labels)
    )#end scale_x_discrete
    gg_hgt = gg_hgt + labs( title = "Height and PFT distribution")
    gg_hgt = gg_hgt + xlab(desc.unit(desc=v_desc          ,unit=untab[[v_unit]]))
    gg_hgt = gg_hgt + ylab(desc.unit(desc="Height classes",unit=untab$m        ))
    gg_hgt = gg_hgt + theme_grey( base_size      = gg_ptsz
                                  , base_family    = "Helvetica"
                                  , base_line_size = 0.5
                                  , base_rect_size = 0.5
    )#end theme_grey
    gg_hgt = gg_hgt + theme( axis.text.x       = 
                                 element_text( size= gg_ptsz
                                               , margin = unit(rep(0.35,times=4),"cm")
                                 )#end element_text
                             , axis.text.y       = 
                                 element_text( size   = gg_ptsz
                                               , margin = unit(rep(0.35,times=4),"cm")
                                 )#end element_text
                             , axis.ticks.length = unit(-0.25,"cm")
                             , legend.position   = "right"
                             , legend.direction  = "vertical"
    )#end theme
    #---~---
    
    
    #---~---
    #   Save plots.
    #---~---
    for (d in sequence(n_device)){
        h_output = paste0("byHeight-",v_var,".",gg_device[d])
        dummy    = ggsave( filename = h_output
                           , plot     = gg_hgt
                           , device   = gg_device[d]
                           , path     = byhgt_path
                           , width    = gg_width *14/11
                           , height   = gg_height
                           , units    = gg_units
                           , dpi      = gg_depth
        )#end ggsave
        
    }#end for (d in sequence(ndevice))
    #---~---
}#end for (v in n_var)
#------------------------------------------------------------------------------------------#

