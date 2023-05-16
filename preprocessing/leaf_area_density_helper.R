library(leafR)
options(stringsAsFactors = FALSE)

calc_leaf_area_density <- function(laz_path) {
    # Calculate LAD from voxelization
    lad_voxels = lad.voxels(laz_path,
                            grain.size = 2)

    lad_profile = lad.profile(lad_voxels)
    return(lad_profile)
}