library(lidR)
options(stringsAsFactors = FALSE)

calc_leaf_area_density <- function(laz_path, dz=0.5, z0=1) {
    # Calculate LAD from voxelization
    laz_data = readLAS(laz_path)
    lad_profile = LAD(laz_data@data$Z,
                      dz=dz,
                      z0=z0)
    return(lad_profile)
}