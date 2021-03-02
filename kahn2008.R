# Read coordinates
ROI <- read.csv("https://raw.githubusercontent.com/zchuri/kahn2008/main/kahn2008_HIP_ROIs.csv")

# Load 'fslr' package
if(!require("fslr")) install.packages("fslr"); library(fslr)

# Read T1w ICBM-152 (aka MNI) 1mm
t1 <- file.path(fsldir(),"data/standard/MNI152_T1_1mm_brain.nii.gz")
t1_nii <- readnii(file.path(fsldir(),"data/standard/MNI152_T1_1mm_brain.nii.gz"))
# Get XYZ volume coordinates from MNI coordinates
xyz_coor <- cbind((ROI[,1]-t1_nii@srow_x[4])*t1_nii@srow_x[1],
                   (ROI[,2]-t1_nii@srow_y[4])*t1_nii@srow_y[2],
                   (ROI[,3]-t1_nii@srow_z[4])*t1_nii@srow_z[3])

# Create spheres with 'fslmaths' with 3mm radius
roi_n <- nrow(xyz_coor)
for(ii in 1:roi_n){
  cat(paste("\nCreating sphere",ii,"\n"))
  # Create and empty NIfTI object with T1 standard dimensions
  coor_nii <- fslmaths(t1,opts="-mul 0")
  # Note that we have to add one voxel to the coordinates
  # cause 'oro.nifti' axis starts in one instead of zero (like FSL axis)
  xyz <- unlist(xyz_coor[ii,])+1
  coor_nii[xyz[1],xyz[2],xyz[3]] <- 1
  # Create sphere (5mm radious)
  fslmaths(coor_nii,outfile = paste0("ROI_",ii),opts="-kernel sphere 3 -fmean -bin")
}

# Merge all volumes
rois <- list.files(pattern = "ROI_")
for(ii in which(nchar(rois)==12)) file.rename(rois[ii],gsub("_","_0",rois[ii]))
rois <- list.files(pattern = "ROI_")
# Lenght of this objects should be 40
# Merge list
fslmerge(rois,"t","kahn_4D")

# Verify if there is any overlap between spheres
sum_nii <- fslmaths("kahn_4D",opts = paste("-Tmean -mul",roi_n))
table(sum_nii)
# As we can see there is an overlap between ROIs (164 voxels)
# So when taking the explicit ROI the volume of interest should be selected in the 4D.

# Now, the atlas can be resize to any other dimension with 'flirt_apply'
# For example to 2mm
t1_2mm <- file.path(fsldir(),"data/standard/MNI152_T1_2mm_brain.nii.gz")
# Flirt Apply
flirt_apply(infile = "kahn_4D", reffile = t1_2mm,
            initmat = paste0(fsldir(),"/etc/flirtsch/ident.mat"),
            outfile = "kahn_4D_2mm", opts = "-interp nearestneighbour")
