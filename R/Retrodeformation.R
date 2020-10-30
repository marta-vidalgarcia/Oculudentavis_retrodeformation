#### 0. Libraries ####

setwd("~/Documents/GITHUB_repos/Oculudentavis_retrodeformation")

library(geomorph)
library(Morpho)
library(rgl)
library(Rvcg)
library(viridis)
source("R/ShapeRotator.R")

#### 1. data import - Specimen 1 (bird-like) ####
# Landmark data
lm_bird <- read.csv("data/Landmarks_bird_FINAL.csv", skip = 1)[,4:6]
lm_bird

# Original mesh
mesh_bird <- file2mesh("data/bird_lowRes.ply")
str(mesh_bird)


#### 1.1. Plot all landmarks and mesh ####
open3d(windowRect = c(20, 30, 800, 800))
shade3d(mesh_bird, color="gray", alpha=0.9)
plot3d(lm_bird, aspect="iso", type="s", size=0.3, col="blue", add=T)
# plot3d(lm_bird[470:478,], aspect="iso", type="s", size=0.3, col="blue", add=T)
# plot3d(lm_bird[479:487,], aspect="iso", type="s", size=0.3, col="green", add=T)
# plot3d(lm_bird[956:964,], aspect="iso", type="s", size=0.3, col="red", add=T)
# plot3d(lm_bird[965:973,], aspect="iso", type="s", size=0.3, col="yellow", add=T)
rgl.close()


#### 1.2. Plot only LMs used in retrodeformation ####
NUM_premaxilla_LEFT <- c(470:550)
NUM_premaxilla_RIGHT <- c(551:631)
NUM_maxilla_LEFT <- c(956:1036)
NUM_maxilla_RIGHT <- c(1037:1117)
NUM_nasal_LEFT <- c(632:712)
NUM_nasal_RIGHT <- c(713:793)

lm_bird_subset <- lm_bird[c(NUM_premaxilla_LEFT, NUM_premaxilla_RIGHT, NUM_maxilla_LEFT,
                                                         NUM_maxilla_RIGHT, NUM_nasal_LEFT, NUM_nasal_RIGHT), ]

open3d(windowRect = c(20, 30, 1800, 800))
shade3d(mesh_bird, color="gray", alpha=0.9)
plot3d(lm_bird_subset, aspect="iso", type="s", size=0.2, col="blue", add=T)
rgl.snapshot("figs/bird_original_LMs_RIGHT.png", top = TRUE) 
rgl.snapshot("figs/bird_original_LMs_LEFT.png", top = TRUE) 
rgl.snapshot("figs/bird_original_LMs_TOP.png", top = TRUE) 
rgl.close()
#### 2. Retrodeformation maxilla & premaxilla - ShapeRotator ####
#### 2.1. Landmark selection, angle, translation to c(0, 0, 0) and rotation to the xz-plane ####
# landmarks
land.a = 1 # distal LM on the nasal ridge
land.b = 46 # last point of the nasal ridge, between the eyes
land.c = 300 # first landmark on the palatine

# angle
angle = 0

# translate
lm_bird_a <- array(data = cbind(as.matrix(lm_bird), as.matrix(lm_bird)), dim = c(dim(lm_bird)[1], dim(lm_bird)[2], 2))

dimnames(lm_bird_a)[[3]] <- c("bird1", "bird2") # duplication of the data as ShapeRotator works with arrays only
str(lm_bird_a)
class(lm_bird_a)

bird_t <- translate(lm_bird_a, land.a)

bird_t[1,,1]

# dataset
data.1 = bird_t

# Rotation
rot.array.1 <- array(data = NA, dim = c(dim(data.1)[1], dim(data.1)[2], dim(data.1)[3]),dimnames = list(NULL, NULL, c(dimnames(data.1)[[3]])));

# This is the angle which land.b will face to land.e
angle.v <-vector.angle(angle)


#Rotating for data.1
for( i in 1:(dim(data.1)[3]) ){
  # Untranslated specimen?
  if ( ! isTRUE(  all.equal((data.1[,,i])[land.a,] , c(0,0,0)) ) ) {
    # warning (sprintf("Landmark A is not at the origin (data.1[,,%d])[%d,] = (%f,%f,%f)?", i, land.a, (data.1[,,i])[land.a,1], (data.1[,,i])[land.a,2], (data.1[,,i])[land.a,3]), immediate. = immediate_on)
  }
  
  # First, we rotate so that land.b rests on the axis (0,1,0).
  rot.array.1[,,i] <- rotate.orientation(data.1[,,i], land.b, c(0,1,0))
  
  #Sanity check: we need to have rot.array.1[,,i][land.b,] = (0, y, 0)
  if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.b,1] , 0)) || ! isTRUE(  all.equal( (rot.array.1[,,i])[land.b,3] , 0)) ) {
    warning (sprintf("First rotation gone wrong in (rot.array.1[,,%d])[%d,] = (%f,%f,%f)?", i, land.b, (rot.array.1[,,i])[land.b,1], (rot.array.1[,,i])[land.b,2], (rot.array.1[,,i])[land.b,3]), immediate. = immediate_on)
  }
  
  # Rotate landmark C to the X-Y plane - we do so by computing the cross product between the unit.landmark.c.pos.proj vector and (-1,0,0),
  # which in the X-Y plane forms an angle given by the angle between these two vectors:
  unit.landmark.c.pos.proj = unit_3D( c( (rot.array.1[,,i])[land.c,1], 0, (rot.array.1[,,i])[land.c,3]))	# Unit vector with landmark.c.pos projected to the X-Z plane.
  rotmat <- rotmat_3D( cross_3D(unit.landmark.c.pos.proj, c(-1,0,0)), angle_3D(unit.landmark.c.pos.proj, c(-1,0,0)) )
  rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
  
  
  # Let us now check that the landmark C is in the right spot. That is, we want C_x < B_x, if not, we rotate by angle pi.
  if ( isTRUE((rot.array.1[,,i])[land.b,1] < (rot.array.1[,,i])[land.c,1]) ) {
    rotmat <- rotmat_3D( c(0,1,0), pi)
    rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
  }
  
  # Sanity check:
  if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.c,3] , 0)) || isTRUE((rot.array.1[,,i])[land.b,1] < (rot.array.1[,,i])[land.c,1])) {
    warning (sprintf("Second rotation gone wrong in (rot.array.1[,,%d])[%d,] = (%f,%f,%f)?", i, land.c, (rot.array.1[,,i])[land.c,1], (rot.array.1[,,i])[land.c,2], (rot.array.1[,,i])[land.c,3]), immediate. = immediate_on)
  }
  
  # Now, we rigidly move this object to sit at an angle given by "angle" to the X-axis.
  rot.array.1[,,i] <- rotate.orientation(rot.array.1[,,i], land.b, angle.v)
}

rotated_bird <- rot.array.1[,,1]

str(rotated_bird)

rotated_bird[1,]
rotated_bird[46,]
rotated_bird[300,] # please note that the ventral side is on the top

write.csv(rotated_bird, "data/rotated_bird.csv")

#### 2.3. RETRODEFORMATION ON THE xz-plane ####
# Patches to deform on the z-axis
NUM_premaxilla_LEFT <- c(470:550)
NUM_premaxilla_RIGHT <- c(551:631)
NUM_maxilla_LEFT <- c(956:1036)
NUM_maxilla_RIGHT <- c(1037:1117)
NUM_nasal_LEFT <- c(632:712)
NUM_nasal_RIGHT <- c(713:793)


# Retrodeformation factor for z
retro_vec <- seq(1.75, 1.5, length.out = 18) # differential deformation as the snout got more deformed than the maxillas closer to the vault

# Retrodeformed patches
### degree of retrodeformation goes in order from the end of the maxilla starting at 1.05 until 2.25 at the start of the premaxilla
# They go in groups of 9 (rows) for premaxilla, and the following 9 for maxilla)
retrodeformed_premax_L <- c(rotated_bird[NUM_premaxilla_LEFT[1:9], 3] * retro_vec[1], rotated_bird[NUM_premaxilla_LEFT[10:18], 3] * retro_vec[2], 
                            rotated_bird[NUM_premaxilla_LEFT[19:27], 3] * retro_vec[3], rotated_bird[NUM_premaxilla_LEFT[28:36], 3] * retro_vec[4],
                            rotated_bird[NUM_premaxilla_LEFT[37:45], 3] * retro_vec[5], rotated_bird[NUM_premaxilla_LEFT[46:54], 3] * retro_vec[6],
                            rotated_bird[NUM_premaxilla_LEFT[55:63], 3] * retro_vec[7], rotated_bird[NUM_premaxilla_LEFT[64:72], 3] * retro_vec[8], 
                            rotated_bird[NUM_premaxilla_LEFT[73:81], 3] * retro_vec[9])

retrodeformed_premax_R <- c(rotated_bird[NUM_premaxilla_RIGHT[1:9], 3] * retro_vec[1], rotated_bird[NUM_premaxilla_RIGHT[10:18], 3] * retro_vec[2], 
                            rotated_bird[NUM_premaxilla_RIGHT[19:27], 3] * retro_vec[3], rotated_bird[NUM_premaxilla_RIGHT[28:36], 3] * retro_vec[4],
                            rotated_bird[NUM_premaxilla_RIGHT[37:45], 3] * retro_vec[5], rotated_bird[NUM_premaxilla_RIGHT[46:54], 3] * retro_vec[6],
                            rotated_bird[NUM_premaxilla_RIGHT[55:63], 3] * retro_vec[7], rotated_bird[NUM_premaxilla_RIGHT[64:72], 3] * retro_vec[8], 
                            rotated_bird[NUM_premaxilla_RIGHT[73:81], 3] * retro_vec[9])

retrodeformed_max_L <- c(rotated_bird[NUM_maxilla_LEFT[1:9], 3] * retro_vec[10], rotated_bird[NUM_maxilla_LEFT[10:18], 3] * retro_vec[11], 
                            rotated_bird[NUM_maxilla_LEFT[19:27], 3] * retro_vec[12], rotated_bird[NUM_maxilla_LEFT[28:36], 3] * retro_vec[13],
                            rotated_bird[NUM_maxilla_LEFT[37:45], 3] * retro_vec[14], rotated_bird[NUM_maxilla_LEFT[46:54], 3] * retro_vec[15],
                            rotated_bird[NUM_maxilla_LEFT[55:63], 3] * retro_vec[16], rotated_bird[NUM_maxilla_LEFT[64:72], 3] * retro_vec[17], 
                            rotated_bird[NUM_maxilla_LEFT[73:81], 3] * retro_vec[18])

retrodeformed_max_R <- c(rotated_bird[NUM_maxilla_RIGHT[1:9], 3] * retro_vec[10], rotated_bird[NUM_maxilla_RIGHT[10:18], 3] * retro_vec[11], 
                            rotated_bird[NUM_maxilla_RIGHT[19:27], 3] * retro_vec[12], rotated_bird[NUM_maxilla_RIGHT[28:36], 3] * retro_vec[13],
                            rotated_bird[NUM_maxilla_RIGHT[37:45], 3] * retro_vec[14], rotated_bird[NUM_maxilla_RIGHT[46:54], 3] * retro_vec[15],
                            rotated_bird[NUM_maxilla_RIGHT[55:63], 3] * retro_vec[16], rotated_bird[NUM_maxilla_RIGHT[64:72], 3] * retro_vec[17], 
                            rotated_bird[NUM_maxilla_RIGHT[73:81], 3] * retro_vec[18])


# retredeformed_nasal_crest <- rotated_bird[1:46, 2] * (1/4)

# Patches back to matrix of all landmarks
retrodeformed_bird <- rotated_bird

retrodeformed_bird[NUM_premaxilla_LEFT, 3] <- retrodeformed_premax_L
retrodeformed_bird[NUM_premaxilla_RIGHT, 3] <- retrodeformed_premax_R
retrodeformed_bird[NUM_maxilla_LEFT, 3] <- retrodeformed_max_L
retrodeformed_bird[NUM_maxilla_RIGHT, 3] <- retrodeformed_max_R

# Save intermediate step landmark data
write.csv(retrodeformed_bird, "data/retrodeformed_bird_6jul2020.csv")


#### 2.4. RETRODEFORMATION plots of retrodeformation on the z axis ####
# Plot landmarks - original (rotated) vs retrodeformed (intermediate step)
open3d(windowRect = c(20, 30, 800, 800))
plot3d(rotated_bird, aspect="iso", type="s", size=0.3, col="blue", add=T)
plot3d(retrodeformed_bird, aspect="iso", type="s", size=0.2, col="green", add=T)


# Save temporary ply meshes
BIRD_rotated <- tps3d(mesh_bird, as.matrix(lm_bird), rotated_bird)
open3d(windowRect = c(20, 30, 800, 800))
shade3d(BIRD_rotated, color="gray", alpha=0.9)
writePLY("data/bird_rotated_lowRes_6july.ply")


### Generate a retrodeformed mesh from the rotated mesh and the landmarks of both
mesh_bird_rot <- file2mesh("data/bird_rotated_lowRes_6july_clean.ply") # Deleted an isolate piece in meshlab
BIRD_retrodeformed <- tps3d(mesh_bird_rot, rotated_bird, retrodeformed_bird)
open3d(windowRect = c(20, 30, 800, 800))
shade3d(BIRD_retrodeformed, color="gray", alpha=0.9)
writePLY("data/bird_retrodeformed_GOOD_6july.ply")

# Plot landmarks - original (rotated) vs retrodeformed - again
open3d(windowRect = c(20, 30, 1200, 900))
shade3d(mesh_bird_rot, color="gray", alpha=0.9)
plot3d(rotated_bird, aspect="iso", type="s", size=0.4, col="green", add=T)
rgl.snapshot("figs/original_LMs_top.png", top = TRUE) 
plot3d(retrodeformed_bird, aspect="iso", type="s", size=0.3, col="blue", add=T)
rgl.snapshot("figs/original_retrodeformed_LMs_top.png", top = TRUE) 

# Plot original (rotated) mesh with both retrodeformed (intermediate step) and original landmarks
open3d(zoom=0.75)
par3d(windowRect= c(0,0,1000,700))
shade3d(mesh_bird_rot, color="gray", alpha=0.9)
plot3d(rotated_bird, aspect="iso", type="s", size=0.28, col="blue", add=T)


#### 2.5. RETRODEFORMATION ON THE xy-PLANE ####
# Retrodeforming the nasal crest

nas_vec <- seq(4.75, 1.05, length.out = 46)
retrodeformed_bird_v3 <- retrodeformed_bird
retrodeformed_bird_v3[1:46, 2] <- retrodeformed_bird[1:46, 2]*nas_vec

# find the points between LM 39 and LM 7

# midpoint between lm7 and lm46
mid_lm23 <- c((retrodeformed_bird[39,1] + retrodeformed_bird_v3[7, 1])/2, (retrodeformed_bird[39,2] + retrodeformed_bird_v3[7, 2])/2,
                                  (retrodeformed_bird[39,3] + retrodeformed_bird_v3[7, 3])/2)

# midpoint lm39 & lm23
mid_lm31 <- c((retrodeformed_bird[39,1] + mid_lm23[1])/2, (retrodeformed_bird[39,2] + mid_lm23[2])/2,
              (retrodeformed_bird[39,3] + mid_lm23[3])/2)

# midpoint lm7 & lm23
mid_lm15 <- c((mid_lm23[1] + retrodeformed_bird_v3[7, 1])/2, (mid_lm23[2] + retrodeformed_bird_v3[7, 2])/2,
              (mid_lm23[3] + retrodeformed_bird_v3[7, 3])/2)

# midpoint lm7 & lm15
mid_lm11 <- c((mid_lm15[1] + retrodeformed_bird_v3[7, 1])/2, (mid_lm15[2] + retrodeformed_bird_v3[7, 2])/2,
              (mid_lm15[3] + retrodeformed_bird_v3[7, 3])/2)

# midpoint lm7 & lm11
mid_lm9 <- c((mid_lm11[1] + retrodeformed_bird_v3[7, 1])/2, (mid_lm11[2] + retrodeformed_bird_v3[7, 2])/2,
              (mid_lm11[3] + retrodeformed_bird_v3[7, 3])/2)

# midpoint lm7 & lm9
mid_lm8 <- c((mid_lm9[1] + retrodeformed_bird_v3[7, 1])/2, (mid_lm9[2] + retrodeformed_bird_v3[7, 2])/2,
             (mid_lm9[3] + retrodeformed_bird_v3[7, 3])/2)

# midpoint lm11 & lm9
mid_lm10 <- c((mid_lm9[1] + mid_lm11[1])/2, (mid_lm9[2] + mid_lm11[2])/2,
             (mid_lm9[3] + mid_lm11[3])/2)

# midpoint lm11 & lm15
mid_lm13 <- c((mid_lm15[1] + mid_lm11[1])/2, (mid_lm15[2] + mid_lm11[2])/2,
              (mid_lm15[3] + mid_lm11[3])/2)

mid_lm14 <- c((mid_lm15[1] + mid_lm13[1])/2, (mid_lm15[2] + mid_lm13[2])/2,
              (mid_lm15[3] + mid_lm13[3])/2)

mid_lm12 <- c((mid_lm11[1] + mid_lm13[1])/2, (mid_lm11[2] + mid_lm13[2])/2,
              (mid_lm11[3] + mid_lm13[3])/2)

# midpoint lm23 & lm15
mid_lm19 <- c((mid_lm15[1] + mid_lm23[1])/2, (mid_lm15[2] + mid_lm23[2])/2,
              (mid_lm15[3] + mid_lm23[3])/2)

mid_lm17 <- c((mid_lm15[1] + mid_lm19[1])/2, (mid_lm15[2] + mid_lm19[2])/2,
              (mid_lm15[3] + mid_lm19[3])/2)

mid_lm16 <- c((mid_lm15[1] + mid_lm17[1])/2, (mid_lm15[2] + mid_lm17[2])/2,
              (mid_lm15[3] + mid_lm17[3])/2)

mid_lm18 <- c((mid_lm17[1] + mid_lm19[1])/2, (mid_lm17[2] + mid_lm19[2])/2,
              (mid_lm17[3] + mid_lm19[3])/2)

mid_lm21 <- c((mid_lm19[1] + mid_lm23[1])/2, (mid_lm19[2] + mid_lm23[2])/2,
              (mid_lm19[3] + mid_lm23[3])/2)

mid_lm20 <- c((mid_lm19[1] + mid_lm21[1])/2, (mid_lm19[2] + mid_lm21[2])/2,
              (mid_lm19[3] + mid_lm21[3])/2)

mid_lm22 <- c((mid_lm21[1] + mid_lm23[1])/2, (mid_lm21[2] + mid_lm23[2])/2,
              (mid_lm21[3] + mid_lm23[3])/2)

# midpoint lm31 & lm23
mid_lm27 <- c((mid_lm31[1] + mid_lm23[1])/2, (mid_lm31[2] + mid_lm23[2])/2,
              (mid_lm31[3] + mid_lm23[3])/2)

mid_lm25 <- c((mid_lm27[1] + mid_lm23[1])/2, (mid_lm27[2] + mid_lm23[2])/2,
              (mid_lm27[3] + mid_lm23[3])/2)

mid_lm24 <- c((mid_lm25[1] + mid_lm23[1])/2, (mid_lm25[2] + mid_lm23[2])/2,
              (mid_lm25[3] + mid_lm23[3])/2)

mid_lm26 <- c((mid_lm27[1] + mid_lm25[1])/2, (mid_lm27[2] + mid_lm25[2])/2,
              (mid_lm27[3] + mid_lm25[3])/2)

mid_lm29 <- c((mid_lm31[1] + mid_lm27[1])/2, (mid_lm31[2] + mid_lm27[2])/2,
              (mid_lm31[3] + mid_lm27[3])/2)

mid_lm30 <- c((mid_lm31[1] + mid_lm29[1])/2, (mid_lm31[2] + mid_lm29[2])/2,
              (mid_lm31[3] + mid_lm29[3])/2)

mid_lm28 <- c((mid_lm27[1] + mid_lm29[1])/2, (mid_lm27[2] + mid_lm29[2])/2,
              (mid_lm27[3] + mid_lm29[3])/2)

# midpoint lm39 & lm31
mid_lm35 <- c((retrodeformed_bird[39,1] + mid_lm31[1])/2, (retrodeformed_bird[39,2] + mid_lm31[2])/2,
              (retrodeformed_bird[39,3] + mid_lm31[3])/2)

mid_lm33 <- c((mid_lm31[1] + mid_lm35[1])/2, (mid_lm31[2] + mid_lm35[2])/2,
              (mid_lm31[3] + mid_lm35[3])/2)

mid_lm34 <- c((mid_lm33[1] + mid_lm35[1])/2, (mid_lm33[2] + mid_lm35[2])/2,
              (mid_lm33[3] + mid_lm35[3])/2)

mid_lm32 <- c((mid_lm31[1] + mid_lm33[1])/2, (mid_lm31[2] + mid_lm33[2])/2,
              (mid_lm31[3] + mid_lm33[3])/2)

mid_lm37 <- c((retrodeformed_bird[39,1] + mid_lm35[1])/2, (retrodeformed_bird[39,2] + mid_lm35[2])/2,
              (retrodeformed_bird[39,3] + mid_lm35[3])/2)

mid_lm36 <- c((mid_lm37[1] + mid_lm35[1])/2, (mid_lm37[2] + mid_lm35[2])/2,
              (mid_lm37[3] + mid_lm35[3])/2)

mid_lm38 <- c((mid_lm37[1] + retrodeformed_bird[39,1])/2, (mid_lm37[2] + retrodeformed_bird[39,2])/2,
              (mid_lm37[3] + retrodeformed_bird[39,3])/2)

### Replace the new lms on the retrodeformed_bird_v3 (intermediate step 3)
retrodeformed_bird_v3[8,] <- mid_lm8
retrodeformed_bird_v3[9,] <- mid_lm9
retrodeformed_bird_v3[10,] <- mid_lm10
retrodeformed_bird_v3[11,] <- mid_lm11
retrodeformed_bird_v3[12,] <- mid_lm12
retrodeformed_bird_v3[13,] <- mid_lm13
retrodeformed_bird_v3[14,] <- mid_lm14
retrodeformed_bird_v3[15,] <- mid_lm15
retrodeformed_bird_v3[16,] <- mid_lm16
retrodeformed_bird_v3[17,] <- mid_lm17
retrodeformed_bird_v3[18,] <- mid_lm18
retrodeformed_bird_v3[19,] <- mid_lm19
retrodeformed_bird_v3[20,] <- mid_lm20
retrodeformed_bird_v3[21,] <- mid_lm21
retrodeformed_bird_v3[22,] <- mid_lm22
retrodeformed_bird_v3[23,] <- mid_lm23
retrodeformed_bird_v3[24,] <- mid_lm24
retrodeformed_bird_v3[25,] <- mid_lm25
retrodeformed_bird_v3[26,] <- mid_lm26
retrodeformed_bird_v3[27,] <- mid_lm27
retrodeformed_bird_v3[28,] <- mid_lm28
retrodeformed_bird_v3[29,] <- mid_lm29
retrodeformed_bird_v3[30,] <- mid_lm30
retrodeformed_bird_v3[31,] <- mid_lm31
retrodeformed_bird_v3[32,] <- mid_lm32
retrodeformed_bird_v3[33,] <- mid_lm33
retrodeformed_bird_v3[34,] <- mid_lm34
retrodeformed_bird_v3[35,] <- mid_lm35
retrodeformed_bird_v3[36,] <- mid_lm36
retrodeformed_bird_v3[37,] <- mid_lm37
retrodeformed_bird_v3[38,] <- mid_lm38

retrodeformed_bird_v3[39:46,] <- retrodeformed_bird[39:46,]

# Smooth out the snout
snout_vec <- seq(1.5, 1, length.out = 20)
retrodeformed_bird_v3[1:20, 2] <- retrodeformed_bird_v3[1:20, 2]/snout_vec

# Save intermediate step 3 landmark data
write.csv(retrodeformed_bird, "data/retrodeformed_bird_v3.csv")

open3d(windowRect = c(20, 30, 1200, 900))
plot3d(rotated_bird, aspect="iso", type="s", size=0.4, col="green", add=T)
plot3d(retrodeformed_bird_v3, aspect="iso", type="s", size=0.3, col="red", add=T)
plot3d(retrodeformed_bird_v3[1:7,], aspect="iso", type="s", size=0.6, col="blue", add=T)
rgl.close()

### Step 4. Retrodeform the most internal semilandmarks on the maxilla

open3d(windowRect = c(20, 30, 1200, 900))
plot3d(rotated_bird, aspect="iso", type="s", size=0.4, col="green", add=T)
plot3d(retrodeformed_bird_v3[c(NUM_maxilla_LEFT, NUM_maxilla_RIGHT), ], aspect="iso", type="s", size=0.3, col="red", add=T)
plot3d(retrodeformed_bird_v3[1:7,], aspect="iso", type="s", size=0.6, col="blue", add=T)

LM_max_L_y <- seq(1:81)[-(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+8)]
LM_max_R_y <- seq(1:81)[-c(1, 10, 19, 28, 37, 46, 55, 64, 73)]

retrodeformed_bird_v4 <- retrodeformed_bird_v3
retrodeformed_bird_v4[c(NUM_maxilla_LEFT[-LM_max_L_y], NUM_maxilla_RIGHT[-LM_max_R_y]), 2] <- retrodeformed_bird_v3[c(NUM_maxilla_LEFT[-LM_max_L_y], NUM_maxilla_RIGHT[-LM_max_R_y]), 2] - 7.5

open3d(windowRect = c(20, 30, 1200, 900))
plot3d(rotated_bird, aspect="iso", type="s", size=0.4, col="green", add=T)
plot3d(retrodeformed_bird_v4, aspect="iso", type="s", size=0.3, col="red", add=T)


# Take out all maxilla, premaxilla (except border) and preforntals in this fourth step of the retrodeformation
premax_L <- seq(1:81)[-c(1, 10, 19, 28, 37, 46, 55, 64, 73)]
premax_R <- seq(1:81)[-(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+8)]

NUM_premaxilla_LEFT
NUM_premaxilla_RIGHT
NUM_maxilla_LEFT
NUM_maxilla_RIGHT
NUM_nasal_LEFT
NUM_nasal_RIGHT

retrodeformed_bird_v3_subset <- retrodeformed_bird_v3[-c(NUM_premaxilla_LEFT[premax_L], NUM_premaxilla_RIGHT[premax_R], NUM_maxilla_LEFT,
                                                         NUM_maxilla_RIGHT, NUM_nasal_LEFT, NUM_nasal_RIGHT), ]

rotated_bird_subset <- rotated_bird[-c(NUM_premaxilla_LEFT[premax_L], NUM_premaxilla_RIGHT[premax_R], NUM_maxilla_LEFT,
                                       NUM_maxilla_RIGHT, NUM_nasal_LEFT, NUM_nasal_RIGHT), ]


open3d(windowRect = c(20, 30, 1200, 900))
BIRD_retrodeformed_v3_subset <- tps3d(mesh_bird_rot, rotated_bird_subset, retrodeformed_bird_v3_subset)
shade3d(BIRD_retrodeformed_v3_subset, color="gray", alpha=0.9)


### Maxilla landmark subset

max_L <- seq(1:81)[-c(1, 9, 10, 18, 19, 27, 28, 36, 37, 45, 46, 54, 55, 63, 64, 72, 73, 81)]
max_R <- seq(1:81)[-c(1, 9, 10, 18, 19, 27, 28, 36, 37, 45, 46, 54, 55, 63, 64, 72, 73, 81)]


retrodeformed_bird_v4_subset <- retrodeformed_bird_v4[-c(NUM_premaxilla_LEFT[premax_L], NUM_premaxilla_RIGHT[premax_R], NUM_maxilla_LEFT[max_L], 
                                                         NUM_maxilla_RIGHT[max_R], NUM_nasal_LEFT, NUM_nasal_RIGHT), ]

rotated_bird_subset <- rotated_bird[-c(NUM_premaxilla_LEFT[premax_L], NUM_premaxilla_RIGHT[premax_R], NUM_maxilla_LEFT[max_L], 
                                       NUM_maxilla_RIGHT[max_R], NUM_nasal_LEFT, NUM_nasal_RIGHT), ]


open3d(windowRect = c(20, 30, 1200, 900))
BIRD_retrodeformed_v4_subset <- tps3d(mesh_bird_rot, rotated_bird_subset, retrodeformed_bird_v4_subset)
shade3d(BIRD_retrodeformed_v4_subset, color="gray", alpha=0.9)


writePLY("data/bird_retrodeformed_FIGURE.ply")


#### 3. Heatmap - original vs retrodeformed mesh differences - Specimen 1 (bird-like) ####

# mesh_bird_retro_v3 <- file2mesh("data/bird_retrodeformed_GOOD_v3.ply")
mesh_bird_retro_v4 <- file2mesh("data/bird_retrodeformed_FIGURE.ply")

open3d(zoom=0.75)
par3d(windowRect= c(0,0,1200,800))
meshDist(mesh_bird_retro_v4, mesh_bird_rot, rampcolors = rev(viridis(n=20)), sign = FALSE)

rgl.snapshot("figs/retrodeformed_bird_FIGURE_top.png", top = TRUE) 
rgl.snapshot("figs/retrodeformed_bird_FIGURE_side_R.png", top = TRUE) 
rgl.snapshot("figs/retrodeformed_bird_FIGURE_side_L.png", top = TRUE) 
rgl.snapshot("figs/retrodeformed_bird_FIGURE_bottom.png", top = TRUE) 
clear3d()
rgl.close()
dev.off() 


#### 4. Import data - specimen 2 (reptile-like) ####
data_skull <- read.csv("data/Landmarks_lizard_FINAL.csv", skip = 1)

lm_reptile <- data_skull[,4:6]
lm_reptile

mesh_reptile <- file2mesh("data/lizard_lowRes.ply")

open3d(windowRect = c(20, 30, 800, 800))
shade3d(mesh_reptile, color="gray", alpha=0.9)

rgl.close()


##### 4.1 Split LM matrix to bones ####

# Curves
Jugal_L <- c(1:25)
Jugal_R <- c(26:65)
Postorbital_L <- c(66:109)
Postorbital_R <- c(110:138)
Vomer_L <- c(139:154)
Vomer_R <- c(155:179)
Vomer_Palatine_L <- c(180:216)
Vomer_Palatine_R <- c(217:244)
Palatine_L1 <- c(245:266)
Palatine_R1 <- c(267:288)
Palatine_L2 <- c(289:307)
Palatine_R2 <- c(308:326)
Orbit_L <- c(327:342)
Orbit_R <- c(343:358)
Premaxilla_L <- c(359:399)
Premaxilla_R <- c(400:440)
Pterygoid_Ecto_Pal_L <- c(441:489)
Pterygoid_Ecto_Pal_R <- c(490:538)
Nasal_crest <- c(539:599)
Frontal_Orbit_L <- c(600:636)
Frontal_Orbit_R <- c(637:673)
Post_Parietal_L <- c(674:695)
Post_Parietal_R <- c(696:717)
Basipterygoid_Basisphenoid_R <- c(718:736)
Basipterygoid_Basisphenoid_L <- c(737:755)
Pterygoid_R <- c(756:774)
Pterygoid_L <- c(775:793)
Postfrontal_L <- c(794:809)
Postfrontal_R<- c(810:829)
Basisphenoid_R <- c(826:832)
Basisphenoid_L <- c(833:839)

# Patches
PATCH_Premaxilla_L <- c(840:920)
PATCH_Premaxilla_R <- c(921:1001)
PATCH_Nasal_L <- c(1002:1082)
PATCH_Nasal_R <- c(1083:1163)
PATCH_Postfrontal_L <- c(1164:1244)
PATCH_Postfrontal_R <- c(1245:1325)
PATCH_Maxilla_L <- c(1326:1406)
PATCH_Maxilla_R <- c(1407:1487)
PATCH_Quadrate_L <- c(1488:1568)
PATCH_Quadrate_R <- c(1731:1811)
PATCH_Supraoccipital_L <- c(1569:1649)
PATCH_Supraoccipital_R <- c(1650:1730)
PATCH_Tabular_L <- c(1812:1892)
PATCH_Tabular_R <- c(1893:1973)


##### 4.2 Plot LMs on mesh ####
open3d(windowRect = c(20, 30, 800, 800))
shade3d(mesh_reptile, color="gray", alpha=0.9)
plot3d(lm_reptile[PATCH_Premaxilla_L, ], aspect="iso", type="s", size=0.6, col="green", add=T)
plot3d(lm_reptile[PATCH_Premaxilla_R, ], aspect="iso", type="s", size=0.6, col="green", add=T)
plot3d(lm_reptile[PATCH_Nasal_L, ], aspect="iso", type="s", size=0.6, col="blue", add=T)
plot3d(lm_reptile[PATCH_Nasal_R, ], aspect="iso", type="s", size=0.6, col="blue", add=T)
plot3d(lm_reptile[PATCH_Maxilla_L, ], aspect="iso", type="s", size=0.6, col="red", add=T)
plot3d(lm_reptile[PATCH_Maxilla_R, ], aspect="iso", type="s", size=0.6, col="red", add=T)
plot3d(lm_reptile[PATCH_Postfrontal_L, ], aspect="iso", type="s", size=0.6, col="yellow", add=T)
plot3d(lm_reptile[PATCH_Postfrontal_R, ], aspect="iso", type="s", size=0.6, col="yellow", add=T)
plot3d(lm_reptile[Orbit_L, ], aspect="iso", type="s", size=0.4, col="black", add=T) 
plot3d(lm_reptile[Orbit_R, ], aspect="iso", type="s", size=0.4, col="black", add=T) 
plot3d(lm_reptile[Frontal_Orbit_L, ], aspect="iso", type="s", size=0.4, col="darkblue", add=T)
plot3d(lm_reptile[Frontal_Orbit_R, ], aspect="iso", type="s", size=0.4, col="darkblue", add=T)

open3d(windowRect = c(20, 30, 800, 800))
shade3d(mesh_reptile, color="gray", alpha=0.9)
plot3d(lm_reptile[Nasal_crest, ], aspect="iso", type="s", size=0.4, col="darkblue", add=T)
plot3d(lm_reptile[Palatine_L2, ], aspect="iso", type="s", size=0.8, col="darkblue", add=T)


open3d(windowRect = c(20, 30, 800, 800))
shade3d(mesh_reptile, color="gray", alpha=0.9)
plot3d(lm_reptile[PATCH_Quadrate_L, ], aspect="iso", type="s", size=0.4, col="darkblue", add=T)
plot3d(lm_reptile[PATCH_Quadrate_R, ], aspect="iso", type="s", size=0.8, col="darkblue", add=T)
plot3d(lm_reptile[PATCH_Supraoccipital_L, ], aspect="iso", type="s", size=0.6, col="green", add=T)
plot3d(lm_reptile[PATCH_Supraoccipital_R, ], aspect="iso", type="s", size=0.6, col="green", add=T)
plot3d(lm_reptile[PATCH_Tabular_L, ], aspect="iso", type="s", size=0.6, col="blue", add=T)
plot3d(lm_reptile[PATCH_Tabular_R, ], aspect="iso", type="s", size=0.6, col="blue", add=T)


##### 4.3 Plot retrodeformation subset LMs on mesh ####

lm_reptile_subset <- lm_reptile[c(Nasal_crest, PATCH_Nasal_L, PATCH_Nasal_R, Orbit_L, Orbit_R, 
   Frontal_Orbit_L, Frontal_Orbit_R, Postfrontal_L, Postfrontal_R, 
   PATCH_Supraoccipital_L, PATCH_Supraoccipital_R, Post_Parietal_L, Post_Parietal_R, Palatine_L2, Palatine_R2,
   PATCH_Premaxilla_L, PATCH_Premaxilla_R, PATCH_Maxilla_L, PATCH_Maxilla_R, PATCH_Nasal_L, PATCH_Nasal_R), ]

open3d(windowRect = c(20, 30, 1800, 800))
shade3d(mesh_reptile, color="gray", alpha=0.9)
plot3d(lm_reptile_subset, aspect="iso", type="s", size=0.4, col="green", add=T)
rgl.snapshot("figs/reptile_original_LMs_RIGHT.png", top = TRUE) 
rgl.snapshot("figs/reptile_original_LMs_LEFT.png", top = TRUE) 
rgl.snapshot("figs/reptile_original_LMs_TOP.png", top = TRUE) 
rgl.close()

#### 5. Rotate the reptile - ShapeRotator ####
#### 5.1. Landmark selection, angle, translation to c(0, 0, 0) and rotation to the xz-plane ####
# landmarks

land.a = 539 # distal LM on the nasal ridge
land.b = 599 # last point of the nasal ridge, between the eyes
land.c = 289 # first landmark on the palatine L2

# angle
angle = 0

# translate
lm_reptile_a <- array(data = cbind(as.matrix(lm_reptile), as.matrix(lm_reptile)), dim = c(dim(lm_reptile)[1], dim(lm_reptile)[2], 2))

dimnames(lm_reptile_a)[[3]] <- c("reptile1", "reptile2") # duplication of the data as ShapeRotator works with arrays only
str(lm_reptile_a)
class(lm_reptile_a)

reptile_t <- translate(lm_reptile_a, land.a)

reptile_t[1,,1]

# dataset
data.1 = reptile_t

# Rotation
rot.array.1 <- array(data = NA, dim = c(dim(data.1)[1], dim(data.1)[2], dim(data.1)[3]),dimnames = list(NULL, NULL, c(dimnames(data.1)[[3]])));

# This is the angle which land.b will face to land.e
angle.v <-vector.angle(angle)


#Rotating for data.1
for( i in 1:(dim(data.1)[3]) ){
  # Untranslated specimen?
  if ( ! isTRUE(  all.equal((data.1[,,i])[land.a,] , c(0,0,0)) ) ) {
    # warning (sprintf("Landmark A is not at the origin (data.1[,,%d])[%d,] = (%f,%f,%f)?", i, land.a, (data.1[,,i])[land.a,1], (data.1[,,i])[land.a,2], (data.1[,,i])[land.a,3]), immediate. = immediate_on)
  }
  
  # First, we rotate so that land.b rests on the axis (0,1,0).
  rot.array.1[,,i] <- rotate.orientation(data.1[,,i], land.b, c(0,1,0))
  
  #Sanity check: we need to have rot.array.1[,,i][land.b,] = (0, y, 0)
  if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.b,1] , 0)) || ! isTRUE(  all.equal( (rot.array.1[,,i])[land.b,3] , 0)) ) {
    warning (sprintf("First rotation gone wrong in (rot.array.1[,,%d])[%d,] = (%f,%f,%f)?", i, land.b, (rot.array.1[,,i])[land.b,1], (rot.array.1[,,i])[land.b,2], (rot.array.1[,,i])[land.b,3]), immediate. = immediate_on)
  }
  
  # Rotate landmark C to the X-Y plane - we do so by computing the cross product between the unit.landmark.c.pos.proj vector and (-1,0,0),
  # which in the X-Y plane forms an angle given by the angle between these two vectors:
  unit.landmark.c.pos.proj = unit_3D( c( (rot.array.1[,,i])[land.c,1], 0, (rot.array.1[,,i])[land.c,3]))	# Unit vector with landmark.c.pos projected to the X-Z plane.
  rotmat <- rotmat_3D( cross_3D(unit.landmark.c.pos.proj, c(-1,0,0)), angle_3D(unit.landmark.c.pos.proj, c(-1,0,0)) )
  rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
  
  
  # Let us now check that the landmark C is in the right spot. That is, we want C_x < B_x, if not, we rotate by angle pi.
  if ( isTRUE((rot.array.1[,,i])[land.b,1] < (rot.array.1[,,i])[land.c,1]) ) {
    rotmat <- rotmat_3D( c(0,1,0), pi)
    rot.array.1[,,i] <- rotveclist_3D ( rotmat, rot.array.1[,,i])
  }
  
  # Sanity check:
  if ( ! isTRUE(  all.equal( (rot.array.1[,,i])[land.c,3] , 0)) || isTRUE((rot.array.1[,,i])[land.b,1] < (rot.array.1[,,i])[land.c,1])) {
    warning (sprintf("Second rotation gone wrong in (rot.array.1[,,%d])[%d,] = (%f,%f,%f)?", i, land.c, (rot.array.1[,,i])[land.c,1], (rot.array.1[,,i])[land.c,2], (rot.array.1[,,i])[land.c,3]), immediate. = immediate_on)
  }
  
  # Now, we rigidly move this object to sit at an angle given by "angle" to the X-axis.
  rot.array.1[,,i] <- rotate.orientation(rot.array.1[,,i], land.b, angle.v)
}

rotated_reptile <- rot.array.1[,,1]

str(rotated_reptile)


rotated_reptile[539,] # land.a = 539 # distal LM on the nasal ridge
rotated_reptile[599,] # land.b = 599 # last point of the nasal ridge, between the eyes
rotated_reptile[289,] # land.c = 289 
# please note that the ventral side is on the top

write.csv(rotated_reptile, "data/rotated_reptile.csv")

rotated_reptile <- read.csv("data/rotated_reptile.csv")[,2:4]
# Save temporary ply meshes
mesh_reptile_rotated <- tps3d(mesh_reptile, as.matrix(lm_reptile), rotated_reptile)
open3d(windowRect = c(20, 30, 800, 800))
shade3d(mesh_reptile_rotated, color="gray", alpha=0.9)
writePLY("data/reptile_rotated_lowRes_10july.ply")

mesh_reptile_rotated <- file2mesh("data/reptile_rotated_lowRes_10july.ply")

#### 6. RETRODEFORMATION - Intermediate step 1 ####
# Plot rotated reptile
NUM_Nasal_LEFT <- c(1002:1082)
NUM_Nasal_RIGHT <- c(1083:1163)
NUM_Nasal_Crest <- c(539:599)
NUM_Orbit_LEFT <- c(327:342)
NUM_Orbit_RIGHT <- c(343:358)
NUM_Frontal_Orbit_LEFT <- c(600:636)
NUM_Frontal_Orbit_RIGHT <- c(637:673)
NUM_Postfrontal_LEFT <- c(1164:1244)
NUM_Postfrontal_RIGHT <- c(1245:1325)
NUM_Supraoccipital_LEFT <- c(1569:1649)
NUM_Supraoccipital_RIGHT <- c(1650:1730)

open3d(windowRect = c(20, 30, 800, 800))
shade3d(mesh_reptile_rotated, color="gray", alpha=0.9)
plot3d(rotated_reptile[NUM_Nasal_Crest, ], aspect="iso", type="s", size=0.3, col="blue", add=T)
plot3d(rotated_reptile[NUM_Nasal_LEFT, ], aspect="iso", type="s", size=0.6, col="yellow", add=T)
plot3d(rotated_reptile[NUM_Nasal_RIGHT, ], aspect="iso", type="s", size=0.6, col="yellow", add=T)
plot3d(rotated_reptile[NUM_Orbit_LEFT, ], aspect="iso", type="s", size=0.4, col="black", add=T) 
plot3d(rotated_reptile[NUM_Orbit_RIGHT, ], aspect="iso", type="s", size=0.4, col="black", add=T) 
plot3d(rotated_reptile[NUM_Frontal_Orbit_LEFT, ], aspect="iso", type="s", size=0.4, col="darkblue", add=T)
plot3d(rotated_reptile[NUM_Frontal_Orbit_RIGHT, ], aspect="iso", type="s", size=0.4, col="darkblue", add=T)
plot3d(rotated_reptile[NUM_Postfrontal_LEFT, ], aspect="iso", type="s", size=0.4, col="red", add=T)
plot3d(rotated_reptile[NUM_Postfrontal_RIGHT, ], aspect="iso", type="s", size=0.4, col="red", add=T)
plot3d(rotated_reptile[Post_Parietal_L, ], aspect="iso", type="s", size=0.4, col="pink", add=T)
plot3d(rotated_reptile[Post_Parietal_R, ], aspect="iso", type="s", size=0.4, col="pink", add=T)

#### 6.1. Retrodeformation Orbits ####
# They skull was flattened dorsi-ventrally and that is quite conspicuous on the orbits
# Curves to retrodeform on the y-axis
# Orbit_L & Orbit_R (plotted in black) - incremental vector
retro_reptile <- rotated_reptile
retro_vec_orbit <- seq(0.1, 1.05, length.out = length(NUM_Orbit_LEFT))

retro_reptile[NUM_Orbit_LEFT, 2] <- retro_reptile[NUM_Orbit_LEFT, 2] - retro_vec_orbit
retro_reptile[NUM_Orbit_RIGHT, 2] <- retro_reptile[NUM_Orbit_RIGHT, 2] - retro_vec_orbit                                                                      

# Nasal_crest - incremental vector half way to the closest point to the vault
retro_vec_nas <- seq(0, 1.45, length.out = (length(NUM_Nasal_Crest) - ceiling(length(NUM_Nasal_Crest)/2) + 1)) 
retro_reptile[NUM_Nasal_Crest[ceiling(length(NUM_Nasal_Crest)/2)]:NUM_Nasal_Crest[length(NUM_Nasal_Crest)], 
              2] <- retro_reptile[NUM_Nasal_Crest[ceiling(length(NUM_Nasal_Crest)/2)]:NUM_Nasal_Crest[length(NUM_Nasal_Crest)], 2] - retro_vec_nas


# Frontal_Orbit_L & Frontal_Orbit_R (plotted in dark blue) - bring it up but the most basal points stay the same
length(NUM_Frontal_Orbit_LEFT)
length(NUM_Frontal_Orbit_RIGHT)

retro_vec_front_orbit <- seq(1.5, 0.75, length.out = length(NUM_Frontal_Orbit_LEFT))

retro_reptile[NUM_Frontal_Orbit_LEFT, 2] <- retro_reptile[NUM_Frontal_Orbit_LEFT, 2] - retro_vec_front_orbit 
retro_reptile[NUM_Frontal_Orbit_RIGHT, 2] <- retro_reptile[NUM_Frontal_Orbit_RIGHT, 2] - retro_vec_front_orbit 

#### 6.2. Retrodeformation Postfrontals ####
# Patches to deform on the y-axis - bring up the whole thing
# PATCH_Postfrontal_L & PATCH_Postfrontal_R (plotted in yellow) - the part closer to the snout 
retro_vec_postfront <- seq(1.15, 0, length.out = 9)

retro_reptile[NUM_Postfrontal_LEFT[1:9], 2] <- retro_reptile[NUM_Postfrontal_LEFT[1:9], 2] - rep(retro_vec_postfront[1], 9)
retro_reptile[NUM_Postfrontal_LEFT[10:18], 2] <- retro_reptile[NUM_Postfrontal_LEFT[10:18], 2] - rep(retro_vec_postfront[2], 9)
retro_reptile[NUM_Postfrontal_LEFT[19:27], 2] <- retro_reptile[NUM_Postfrontal_LEFT[19:27], 2] - rep(retro_vec_postfront[3], 9)
retro_reptile[NUM_Postfrontal_LEFT[28:36], 2] <- retro_reptile[NUM_Postfrontal_LEFT[28:36], 2] - rep(retro_vec_postfront[4], 9)
retro_reptile[NUM_Postfrontal_LEFT[37:45], 2] <- retro_reptile[NUM_Postfrontal_LEFT[37:45], 2] - rep(retro_vec_postfront[5], 9)
retro_reptile[NUM_Postfrontal_LEFT[46:54], 2] <- retro_reptile[NUM_Postfrontal_LEFT[46:54], 2] - rep(retro_vec_postfront[6], 9)
retro_reptile[NUM_Postfrontal_LEFT[55:63], 2] <- retro_reptile[NUM_Postfrontal_LEFT[55:63], 2] - rep(retro_vec_postfront[7], 9)
retro_reptile[NUM_Postfrontal_LEFT[64:72], 2] <- retro_reptile[NUM_Postfrontal_LEFT[64:72], 2] - rep(retro_vec_postfront[8], 9)
retro_reptile[NUM_Postfrontal_LEFT[73:81], 2] <- retro_reptile[NUM_Postfrontal_LEFT[73:81], 2] - rep(retro_vec_postfront[9], 9)

retro_reptile[NUM_Postfrontal_RIGHT[1:9], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[1:9], 2] - rep(retro_vec_postfront[1], 9)
retro_reptile[NUM_Postfrontal_RIGHT[10:18], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[10:18], 2] - rep(retro_vec_postfront[2], 9)
retro_reptile[NUM_Postfrontal_RIGHT[19:27], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[19:27], 2] - rep(retro_vec_postfront[3], 9)
retro_reptile[NUM_Postfrontal_RIGHT[28:36], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[28:36], 2] - rep(retro_vec_postfront[4], 9)
retro_reptile[NUM_Postfrontal_RIGHT[37:45], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[37:45], 2] - rep(retro_vec_postfront[5], 9)
retro_reptile[NUM_Postfrontal_RIGHT[46:54], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[46:54], 2] - rep(retro_vec_postfront[6], 9)
retro_reptile[NUM_Postfrontal_RIGHT[55:63], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[55:63], 2] - rep(retro_vec_postfront[7], 9)
retro_reptile[NUM_Postfrontal_RIGHT[64:72], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[64:72], 2] - rep(retro_vec_postfront[8], 9)
retro_reptile[NUM_Postfrontal_RIGHT[73:81], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[73:81], 2] - rep(retro_vec_postfront[9], 9)

### Now bring only the center of the postfrontals up
# Retrodeformation factor for y
retro_vec_postfront2_left <- seq(0, 0.45, length.out = 9) # differential deformation as the top of the postfrontals has been more flattened than the borders
retro_vec_postfront2_right <- seq(0.45, 0, length.out = 9)

### Left
retro_reptile[NUM_Postfrontal_LEFT[1:9], 2] <- retro_reptile[NUM_Postfrontal_LEFT[1:9], 2] - retro_vec_postfront2_left
retro_reptile[NUM_Postfrontal_LEFT[10:18], 2] <- retro_reptile[NUM_Postfrontal_LEFT[10:18], 2] - retro_vec_postfront2_left
retro_reptile[NUM_Postfrontal_LEFT[19:27], 2] <- retro_reptile[NUM_Postfrontal_LEFT[19:27], 2] - retro_vec_postfront2_left
retro_reptile[NUM_Postfrontal_LEFT[28:36], 2] <- retro_reptile[NUM_Postfrontal_LEFT[28:36], 2] - retro_vec_postfront2_left
retro_reptile[NUM_Postfrontal_LEFT[37:45], 2] <- retro_reptile[NUM_Postfrontal_LEFT[37:45], 2] - retro_vec_postfront2_left
retro_reptile[NUM_Postfrontal_LEFT[46:54], 2] <- retro_reptile[NUM_Postfrontal_LEFT[46:54], 2] - retro_vec_postfront2_left
retro_reptile[NUM_Postfrontal_LEFT[55:63], 2] <- retro_reptile[NUM_Postfrontal_LEFT[55:63], 2] - retro_vec_postfront2_left
retro_reptile[NUM_Postfrontal_LEFT[64:72], 2] <- retro_reptile[NUM_Postfrontal_LEFT[64:72], 2] - retro_vec_postfront2_left
retro_reptile[NUM_Postfrontal_LEFT[73:81], 2] <- retro_reptile[NUM_Postfrontal_LEFT[73:81], 2] - retro_vec_postfront2_left

### Right
retro_reptile[NUM_Postfrontal_RIGHT[1:9], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[1:9], 2] - retro_vec_postfront2_right
retro_reptile[NUM_Postfrontal_RIGHT[10:18], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[10:18], 2] - retro_vec_postfront2_right
retro_reptile[NUM_Postfrontal_RIGHT[19:27], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[19:27], 2] - retro_vec_postfront2_right
retro_reptile[NUM_Postfrontal_RIGHT[28:36], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[28:36], 2] - retro_vec_postfront2_right
retro_reptile[NUM_Postfrontal_RIGHT[37:45], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[37:45], 2] - retro_vec_postfront2_right
retro_reptile[NUM_Postfrontal_RIGHT[46:54], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[46:54], 2] - retro_vec_postfront2_right
retro_reptile[NUM_Postfrontal_RIGHT[55:63], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[55:63], 2] - retro_vec_postfront2_right
retro_reptile[NUM_Postfrontal_RIGHT[64:72], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[64:72], 2] - retro_vec_postfront2_right
retro_reptile[NUM_Postfrontal_RIGHT[73:81], 2] <- retro_reptile[NUM_Postfrontal_RIGHT[73:81], 2] - retro_vec_postfront2_right


#### 6.3. Retrodeformation Nasals ####
# PATCH_Nasal_L & PATCH_Nasal_R (plotted in blue) - bring the whole thing up differentially. The most frontal part does not come up that much.
# It has been more deformed on the first region
retro_vec_nasal_patch <- seq(0, 1.35, length.out = 9)

retro_reptile[NUM_Nasal_LEFT[1:9], 2] <- retro_reptile[NUM_Nasal_LEFT[1:9], 2] - rep(retro_vec_nasal_patch[1], 9)
retro_reptile[NUM_Nasal_LEFT[10:18], 2] <- retro_reptile[NUM_Nasal_LEFT[10:18], 2] - rep(retro_vec_nasal_patch[2], 9)
retro_reptile[NUM_Nasal_LEFT[19:27], 2] <- retro_reptile[NUM_Nasal_LEFT[19:27], 2] - rep(retro_vec_nasal_patch[3], 9)
retro_reptile[NUM_Nasal_LEFT[28:36], 2] <- retro_reptile[NUM_Nasal_LEFT[28:36], 2] - rep(retro_vec_nasal_patch[4], 9)
retro_reptile[NUM_Nasal_LEFT[37:45], 2] <- retro_reptile[NUM_Nasal_LEFT[37:45], 2] - rep(retro_vec_nasal_patch[5], 9)
retro_reptile[NUM_Nasal_LEFT[46:54], 2] <- retro_reptile[NUM_Nasal_LEFT[46:54], 2] - rep(retro_vec_nasal_patch[6], 9)
retro_reptile[NUM_Nasal_LEFT[55:63], 2] <- retro_reptile[NUM_Nasal_LEFT[55:63], 2] - rep(retro_vec_nasal_patch[7], 9)
retro_reptile[NUM_Nasal_LEFT[64:72], 2] <- retro_reptile[NUM_Nasal_LEFT[64:72], 2] - rep(retro_vec_nasal_patch[8], 9)
retro_reptile[NUM_Nasal_LEFT[73:81], 2] <- retro_reptile[NUM_Nasal_LEFT[73:81], 2] - rep(retro_vec_nasal_patch[9], 9)

retro_reptile[NUM_Nasal_RIGHT[1:9], 2] <- retro_reptile[NUM_Nasal_RIGHT[1:9], 2] - rep(retro_vec_nasal_patch[1], 9)
retro_reptile[NUM_Nasal_RIGHT[10:18], 2] <- retro_reptile[NUM_Nasal_RIGHT[10:18], 2] - rep(retro_vec_nasal_patch[2], 9)
retro_reptile[NUM_Nasal_RIGHT[19:27], 2] <- retro_reptile[NUM_Nasal_RIGHT[19:27], 2] - rep(retro_vec_nasal_patch[3], 9)
retro_reptile[NUM_Nasal_RIGHT[28:36], 2] <- retro_reptile[NUM_Nasal_RIGHT[28:36], 2] - rep(retro_vec_nasal_patch[4], 9)
retro_reptile[NUM_Nasal_RIGHT[37:45], 2] <- retro_reptile[NUM_Nasal_RIGHT[37:45], 2] - rep(retro_vec_nasal_patch[5], 9)
retro_reptile[NUM_Nasal_RIGHT[46:54], 2] <- retro_reptile[NUM_Nasal_RIGHT[46:54], 2] - rep(retro_vec_nasal_patch[6], 9)
retro_reptile[NUM_Nasal_RIGHT[55:63], 2] <- retro_reptile[NUM_Nasal_RIGHT[55:63], 2] - rep(retro_vec_nasal_patch[7], 9)
retro_reptile[NUM_Nasal_RIGHT[64:72], 2] <- retro_reptile[NUM_Nasal_RIGHT[64:72], 2] - rep(retro_vec_nasal_patch[8], 9)
retro_reptile[NUM_Nasal_RIGHT[73:81], 2] <- retro_reptile[NUM_Nasal_RIGHT[73:81], 2] - rep(retro_vec_nasal_patch[9], 9)

#### 6.4. Retrodeformation SUPRAOCCIPITAL ####
NUM_Supraoccipital_LEFT <- c(1569:1649)
NUM_Supraoccipital_RIGHT <- c(1650:1730)

retro_vec_postpar <- seq(1.05, 0, length.out = 9)

# Left
retro_reptile[NUM_Supraoccipital_LEFT[1:9], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[1:9], 2] - rep(retro_vec_postpar[1], 9)
retro_reptile[NUM_Supraoccipital_LEFT[10:18], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[10:18], 2] - rep(retro_vec_postpar[2], 9)
retro_reptile[NUM_Supraoccipital_LEFT[19:27], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[19:27], 2] - rep(retro_vec_postpar[3], 9)
retro_reptile[NUM_Supraoccipital_LEFT[28:36], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[28:36], 2] - rep(retro_vec_postpar[4], 9)
retro_reptile[NUM_Supraoccipital_LEFT[37:45], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[37:45], 2] - rep(retro_vec_postpar[5], 9)
retro_reptile[NUM_Supraoccipital_LEFT[46:54], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[46:54], 2] - rep(retro_vec_postpar[6], 9)
retro_reptile[NUM_Supraoccipital_LEFT[55:63], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[55:63], 2] - rep(retro_vec_postpar[7], 9)
retro_reptile[NUM_Supraoccipital_LEFT[64:72], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[64:72], 2] - rep(retro_vec_postpar[8], 9)
retro_reptile[NUM_Supraoccipital_LEFT[73:81], 2] <- retro_reptile[NUM_Supraoccipital_LEFT[73:81], 2] - rep(retro_vec_postpar[9], 9)

# Right
retro_reptile[NUM_Supraoccipital_RIGHT[1:9], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[1:9], 2] - rep(retro_vec_postpar[1], 9)
retro_reptile[NUM_Supraoccipital_RIGHT[10:18], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[10:18], 2] - rep(retro_vec_postpar[2], 9)
retro_reptile[NUM_Supraoccipital_RIGHT[19:27], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[19:27], 2] - rep(retro_vec_postpar[3], 9)
retro_reptile[NUM_Supraoccipital_RIGHT[28:36], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[28:36], 2] - rep(retro_vec_postpar[4], 9)
retro_reptile[NUM_Supraoccipital_RIGHT[37:45], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[37:45], 2] - rep(retro_vec_postpar[5], 9)
retro_reptile[NUM_Supraoccipital_RIGHT[46:54], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[46:54], 2] - rep(retro_vec_postpar[6], 9)
retro_reptile[NUM_Supraoccipital_RIGHT[55:63], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[55:63], 2] - rep(retro_vec_postpar[7], 9)
retro_reptile[NUM_Supraoccipital_RIGHT[64:72], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[64:72], 2] - rep(retro_vec_postpar[8], 9)
retro_reptile[NUM_Supraoccipital_RIGHT[73:81], 2] <- retro_reptile[NUM_Supraoccipital_RIGHT[73:81], 2] - rep(retro_vec_postpar[9], 9)

### Bring the top of the supraoccipitals closer to the parietals (so that they can be tucked in) - x axis (1.5% closer)
retro_vec_postpar2 <- seq(0.98, 1, length.out = 9)

# Left
retro_reptile[NUM_Supraoccipital_LEFT[c(1, 10, 19, 28, 37, 46, 55, 64, 73)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[c(1, 10, 19, 28, 37, 46, 55, 64, 73)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+1)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+1)], 1]* retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+2)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+2)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+3)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+3)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+4)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+4)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+5)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+5)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+6)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+6)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+7)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+7)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+8)], 1] <- retro_reptile[NUM_Supraoccipital_LEFT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+8)], 1] * retro_vec_postpar2

# Right
retro_reptile[NUM_Supraoccipital_RIGHT[c(1, 10, 19, 28, 37, 46, 55, 64, 73)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[c(1, 10, 19, 28, 37, 46, 55, 64, 73)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+1)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+1)], 1]* retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+2)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+2)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+3)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+3)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+4)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+4)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+5)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+5)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+6)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+6)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+7)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+7)], 1] * retro_vec_postpar2
retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+8)], 1] <- retro_reptile[NUM_Supraoccipital_RIGHT[(c(1, 10, 19, 28, 37, 46, 55, 64, 73)+8)], 1] * retro_vec_postpar2

#### 6.5. POST PARIETALS retrodeformed in an angle ####
# x-axis 
retro_vec_post_par <- seq(0.75, 0.25, length.out = 22)
retro_reptile[Post_Parietal_L, 1] <- retro_reptile[Post_Parietal_L, 1] - retro_vec_post_par
retro_reptile[Post_Parietal_R, 1] <- retro_reptile[Post_Parietal_R, 1] - retro_vec_post_par

# y-axis
retro_reptile[Post_Parietal_L, 2] <- retro_reptile[Post_Parietal_L, 2] - 0.5
retro_reptile[Post_Parietal_R, 2] <- retro_reptile[Post_Parietal_R, 2] - 0.5


#### 6.6. PLOT RETRODEFORMED REPTILE ####

open3d(windowRect = c(20, 30, 800, 800))
shade3d(mesh_reptile_rotated, color="gray", alpha=0.9)
plot3d(retro_reptile[NUM_Nasal_Crest, ], aspect="iso", type="s", size=0.3, col="blue", add=T)
plot3d(retro_reptile[NUM_Nasal_LEFT, ], aspect="iso", type="s", size=0.6, col="yellow", add=T)
plot3d(retro_reptile[NUM_Nasal_RIGHT, ], aspect="iso", type="s", size=0.6, col="yellow", add=T)
plot3d(retro_reptile[NUM_Orbit_LEFT, ], aspect="iso", type="s", size=0.4, col="black", add=T) 
plot3d(retro_reptile[NUM_Orbit_RIGHT, ], aspect="iso", type="s", size=0.4, col="black", add=T) 
plot3d(retro_reptile[NUM_Frontal_Orbit_LEFT, ], aspect="iso", type="s", size=0.4, col="darkblue", add=T)
plot3d(retro_reptile[NUM_Frontal_Orbit_RIGHT, ], aspect="iso", type="s", size=0.4, col="darkblue", add=T)
plot3d(retro_reptile[NUM_Postfrontal_LEFT, ], aspect="iso", type="s", size=0.4, col="red", add=T)
plot3d(retro_reptile[NUM_Postfrontal_RIGHT, ], aspect="iso", type="s", size=0.4, col="red", add=T)
plot3d(retro_reptile[NUM_Supraoccipital_LEFT, ], aspect="iso", type="s", size=0.4, col="pink", add=T)
plot3d(retro_reptile[NUM_Supraoccipital_RIGHT, ], aspect="iso", type="s", size=0.4, col="pink", add=T)
plot3d(retro_reptile[Post_Parietal_L, ], aspect="iso", type="s", size=0.6, col="black", add=T)
plot3d(retro_reptile[Post_Parietal_R, ], aspect="iso", type="s", size=0.6, col="black", add=T)
plot3d(retro_reptile[Palatine_L2, ], aspect="iso", type="s", size=0.4, col="green", add=T) #Palatine_L2 <- c(289:307)
plot3d(retro_reptile[Palatine_R2, ], aspect="iso", type="s", size=0.4, col="green", add=T) #Palatine_R2 <- c(308:326)


open3d(windowRect = c(20, 30, 1200, 900))
REPTILE_retrodeformed_v3 <- tps3d(mesh_reptile_rotated, as.matrix(rotated_reptile[-NUM_Frontal_Orbit_LEFT, ]), as.matrix(retro_reptile[-NUM_Frontal_Orbit_LEFT, ]))
shade3d(REPTILE_retrodeformed_v3, color="gray", alpha=0.9)
# plot3d(retrodeformed_bird_v3_subset, aspect="iso", type="s", size=0.3, col="blue", add=T)

writePLY("data/reptile_retrodeformed_FIGURE.ply")
REPTILE_FIGURE_clean <- file2mesh("data/reptile_retrodeformed_FIGURE_clean.ply") # took an isolated group of faces out in meshlab
write.csv(retro_reptile, "data/retrodeformed_reptile_v3.csv", row.names = F)


open3d(zoom=0.75)
par3d(windowRect= c(0,0,1200,800))
meshDist(REPTILE_retrodeformed_FIGURE_clean, mesh_reptile_rotated, rampcolors = rev(viridis(n=20)), sign = FALSE)

### SAVE SNAPSHOTS
rgl.snapshot("figs/retrodeformed_reptile_FIGURE_top.png", top = TRUE) 
rgl.snapshot("figs/retrodeformed_reptile_FIGURE_side_R.png", top = TRUE) 
rgl.snapshot("figs/retrodeformed_reptile_FIGURE_side_L.png", top = TRUE) 
rgl.snapshot("figs/retrodeformed_reptile_FIGURE_bottom.png", top = TRUE) 
clear3d()
rgl.close()
dev.off() 
