#### 0. ShapeRotator ####
# We have included here the internal functions of the R package ShapeRotator (Vidal-Garcia et al, 2018) for review purposes
# Obtain the length of a vector u
norm_3D <- function (u)
{
  return( sqrt(u[1]^2 + u[2]^2 + u[3]^2));}

# Return the unit renormalisation of the vector from v
unit_3D <- function (v)
{
  return ( v/norm_3D(v) );
}

# Compute the dot product between vectors u and v
dot_3D <- function (u,v){
  d = u[1] * v[1] + u[2] * v[2] + u[3] * v[3];
  return (d);
}

# Take the cross product between the vectors u and v
cross_3D <- function (u,v)
{
  w = c( u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1] );
  return (w);
}

# Compute the angle between two vectors
angle_3D <- function (u,v)
{
  a = acos( dot_3D(u,v)/ (norm_3D(u) * norm_3D(v)) );
  return (a);
}

# Return M * v, the action of the matrix M on a vector v.
mmult_3D <- function (M, v)
{
  w = c( M[1,1] * v[1] + M [1, 2] * v[2] + M [1,3] * v[3],
         M [2,1] * v[1] + M [2,2] * v[2] + M[2,3] * v[3],
         M [3,1] * v[1] + M [3,2] * v[2] + M[3,3] * v[3]);
  return (w);
}

# Return a rotation matrix given an axis v and an angle t:
rotmat_3D <- function (v, t)
{
  # Create a 3 x 3 matrix with 0's.
  R = matrix (0, 3, 3);
  # normalise v to u, because the rotation matrix needs a unit vector
  u = unit_3D (v);
  # The first entry of the rotation matrix is: 
  R[1,] = c (cos(t) + u[1]^2 *(1 - cos(t)), u[1]* u[2] * (1 - cos(t)) - u[3] * sin(t), u[1]* u[3]* (1 - cos(t)) + u[2] *sin(t) );
  R[2,] = c (u[2]* u[1] * (1- cos(t)) + u[3] * sin(t), cos(t) + u[2]^2 *(1 - cos(t)), u[2]*u[3] * (1 - cos(t)) - u[1] *sin(t) );
  R[3,] = c ( u[3]* u[1] * (1 - cos(t)) - u[2] * sin(t), u[3]* u[2]* (1 - cos(t)) + u[1] * sin(t), cos(t) + u[3]^2 * (1 - cos(t)) );
  return (R);
}

# Given a rotation matrix and a list of vectors, return the accompanying
# list of rotated vectors.

rotveclist_3D <- function (R, vlist)
{
  # Create an empty vector list for the rotation
  rvlist <-matrix(NA, nrow = dim(vlist)[1], ncol = 3) # Define size of the matrix to fill in the loop
  # Run through and perform the rotatio 
  for( i in 1: dim(vlist) [1] ){
    rvlist [i, ] = mmult_3D (R, vlist [i, ]);
  }
  return (rvlist);
}

# Compute the angle between [-\pi, pi) given ajacent and opposite sides
angle_tan <- function(adj, opp)
{
  # When the vector is at pi/2 or -pi/2, can't feed a vector that's zero though.
  if ( isTRUE (adj == 0) ) {
    if ( isTRUE (opp == 0) ) {
      angle <- NaN		# Zero vector, no angle.
    } else if ( isTRUE (opp < 0 ) ) {
      angle <- - pi/2		#
    } else {
      angle <- pi/2
    }
    
    return (angle)
  }
  
  # Outside of the pi/2 region.
  angle <- atan (opp/adj)
  
  return (angle)
}


# This returns the angle of the vector in the Y-Z plane. A rotation around the X axis of this angle
# brings the Z component to zero.
X_to_Z_rot_angle_3D = function (v)
{
  return (angle_tan (v[2],v[3]) )
}

# The angle in the X-Z plane. A rotation around the Y axis of this angle brings the Y component to zero.
Y_to_Z_rot_angle_3D = function (v)
{
  return (angle_tan (v[1],v[3]) )
}

# function vector.angle (vector from angle)
vector.angle <- function (angle)
{
  theta <- angle* pi/180;
  return( c(cos(theta), sin(theta), 0) ); # This simply produces a *unit* vector to the specified angle theta (in radians)
}

translate <- function (T, landmark, origin = c(0,0,0))
{
  if(!inherits(T, "array")){
    stop("dataset must be of class 'array'")
  }
  translist <- array(data = NA, dim = c(dim(T)[1], dim(T)[2], dim(T)[3]),dimnames = list(NULL, NULL, c(dimnames(T)[[3]])));
  for( i in 1:dim(T)[3] ){
    for (j in 1: dim(T)[1]) {
      translist [j,,i] <- T [j,,i] - T[landmark,,i] + origin
    }
  }
  return (translist);
}

# Function that rotates an array, given a landmark and a vector (to which the rotation should take place).
# The resulting array will be such that the landmark arr[land,] will actually be on the same line, including orientation
# as the vector vec.
rotate.orientation <- function(arr, land, vec)
{
  # Create an array of the right dimensions
  rot <- matrix(NA, nrow = dim(arr)[1], ncol = 3)
  # Lets normalise and use unit vectors (vec  can be non-unital)
  uvec <- unit_3D(vec)
  # Rotation axis
  axis <- cross_3D(arr[land,], uvec)
  # ROtation matrix
  rotmat <- rotmat_3D(axis, angle_3D(arr[land,], uvec))
  # Rotate and return
  rot <- rotveclist_3D(rotmat, arr)
  
  # Check if rotmat[land,] is in the same orentation as vec.
  # Logic: svec which is the rotated vector should be the same as uvec if it is in the right orientation.
  # Unfortunately, due to numerical error, we can't just equality, which is equivalent to asking that ||svec - uvec|| is close to zero.
  # However, the numerically stable solution is to check that || svec - uvec || > || svec + uvec ||, which means that you have the wrong orientation.
  # In that case, we rotate again on the axis given by the variable "axis", but this time, we go around \pi radians.
  svec <- unit_3D(rot[land,])
  if ( isTRUE (norm_3D(svec - uvec) > norm_3D(svec + uvec)) ) {
    rotmat <- rotmat_3D(axis, pi)
    rot <- rotveclist_3D(rotmat, rot)
  }
  
  return(rot)
}

