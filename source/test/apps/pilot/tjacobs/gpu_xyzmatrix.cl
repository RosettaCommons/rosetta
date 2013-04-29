struct xyzmatrix
{
  float xx_, xy_, xz_;
  float yx_, yy_, yz_;
  float zx_, zy_, zz_;
};

struct xyzvector
{
  float x_, y_, z_;
};

/// Modify matrix 1, m1, by performing a right-multiply into it of m2
void
right_multiply_by_private(
  __private struct xyzmatrix * m1,
  __private struct xyzmatrix * m2
)
{
  float x, y, z; // Temporaries

  // First row
  x = ( m1->xx_ * m2->xx_ ) + ( m1->xy_ * m2->yx_ ) + ( m1->xz_ * m2->zx_ );
  y = ( m1->xx_ * m2->xy_ ) + ( m1->xy_ * m2->yy_ ) + ( m1->xz_ * m2->zy_ );
  z = ( m1->xx_ * m2->xz_ ) + ( m1->xy_ * m2->yz_ ) + ( m1->xz_ * m2->zz_ );
  m1->xx_ = x; m1->xy_ = y; m1->xz_ = z;

  // Second row
  x = ( m1->yx_ * m2->xx_ ) + ( m1->yy_ * m2->yx_ ) + ( m1->yz_ * m2->zx_ );
  y = ( m1->yx_ * m2->xy_ ) + ( m1->yy_ * m2->yy_ ) + ( m1->yz_ * m2->zy_ );
  z = ( m1->yx_ * m2->xz_ ) + ( m1->yy_ * m2->yz_ ) + ( m1->yz_ * m2->zz_ );
  m1->yx_ = x; m1->yy_ = y; m1->yz_ = z;

  // Third row
  x = ( m1->zx_ * m2->xx_ ) + ( m1->zy_ * m2->yx_ ) + ( m1->zz_ * m2->zx_ );
  y = ( m1->zx_ * m2->xy_ ) + ( m1->zy_ * m2->yy_ ) + ( m1->zz_ * m2->zy_ );
  z = ( m1->zx_ * m2->xz_ ) + ( m1->zy_ * m2->yz_ ) + ( m1->zz_ * m2->zz_ );
  m1->zx_ = x; m1->zy_ = y; m1->zz_ = z;
}

/// Modify matrix 1, m1, by performing a right-multiply into it of m2-transposed
void
right_multiply_by_transpose_private(
  __private struct xyzmatrix * m1,
  __private struct xyzmatrix * m2
)
{
  float x, y, z; // Temporaries

  // First row
  x = ( m1->xx_ * m2->xx_ ) + ( m1->xy_ * m2->xy_ ) + ( m1->xz_ * m2->xz_ );
  y = ( m1->xx_ * m2->yx_ ) + ( m1->xy_ * m2->yy_ ) + ( m1->xz_ * m2->yz_ );
  z = ( m1->xx_ * m2->zx_ ) + ( m1->xy_ * m2->zy_ ) + ( m1->xz_ * m2->zz_ );
  m1->xx_ = x; m1->xy_ = y; m1->xz_ = z;

  // Second row
  x = ( m1->yx_ * m2->xx_ ) + ( m1->yy_ * m2->xy_ ) + ( m1->yz_ * m2->xz_ );
  y = ( m1->yx_ * m2->yx_ ) + ( m1->yy_ * m2->yy_ ) + ( m1->yz_ * m2->yz_ );
  z = ( m1->yx_ * m2->zx_ ) + ( m1->yy_ * m2->zy_ ) + ( m1->yz_ * m2->zz_ );
  m1->yx_ = x; m1->yy_ = y; m1->yz_ = z;

  // Third row
  x = ( m1->zx_ * m2->xx_ ) + ( m1->zy_ * m2->xy_ ) + ( m1->zz_ * m2->xz_ );
  y = ( m1->zx_ * m2->yx_ ) + ( m1->zy_ * m2->yy_ ) + ( m1->zz_ * m2->yz_ );
  z = ( m1->zx_ * m2->zx_ ) + ( m1->zy_ * m2->zy_ ) + ( m1->zz_ * m2->zz_ );
  m1->zx_ = x; m1->zy_ = y; m1->zz_ = z;

}


/// Modify matrix 1, m1, by performing a left-multiply of matrix m2 into it.
void
left_multiply_by_private(
  __private struct xyzmatrix * m1,
  __private struct xyzmatrix * m2
)
{
  float x, y, z; // Temporaries

  // First column
  x = ( m2->xx_ * m1->xx_ ) + ( m2->xy_ * m1->yx_ ) + ( m2->xz_ * m1->zx_ );
  y = ( m2->yx_ * m1->xx_ ) + ( m2->yy_ * m1->yx_ ) + ( m2->yz_ * m1->zx_ );
  z = ( m2->zx_ * m1->xx_ ) + ( m2->zy_ * m1->yx_ ) + ( m2->zz_ * m1->zx_ );
  m1->xx_ = x; m1->yx_ = y; m1->zx_ = z;

  // Second column
  x = ( m2->xx_ * m1->xy_ ) + ( m2->xy_ * m1->yy_ ) + ( m2->xz_ * m1->zy_ );
  y = ( m2->yx_ * m1->xy_ ) + ( m2->yy_ * m1->yy_ ) + ( m2->yz_ * m1->zy_ );
  z = ( m2->zx_ * m1->xy_ ) + ( m2->zy_ * m1->yy_ ) + ( m2->zz_ * m1->zy_ );
  m1->xy_ = x; m1->yy_ = y; m1->zy_ = z;

  // Third column
  x = ( m2->xx_ * m1->xz_ ) + ( m2->xy_ * m1->yz_ ) + ( m2->xz_ * m1->zz_ );
  y = ( m2->yx_ * m1->xz_ ) + ( m2->yy_ * m1->yz_ ) + ( m2->yz_ * m1->zz_ );
  z = ( m2->zx_ * m1->xz_ ) + ( m2->zy_ * m1->yz_ ) + ( m2->zz_ * m1->zz_ );
  m1->xz_ = x; m1->yz_ = y; m1->zz_ = z;

}

/// Modify matrix 1, m1, by performing a left-multiply into it of m2-transposed
void
left_multiply_by_transpose_private(
  __private struct xyzmatrix * m1,
  __private struct xyzmatrix * m2
)
{
  float x, y, z; // Temporaries

  // First column
  x = ( m2->xx_ * m1->xx_ ) + ( m2->yx_ * m1->yx_ ) + ( m2->zx_ * m1->zx_ );
  y = ( m2->xy_ * m1->xx_ ) + ( m2->yy_ * m1->yx_ ) + ( m2->zy_ * m1->zx_ );
  z = ( m2->xz_ * m1->xx_ ) + ( m2->yz_ * m1->yx_ ) + ( m2->zz_ * m1->zx_ );
  m1->xx_ = x; m1->yx_ = y; m1->zx_ = z;

  // Second column
  x = ( m2->xx_ * m1->xy_ ) + ( m2->yx_ * m1->yy_ ) + ( m2->zx_ * m1->zy_ );
  y = ( m2->xy_ * m1->xy_ ) + ( m2->yy_ * m1->yy_ ) + ( m2->zy_ * m1->zy_ );
  z = ( m2->xz_ * m1->xy_ ) + ( m2->yz_ * m1->yy_ ) + ( m2->zz_ * m1->zy_ );
  m1->xy_ = x; m1->yy_ = y; m1->zy_ = z;

  // Third column
  x = ( m2->xx_ * m1->xz_ ) + ( m2->yx_ * m1->yz_ ) + ( m2->zx_ * m1->zz_ );
  y = ( m2->xy_ * m1->xz_ ) + ( m2->yy_ * m1->yz_ ) + ( m2->zy_ * m1->zz_ );
  z = ( m2->xz_ * m1->xz_ ) + ( m2->yz_ * m1->yz_ ) + ( m2->zz_ * m1->zz_ );
  m1->xz_ = x; m1->yz_ = y; m1->zz_ = z;

}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
void set_xyzmatrix_value_private(
  __private struct xyzmatrix * m,
  int row,
  int column,
  float value
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    m->xx_ = value; break;
  case 1:
    m->xy_ = value; break;
  case 2:
    m->xz_ = value; break;
  case 3 :
    m->yx_ = value; break;
  case 4:
    m->yy_ = value; break;
  case 5:
    m->yz_ = value; break;
  case 6 :
    m->zx_ = value; break;
  case 7:
    m->zy_ = value; break;
  case 8:
    m->zz_ = value; break;
  default :
    // ERROR!
    // The GPU can't exactly throw an error, so
    // try to leave a message in the input xyzmatrix
    // saying something went wrong.
    m->xx_ = m->xy_ = m->xz_ = 0;
    m->yx_ = m->yy_ = m->yz_ = 0;
    m->zx_ = m->zy_ = m->zz_ = 0;
    break;
  }
}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
float get_xyzmatrix_value_private(
  __private struct xyzmatrix * m,
  int row,
  int column
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    return m->xx_;
  case 1:
    return m->xy_;
  case 2:
    return m->xz_;
  case 3 :
    return m->yx_;
  case 4:
    return m->yy_;
  case 5:
    return m->yz_;
  case 6 :
    return m->zx_;
  case 7:
    return m->zy_;
  case 8:
    return m->zz_;
  }
  return 1/0; // error rowmajor > 8
}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
void set_xyzmatrix_value_global(
  __global struct xyzmatrix * m,
  int row,
  int column,
  float value
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    m->xx_ = value; break;
  case 1:
    m->xy_ = value; break;
  case 2:
    m->xz_ = value; break;
  case 3 :
    m->yx_ = value; break;
  case 4:
    m->yy_ = value; break;
  case 5:
    m->yz_ = value; break;
  case 6 :
    m->zx_ = value; break;
  case 7:
    m->zy_ = value; break;
  case 8:
    m->zz_ = value; break;
  default :
    // ERROR!
    // The GPU can't exactly throw an error, so
    // try to leave a message in the input xyzmatrix
    // saying something went wrong.
    m->xx_ = m->xy_ = m->xz_ = 0;
    m->yx_ = m->yy_ = m->yz_ = 0;
    m->zx_ = m->zy_ = m->zz_ = 0;
    break;
  }
}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
float get_xyzmatrix_value_global(
  __global struct xyzmatrix * m,
  int row,
  int column
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    return m->xx_;
  case 1:
    return m->xy_;
  case 2:
    return m->xz_;
  case 3 :
    return m->yx_;
  case 4:
    return m->yy_;
  case 5:
    return m->yz_;
  case 6 :
    return m->zx_;
  case 7:
    return m->zy_;
  case 8:
    return m->zz_;
  }
  return 1/0; // error rowmajor > 8
}

/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
void set_xyzmatrix_value_local(
  __local struct xyzmatrix * m,
  int row,
  int column,
  float value
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    m->xx_ = value; break;
  case 1:
    m->xy_ = value; break;
  case 2:
    m->xz_ = value; break;
  case 3 :
    m->yx_ = value; break;
  case 4:
    m->yy_ = value; break;
  case 5:
    m->yz_ = value; break;
  case 6 :
    m->zx_ = value; break;
  case 7:
    m->zy_ = value; break;
  case 8:
    m->zz_ = value; break;
  default :
    // ERROR!
    // The GPU can't exactly throw an error, so
    // try to leave a message in the input xyzmatrix
    // saying something went wrong.
    m->xx_ = m->xy_ = m->xz_ = 0;
    m->yx_ = m->yy_ = m->yz_ = 0;
    m->zx_ = m->zy_ = m->zz_ = 0;
    break;
  }
}


/// Set the value for one of the 9 elements in a 3x3 xyzmatrix given a row
/// index / column index pair, indexing from 0.
/// The first row is:     xx, xy, xz;
/// the second row is:    yx, yy, yz;
/// and the third row is: zx, zy, zz;
float get_xyzmatrix_value_local(
  __local struct xyzmatrix * m,
  int row,
  int column
)
{
  int rowmajor = row*3 + column;
  switch ( rowmajor ) {
  case 0 :
    return m->xx_;
  case 1:
    return m->xy_;
  case 2:
    return m->xz_;
  case 3 :
    return m->yx_;
  case 4:
    return m->yy_;
  case 5:
    return m->yz_;
  case 6 :
    return m->zx_;
  case 7:
    return m->zy_;
  case 8:
    return m->zz_;
  }
  return 1/0; // error rowmajor > 8
}

/// set a value in an xyzvector, indexing from 0
void
set_xyzvector_value_private(
  __private struct xyzvector * v,
  int index,
  float value
)
{
  switch ( index ) {
  case 0 :
    v->x_ = value; break;
  case 1:
    v->y_ = value; break;
  case 2:
    v->z_ = value; break;
  } 
}

/// get a value out of an xyzvector, indexing from 0
float
get_xyzvector_value_private(
  __private struct xyzvector * v,
  int index
)
{
  switch ( index ) {
  case 0 :
    return v->x_;
  case 1 :
    return v->y_;
  case 2 :
    return v->z_;
  }
}


/// set a value in an xyzvector, indexing from 0
void
set_xyzvector_value_local(
  __local struct xyzvector * v,
  int index,
  float value
)
{
  switch ( index ) {
  case 0 :
    v->x_ = value; break;
  case 1:
    v->y_ = value; break;
  case 2:
    v->z_ = value; break;
  } 
}

/// get a value out of an xyzvector, indexing from 0
float
get_xyzvector_value_local(
  __local struct xyzvector * v,
  int index
)
{
  switch ( index ) {
  case 0 :
    return v->x_;
  case 1 :
    return v->y_;
  case 2 :
    return v->z_;
  }
}


/// set a value in an xyzvector, indexing from 0
void
set_xyzvector_value_global(
  __global struct xyzvector * v,
  int index,
  float value
)
{
  switch ( index ) {
  case 0 :
    v->x_ = value; break;
  case 1:
    v->y_ = value; break;
  case 2:
    v->z_ = value; break;
  } 
}

/// get a value out of an xyzvector, indexing from 0
float
get_xyzvector_value_global(
  __global struct xyzvector * v,
  int index
)
{
  switch ( index ) {
  case 0 :
    return v->x_;
  case 1 :
    return v->y_;
  case 2 :
    return v->z_;
  }
}

void
set_xyzmatrix_to_identity_private(
  __private struct xyzmatrix * m
)
{
  m->xx_ = 1.0f; m->xy_ = 0.0f; m->xz_ = 0.0f;
  m->yx_ = 0.0f; m->yy_ = 1.0f; m->yz_ = 0.0f;
  m->zx_ = 0.0f; m->zy_ = 0.0f; m->zz_ = 1.0f;
}

/// swap the off-diagonal entries
void
xyzmatrix_transpose_private(
  __private struct xyzmatrix * m
)
{
  float temp;
  temp = m->xy_;
  m->xy_ = m->yx_;
  m->yx_ = temp;

  temp = m->xz_;
  m->xz_ = m->zx_;
  m->zx_ = temp;

  temp = m->yz_;
  m->yz_ = m->zy_;
  m->zy_ = temp;
}


/// compute m*v
struct xyzvector
xyzmatrix_xyzvector_multiply_private(
  __private struct xyzmatrix * m,
  __private struct xyzvector * v
)
{
  struct xyzvector p;
  p.x_ = m->xx_ * v->x_ + m->xy_ * v->y_ + m->xz_ * v->z_;
  p.y_ = m->yx_ * v->x_ + m->yy_ * v->y_ + m->yz_ * v->z_;
  p.z_ = m->zx_ * v->x_ + m->zy_ * v->y_ + m->zz_ * v->z_;
  return p;
}

float
xyzvector_dot_product_private(
  __private struct xyzvector * v1,
  __private struct xyzvector * v2
)
{
  return v1->x_*v2->x_ + v1->y_*v2->y_ + v1->z_*v2->z_;
}

// compute the square magnitude of an xyzvector
float
xyzvector_square_magnitude_private(
  __private struct xyzvector * v
)
{
  return v->x_*v->x_ + v->y_*v->y_ + v->z_*v->z_;
}

/// compute the square distance between two points, held in xyzvectors
float
xyzvector_square_distance_private(
  __private struct xyzvector * v1,
  __private struct xyzvector * v2
)
{
  float sum2 = 0;
  float val;
  val = v1->x_ - v2->x_;
  sum2 += val*val;
  val = v1->y_ - v2->y_;
  sum2 += val*val;
  val = v1->z_ - v2->z_;
  sum2 += val*val;
  return sum2;
}

void
jacobi_rotation(
  __private struct xyzmatrix * m,
  int i,
  int j,
  __private struct xyzmatrix * r
)
{
  // precondition: 0 <= i < 3
  // precondition: 0 <= j < 3
  // precondition: i != j

  float const tau = ( get_xyzmatrix_value_private( m,i,i) - get_xyzmatrix_value_private(m,j,j) ) / ( 2 * get_xyzmatrix_value_private(m,i,j) );
  float const t = ( tau < 0.f ? -1.f : 1.f ) / ( fabs( tau ) + sqrt( 1.f + ( tau * tau ) ) );

  float const c = 1.f / sqrt( 1.f + ( t * t ) );
  float const s = c * t;

  set_xyzmatrix_to_identity_private( r );
  set_xyzmatrix_value_private(r,i,i, c); set_xyzmatrix_value_private(r,i,j,-s);
  set_xyzmatrix_value_private(r,j,i, s); set_xyzmatrix_value_private(r,j,j, c);
}

struct xyzvector
eigenvector_jacobi(
  __private struct xyzmatrix * a,
  float tol,
  __private struct xyzmatrix * J
)
{
  // Copy matrix as it will be modified by the algorithm
  struct xyzmatrix m = *a;

  // Initialize the off-diagonal, upper-triangular sum of squares
  float off = ( m.xy_ * m.xy_ ) + ( m.xz_ * m.xz_ ) + ( m.yz_ * m.yz_ );

  set_xyzmatrix_to_identity_private( J );
  int i, j, n_iterations = 0;
  struct xyzmatrix r;
  while ( off > tol ) {
    ++n_iterations;

    // Ensure number of iterations does not exceed 50:
    // otherwise, re-evaluate stopping criterion
    // assert( n_iterations <= 50 );
    if ( n_iterations == 50 ) {
      break;
    }

    // Determine index of upper-triangular element that will be zeroed out
    if ( fabs( m.xy_ ) >= fabs( m.xz_ ) ) {
      if ( fabs( m.xy_ ) >= fabs( m.yz_ ) ) {
	i = 0; j = 1;
      } else {
	i = 1; j = 2;
      }
    } else if ( fabs( m.xz_ ) >= fabs( m.yz_ ) ) {
      i = 0; j = 2;
    } else {
      i = 1; j = 2;
    }

    // After four iterations, skip the rotation if the off-diagonal element is small
    float const ij_scaled = fabs( 100.0f * get_xyzmatrix_value_private(&m,i,j) );
    if ( ( n_iterations > 4 )
        && fabs( get_xyzmatrix_value_private(&m,i,i) ) + ij_scaled == fabs( get_xyzmatrix_value_private(&m,i,i) )
        && fabs( get_xyzmatrix_value_private(&m,j,j) ) + ij_scaled == fabs( get_xyzmatrix_value_private(&m,j,j) ) ) {
      set_xyzmatrix_value_private( &m, i, j, 0.f );
      set_xyzmatrix_value_private( &m, j, i, 0.f );
    } else {
      // Compute the rotation matrix
      jacobi_rotation( &m, i, j, &r );

      // Zero out the i,j and j,i elements
      right_multiply_by_private( &m, &r );
      left_multiply_by_transpose_private( &m, &r );

      // Accumulate the rotation transformations to form the matrix of eigenvectors
      right_multiply_by_private( J, &r );
    }

    // Recalculate the off-diagonal, upper-triangular sum of squares
    off = ( m.xy_ * m.xy_ ) + ( m.xz_ * m.xz_ ) + ( m.yz_ * m.yz_ );
  }

  struct xyzvector eigenvals;
  eigenvals.x_ = m.xx_;
  eigenvals.y_ = m.yy_;
  eigenvals.z_ = m.zz_;
  return eigenvals;

}


/// replaces the third eigenvector by taking cross product of
/// of the first two eigenvectors
void
fixEigenvector(
  __private struct xyzmatrix * m
)
{
  m->xz_ = m->yx_*m->zy_ - m->zx_*m->yy_;
  m->yz_ = m->zx_*m->xy_ - m->xx_*m->zy_;
  m->zz_ = m->xx_*m->yy_ - m->yx_*m->xy_;
  float norm = sqrt( 1 / (m->xz_*m->xz_ + m->yz_*m->yz_ + m->zz_*m->zz_) );
  m->xz_ *= norm;
  m->yz_ *= norm;
  m->zz_ *= norm;
}

/// Transform the coordinates to align their center of mass with the origin
/// and return the "moment of inertia" -- the sum of the square-distance
/// of those coordinates from the origin -- of those transformed coordinates
float
transform_center_of_mass_to_origin(
  __global struct xyzvector * XX,
  __private struct xyzvector * center_of_mass,
  int Npoints
)
{

  float sumx, sumy, sumz = 0.0;
  for ( int ii = 0; ii < Npoints; ++ii ) {
    sumx += XX[ii].x_;
    sumy += XX[ii].y_;
    sumz += XX[ii].z_;
  }
  sumx /= Npoints;
  sumy /= Npoints;
  sumz /= Npoints;
  center_of_mass->x_ = sumx;
  center_of_mass->y_ = sumy;
  center_of_mass->z_ = sumz;
  float moment_of_inertia = 0;
  for ( int ii = 0; ii < Npoints; ++ii ) {
    struct xyzvector iicoord = XX[ ii ];
    iicoord.x_ -= sumx;
    iicoord.y_ -= sumy;
    iicoord.z_ -= sumz;
    moment_of_inertia += iicoord.x_*iicoord.x_ + iicoord.y_*iicoord.y_ + iicoord.z_*iicoord.z_;
    XX[ii] = iicoord;
  }
  return moment_of_inertia;
}

/// RMS superposition routine taken from rms.cc, and adapted to use the gpu
/// data structures / functions / 0-based indexing that has
/// so far been written
void
findUU(
  __global struct xyzvector * XX,
  __global struct xyzvector * YY,
  __global float * WW,
  int Npoints,
  __private struct xyzmatrix * UU,
  __private float * sigma3//,
  //__global struct xyzmatrix * m_moment_out
)
{
  int sort[3];
  struct xyzmatrix eVec;
  struct xyzmatrix bb;
  struct xyzvector w_w;
  struct xyzmatrix m_moment;
  struct xyzmatrix rr_moment;
  float temp1;
  float temp2;
  float temp3;
  struct xyzvector Ra;

  if ( Npoints < 1 ) {
    // return identity rotation matrix to moron
    set_xyzmatrix_to_identity_private( UU );
    *sigma3 = 0.0;
    return;
  }

  // align center of mass to origin
  // apl -- reorder these loops to get better memory access patterns
  for ( int k = 0; k < 3; ++k ) {
    temp1 = 0.0;
    temp2 = 0.0;
    temp3 = 0.0;
    for ( int j = 0; j < Npoints; ++j ) {
      float jweight = WW[j];
      temp1 += get_xyzvector_value_global( &XX[j], k ) * jweight;
      temp2 += get_xyzvector_value_global( &YY[j], k ) * jweight;
      temp3 += jweight;
    }
    if (temp3 > 0.001) temp1 /= temp3;
    if (temp3 > 0.001) temp2 /= temp3;

    for ( int j = 0; j < Npoints; ++j ) {
      set_xyzvector_value_global( &XX[j], k, get_xyzvector_value_global( &XX[j], k ) - temp1 );
      set_xyzvector_value_global( &YY[j], k, get_xyzvector_value_global( &YY[j], k ) - temp2 );
    }
  }

  // Make cross moments matrix   INCLUDE THE WEIGHTS HERE
  for ( int k = 0; k < 3; ++k ) {
    for ( int j = 0; j < 3; ++j ) {
      temp1 = 0.0;
      for ( int i = 0; i < Npoints; ++i ) {
	temp1 += WW[i] * get_xyzvector_value_global(&YY[i],k) * get_xyzvector_value_global(&XX[i],j);
      }
      set_xyzmatrix_value_private( &m_moment, k, j, temp1 );
    }
  }

  // Multiply CROSS MOMENTS by transpose
  //BlankMatrixMult(m_moment,3,3,1,m_moment,3,0,rr_moment);
  rr_moment = m_moment;
  left_multiply_by_transpose_private( & rr_moment, & m_moment );

  // Copy to/from xyzMatrix/xyzVector since rest of functions use FArrays
  //xyzMatrix< numeric::Real > xyz_rr_moment( xyzMatrix< numeric::Real >::cols( &rr_moment( 1,1 ) ) );
  //xyzVector< numeric::Real > xyz_w_w;
  //xyzMatrix< numeric::Real > xyz_eVec;

  // Find eigenvalues, eigenvectors of symmetric matrix rr_moment
  w_w = eigenvector_jacobi( & rr_moment, 1E-9, & eVec );

  // Copy eigenvalues/vectors back to FArray
  //for ( int i = 1; i <= 3; ++i ) {
  //  w_w( i ) = xyz_w_w( i );
  //  for ( int j = 1; j <= 3; ++j ) {
  //    eVec( i, j ) = xyz_eVec( i, j );
  //  }
  //}

  // explicitly coded 3 level index sort using eigenvalues
  for ( int i = 0; i < 3; ++i ) {
    sort[ i ] = i;
  }

  if ( get_xyzvector_value_private( &w_w,0) < get_xyzvector_value_private(&w_w,1) ) {
    sort[ 1 ] = 0;
    sort[ 0 ] = 1;
  }

  if ( get_xyzvector_value_private( &w_w, sort[ 1 ] ) < get_xyzvector_value_private( &w_w, 2) ) {
    sort[ 2 ] = sort[ 1 ];
    sort[ 1 ] = 2;

    if ( get_xyzvector_value_private( &w_w, sort[ 0 ] ) < get_xyzvector_value_private( &w_w, 2) ) {
      sort[ 1 ] = sort[ 0 ];
      sort[ 0 ] = 2;
    }
  }
  // sort is now an index to order of eigen values

  if ( get_xyzvector_value_private( &w_w, sort[1] ) == 0.0 ) { // holy smokes, two eigen values are zeros
    // return identity rotation matrix to moron
    set_xyzmatrix_to_identity_private( UU );
    if ( get_xyzvector_value_private( &w_w, sort[0] ) < 0.0f ) {
      set_xyzvector_value_private( &w_w, sort[0], fabs( get_xyzvector_value_private( &w_w, sort[0]) ) );
    }
    *sigma3 = sqrt( get_xyzvector_value_private( &w_w, sort[0] ));

    return; // make like a prom dress and slip off
  }

  // sort eigen values
  temp1 = get_xyzvector_value_private( &w_w, sort[0] );
  temp2 = get_xyzvector_value_private( &w_w, sort[1] );
  set_xyzvector_value_private( &w_w, 2, get_xyzvector_value_private( &w_w, sort[2] ) );
  set_xyzvector_value_private( &w_w, 1, temp2 );
  set_xyzvector_value_private( &w_w, 0, temp1 );
  // sort first two eigen vectors (dont care about third)
  for ( int i = 0; i < 3; ++i ) {
    temp1 = get_xyzmatrix_value_private( &eVec, i, sort[0] );
    temp2 = get_xyzmatrix_value_private( &eVec, i, sort[1] );
    set_xyzmatrix_value_private( &eVec, i, 0, temp1 );
    set_xyzmatrix_value_private( &eVec, i, 1, temp2 );
  }

  // april 20: the fix not only fixes bad eigen vectors but solves a problem of
  // forcing a right-handed coordinate system

  fixEigenvector( &eVec);

  // at this point we now have three good eigenvectors in a right hand
  // coordinate system.

  // make bb basis vectors   = moments*eVec

  //BlankMatrixMult(m_moment,3,3,0,eVec,3,0,bb);
  bb = eVec;
  left_multiply_by_private( &bb, &m_moment );

  //     std::cerr << "m_moment" << std::endl;
  // squirrel away a free copy of the third eigenvector before normalization/fix
  //for ( int j = 1; j <= 3; ++j ) {
  //  Ra(j) = bb(j,3);
  //}
  Ra.x_ = bb.xz_;
  Ra.y_ = bb.yz_;
  Ra.z_ = bb.zz_;

  //m_moment_out->xx_ = Ra.x_;
  //m_moment_out->xy_ = Ra.y_;
  //m_moment_out->xz_ = Ra.z_;

  // normalize first two bb-basis vectors
  // dont care about third since were going to replace it with b1xb2
  // this also avoids problem of possible zero third eigen value
  for ( int j = 0; j < 2; ++j ) {
    temp1 = 1.0/sqrt(get_xyzvector_value_private( &w_w, j )); // zero checked for above
    for ( int k = 0; k < 3; ++k ) { // x,y,z
      set_xyzmatrix_value_private( &bb, k, j, get_xyzmatrix_value_private( &bb, k, j ) * temp1 );
    }
  }


  //  fix things so that bb eigenvecs are right handed
  fixEigenvector(&bb); // need to fix this one too
  //*m_moment_out = bb;

  // find  product of eVec and bb matrices
  // BlankMatrixMult(eVec,3,3,0,bb,3,1,UU);
  // result is returned in UU.
  //*UU = bb;
  //xyzmatrix_transpose_private( UU );
  //left_multiply_by_private( UU, & eVec );
  *UU = eVec;
  right_multiply_by_transpose_private( UU, & bb );

  //*m_moment_out = *UU;

  // and lastly determine a value used in another function to compute the rms
  float private_sigma3 = 0.0;
  for ( int j = 0; j < 3; ++j ) {
    private_sigma3 += get_xyzmatrix_value_private(&bb,j,2) * get_xyzvector_value_private(&Ra,j);
  }

  //m_moment_out->xx_ = private_sigma3;

  //cems the abs() fixes some round off error situations where the w_w values are
  //cems very small and accidentally negative.  (theoretically they are positive,
  //cems but in practice round off error makes them negative)
  if ( private_sigma3 < 0.0 ) {
    *sigma3 = sqrt(fabs( get_xyzvector_value_private( &w_w, 0)) ) +
      sqrt(fabs(get_xyzvector_value_private( &w_w, 1 )) ) -
      sqrt(fabs(get_xyzvector_value_private( &w_w, 2 )) );
  } else {

    *sigma3 = sqrt(fabs( get_xyzvector_value_private( &w_w, 0)) ) +
      sqrt(fabs(get_xyzvector_value_private( &w_w, 1 )) ) +
      sqrt(fabs(get_xyzvector_value_private( &w_w, 2 )) );

  }


}


/// Variant on findUU where the coordinates XX and YY are already assumed to have been
/// translated to the origin, and the weight vector WW is assumed to be one everywhere.
void
findUU_no_translation_to_origin(
  __global struct xyzvector * XX,
  __global struct xyzvector * YY,
  int Npoints,
  __private struct xyzmatrix * UU,
  __private float * sigma3//,
  //__global struct xyzmatrix * m_moment_out
)
{
  int sort[3];
  struct xyzmatrix eVec;
  struct xyzmatrix bb;
  struct xyzvector w_w;
  struct xyzmatrix m_moment;
  struct xyzmatrix rr_moment;
  float temp1;
  float temp2;
  float temp3;
  struct xyzvector Ra;

  if ( Npoints < 1 ) {
    // return identity rotation matrix to moron
    set_xyzmatrix_to_identity_private( UU );
    *sigma3 = 0.0;
    return;
  }

  // Make cross moments matrix
  for ( int k = 0; k < 3; ++k ) {
    for ( int j = 0; j < 3; ++j ) {
      temp1 = 0.0;
      for ( int i = 0; i < Npoints; ++i ) {
	temp1 += get_xyzvector_value_global(&YY[i],k) * get_xyzvector_value_global(&XX[i],j);
      }
      set_xyzmatrix_value_private( &m_moment, k, j, temp1 );
    }
  }

  // Multiply CROSS MOMENTS by transpose
  //BlankMatrixMult(m_moment,3,3,1,m_moment,3,0,rr_moment);
  rr_moment = m_moment;
  left_multiply_by_transpose_private( & rr_moment, & m_moment );

  // Copy to/from xyzMatrix/xyzVector since rest of functions use FArrays
  //xyzMatrix< numeric::Real > xyz_rr_moment( xyzMatrix< numeric::Real >::cols( &rr_moment( 1,1 ) ) );
  //xyzVector< numeric::Real > xyz_w_w;
  //xyzMatrix< numeric::Real > xyz_eVec;

  // Find eigenvalues, eigenvectors of symmetric matrix rr_moment
  w_w = eigenvector_jacobi( & rr_moment, 1E-9, & eVec );
  //if ( m_moment_out ) { m_moment_out->xx_ = w_w.x_; m_moment_out->xy_ = w_w.y_; m_moment_out->xz_ = w_w.z_; }

  // Copy eigenvalues/vectors back to FArray
  //for ( int i = 1; i <= 3; ++i ) {
  //  w_w( i ) = xyz_w_w( i );
  //  for ( int j = 1; j <= 3; ++j ) {
  //    eVec( i, j ) = xyz_eVec( i, j );
  //  }
  //}

  // explicitly coded 3 level index sort using eigenvalues
  for ( int i = 0; i < 3; ++i ) {
    sort[ i ] = i;
  }

  if ( get_xyzvector_value_private( &w_w,0) < get_xyzvector_value_private(&w_w,1) ) {
    sort[ 1 ] = 0;
    sort[ 0 ] = 1;
  }

  if ( get_xyzvector_value_private( &w_w, sort[ 1 ] ) < get_xyzvector_value_private( &w_w, 2) ) {
    sort[ 2 ] = sort[ 1 ];
    sort[ 1 ] = 2;

    if ( get_xyzvector_value_private( &w_w, sort[ 0 ] ) < get_xyzvector_value_private( &w_w, 2) ) {
      sort[ 1 ] = sort[ 0 ];
      sort[ 0 ] = 2;
    }
  }
  // sort is now an index to order of eigen values
  //if ( m_moment_out ) { m_moment_out->yx_ = sort[0]; m_moment_out->yy_ = sort[1]; m_moment_out->yz_ = sort[2]; }

  if ( get_xyzvector_value_private( &w_w, sort[1] ) == 0.0 ) { // holy smokes, two eigen values are zeros
    // return identity rotation matrix to moron
    set_xyzmatrix_to_identity_private( UU );
    if ( get_xyzvector_value_private( &w_w, sort[0] ) < 0.0f ) {
      set_xyzvector_value_private( &w_w, sort[0], fabs( get_xyzvector_value_private( &w_w, sort[0]) ) );
    }
    *sigma3 = sqrt( get_xyzvector_value_private( &w_w, sort[0] ));

    return; // make like a prom dress and slip off
  }

  /// before sorting the eigenvectors
  //if ( m_moment_out ) { *m_moment_out = eVec; }

  // sort eigen values
  temp1 = get_xyzvector_value_private( &w_w, sort[0] );
  temp2 = get_xyzvector_value_private( &w_w, sort[1] );
  set_xyzvector_value_private( &w_w, 2, get_xyzvector_value_private( &w_w, sort[2] ) );
  set_xyzvector_value_private( &w_w, 1, temp2 );
  set_xyzvector_value_private( &w_w, 0, temp1 );
  // sort first two eigen vectors (dont care about third)
  for ( int i = 0; i < 3; ++i ) {
    temp1 = get_xyzmatrix_value_private( &eVec, i, sort[0] );
    temp2 = get_xyzmatrix_value_private( &eVec, i, sort[1] );
    set_xyzmatrix_value_private( &eVec, i, 0, temp1 );
    set_xyzmatrix_value_private( &eVec, i, 1, temp2 );
  }

  /// after sorting the eigenvectors
  //if ( m_moment_out ) { *m_moment_out = eVec; }

  // april 20: the fix not only fixes bad eigen vectors but solves a problem of
  // forcing a right-handed coordinate system

  fixEigenvector( &eVec);

  /// after "fixing" the sorted eigenvectors
  //if ( m_moment_out ) { *m_moment_out = eVec; }


  // at this point we now have three good eigenvectors in a right hand
  // coordinate system.

  // make bb basis vectors   = moments*eVec

  //BlankMatrixMult(m_moment,3,3,0,eVec,3,0,bb);
  bb = eVec;
  left_multiply_by_private( &bb, &m_moment );

  //     std::cerr << "m_moment" << std::endl;
  // squirrel away a free copy of the third eigenvector before normalization/fix
  //for ( int j = 1; j <= 3; ++j ) {
  //  Ra(j) = bb(j,3);
  //}
  Ra.x_ = bb.xz_;
  Ra.y_ = bb.yz_;
  Ra.z_ = bb.zz_;

  //if ( m_moment_out ) { m_moment_out->xx_ = Ra.x_; m_moment_out->xy_ = Ra.y_; m_moment_out->xz_ = Ra.z_; }
  //if ( m_moment_out ) { *m_moment_out = bb; }

  // normalize first two bb-basis vectors
  // dont care about third since were going to replace it with b1xb2
  // this also avoids problem of possible zero third eigen value
  for ( int j = 0; j < 2; ++j ) {
    temp1 = 1.0/sqrt(get_xyzvector_value_private( &w_w, j )); // zero checked for above
    for ( int k = 0; k < 3; ++k ) { // x,y,z
      set_xyzmatrix_value_private( &bb, k, j, get_xyzmatrix_value_private( &bb, k, j ) * temp1 );
    }
  }

  //if ( m_moment_out ) *m_moment_out = bb;

  //  fix things so that bb eigenvecs are right handed
  fixEigenvector(&bb); // need to fix this one too
  //if ( m_moment_out ) *m_moment_out = bb;

  // find  product of eVec and bb matrices
  // BlankMatrixMult(eVec,3,3,0,bb,3,1,UU);
  // result is returned in UU.
  //*UU = bb;
  //xyzmatrix_transpose_private( UU );
  //left_multiply_by_private( UU, & eVec );
  *UU = eVec;
  right_multiply_by_transpose_private( UU, & bb );

  // and lastly determine a value used in another function to compute the rms
  float private_sigma3 = 0.0;
  for ( int j = 0; j < 3; ++j ) {
    private_sigma3 += get_xyzmatrix_value_private(&bb,j,2) * get_xyzvector_value_private(&Ra,j);
  }

  //if ( m_moment_out ) { m_moment_out->yx_ = get_xyzmatrix_value_private(&bb,2,0); m_moment_out->yy_ = get_xyzmatrix_value_private(&bb,2,1); m_moment_out->yz_ = get_xyzmatrix_value_private(&bb,2,2); }

  //if ( m_moment_out ) m_moment_out->xx_ = private_sigma3;

  //cems the abs() fixes some round off error situations where the w_w values are
  //cems very small and accidentally negative.  (theoretically they are positive,
  //cems but in practice round off error makes them negative)
  if ( private_sigma3 < 0.0 ) {
    *sigma3 = sqrt(fabs( get_xyzvector_value_private( &w_w, 0)) ) +
      sqrt(fabs(get_xyzvector_value_private( &w_w, 1 )) ) -
      sqrt(fabs(get_xyzvector_value_private( &w_w, 2 )) );
  } else {

    *sigma3 = sqrt(fabs( get_xyzvector_value_private( &w_w, 0)) ) +
      sqrt(fabs(get_xyzvector_value_private( &w_w, 1 )) ) +
      sqrt(fabs(get_xyzvector_value_private( &w_w, 2 )) );

  }


}

/// The sigma3 value returned by a call to findUU and the
/// moments of inertia of the other two sets of coordinates are
/// all you need in order to calculate the RMSD of two sets of points
float
calculate_rmsd_fast(
  __global float * node1_coord_at_origin_moment_of_inertia,
  __global float * node2_coord_at_origin_moment_of_inertia,
  int node_npoints,
  float sigma3
)
{
  float rms_sum = *node1_coord_at_origin_moment_of_inertia + *node2_coord_at_origin_moment_of_inertia - 2*sigma3;
  return sqrt( fabs( rms_sum / node_npoints ));
}

float
calculate_helix_clash_score( 
  __global struct xyzvector * node1_other_coords,
  int node1_n_other_coords,
  __global int * node1_center_coords,
  int node1_ncenters,
  __global struct xyzvector * node2_other_coords,
  int node2_n_other_coords,
  __global int * node2_center_coords,
  int node2_ncenters,
  __private struct xyzmatrix * UU//,
  //__global float * all_square_distances
)
{
  float closest_contact = 12345;

  //int count_square_distances = -1;
  for ( int ii = 0; ii < node1_ncenters; ++ii ) {
    struct xyzvector iicoord = node1_other_coords[ node1_center_coords[ ii ]];
    for ( int jj = 0; jj < node2_n_other_coords; ++jj ) {
      struct xyzvector jjcoord = node2_other_coords[ jj ];
      struct xyzvector jjrotated = xyzmatrix_xyzvector_multiply_private( UU, &jjcoord );
      float dist2 = xyzvector_square_distance_private( &iicoord, & jjrotated );
      closest_contact = dist2 < closest_contact ? dist2 : closest_contact;
      //all_square_distances[ ++count_square_distances ] = dist2;
    }
  }

  for ( int ii = 0; ii < node2_ncenters; ++ii ) {
    struct xyzvector iicoord = node2_other_coords[ node2_center_coords[ ii ]];
    struct xyzvector iirotated = xyzmatrix_xyzvector_multiply_private( UU, &iicoord );
    for ( int jj = 0; jj < node1_n_other_coords; ++jj ) {
      struct xyzvector jjcoord = node1_other_coords[ jj ];
      float dist2 = xyzvector_square_distance_private( & iirotated, & jjcoord );
      closest_contact = dist2 < closest_contact ? dist2 : closest_contact;
      //all_square_distances[ ++count_square_distances ] = dist2;
    }
  }

  // TEMP DEBUG
  //for ( int ii = 0; ii < node2_n_other_coords; ++ii ) {
  //  struct xyzvector iicoord = node2_other_coords[ ii ];
  //  node2_other_coords[ ii ] = xyzmatrix_xyzvector_multiply_private( UU, &iicoord );
  //}
  return closest_contact;
}

void
calculate_rmsd_and_clash_score_for_bundles(
  __global struct xyzvector * node1_coords_at_origin,
  __global float * node1_coords_at_origin_moment_of_inertia,
  __global struct xyzvector * node2_coords_at_origin,
  __global float * node2_coords_at_origin_moment_of_inertia,
  int node_npoints,
  __global struct xyzvector * node1_other_coords,
  int node1_n_other_coords,
  __global int * node1_center_coords,
  int node1_ncenters,
  __global struct xyzvector * node2_other_coords,
  int node2_n_other_coords,
  __global int * node2_center_coords,
  int node2_ncenters,
  __private float * rms,
  __private float * clash//,
  //__global float * all_square_distances // <-- either debug collision calculation
  //__global struct xyzmatrix * UUout, /// <--- or use the next three to debug the findUU calculation
  //__global float * sigma3out ,
  //__global struct xyzmatrix * m_moment_out
)
{
  struct xyzmatrix UU;
  float sigma3;
  findUU_no_translation_to_origin(
    node1_coords_at_origin,
    node2_coords_at_origin,
    node_npoints, &UU, &sigma3//, m_moment_out
  );

  //*UUout = UU;
  //*sigma3out = sigma3;

  *rms = calculate_rmsd_fast(
    node1_coords_at_origin_moment_of_inertia,
    node2_coords_at_origin_moment_of_inertia,
    node_npoints, sigma3 );

  *clash = calculate_helix_clash_score(
    node1_other_coords, node1_n_other_coords, node1_center_coords, node1_ncenters,
    node2_other_coords, node2_n_other_coords, node2_center_coords, node2_ncenters,
    &UU//, all_square_distances
  );
}

__kernel
void
compute_rmsd_and_clash_scores(
  int block_size,
  int n_nodes_1,
  int n_nodes_2,
  int n_atoms_in_rms_calc,
  __global struct xyzvector * rmsd_coords_1,
  __global struct xyzvector * rmsd_coords_2,
  __global float * rmsd_coords_moment_of_inertia_1,
  __global float * rmsd_coords_moment_of_inertia_2,
  __global struct xyzvector * col_coords_1,
  __global struct xyzvector * col_coords_2,
  __global int * col_coords_offsets_1,
  __global int * col_coords_offsets_2,
  __global int * n_coords_for_col_calc_1,
  __global int * n_coords_for_col_calc_2,
  __global int * comparison_coord_ind_list_1,
  __global int * comparison_coord_ind_list_2,
  __global int * comparison_coord_ind_offset_1,
  __global int * comparison_coord_ind_offset_2,
  __global int * n_comparison_coords_1,
  __global int * n_comparison_coords_2,
  __global float * rmsd_table,
  __global float * collision_table,
  __global unsigned char * calc_rmsd_table//,
  //__global float * all_square_distances
  //__global struct xyzmatrix * UUout,
  //__global float * sigma3out,
  //__global struct xyzmatrix * m_moment_out
)
{
  int n1 = get_global_id(0);
  int n2 = get_global_id(1);

  if ( n1 < n_nodes_1 && n2 < n_nodes_2 && calc_rmsd_table[ n1*block_size + n2 ] == 1 ) {
    float rms, collision;
    //printf( "%x %x :: %x %x\n", rmsd_coords_1, rmsd_coords_2, &rmsd_coords_1[ n1 * n_atoms_in_rms_calc ], &rmsd_coords_2[ n2 * n_atoms_in_rms_calc ] );
    calculate_rmsd_and_clash_score_for_bundles(
      &rmsd_coords_1[ n1 * n_atoms_in_rms_calc ], &rmsd_coords_moment_of_inertia_1[ n1 ],
      &rmsd_coords_2[ n2 * n_atoms_in_rms_calc ], &rmsd_coords_moment_of_inertia_2[ n2 ],
      n_atoms_in_rms_calc,
      &col_coords_1[ col_coords_offsets_1[ n1 ] ], n_coords_for_col_calc_1[ n1 ],
      &comparison_coord_ind_list_1[ comparison_coord_ind_offset_1[ n1 ] ], n_comparison_coords_1[ n1 ], 
      &col_coords_2[ col_coords_offsets_2[ n2 ] ], n_coords_for_col_calc_2[ n2 ],
      &comparison_coord_ind_list_2[ comparison_coord_ind_offset_2[ n2 ] ], n_comparison_coords_2[ n2 ],
      &rms, &collision
      //, UUout, sigma3out, m_moment_out // <-- debugging variables
      //, all_square_distances
    );

    rmsd_table[ n1*block_size + n2 ] = rms;
    collision_table[ n1*block_size + n2 ] = collision;
  }
}

/// find the center of mass of a list of coordinates used in a future RMS calculation
/// and translate them so that their center of mass is at the origin.  Next compute
/// the moment of inertia for those translated points.  Finally, translate a second
/// set of coordinates using the same translation vector as for the first set of
/// coordinates.
__kernel
void
translate_rmscoords_to_origin(
  int n_nodes,
  int n_rms_coords_per_node,
  __global struct xyzvector * rms_coords,
  __global float * node_moments_of_inertia,
  __global struct xyzvector * collision_coords,
  __global int * col_coords_offsets,
  __global int * n_coords_for_col_calc
)
{
  if ( get_global_id(0) < n_nodes ) {
    int node_id = get_global_id(0);
    struct xyzvector rms_coords_center_of_mass;
    float moment_of_inertia = transform_center_of_mass_to_origin(
      & rms_coords[ node_id * n_rms_coords_per_node ],
      & rms_coords_center_of_mass,
      n_rms_coords_per_node );
    int node_num_collision_coords = n_coords_for_col_calc[ node_id ];
    int offset = col_coords_offsets[ node_id ];
    for ( int ii = 0; ii < node_num_collision_coords; ++ii ) {
      struct xyzvector iicolcoord = collision_coords[ ii+offset ];
      iicolcoord.x_ -= rms_coords_center_of_mass.x_;
      iicolcoord.y_ -= rms_coords_center_of_mass.y_;
      iicolcoord.z_ -= rms_coords_center_of_mass.z_;
      collision_coords[ ii+offset ] = iicolcoord;
    }
    node_moments_of_inertia[ node_id ] = moment_of_inertia;
  }
}

__kernel
void
test_right_multiply_by_private(
  __global struct xyzmatrix * m1_in,
  __global struct xyzmatrix * m2_in
)
{
  struct xyzmatrix m1, m2;
  if ( get_global_id(0) == 0 ) {
    m1 = *m1_in;
    m2 = *m2_in;
    right_multiply_by_private( &m1, &m2 );
    *m1_in = m1;
  }
}

__kernel
void
test_right_multiply_by_transpose_private(
  __global struct xyzmatrix * m1_in,
  __global struct xyzmatrix * m2_in
)
{
  struct xyzmatrix m1, m2;
  if ( get_global_id(0) == 0 ) {
    m1 = *m1_in;
    m2 = *m2_in;
    right_multiply_by_transpose_private( &m1, &m2 );
    *m1_in = m1;
  }
}

__kernel
void
test_left_multiply_by_private(
  __global struct xyzmatrix * m1_in,
  __global struct xyzmatrix * m2_in
)
{
  struct xyzmatrix m1, m2;
  if ( get_global_id(0) == 0 ) {
    m1 = *m1_in;
    m2 = *m2_in;
    left_multiply_by_private( &m1, &m2 );
    *m1_in = m1;
  }
}
__kernel
void
test_left_multiply_by_transpose_private(
  __global struct xyzmatrix * m1_in,
  __global struct xyzmatrix * m2_in
)
{
  struct xyzmatrix m1, m2;
  if ( get_global_id(0) == 0 ) {
    m1 = *m1_in;
    m2 = *m2_in;
    left_multiply_by_transpose_private( &m1, &m2 );
    *m1_in = m1;
  }
}

__kernel
void
test_transpose_matrix_private(
  __global struct xyzmatrix * m
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzmatrix mprivate = *m;
    xyzmatrix_transpose_private( &mprivate );
    *m = mprivate;
  }
}

__kernel
void
test_xyzmatrix_xyzvector_multiply_private(
  __global struct xyzmatrix * m,
  __global struct xyzvector * v
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzmatrix mp = *m;
    struct xyzvector vp = *v;
    struct xyzvector vres = xyzmatrix_xyzvector_multiply_private( &mp, &vp );
    *v = vres;
  }
}

/// given an array of 9 xyz matrices,
/// sets one value in each of them to the given
/// input value, in row major order
__kernel
void
test_set_xyzmatrix_value_private(
  __global struct xyzmatrix * marray,
  float value
)
{
  if ( get_global_id(0) == 0 ) {
    int count = 0;
    for ( int ii = 0; ii < 3; ++ii ) {
      for ( int jj = 0; jj < 3; ++jj ) {
	struct xyzmatrix mprivate;
	mprivate = marray[ count ];
	set_xyzmatrix_value_private( & mprivate, ii, jj, value );
	marray[count] = mprivate;
	++count;
      }
    }
  }
}

/// Arguments: an xyzmatrix, an array of 9 floats in which to store the contents
/// of the matrix extracted in row-major order
__kernel
void
test_get_xyzmatrix_value_private(
  __global struct xyzmatrix * m,
  __global float * return_values
)
{
  if ( get_global_id(0) == 0 ) {
    int count = 0;
    struct xyzmatrix mprivate = *m;
    for ( int ii = 0; ii < 3; ++ii ) {
      for ( int jj = 0; jj < 3; ++jj ) {
	return_values[ count ] = get_xyzmatrix_value_private( & mprivate, ii, jj );
	++count;
      }
    }
  }
}

/// given an array of 3 xyz vectors,
/// sets one value in each of them to the given
/// input value
__kernel
void
test_set_xyzvector_value_private(
  __global struct xyzvector * varray,
  float value
)
{
  if ( get_global_id(0) == 0 ) {
    for ( int ii = 0; ii < 3; ++ii ) {
      struct xyzvector vprivate;
      vprivate = varray[ ii ];
      set_xyzvector_value_private( & vprivate, ii, value );
      varray[ii] = vprivate;
    }
  }
}

/// Arguments: an xyzvector, an array of 3 floats in which to store the contents
/// of the vector
__kernel
void
test_get_xyzvector_value_private(
  __global struct xyzvector * v,
  __global float * return_values
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzvector vprivate = *v;
    for ( int ii = 0; ii < 3; ++ii ) {
      return_values[ ii ] = get_xyzvector_value_private( & vprivate, ii );
    }
  }
}

/// Test that set_xyzmatrix_to_identity_private accurately creates an identity matrix
__kernel
void
test_set_xyzmatrix_to_identity_private(
  __global struct xyzmatrix * m
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzmatrix mprivate = *m;
    set_xyzmatrix_to_identity_private( &mprivate );
    *m = mprivate;
  }
}

__kernel
void
test_xyzvector_square_magnitude(
  __global struct xyzvector * v,
  __global float * square_mag
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzvector vp = *v;
    *square_mag = xyzvector_square_magnitude_private( &vp );
  }
}

__kernel
void
test_xyzvector_square_distance(
  __global struct xyzvector * v1,
  __global struct xyzvector * v2,
  __global float * square_dist
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzvector v1p = *v1;
    struct xyzvector v2p = *v2;
    *square_dist = xyzvector_square_distance_private( &v1p, &v2p );
  }
}

__kernel
void
test_xyzvector_dot_product(
  __global struct xyzvector * v1,
  __global struct xyzvector * v2,
  __global float * dot_product
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzvector v1p = *v1;
    struct xyzvector v2p = *v2;
    *dot_product = xyzvector_dot_product_private( &v1p, &v2p );
  }
}

/// Test that performing a jacobi rotation works:
/// Arguments: in_matrix -- one starting matrix
/// Arguments: out_matrices -- six output rotation matrices for
/// all combos satisfying the jacobi_rotation preconditions
__kernel
void
test_jacobi_xyzmatrix(
  __global struct xyzmatrix * in_matrix,
  __global struct xyzmatrix * out_matrices
)
{
  if ( get_global_id(0) == 0 ) {
    int count=0;
    struct xyzmatrix mprivate = *in_matrix;
    for ( int ii = 0; ii < 3; ++ii ) {
      for ( int jj = 0; jj < 3; ++jj ) {
	if ( ii == jj ) continue;
	struct xyzmatrix r;
	jacobi_rotation( &mprivate, ii, jj, &r );
	out_matrices[ count ] = r;
	++count;
      }
    }
  }
}
  
/// Test eigenvect_jacobi
__kernel
void
test_eigenvector_jacobi(
  __global struct xyzmatrix * in_matrix,
  __global struct xyzmatrix * eigenvectors,
  __global struct xyzvector * eigenvalues,
  float tolerance
)
{
  if ( get_global_id( 0 ) == 0 ) {
    struct xyzvector evals;
    struct xyzmatrix m = *in_matrix;
    struct xyzmatrix J;
    evals = eigenvector_jacobi( &m, tolerance, &J );
    *eigenvectors = J;
    *eigenvalues = evals;
  }
}

__kernel
void
test_findUU(
  __global struct xyzvector * XX,
  __global struct xyzvector * YY,
  __global float * WW,
  int Npoints,
  __global struct xyzmatrix * UUout,
  __global float * sigma3out//,
  //__global struct xyzmatrix * m_moment_out
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzmatrix UU;
    float sigma3;
    findUU( XX, YY, WW, Npoints, & UU, & sigma3 ); //, m_moment_out );
    *UUout = UU;
    *sigma3out = sigma3;
  }
}

__kernel
void
test_translate_coords_to_origin(
  __global struct xyzvector * XX,
  __global float * moment_of_inertia,
  int npoints
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzvector coords_center_of_mass;
    float moi = transform_center_of_mass_to_origin( XX, & coords_center_of_mass, npoints );
    *moment_of_inertia = moi;
  }
}

__kernel
void
test_findUU_no_translation_to_origin(
  __global struct xyzvector * XX, // <-- must have already been translated so that the center of mass is at the origin
  __global struct xyzvector * YY, // <-- must have already been translated so that the center of mass is at the origin
  int npoints,
  __global struct xyzmatrix * UUout,
  __global float * sigma3out//,
  //__global struct xyzmatrix * m_moment_out
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzmatrix UU;
    float sigma3;
    findUU_no_translation_to_origin( XX, YY, npoints, & UU, &sigma3 );//, 0 );
    *UUout = UU;
    *sigma3out = sigma3;
  }
}

__kernel
void
test_calculate_rmsd_fast_no_translation_to_origin(
  __global struct xyzvector * XX,
  __global struct xyzvector * YY,
  int npoints,
  __global float * xx_moi,
  __global float * yy_moi,
  __global float * rms_computed
)
{
  if ( get_global_id(0) == 0 ) {
    struct xyzmatrix UU;
    float sigma3;
    findUU_no_translation_to_origin( XX, YY, npoints, &UU, &sigma3 ); // , 0 );
    *rms_computed = calculate_rmsd_fast( xx_moi, yy_moi, npoints, sigma3 );
  }
}
