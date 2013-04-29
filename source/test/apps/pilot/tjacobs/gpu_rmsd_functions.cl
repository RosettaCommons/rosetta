//struct xyzmatrix;
//struct xyzvector;

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

void
sort_three_values(
  __private struct xyzvector * v, // three values held in an xyzvector in the x_, y_, and z_ positions
  __private int * sort // store the indices of the sorted values in this array
)
{

  // explicitly coded 3 level index sort using eigenvalues
  for ( int i = 0; i < 3; ++i ) {
    sort[ i ] = i;
  }

  if ( get_xyzvector_value_private( v,0) < get_xyzvector_value_private( v,1) ) {
    sort[ 1 ] = 0;
    sort[ 0 ] = 1;
  }

  if ( get_xyzvector_value_private( v, sort[ 1 ] ) < get_xyzvector_value_private( v, 2) ) {
    sort[ 2 ] = sort[ 1 ];
    sort[ 1 ] = 2;

    if ( get_xyzvector_value_private( v, sort[ 0 ] ) < get_xyzvector_value_private( v, 2) ) {
      sort[ 1 ] = sort[ 0 ];
      sort[ 0 ] = 2;
    }
  }
}

void
sort_eigenvectors_and_eigenvalues(
  __private int * sort, // indices of the eigenvalues, stored largest to smallest
  __private struct xyzvector * eVal, // eigenvalues
  __private struct xyzmatrix * eVec // eigenvectors
)
{
  // sort eigen values
  float temp1 = get_xyzvector_value_private( eVal, sort[0] );
  float temp2 = get_xyzvector_value_private( eVal, sort[1] );
  set_xyzvector_value_private( eVal, 2, get_xyzvector_value_private( eVal, sort[2] ) );
  set_xyzvector_value_private( eVal, 1, temp2 );
  set_xyzvector_value_private( eVal, 0, temp1 );
  // sort first two eigen vectors (dont care about third)
  for ( int i = 0; i < 3; ++i ) {
    temp1 = get_xyzmatrix_value_private( eVec, i, sort[0] );
    temp2 = get_xyzmatrix_value_private( eVec, i, sort[1] );
    set_xyzmatrix_value_private( eVec, i, 0, temp1 );
    set_xyzmatrix_value_private( eVec, i, 1, temp2 );
  }
}

void
make_cross_moments_matrix_unweighted(
  __global struct xyzvector * XX,
  __global struct xyzvector * YY,
  __private struct xyzmatrix * m_moment,
  int Npoints
)
{
  for ( int ii = 0; ii < 3; ++ii ) for ( int jj = 0; jj < 3; ++jj ) set_xyzmatrix_value_private( m_moment, ii, jj, 0.0 );

  for ( int jj = 0; jj < Npoints; ++jj ) {
    struct xyzvector xx = XX[jj];
    struct xyzvector yy = YY[jj];
    for ( int j = 0; j < 3; ++j ) {
      for ( int k = 0; k < 3; ++k ) {
	set_xyzmatrix_value_private( m_moment, k, j,
	  get_xyzmatrix_value_private(m_moment, k, j ) +
	  get_xyzvector_value_private( &yy, k ) *
	  get_xyzvector_value_private( &xx, j ) );
      }
    }
  }
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

  rr_moment = m_moment;
  left_multiply_by_transpose_private( & rr_moment, & m_moment );

  // Find eigenvalues, eigenvectors of symmetric matrix rr_moment
  w_w = eigenvector_jacobi( & rr_moment, 1E-9, & eVec );
  sort_three_values( &w_w, sort );
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
  
  // sort the eigenvalues, and also sort the first two eigenvectors
  sort_eigenvectors_and_eigenvalues( sort, &w_w, &eVec );

  // april 20: the fix not only fixes bad eigen vectors but solves a problem of
  // forcing a right-handed coordinate system
  fixEigenvector( &eVec);

  // at this point we now have three good eigenvectors in a right hand
  // coordinate system.

  // make bb basis vectors   = moments*eVec
  bb = eVec;
  left_multiply_by_private( &bb, &m_moment );

  //     std::cerr << "m_moment" << std::endl;
  // squirrel away a free copy of the third eigenvector before normalization/fix
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

  make_cross_moments_matrix_unweighted( XX, YY, &m_moment, Npoints );

  // Multiply CROSS MOMENTS by transpose
  rr_moment = m_moment;
  left_multiply_by_transpose_private( & rr_moment, & m_moment );

  // Find eigenvalues, eigenvectors of symmetric matrix rr_moment
  w_w = eigenvector_jacobi( & rr_moment, 1E-9, & eVec );
  //if ( m_moment_out ) { m_moment_out->xx_ = w_w.x_; m_moment_out->xy_ = w_w.y_; m_moment_out->xz_ = w_w.z_; }

  sort_three_values( &w_w, sort );
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

  // sort the eigenvalues, and also sort the first two eigenvectors
  sort_eigenvectors_and_eigenvalues( sort, &w_w, &eVec );

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
  bb = eVec;
  left_multiply_by_private( &bb, &m_moment );

  //     std::cerr << "m_moment" << std::endl;
  // squirrel away a free copy of the third eigenvector before normalization/fix
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

