//struct xyzmatrix;
//struct xyzvector;

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
