//struct xyzmatrix;
//struct xyzvector;

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
