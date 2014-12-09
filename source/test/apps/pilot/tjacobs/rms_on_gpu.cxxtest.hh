// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/apps/pilot/tjacobs/rms_on_gpu.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

#include <basic/gpu/GPU.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

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

class rms_on_gpu_Tests : public CxxTest::TestSuite
{
public:
	void setUp() {
		core_init();
	}

	std::vector< std::string >
	gpu_rmsd_test_programs() {
		std::vector< std::string > programs;
		programs.push_back( "/home/andrew/scr/GIT/rosetta/rosetta_source/test/apps/pilot/tjacobs/gpu_xyzfunctions.cl" );
		programs.push_back( "/home/andrew/scr/GIT/rosetta/rosetta_source/test/apps/pilot/tjacobs/gpu_rmsd_functions.cl" );
		programs.push_back( "/home/andrew/scr/GIT/rosetta/rosetta_source/test/apps/pilot/tjacobs/gpu_helical_bundle_rms_and_clash_calculations.cl" );
		programs.push_back( "/home/andrew/scr/GIT/rosetta/rosetta_source/test/apps/pilot/tjacobs/gpu_xyzfunctions_tests.cl" );
		programs.push_back( "/home/andrew/scr/GIT/rosetta/rosetta_source/test/apps/pilot/tjacobs/gpu_rmsd_functions_tests.cl" );
		return programs;
	}

	void test_right_multiply_by() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_right_mult_private = gpu.BuildKernel( "test_right_multiply_by_private" );
		if ( ! test_right_mult_private ) {
			std::cerr << "Failed to build test_right_multiply_by_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		xyzmatrix m2;
		m2.xx_ = 4.1; m2.xy_ = 5.1; m2.xz_ = 6.1;
		m2.yx_ = 4.2; m2.yy_ = 5.2; m2.yz_ = 6.2;
		m2.zx_ = 4.3; m2.zy_ = 5.3; m2.zz_ = 6.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix), &m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}
		cl_mem d_m2 = gpu.AllocateMemory( sizeof( xyzmatrix), &m2, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m2 ) {
			std::cerr << "Failed to allocate d_m2" << std::endl;
			TS_ASSERT( d_m2 );
			return;
		}

		//gpu.WriteData( d_m1, &m1, sizeof( xyzmatrix ), CL_TRUE );
		//gpu.WriteData( d_m2, &m2, sizeof( xyzmatrix ), CL_TRUE );

		if ( ! gpu.ExecuteKernel( "test_right_multiply_by_private", 32, 32, 32, GPU_DEVMEM, d_m1, GPU_DEVMEM, d_m2, NULL ) ) {
			std::cerr << "Failed to execute test_right_multiply_by_private" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix mresult;
		mresult.xx_ = mresult.xy_ = mresult.xz_ = -1.0;
		mresult.yx_ = mresult.yy_ = mresult.yz_ = -1.0;
		mresult.zx_ = mresult.zy_ = mresult.zz_ = -1.0;
		gpu.ReadData( &mresult, d_m1, sizeof( xyzmatrix ), CL_TRUE );


		numeric::xyzMatrix< float > nm1;
		nm1.xx( 1.1 ); nm1.xy( 2.1 ); nm1.xz( 3.1 );
		nm1.yx( 1.2 ); nm1.yy( 2.2 ); nm1.yz( 3.2 );
		nm1.zx( 1.3 ); nm1.zy( 2.3 ); nm1.zz( 3.3 );

		numeric::xyzMatrix< float > nm2;
		nm2.xx( 4.1 ); nm2.xy( 5.1 ); nm2.xz( 6.1 );
		nm2.yx( 4.2 ); nm2.yy( 5.2 ); nm2.yz( 6.2 );
		nm2.zx( 4.3 ); nm2.zy( 5.3 ); nm2.zz( 6.3 );

		nm1.right_multiply_by( nm2 );

		TS_ASSERT_DELTA( mresult.xx_, nm1.xx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.xy_, nm1.xy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.xz_, nm1.xz(), 1e-5 );

		TS_ASSERT_DELTA( mresult.yx_, nm1.yx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.yy_, nm1.yy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.yz_, nm1.yz(), 1e-5 );

		TS_ASSERT_DELTA( mresult.zx_, nm1.zx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.zy_, nm1.zy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.zz_, nm1.zz(), 1e-5 );

#endif
	}

	void test_right_multiply_by_transpose() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_right_mult_private = gpu.BuildKernel( "test_right_multiply_by_transpose_private" );
		if ( ! test_right_mult_private ) {
			std::cerr << "Failed to build test_right_multiply_by_transpose_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		xyzmatrix m2;
		m2.xx_ = 4.1; m2.xy_ = 5.1; m2.xz_ = 6.1;
		m2.yx_ = 4.2; m2.yy_ = 5.2; m2.yz_ = 6.2;
		m2.zx_ = 4.3; m2.zy_ = 5.3; m2.zz_ = 6.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix), &m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}
		cl_mem d_m2 = gpu.AllocateMemory( sizeof( xyzmatrix), &m2, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m2 ) {
			std::cerr << "Failed to allocate d_m2" << std::endl;
			TS_ASSERT( d_m2 );
			return;
		}

		//gpu.WriteData( d_m1, &m1, sizeof( xyzmatrix ), CL_TRUE );
		//gpu.WriteData( d_m2, &m2, sizeof( xyzmatrix ), CL_TRUE );

		if ( ! gpu.ExecuteKernel( "test_right_multiply_by_transpose_private", 32, 32, 32, GPU_DEVMEM, d_m1, GPU_DEVMEM, d_m2, NULL ) ) {
			std::cerr << "Failed to execute test_right_multiply_by_transpose_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix mresult;
		mresult.xx_ = mresult.xy_ = mresult.xz_ = -1.0;
		mresult.yx_ = mresult.yy_ = mresult.yz_ = -1.0;
		mresult.zx_ = mresult.zy_ = mresult.zz_ = -1.0;
		gpu.ReadData( &mresult, d_m1, sizeof( xyzmatrix ), CL_TRUE );


		numeric::xyzMatrix< float > nm1;
		nm1.xx( 1.1 ); nm1.xy( 2.1 ); nm1.xz( 3.1 );
		nm1.yx( 1.2 ); nm1.yy( 2.2 ); nm1.yz( 3.2 );
		nm1.zx( 1.3 ); nm1.zy( 2.3 ); nm1.zz( 3.3 );

		numeric::xyzMatrix< float > nm2;
		nm2.xx( 4.1 ); nm2.xy( 5.1 ); nm2.xz( 6.1 );
		nm2.yx( 4.2 ); nm2.yy( 5.2 ); nm2.yz( 6.2 );
		nm2.zx( 4.3 ); nm2.zy( 5.3 ); nm2.zz( 6.3 );

		nm1.right_multiply_by_transpose( nm2 );

		TS_ASSERT_DELTA( mresult.xx_, nm1.xx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.xy_, nm1.xy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.xz_, nm1.xz(), 1e-5 );

		TS_ASSERT_DELTA( mresult.yx_, nm1.yx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.yy_, nm1.yy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.yz_, nm1.yz(), 1e-5 );

		TS_ASSERT_DELTA( mresult.zx_, nm1.zx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.zy_, nm1.zy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.zz_, nm1.zz(), 1e-5 );

#endif
	}

	void test_left_multiply_by() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_left_mult = gpu.BuildKernel( "test_left_multiply_by_private" );
		if ( ! test_left_mult ) {
			std::cerr << "Failed to build test_left_multiply_by_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		xyzmatrix m2;
		m2.xx_ = 4.1; m2.xy_ = 5.1; m2.xz_ = 6.1;
		m2.yx_ = 4.2; m2.yy_ = 5.2; m2.yz_ = 6.2;
		m2.zx_ = 4.3; m2.zy_ = 5.3; m2.zz_ = 6.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix), &m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}
		cl_mem d_m2 = gpu.AllocateMemory( sizeof( xyzmatrix), &m2, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m2 ) {
			std::cerr << "Failed to allocate d_m2" << std::endl;
			TS_ASSERT( d_m2 );
			return;
		}

		//gpu.WriteData( d_m1, &m1, sizeof( xyzmatrix ), CL_TRUE );
		//gpu.WriteData( d_m2, &m2, sizeof( xyzmatrix ), CL_TRUE );

		if ( ! gpu.ExecuteKernel( "test_left_multiply_by_private", 32, 32, 32, GPU_DEVMEM, d_m1, GPU_DEVMEM, d_m2, NULL ) ) {
			std::cerr << "Failed to execute test_left_multiply_by_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix mresult;
		mresult.xx_ = mresult.xy_ = mresult.xz_ = -1.0;
		mresult.yx_ = mresult.yy_ = mresult.yz_ = -1.0;
		mresult.zx_ = mresult.zy_ = mresult.zz_ = -1.0;
		gpu.ReadData( &mresult, d_m1, sizeof( xyzmatrix ), CL_TRUE );


		numeric::xyzMatrix< float > nm1;
		nm1.xx( 1.1 ); nm1.xy( 2.1 ); nm1.xz( 3.1 );
		nm1.yx( 1.2 ); nm1.yy( 2.2 ); nm1.yz( 3.2 );
		nm1.zx( 1.3 ); nm1.zy( 2.3 ); nm1.zz( 3.3 );

		numeric::xyzMatrix< float > nm2;
		nm2.xx( 4.1 ); nm2.xy( 5.1 ); nm2.xz( 6.1 );
		nm2.yx( 4.2 ); nm2.yy( 5.2 ); nm2.yz( 6.2 );
		nm2.zx( 4.3 ); nm2.zy( 5.3 ); nm2.zz( 6.3 );

		nm1.left_multiply_by( nm2 );

		TS_ASSERT_DELTA( mresult.xx_, nm1.xx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.xy_, nm1.xy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.xz_, nm1.xz(), 1e-5 );

		TS_ASSERT_DELTA( mresult.yx_, nm1.yx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.yy_, nm1.yy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.yz_, nm1.yz(), 1e-5 );

		TS_ASSERT_DELTA( mresult.zx_, nm1.zx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.zy_, nm1.zy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.zz_, nm1.zz(), 1e-5 );

#endif
	}

	void test_left_multiply_by_transpose() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_left_mult_by_transpose_private = gpu.BuildKernel( "test_left_multiply_by_transpose_private" );
		if ( ! test_left_mult_by_transpose_private ) {
			std::cerr << "Failed to build test_left_multiply_by_transpose_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		//cl_kernel test_left_mult_by_transpose_local = gpu.BuildKernel( "test_left_multiply_by_transpose_local" );
		//if ( ! test_left_mult_by_transpose_local ) {
		//	std::cerr << "Failed to build test_left_multiply_by_transpose_local kernel" << std::endl;
		//	TS_ASSERT( false );
		//	return;
		//}

		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		xyzmatrix m2;
		m2.xx_ = 4.1; m2.xy_ = 5.1; m2.xz_ = 6.1;
		m2.yx_ = 4.2; m2.yy_ = 5.2; m2.yz_ = 6.2;
		m2.zx_ = 4.3; m2.zy_ = 5.3; m2.zz_ = 6.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix), &m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}
		cl_mem d_m2 = gpu.AllocateMemory( sizeof( xyzmatrix), &m2, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m2 ) {
			std::cerr << "Failed to allocate d_m2" << std::endl;
			TS_ASSERT( d_m2 );
			return;
		}

		//gpu.WriteData( d_m1, &m1, sizeof( xyzmatrix ), CL_TRUE );
		//gpu.WriteData( d_m2, &m2, sizeof( xyzmatrix ), CL_TRUE );

		if ( ! gpu.ExecuteKernel( "test_left_multiply_by_transpose_private", 32, 32, 32, GPU_DEVMEM, d_m1, GPU_DEVMEM, d_m2, NULL ) ) {
			std::cerr << "Failed to execute test_left_multiply_by_transpose_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix mresult;
		mresult.xx_ = mresult.xy_ = mresult.xz_ = -1.0;
		mresult.yx_ = mresult.yy_ = mresult.yz_ = -1.0;
		mresult.zx_ = mresult.zy_ = mresult.zz_ = -1.0;
		gpu.ReadData( &mresult, d_m1, sizeof( xyzmatrix ), CL_TRUE );


		numeric::xyzMatrix< float > nm1;
		nm1.xx( 1.1 ); nm1.xy( 2.1 ); nm1.xz( 3.1 );
		nm1.yx( 1.2 ); nm1.yy( 2.2 ); nm1.yz( 3.2 );
		nm1.zx( 1.3 ); nm1.zy( 2.3 ); nm1.zz( 3.3 );

		numeric::xyzMatrix< float > nm2;
		nm2.xx( 4.1 ); nm2.xy( 5.1 ); nm2.xz( 6.1 );
		nm2.yx( 4.2 ); nm2.yy( 5.2 ); nm2.yz( 6.2 );
		nm2.zx( 4.3 ); nm2.zy( 5.3 ); nm2.zz( 6.3 );

		nm1.left_multiply_by_transpose( nm2 );

		TS_ASSERT_DELTA( mresult.xx_, nm1.xx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.xy_, nm1.xy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.xz_, nm1.xz(), 1e-5 );

		TS_ASSERT_DELTA( mresult.yx_, nm1.yx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.yy_, nm1.yy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.yz_, nm1.yz(), 1e-5 );

		TS_ASSERT_DELTA( mresult.zx_, nm1.zx(), 1e-5 );
		TS_ASSERT_DELTA( mresult.zy_, nm1.zy(), 1e-5 );
		TS_ASSERT_DELTA( mresult.zz_, nm1.zz(), 1e-5 );

#endif
	}

	void test_transpose_matrix() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_transpose_matrix_private = gpu.BuildKernel( "test_transpose_matrix_private" );
		if ( ! test_transpose_matrix_private ) {
			std::cerr << "Failed to build test_transpose_matrix_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix), &m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_transpose_matrix_private", 32, 32, 32, GPU_DEVMEM, d_m1, NULL ) ) {
			std::cerr << "Failed to execute test_transpose_matrix_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix mresult;
		mresult.xx_ = mresult.xy_ = mresult.xz_ = -1.0;
		mresult.yx_ = mresult.yy_ = mresult.yz_ = -1.0;
		mresult.zx_ = mresult.zy_ = mresult.zz_ = -1.0;
		gpu.ReadData( &mresult, d_m1, sizeof( xyzmatrix ), CL_TRUE );


		TS_ASSERT_DELTA( mresult.xx_, m1.xx_, 1e-5 );
		TS_ASSERT_DELTA( mresult.xy_, m1.yx_, 1e-5 );
		TS_ASSERT_DELTA( mresult.xz_, m1.zx_, 1e-5 );

		TS_ASSERT_DELTA( mresult.yx_, m1.xy_, 1e-5 );
		TS_ASSERT_DELTA( mresult.yy_, m1.yy_, 1e-5 );
		TS_ASSERT_DELTA( mresult.yz_, m1.zy_, 1e-5 );

		TS_ASSERT_DELTA( mresult.zx_, m1.xz_, 1e-5 );
		TS_ASSERT_DELTA( mresult.zy_, m1.yz_, 1e-5 );
		TS_ASSERT_DELTA( mresult.zz_, m1.zz_, 1e-5 );

#endif
	}

	void test_xyzmatrix_xyzvector_multiply() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_xyzmatrix_xyzvector_multiply_private = gpu.BuildKernel( "test_xyzmatrix_xyzvector_multiply_private" );
		if ( ! test_xyzmatrix_xyzvector_multiply_private ) {
			std::cerr << "Failed to build test_xyzmatrix_xyzvector_multiply_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		xyzvector v1;
		v1.x_ = 3.4;
		v1.y_ = -2.1;
		v1.z_ = 10.2;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix ), &m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}

		cl_mem d_v1 = gpu.AllocateMemory( sizeof( xyzvector ), &v1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_v1 ) {
			std::cerr << "Failed to allocate d_v1" << std::endl;
			TS_ASSERT( d_v1 );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_xyzmatrix_xyzvector_multiply_private", 32, 32, 32, GPU_DEVMEM, d_m1, GPU_DEVMEM, d_v1, NULL ) ) {
			std::cerr << "Failed to execute test_xyzmatrix_xyzvector_multiply_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzvector vresult;
		gpu.ReadData( &vresult, d_v1, sizeof( xyzvector ), CL_TRUE );

		numeric::xyzVector< float > nv1; nv1.x() = v1.x_; nv1.y() = v1.y_; nv1.z() = v1.z_;
		numeric::xyzMatrix< float > nm1;
		nm1.xx( 1.1 ); nm1.xy( 2.1 ); nm1.xz( 3.1 );
		nm1.yx( 1.2 ); nm1.yy( 2.2 ); nm1.yz( 3.2 );
		nm1.zx( 1.3 ); nm1.zy( 2.3 ); nm1.zz( 3.3 );

		numeric::xyzVector< float > nvp = nm1 * nv1;

		TS_ASSERT_DELTA( vresult.x_, nvp.x(), 1e-5 );
		TS_ASSERT_DELTA( vresult.y_, nvp.y(), 1e-5 );
		TS_ASSERT_DELTA( vresult.z_, nvp.z(), 1e-5 );

#endif
	}

	void test_set_xyzmatrix_value() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_set_xyzmatrix_value_private = gpu.BuildKernel( "test_set_xyzmatrix_value_private" );
		if ( ! test_set_xyzmatrix_value_private ) {
			std::cerr << "Failed to build test_set_xyzmatrix_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		//cl_kernel test_left_mult_by_transpose_local = gpu.BuildKernel( "test_left_multiply_by_transpose_local" );
		//if ( ! test_left_mult_by_transpose_local ) {
		//	std::cerr << "Failed to build test_left_multiply_by_transpose_local kernel" << std::endl;
		//	TS_ASSERT( false );
		//	return;
		//}

		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		utility::vector1< xyzmatrix > marray( 9 );
		for ( core::Size ii = 1; ii <= 9; ++ii ) marray[ ii ] = m1;

		cl_mem d_marray = gpu.AllocateMemory( 9 * sizeof( xyzmatrix ), & marray[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_marray ) {
			std::cerr << "Failed to allocate d_marray" << std::endl;
			TS_ASSERT( d_marray );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_set_xyzmatrix_value_private", 32, 32, 32, GPU_DEVMEM, d_marray, GPU_FLOAT, -1234.5, NULL ) ) {
			std::cerr << "Failed to execute test_set_xyzmatrix_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< xyzmatrix > mresult( 9 );
		for ( core::Size ii = 1; ii <= 9; ++ii ) {
			mresult[ii].xx_ = mresult[ii].xy_ = mresult[ii].xz_ = -1.0;
			mresult[ii].yx_ = mresult[ii].yy_ = mresult[ii].yz_ = -1.0;
			mresult[ii].zx_ = mresult[ii].zy_ = mresult[ii].zz_ = -1.0;
		}
		gpu.ReadData( &mresult[1], d_marray, 9*sizeof( xyzmatrix ), CL_TRUE );

		for ( core::Size ii = 1; ii <= 9; ++ii ) {
			if ( ii != 1 ) {TS_ASSERT_DELTA( mresult[ ii ].xx_, m1.xx_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].xx_, -1234.5, 1e-5 );}
			if ( ii != 2 ) {TS_ASSERT_DELTA( mresult[ ii ].xy_, m1.xy_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].xy_, -1234.5, 1e-5 );}
			if ( ii != 3 ) {TS_ASSERT_DELTA( mresult[ ii ].xz_, m1.xz_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].xz_, -1234.5, 1e-5 );}
			if ( ii != 4 ) {TS_ASSERT_DELTA( mresult[ ii ].yx_, m1.yx_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].yx_, -1234.5, 1e-5 );}
			if ( ii != 5 ) {TS_ASSERT_DELTA( mresult[ ii ].yy_, m1.yy_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].yy_, -1234.5, 1e-5 );}
			if ( ii != 6 ) {TS_ASSERT_DELTA( mresult[ ii ].yz_, m1.yz_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].yz_, -1234.5, 1e-5 );}
			if ( ii != 7 ) {TS_ASSERT_DELTA( mresult[ ii ].zx_, m1.zx_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].zx_, -1234.5, 1e-5 );}
			if ( ii != 8 ) {TS_ASSERT_DELTA( mresult[ ii ].zy_, m1.zy_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].zy_, -1234.5, 1e-5 );}
			if ( ii != 9 ) {TS_ASSERT_DELTA( mresult[ ii ].zz_, m1.zz_, 1e-5 );} else {TS_ASSERT_DELTA( mresult[ ii ].zz_, -1234.5, 1e-5 );}
		}
#endif
	}

	void test_get_xyzmatrix_value() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_get_xyzmatrix_value_private = gpu.BuildKernel( "test_get_xyzmatrix_value_private" );
		if ( ! test_get_xyzmatrix_value_private ) {
			std::cerr << "Failed to build test_get_xyzmatrix_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		//cl_kernel test_left_mult_by_transpose_local = gpu.BuildKernel( "test_left_multiply_by_transpose_local" );
		//if ( ! test_left_mult_by_transpose_local ) {
		//	std::cerr << "Failed to build test_left_multiply_by_transpose_local kernel" << std::endl;
		//	TS_ASSERT( false );
		//	return;
		//}

		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix ), & m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}
		utility::vector1< float > vals_init(9,-1.0f);
		cl_mem d_vals = gpu.AllocateMemory( 9 * sizeof( float ), &vals_init[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_vals ) {
			std::cerr << "Failed to allocate d_vals" << std::endl;
			TS_ASSERT( d_vals );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_get_xyzmatrix_value_private", 32, 32, 32, GPU_DEVMEM, d_m1, GPU_DEVMEM, d_vals, NULL ) ) {
			std::cerr << "Failed to execute test_get_xyzmatrix_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< float > valresults( 9, -2.0f );
		gpu.ReadData( &valresults[1], d_vals, 9*sizeof( float ), CL_TRUE );

		for ( core::Size ii = 1; ii <= 9; ++ii ) {
			if ( ii == 1 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.xx_, 1e-5 );}
			if ( ii == 2 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.xy_, 1e-5 );}
			if ( ii == 3 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.xz_, 1e-5 );}
			if ( ii == 4 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.yx_, 1e-5 );}
			if ( ii == 5 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.yy_, 1e-5 );}
			if ( ii == 6 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.yz_, 1e-5 );}
			if ( ii == 7 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.zx_, 1e-5 );}
			if ( ii == 8 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.zy_, 1e-5 );}
			if ( ii == 9 ) {TS_ASSERT_DELTA( valresults[ ii ], m1.zz_, 1e-5 );}
		}

#endif
	}

	void test_set_xyzvector_value() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_set_xyzvector_value_private = gpu.BuildKernel( "test_set_xyzvector_value_private" );
		if ( ! test_set_xyzvector_value_private ) {
			std::cerr << "Failed to build test_set_xyzvector_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzvector v1;
		v1.x_ = 3.4;
		v1.y_ = -2.1;
		v1.z_ = 10.2;

		utility::vector1< xyzvector > varray( 3 );
		for ( core::Size ii = 1; ii <= 3; ++ii ) varray[ ii ] = v1;

		cl_mem d_varray = gpu.AllocateMemory( 3 * sizeof( xyzvector ), & varray[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_varray ) {
			std::cerr << "Failed to allocate d_varray" << std::endl;
			TS_ASSERT( d_varray );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_set_xyzvector_value_private", 32, 32, 32, GPU_DEVMEM, d_varray, GPU_FLOAT, -1234.5, NULL ) ) {
			std::cerr << "Failed to execute test_set_xyzvector_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< xyzvector > vresult( 3 );
		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			vresult[ii].x_ = vresult[ii].y_ = vresult[ii].z_ = -1.0;
		}
		gpu.ReadData( &vresult[1], d_varray, 3*sizeof( xyzvector ), CL_TRUE );

		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			if ( ii != 1 ) {TS_ASSERT_DELTA( vresult[ ii ].x_, v1.x_, 1e-5 );} else {TS_ASSERT_DELTA( vresult[ ii ].x_, -1234.5, 1e-5 );}
			if ( ii != 2 ) {TS_ASSERT_DELTA( vresult[ ii ].y_, v1.y_, 1e-5 );} else {TS_ASSERT_DELTA( vresult[ ii ].y_, -1234.5, 1e-5 );}
			if ( ii != 3 ) {TS_ASSERT_DELTA( vresult[ ii ].z_, v1.z_, 1e-5 );} else {TS_ASSERT_DELTA( vresult[ ii ].z_, -1234.5, 1e-5 );}
		}
#endif
	}

	void test_get_xyzvector_value() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_get_xyzvector_value_private = gpu.BuildKernel( "test_get_xyzvector_value_private" );
		if ( ! test_get_xyzvector_value_private ) {
			std::cerr << "Failed to build test_get_xyzvector_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzvector v1;
		v1.x_ = 3.4;
		v1.y_ = -2.1;
		v1.z_ = 10.2;

		cl_mem d_v1 = gpu.AllocateMemory( sizeof( xyzvector ), & v1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_v1 ) {
			std::cerr << "Failed to allocate d_v1" << std::endl;
			TS_ASSERT( d_v1 );
			return;
		}

		utility::vector1< float > vals_init(3,-1.0f);
		cl_mem d_vals = gpu.AllocateMemory( 3 * sizeof( float ), &vals_init[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_vals ) {
			std::cerr << "Failed to allocate d_vals" << std::endl;
			TS_ASSERT( d_vals );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_get_xyzvector_value_private", 32, 32, 32, GPU_DEVMEM, d_v1, GPU_DEVMEM, d_vals, NULL ) ) {
			std::cerr << "Failed to execute test_get_xyzvector_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< float > valresults( 3, -2.0f );
		gpu.ReadData( &valresults[1], d_vals, 3*sizeof( float ), CL_TRUE );

		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			if ( ii == 1 ) {TS_ASSERT_DELTA( valresults[ ii ], v1.x_, 1e-5 );}
			if ( ii == 2 ) {TS_ASSERT_DELTA( valresults[ ii ], v1.y_, 1e-5 );}
			if ( ii == 3 ) {TS_ASSERT_DELTA( valresults[ ii ], v1.z_, 1e-5 );}
		}

#endif
	}

	void test_xyzvector_square_magnitude() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_xyzvector_square_magnitude = gpu.BuildKernel( "test_xyzvector_square_magnitude" );
		if ( ! test_xyzvector_square_magnitude ) {
			std::cerr << "Failed to build test_xyzvector_square_magnitude kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzvector v1;
		v1.x_ = 3.4;
		v1.y_ = -2.1;
		v1.z_ = 10.2;

		cl_mem d_v1 = gpu.AllocateMemory( sizeof( xyzvector ), & v1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_v1 ) {
			std::cerr << "Failed to allocate d_v1" << std::endl;
			TS_ASSERT( d_v1 );
			return;
		}

		cl_mem d_val = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_val ) {
			std::cerr << "Failed to allocate d_val" << std::endl;
			TS_ASSERT( d_val );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_xyzvector_square_magnitude", 32, 32, 32, GPU_DEVMEM, d_v1, GPU_DEVMEM, d_val, NULL ) ) {
			std::cerr << "Failed to execute test_get_xyzvector_value_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		float valresult;
		gpu.ReadData( &valresult, d_val, sizeof( float ), CL_TRUE );

		numeric::xyzVector< float > nv1; nv1.x() = v1.x_; nv1.y() = v1.y_; nv1.z() = v1.z_;
		TS_ASSERT_DELTA( valresult, nv1.length_squared(), 1e-5 );

#endif
	}

	void test_xyzvector_square_distance() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_xyzvector_square_distance  = gpu.BuildKernel( "test_xyzvector_square_distance" );
		if ( ! test_xyzvector_square_distance ) {
			std::cerr << "Failed to build test_xyzvector_square_distance kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzvector v1;
		v1.x_ = 3.4;
		v1.y_ = -2.1;
		v1.z_ = 10.2;

		xyzvector v2;
		v2.x_ = 2.4;
		v2.y_ = -5.1;
		v2.z_ = -0.2;

		cl_mem d_v1 = gpu.AllocateMemory( sizeof( xyzvector ), & v1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_v1 ) {
			std::cerr << "Failed to allocate d_v1" << std::endl;
			TS_ASSERT( d_v1 );
			return;
		}

		cl_mem d_v2 = gpu.AllocateMemory( sizeof( xyzvector ), & v2, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_v2 ) {
			std::cerr << "Failed to allocate d_v2" << std::endl;
			TS_ASSERT( d_v2 );
			return;
		}

		cl_mem d_val = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_val ) {
			std::cerr << "Failed to allocate d_val" << std::endl;
			TS_ASSERT( d_val );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_xyzvector_square_distance", 32, 32, 32, GPU_DEVMEM, d_v1, GPU_DEVMEM, d_v2, GPU_DEVMEM, d_val, NULL ) ) {
			std::cerr << "Failed to execute test_xyzvector_square_distance kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		float valresult;
		gpu.ReadData( &valresult, d_val, sizeof( float ), CL_TRUE );

		numeric::xyzVector< float > nv1; nv1.x() = v1.x_; nv1.y() = v1.y_; nv1.z() = v1.z_;
		numeric::xyzVector< float > nv2; nv2.x() = v2.x_; nv2.y() = v2.y_; nv2.z() = v2.z_;
		TS_ASSERT_DELTA( valresult, nv1.distance_squared(nv2), 1e-5 );

#endif
	}

	void test_xyzvector_dot_product() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_xyzvector_dot_product  = gpu.BuildKernel( "test_xyzvector_dot_product" );
		if ( ! test_xyzvector_dot_product ) {
			std::cerr << "Failed to build test_xyzvector_dot_product kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzvector v1;
		v1.x_ = 3.4;
		v1.y_ = -2.1;
		v1.z_ = 10.2;

		xyzvector v2;
		v2.x_ = 2.4;
		v2.y_ = -5.1;
		v2.z_ = -0.2;

		cl_mem d_v1 = gpu.AllocateMemory( sizeof( xyzvector ), & v1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_v1 ) {
			std::cerr << "Failed to allocate d_v1" << std::endl;
			TS_ASSERT( d_v1 );
			return;
		}

		cl_mem d_v2 = gpu.AllocateMemory( sizeof( xyzvector ), & v2, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_v2 ) {
			std::cerr << "Failed to allocate d_v2" << std::endl;
			TS_ASSERT( d_v2 );
			return;
		}

		cl_mem d_val = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_val ) {
			std::cerr << "Failed to allocate d_val" << std::endl;
			TS_ASSERT( d_val );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_xyzvector_dot_product", 32, 32, 32, GPU_DEVMEM, d_v1, GPU_DEVMEM, d_v2, GPU_DEVMEM, d_val, NULL ) ) {
			std::cerr << "Failed to execute test_xyzvector_dot_product kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		float valresult;
		gpu.ReadData( &valresult, d_val, sizeof( float ), CL_TRUE );

		numeric::xyzVector< float > nv1; nv1.x() = v1.x_; nv1.y() = v1.y_; nv1.z() = v1.z_;
		numeric::xyzVector< float > nv2; nv2.x() = v2.x_; nv2.y() = v2.y_; nv2.z() = v2.z_;
		TS_ASSERT_DELTA( valresult, nv1.dot(nv2), 1e-5 );

#endif
	}

	void test_set_xyzmatrix_to_identity() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_set_xyzmatrix_to_identity_private = gpu.BuildKernel( "test_set_xyzmatrix_to_identity_private" );
		if ( ! test_set_xyzmatrix_to_identity_private ) {
			std::cerr << "Failed to build test_set_xyzmatrix_to_identity_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		//cl_kernel test_left_mult_by_transpose_local = gpu.BuildKernel( "test_left_multiply_by_transpose_local" );
		//if ( ! test_left_mult_by_transpose_local ) {
		//	std::cerr << "Failed to build test_left_multiply_by_transpose_local kernel" << std::endl;
		//	TS_ASSERT( false );
		//	return;
		//}

		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix ), & m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_set_xyzmatrix_to_identity_private", 32, 32, 32, GPU_DEVMEM, d_m1, NULL ) ) {
			std::cerr << "Failed to execute test_set_xyzmatrix_to_identity_private kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix mresult;
		gpu.ReadData( &mresult, d_m1, sizeof( xyzmatrix ), CL_TRUE );

		TS_ASSERT_DELTA( mresult.xx_, 1.0f, 1e-5 );
		TS_ASSERT_DELTA( mresult.xy_, 0.0f, 1e-5 );
		TS_ASSERT_DELTA( mresult.xz_, 0.0f, 1e-5 );
		TS_ASSERT_DELTA( mresult.yx_, 0.0f, 1e-5 );
		TS_ASSERT_DELTA( mresult.yy_, 1.0f, 1e-5 );
		TS_ASSERT_DELTA( mresult.yz_, 0.0f, 1e-5 );
		TS_ASSERT_DELTA( mresult.zx_, 0.0f, 1e-5 );
		TS_ASSERT_DELTA( mresult.zy_, 0.0f, 1e-5 );
		TS_ASSERT_DELTA( mresult.zz_, 1.0f, 1e-5 );

#endif
	}

	void test_jacobi_rotation() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_jacobi_xyzmatrix = gpu.BuildKernel( "test_jacobi_xyzmatrix" );
		if ( ! test_jacobi_xyzmatrix ) {
			std::cerr << "Failed to build test_jacobi_xyzmatrix kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		//cl_kernel test_left_mult_by_transpose_local = gpu.BuildKernel( "test_left_multiply_by_transpose_local" );
		//if ( ! test_left_mult_by_transpose_local ) {
		//	std::cerr << "Failed to build test_left_multiply_by_transpose_local kernel" << std::endl;
		//	TS_ASSERT( false );
		//	return;
		//}

		xyzmatrix m1;
		m1.xx_ = 1.1; m1.xy_ = 2.1; m1.xz_ = 3.1;
		m1.yx_ = 1.2; m1.yy_ = 2.2; m1.yz_ = 3.2;
		m1.zx_ = 1.3; m1.zy_ = 2.3; m1.zz_ = 3.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix ), & m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}

		utility::vector1< xyzmatrix > marray( 6 );
		for ( core::Size ii = 1; ii <= 6; ++ii ) marray[ ii ] = m1;

		cl_mem d_marray = gpu.AllocateMemory( 6 * sizeof( xyzmatrix ), & marray[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_marray ) {
			std::cerr << "Failed to allocate d_marray" << std::endl;
			TS_ASSERT( d_marray );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_jacobi_xyzmatrix", 32, 32, 32, GPU_DEVMEM, d_m1, GPU_DEVMEM, d_marray, NULL ) ) {
			std::cerr << "Failed to execute test_jacobi_xyzmatrix kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< xyzmatrix > mresult( 6 );
		for ( core::Size ii = 1; ii <= 6; ++ii ) {
			mresult[ii].xx_ = mresult[ii].xy_ = mresult[ii].xz_ = -1.0;
			mresult[ii].yx_ = mresult[ii].yy_ = mresult[ii].yz_ = -1.0;
			mresult[ii].zx_ = mresult[ii].zy_ = mresult[ii].zz_ = -1.0;
		}
		gpu.ReadData( &mresult[1], d_marray, 6*sizeof( xyzmatrix ), CL_TRUE );

		core::Size count = 1;
		numeric::xyzMatrix< float > m, r;
		m.xx() = m1.xx_; m.xy() = m1.xy_; m.xz() = m1.xz_;
		m.yx() = m1.yx_; m.yy() = m1.yy_; m.yz() = m1.yz_;
		m.zx() = m1.zx_; m.zy() = m1.zy_; m.zz() = m1.zz_;

		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			for ( core::Size jj = 1; jj <= 3; ++jj ) {
				if ( ii == jj ) continue;
				numeric::jacobi_rotation( m, ii, jj, r );
				TS_ASSERT_DELTA( mresult[ count ].xx_, r.xx(), 1e-5 );
				TS_ASSERT_DELTA( mresult[ count ].xy_, r.xy(), 1e-5 );
				TS_ASSERT_DELTA( mresult[ count ].xz_, r.xz(), 1e-5 );
				TS_ASSERT_DELTA( mresult[ count ].yx_, r.yx(), 1e-5 );
				TS_ASSERT_DELTA( mresult[ count ].yy_, r.yy(), 1e-5 );
				TS_ASSERT_DELTA( mresult[ count ].yz_, r.yz(), 1e-5 );
				TS_ASSERT_DELTA( mresult[ count ].zx_, r.zx(), 1e-5 );
				TS_ASSERT_DELTA( mresult[ count ].zy_, r.zy(), 1e-5 );
				TS_ASSERT_DELTA( mresult[ count ].zz_, r.zz(), 1e-5 );
				++count;
			}
		}
#endif
	}

	void test_eigenvector_jacobi() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_eigenvector_jacobi = gpu.BuildKernel( "test_eigenvector_jacobi" );
		if ( ! test_eigenvector_jacobi ) {
			std::cerr << "Failed to build test_eigenvector_jacobi kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		//cl_kernel test_left_mult_by_transpose_local = gpu.BuildKernel( "test_left_multiply_by_transpose_local" );
		//if ( ! test_left_mult_by_transpose_local ) {
		//	std::cerr << "Failed to build test_left_multiply_by_transpose_local kernel" << std::endl;
		//	TS_ASSERT( false );
		//	return;
		//}

		// Need a symmetric matrix to give to the eigenvector_jacobi function
		xyzmatrix m1;
		m1.xx_ =  1.1; m1.xy_ = -2.1; m1.xz_ =  3.1;
		m1.yx_ = -2.1; m1.yy_ =  2.2; m1.yz_ =  3.3;
		m1.zx_ =  3.1; m1.zy_ =  3.3; m1.zz_ = -3.3;

		cl_mem d_m1 = gpu.AllocateMemory( sizeof( xyzmatrix ), & m1, CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_m1 ) {
			std::cerr << "Failed to allocate d_m1" << std::endl;
			TS_ASSERT( d_m1 );
			return;
		}

		cl_mem d_eigenvectors = gpu.AllocateMemory( sizeof( xyzmatrix ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_eigenvectors ) {
			std::cerr << "Failed to allocate d_eigenvectors" << std::endl;
			TS_ASSERT( d_eigenvectors );
			return;
		}

		cl_mem d_eigenvalues = gpu.AllocateMemory( sizeof( xyzvector ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_eigenvalues ) {
			std::cerr << "Failed to allocate d_eigenvalues" << std::endl;
			TS_ASSERT( d_eigenvalues );
			return;
		}

		if ( ! gpu.ExecuteKernel( "test_eigenvector_jacobi", 32, 32, 32, GPU_DEVMEM, d_m1, GPU_DEVMEM, d_eigenvectors, GPU_DEVMEM, d_eigenvalues, GPU_FLOAT, 0.01, NULL ) ) {
			std::cerr << "Failed to execute test_eigenvector_jacobi kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix result_eigenvectors;
		result_eigenvectors.xx_ = result_eigenvectors.xy_ = result_eigenvectors.xz_ = -1.0;
		result_eigenvectors.yx_ = result_eigenvectors.yy_ = result_eigenvectors.yz_ = -1.0;
		result_eigenvectors.zx_ = result_eigenvectors.zy_ = result_eigenvectors.zz_ = -1.0;

		gpu.ReadData( &result_eigenvectors, d_eigenvectors, sizeof( xyzmatrix ), CL_TRUE );

		xyzvector result_eigenvalues;
		gpu.ReadData( &result_eigenvalues, d_eigenvalues, sizeof( xyzvector ), CL_TRUE );

		numeric::xyzMatrix< float > m, J;
		m.xx() = m1.xx_; m.xy() = m1.xy_; m.xz() = m1.xz_;
		m.yx() = m1.yx_; m.yy() = m1.yy_; m.yz() = m1.yz_;
		m.zx() = m1.zx_; m.zy() = m1.zy_; m.zz() = m1.zz_;

		numeric::xyzVector< float > evals = numeric::eigenvector_jacobi( m, 0.01f, J );
		TS_ASSERT_DELTA( result_eigenvectors.xx_, J.xx(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvectors.xy_, J.xy(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvectors.xz_, J.xz(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvectors.yx_, J.yx(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvectors.yy_, J.yy(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvectors.yz_, J.yz(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvectors.zx_, J.zx(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvectors.zy_, J.zy(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvectors.zz_, J.zz(), 1e-5 );

		TS_ASSERT_DELTA( result_eigenvalues.x_, evals.x(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvalues.y_, evals.y(), 1e-5 );
		TS_ASSERT_DELTA( result_eigenvalues.z_, evals.z(), 1e-5 );
#endif
	}

	void my_findUU(
		ObjexxFCL::FArray2< numeric::Real > & XX,
		ObjexxFCL::FArray2< numeric::Real > & YY,
		ObjexxFCL::FArray1< numeric::Real > const & WW,
		int Npoints,
		ObjexxFCL::FArray2< numeric::Real > & UU,
		numeric::Real & sigma3
	)
	{
		using namespace numeric::model_quality;
		using namespace ObjexxFCL;
		using numeric::xyzMatrix;
		using numeric::xyzVector;

		ObjexxFCL::FArray1D_int sort( 3 );
		FArray2D< numeric::Real > eVec( 3, 3 );
		FArray2D< numeric::Real > bb( 3, 3 );
		FArray1D< numeric::Real > w_w( 3 );
		FArray2D< numeric::Real > m_moment( 3, 3 );
		FArray2D< numeric::Real > rr_moment( 3, 3 );
		numeric::Real temp1;
		numeric::Real temp2;
		numeric::Real temp3;
		FArray1D< numeric::Real > Ra( 3 );

		if ( Npoints < 1 ) {

			// return identity rotation matrix to moron
			for ( int i = 1; i <= 3; ++i ) {
				for ( int k = 1; k <= 3; ++k ) {
					UU(i,k) = 0.0;
					if ( i == k ) UU(i,i) = 1.0;
				}
			}
			sigma3 = 0.0;
			return;
		}

		// align center of mass to origin
		for ( int k = 1; k <= 3; ++k ) {
			temp1 = 0.0;
			temp2 = 0.0;
			temp3 = 0.0;
			for ( int j = 1; j <= Npoints; ++j ) {
				temp1 += XX(k,j) * WW(j);
				temp2 += YY(k,j) * WW(j);
				temp3 += WW(j);
			}
			if (temp3 > 0.001) temp1 /= temp3;
			if (temp3 > 0.001) temp2 /= temp3;

			for ( int j = 1; j <= Npoints; ++j ) {
				XX(k,j) -= temp1;
				YY(k,j) -= temp2;
			}
		}

		// Make cross moments matrix   INCLUDE THE WEIGHTS HERE
		for ( int k = 1; k <= 3; ++k ) {
			for ( int j = 1; j <= 3; ++j ) {
				temp1 = 0.0;
				for ( int i = 1; i <= Npoints; ++i ) {
					temp1 += WW(i) * YY(k,i) * XX(j,i);
				}
				m_moment(k,j) = temp1;
			}
		}

		//std::cout << "m_moment:\n";
		//std::cout << m_moment(1,1) << " " << m_moment(1,2) << " " << m_moment(1,3) << "\n";
		//std::cout << m_moment(2,1) << " " << m_moment(2,2) << " " << m_moment(2,3) << "\n";
		//std::cout << m_moment(3,1) << " " << m_moment(3,2) << " " << m_moment(3,3) << std::endl;

		// Multiply CROSS MOMENTS by transpose
		BlankMatrixMult(m_moment,3,3,1,m_moment,3,0,rr_moment);

		//std::cout << "rr_moment:\n";
		//std::cout << rr_moment(1,1) << " " << rr_moment(1,2) << " " << rr_moment(1,3) << "\n";
		//std::cout << rr_moment(2,1) << " " << rr_moment(2,2) << " " << rr_moment(2,3) << "\n";
		//std::cout << rr_moment(3,1) << " " << rr_moment(3,2) << " " << rr_moment(3,3) << std::endl;

		// Copy to/from xyzMatrix/xyzVector since rest of functions use FArrays
		xyzMatrix< numeric::Real > xyz_rr_moment( xyzMatrix< numeric::Real >::cols( &rr_moment( 1,1 ) ) );
		xyzVector< numeric::Real > xyz_w_w;
		xyzMatrix< numeric::Real > xyz_eVec;

		// Find eigenvalues, eigenvectors of symmetric matrix rr_moment
		xyz_w_w = eigenvector_jacobi( xyz_rr_moment, (numeric::Real) 1E-9, xyz_eVec );

		// Copy eigenvalues/vectors back to FArray
		for ( int i = 1; i <= 3; ++i ) {
			w_w( i ) = xyz_w_w( i );
			for ( int j = 1; j <= 3; ++j ) {
				eVec( i, j ) = xyz_eVec( i, j );
			}
		}

		std::cout << "w_w:\n";
		std::cout << w_w(1) << " " << w_w(2) << " " << w_w(3) << std::endl;

		std::cout << "eVec unsorted:\n";
		std::cout << eVec(1,1) << " " << eVec(1,2) << " " << eVec(1,3) << "\n";
		std::cout << eVec(2,1) << " " << eVec(2,2) << " " << eVec(2,3) << "\n";
		std::cout << eVec(3,1) << " " << eVec(3,2) << " " << eVec(3,3) << std::endl;


		// explicitly coded 3 level index sort using eigenvalues
		for ( int i = 1; i <= 3; ++i ) {
			sort(i) = i;
		}

		if ( w_w(1) < w_w(2) ) {
			sort(2) = 1;
			sort(1) = 2;
		}

		if ( w_w(sort(2)) < w_w(3) ) {
			sort(3) = sort(2);
			sort(2) = 3;

			if ( w_w(sort(1)) < w_w(3) ) {
				sort(2) = sort(1);
				sort(1) = 3;
			}
		}

		std::cout << "sort:\n";
		std::cout << sort(1) << " " << sort( 2 ) << " " << sort(3) << std::endl;

		// sort is now an index to order of eigen values

		if ( w_w(sort(2)) == 0.0 ) { // holy smokes, two eigen values are zeros
			// return identity rotation matrix to moron
			for ( int i = 1; i <= 3; ++i ) {
				for ( int k = 1; k <= 3; ++k ) {
					UU(i,k) = 0.0;
				}
				UU(i,i) = 1.0;
			}
			if ( w_w(sort(1)) < 0.0 ) {
				w_w(sort(1)) = std::abs(w_w(sort(1)));
			}
			sigma3 = std::sqrt(w_w(sort(1)));

			return; // make like a prom dress and slip off
		}

		// sort eigen values
		temp1 = w_w(sort(1));
		temp2 = w_w(sort(2));
		w_w(3) = w_w(sort(3));
		w_w(2) = temp2;
		w_w(1) = temp1;
		// sort first two eigen vectors (dont care about third)
		for ( int i = 1; i <= 3; ++i ) {
			temp1 = eVec(i,sort(1));
			temp2 = eVec(i,sort(2));
			eVec(i,1) = temp1;
			eVec(i,2) = temp2;
		}

		std::cout << "sorted eVec:\n";
		std::cout << eVec(1,1) << " " << eVec(1,2) << " " << eVec(1,3) << "\n";
		std::cout << eVec(2,1) << " " << eVec(2,2) << " " << eVec(2,3) << "\n";
		std::cout << eVec(3,1) << " " << eVec(3,2) << " " << eVec(3,3) << std::endl;


		// april 20: the fix not only fixes bad eigen vectors but solves a problem of
		// forcing a right-handed coordinate system

		fixEigenvector(eVec);

		std::cout << "fixed eVec:\n";
		std::cout << eVec(1,1) << " " << eVec(1,2) << " " << eVec(1,3) << "\n";
		std::cout << eVec(2,1) << " " << eVec(2,2) << " " << eVec(2,3) << "\n";
		std::cout << eVec(3,1) << " " << eVec(3,2) << " " << eVec(3,3) << std::endl;


		// at this point we now have three good eigenvectors in a right hand
		// coordinate system.

		// make bb basis vectors   = moments*eVec

		BlankMatrixMult(m_moment,3,3,0,eVec,3,0,bb);
		//     std::cerr << "m_moment" << std::endl;
		// squirrel away a free copy of the third eigenvector before normalization/fix
		for ( int j = 1; j <= 3; ++j ) {
			Ra(j) = bb(j,3);
		}

		std::cout << "Ra: " << Ra(1) << " " << Ra(2) << " " << Ra(3) << std::endl;

		std::cout << "bb:\n";
		std::cout << bb(1,1) << " " << bb(1,2) << " " << bb(1,3) << "\n";
		std::cout << bb(2,1) << " " << bb(2,2) << " " << bb(2,3) << "\n";
		std::cout << bb(3,1) << " " << bb(3,2) << " " << bb(3,3) << std::endl;


		// normalize first two bb-basis vectors
		// dont care about third since were going to replace it with b1xb2
		// this also avoids problem of possible zero third eigen value
		for ( int j = 1; j <= 2; ++j ) {
			temp1 = 1.0/std::sqrt(w_w(j)); // zero checked for above
			for ( int k = 1; k <= 3; ++k ) { // x,y,z
				bb(k,j) *= temp1;
			}
		}

		std::cout << "normalize bb:\n";
		std::cout << bb(1,1) << " " << bb(1,2) << " " << bb(1,3) << "\n";
		std::cout << bb(2,1) << " " << bb(2,2) << " " << bb(2,3) << "\n";
		std::cout << bb(3,1) << " " << bb(3,2) << " " << bb(3,3) << std::endl;


		//  fix things so that bb eigenvecs are right handed

		fixEigenvector(bb); // need to fix this one too

		std::cout << "fixed bb:\n";
		std::cout << bb(1,1) << " " << bb(1,2) << " " << bb(1,3) << "\n";
		std::cout << bb(2,1) << " " << bb(2,2) << " " << bb(2,3) << "\n";
		std::cout << bb(3,1) << " " << bb(3,2) << " " << bb(3,3) << std::endl;

		// find  product of eVec and bb matrices

		BlankMatrixMult(eVec,3,3,0,bb,3,1,UU);
		// result is returned in UU.

		std::cout << "UU:\n";
		std::cout << UU(1,1) << " " << UU(1,2) << " " << UU(1,3) << "\n";
		std::cout << UU(2,1) << " " << UU(2,2) << " " << UU(2,3) << "\n";
		std::cout << UU(3,1) << " " << UU(3,2) << " " << UU(3,3) << std::endl;

		// and lastly determine a value used in another function to compute the rms
		sigma3 = 0.0;
		for ( int j = 1; j <= 3; ++j ) {
			sigma3 += bb(j,3)*Ra(j);
		}
		std::cout << "sigma3 before fix: " << sigma3 << std::endl;

		//cems the abs() fixes some round off error situations where the w_w values are
		//cems very small and accidentally negative.  (theoretically they are positive,
		//cems but in practice round off error makes them negative)
		if ( sigma3 < 0.0 ) {
			sigma3 = std::sqrt(std::abs(w_w(1))) + std::sqrt(std::abs(w_w(2))) -
				std::sqrt(std::abs(w_w(3)));
		} else {
			sigma3 = std::sqrt(std::abs(w_w(1))) + std::sqrt(std::abs(w_w(2))) +
				std::sqrt(std::abs(w_w(3)));
		}
	}



	void test_findUU() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_findUU = gpu.BuildKernel( "test_findUU" );
		if ( ! test_findUU ) {
			std::cerr << "Failed to build test_findUU kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		// Need a bunch of coordinates to give to findUU
		utility::vector1< xyzvector > xx( 4 );
		utility::vector1< xyzvector > yy( 4 );
		xx[1].x_ = 1.0; xx[1].y_ = 0.0; xx[1].z_ = 0.0;
		xx[2].x_ = 0.0; xx[2].y_ = 1.3; xx[2].z_ = 0.0;
		xx[3].x_ = 0.0; xx[3].y_ = 0.0; xx[3].z_ = 1.0;
		xx[4].x_ = 2.0; xx[4].y_ = 2.0; xx[4].z_ = 2.0;

		yy[1].x_ = 1.0; yy[1].y_ = 0.0; yy[1].z_ = 0.0;
		yy[2].x_ = 0.0; yy[2].y_ = 0.0; yy[2].z_ = 1.4;
		yy[3].x_ = 0.0; yy[3].y_ =-1.0; yy[3].z_ = 0.0;
		yy[4].x_ = 2.1; yy[4].y_ =-1.9; yy[4].z_ = 2.2;

		utility::vector1< float > ww( 4, 1.0 );

		cl_mem d_xx = gpu.AllocateMemory( sizeof( xyzvector ) * 4, & xx[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_xx ) {
			std::cerr << "Failed to allocate d_xx" << std::endl;
			TS_ASSERT( d_xx );
			return;
		}

		cl_mem d_yy = gpu.AllocateMemory( sizeof( xyzvector ) * 4, & yy[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_yy ) {
			std::cerr << "Failed to allocate d_yy" << std::endl;
			TS_ASSERT( d_yy );
			return;
		}

		cl_mem d_ww = gpu.AllocateMemory( sizeof( float ) * 4, & ww[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_ww ) {
			std::cerr << "Failed to allocate d_ww" << std::endl;
			TS_ASSERT( d_ww );
			return;
		}

		cl_mem d_uu = gpu.AllocateMemory( sizeof( xyzmatrix ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_uu ) {
			std::cerr << "Failed to allocate d_uu" << std::endl;
			TS_ASSERT( d_uu );
			return;
		}

		cl_mem d_sigma3 = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_sigma3 ) {
			std::cerr << "Failed to allocate d_sigma3" << std::endl;
			TS_ASSERT( d_sigma3 );
			return;
		}

		cl_mem d_m_moment = gpu.AllocateMemory( sizeof( xyzmatrix ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_m_moment ) {
			std::cerr << "Failed to allocate d_m_moment" << std::endl;
			TS_ASSERT( d_m_moment );
			return;
		}

		bool kernel_success = gpu.ExecuteKernel( "test_findUU", 32, 32, 32,
			GPU_DEVMEM, d_xx,
			GPU_DEVMEM, d_yy,
			GPU_DEVMEM, d_ww,
			GPU_INT, 4,
			GPU_DEVMEM, d_uu,
			GPU_DEVMEM, d_sigma3,
			//GPU_DEVMEM, d_m_moment,
			NULL );

		if ( ! kernel_success ) {
			std::cerr << "Failed to execute test_findUU kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		xyzmatrix result_uu;
		result_uu.xx_ = result_uu.xy_ = result_uu.xz_ = -1.0;
		result_uu.yx_ = result_uu.yy_ = result_uu.yz_ = -1.0;
		result_uu.zx_ = result_uu.zy_ = result_uu.zz_ = -1.0;

		gpu.ReadData( &result_uu, d_uu, sizeof( xyzmatrix ), CL_TRUE );

		float result_sigma3;
		gpu.ReadData( &result_sigma3, d_sigma3, sizeof( float ), CL_TRUE );

		//struct xyzmatrix result_m_moment;
		//gpu.ReadData( &result_m_moment, d_m_moment, sizeof( xyzmatrix ), CL_TRUE );
		//std::cout << "m_moment from GPU\n";
		//std::cout << result_m_moment.xx_ << " " << result_m_moment.xy_ << " " << result_m_moment.xz_ << "\n";
		//std::cout << result_m_moment.yx_ << " " << result_m_moment.yy_ << " " << result_m_moment.yz_ << "\n";
		//std::cout << result_m_moment.zx_ << " " << result_m_moment.zy_ << " " << result_m_moment.zz_ << std::endl;;

		ObjexxFCL::FArray2D< core::Real > XX( 3, 4 ), YY( 3, 4 );
		ObjexxFCL::FArray1D< core::Real > WW( 4, 1.0f );
		ObjexxFCL::FArray2D< core::Real > UU( 3, 3 );
		for ( core::Size ii = 1; ii <= 4; ++ii ) {
			XX(1,ii) = xx[ii].x_;
			XX(2,ii) = xx[ii].y_;
			XX(3,ii) = xx[ii].z_;
			YY(1,ii) = yy[ii].x_;
			YY(2,ii) = yy[ii].y_;
			YY(3,ii) = yy[ii].z_;
		}

		core::Real sigma3;
		my_findUU( XX, YY, WW, 4, UU, sigma3 );

		TS_ASSERT_DELTA( result_uu.xx_, UU(1,1), 1e-5 );
		TS_ASSERT_DELTA( result_uu.xy_, UU(1,2), 1e-5 );
		TS_ASSERT_DELTA( result_uu.xz_, UU(1,3), 1e-5 );
		TS_ASSERT_DELTA( result_uu.yx_, UU(2,1), 1e-5 );
		TS_ASSERT_DELTA( result_uu.yy_, UU(2,2), 1e-5 );
		TS_ASSERT_DELTA( result_uu.yz_, UU(2,3), 1e-5 );
		TS_ASSERT_DELTA( result_uu.zx_, UU(3,1), 1e-5 );
		TS_ASSERT_DELTA( result_uu.zy_, UU(3,2), 1e-5 );
		TS_ASSERT_DELTA( result_uu.zz_, UU(3,3), 1e-5 );

		TS_ASSERT_DELTA( result_sigma3, sigma3, 1e-5 );

#endif
	}


	void test_transform_center_of_mass_to_origin_and_findUU_no_translation_to_origin() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_translate_coords_to_origin = gpu.BuildKernel( "test_translate_coords_to_origin" );
		if ( ! test_translate_coords_to_origin ) {
			std::cerr << "Failed to build test_translate_coords_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_findUU_no_translation_to_origin = gpu.BuildKernel( "test_findUU_no_translation_to_origin" );
		if ( ! test_findUU_no_translation_to_origin ) {
			std::cerr << "Failed to build test_findUU_no_translation_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		// Need a bunch of coordinates to give to findUU
		utility::vector1< xyzvector > xx( 4 );
		utility::vector1< xyzvector > yy( 4 );
		xx[1].x_ = 1.0; xx[1].y_ = 0.0; xx[1].z_ = 0.0;
		xx[2].x_ = 0.0; xx[2].y_ = 1.0; xx[2].z_ = 0.0;
		xx[3].x_ = 0.0; xx[3].y_ = 0.0; xx[3].z_ = 1.0;
		xx[4].x_ = 2.0; xx[4].y_ = 2.0; xx[4].z_ = 2.0;

		yy[1].x_ = 1.0; yy[1].y_ = 0.0; yy[1].z_ = 0.0;
		yy[2].x_ = 0.0; yy[2].y_ = 0.0; yy[2].z_ = 1.0;
		yy[3].x_ = 0.0; yy[3].y_ =-1.0; yy[3].z_ = 0.0;
		yy[4].x_ = 2.1; yy[4].y_ =-1.9; yy[4].z_ = 2.2;

		utility::vector1< float > ww( 4, 1.0 );

		cl_mem d_xx = gpu.AllocateMemory( sizeof( xyzvector ) * 4, & xx[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_xx ) {
			std::cerr << "Failed to allocate d_xx" << std::endl;
			TS_ASSERT( d_xx );
			return;
		}

		cl_mem d_yy = gpu.AllocateMemory( sizeof( xyzvector ) * 4, & yy[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_yy ) {
			std::cerr << "Failed to allocate d_yy" << std::endl;
			TS_ASSERT( d_yy );
			return;
		}

		cl_mem d_xxmoi = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_xxmoi ) {
			std::cerr << "Failed to allocate d_xxmoi" << std::endl;
			TS_ASSERT( d_xxmoi );
			return;
		}

		cl_mem d_yymoi = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_yymoi ) {
			std::cerr << "Failed to allocate d_yymoi" << std::endl;
			TS_ASSERT( d_yymoi );
			return;
		}

		cl_mem d_uu = gpu.AllocateMemory( sizeof( xyzmatrix ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_uu ) {
			std::cerr << "Failed to allocate d_uu" << std::endl;
			TS_ASSERT( d_uu );
			return;
		}

		cl_mem d_sigma3 = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_sigma3 ) {
			std::cerr << "Failed to allocate d_sigma3" << std::endl;
			TS_ASSERT( d_sigma3 );
			return;
		}

		//cl_mem d_m_moment = gpu.AllocateMemory( sizeof( xyzmatrix ), NULL, CL_MEM_READ_WRITE );
		//if ( ! d_m_moment ) {
		//	std::cerr << "Failed to allocate d_m_moment" << std::endl;
		//	TS_ASSERT( d_m_moment );
		//	return;
		//}

		bool kernel_success;
		kernel_success = gpu.ExecuteKernel( "test_translate_coords_to_origin", 32, 32, 32,
			GPU_DEVMEM, d_xx,
			GPU_DEVMEM, d_xxmoi,
			GPU_INT, 4,
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute test_translate_coords_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		kernel_success = gpu.ExecuteKernel( "test_translate_coords_to_origin", 32, 32, 32,
			GPU_DEVMEM, d_yy,
			GPU_DEVMEM, d_yymoi,
			GPU_INT, 4,
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute test_translate_coords_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		kernel_success = gpu.ExecuteKernel( "test_findUU_no_translation_to_origin", 32, 32, 32,
			GPU_DEVMEM, d_xx,
			GPU_DEVMEM, d_yy,
			GPU_INT, 4,
			GPU_DEVMEM, d_uu,
			GPU_DEVMEM, d_sigma3,
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute test_findUU_no_translation_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< xyzvector > result_xx( 4 ), result_yy( 4 );
		gpu.ReadData( & result_xx[1], d_xx, sizeof( xyzvector ) * 4, CL_TRUE );
		gpu.ReadData( & result_yy[1], d_yy, sizeof( xyzvector ) * 4, CL_TRUE );

		xyzmatrix result_uu;
		result_uu.xx_ = result_uu.xy_ = result_uu.xz_ = -1.0;
		result_uu.yx_ = result_uu.yy_ = result_uu.yz_ = -1.0;
		result_uu.zx_ = result_uu.zy_ = result_uu.zz_ = -1.0;

		gpu.ReadData( &result_uu, d_uu, sizeof( xyzmatrix ), CL_TRUE );

		float result_sigma3;
		gpu.ReadData( &result_sigma3, d_sigma3, sizeof( float ), CL_TRUE );

		//struct xyzmatrix result_m_moment;
		//gpu.ReadData( &result_m_moment, d_m_moment, sizeof( xyzmatrix ), CL_TRUE );
		//std::cout << "m_moment from GPU\n";
		//std::cout << result_m_moment.xx_ << " " << result_m_moment.xy_ << " " << result_m_moment.xz_ << "\n";
		//std::cout << result_m_moment.yx_ << " " << result_m_moment.yy_ << " " << result_m_moment.yz_ << "\n";
		//std::cout << result_m_moment.zx_ << " " << result_m_moment.zy_ << " " << result_m_moment.zz_ << std::endl;;

		ObjexxFCL::FArray2D< core::Real > XX( 3, 4 ), YY( 3, 4 );
		ObjexxFCL::FArray1D< core::Real > WW( 4, 1.0f );
		ObjexxFCL::FArray2D< core::Real > UU( 3, 3 );
		for ( core::Size ii = 1; ii <= 4; ++ii ) {
			XX(1,ii) = xx[ii].x_;
			XX(2,ii) = xx[ii].y_;
			XX(3,ii) = xx[ii].z_;
			YY(1,ii) = yy[ii].x_;
			YY(2,ii) = yy[ii].y_;
			YY(3,ii) = yy[ii].z_;
		}

		core::Real sigma3;
		numeric::model_quality::findUU( XX, YY, WW, 4, UU, sigma3 );

		for ( core::Size ii = 1; ii <= 4; ++ii ) {
			TS_ASSERT_DELTA( result_xx[ii].x_, XX(1,ii), 1e-5 );
			TS_ASSERT_DELTA( result_xx[ii].y_, XX(2,ii), 1e-5 );
			TS_ASSERT_DELTA( result_xx[ii].z_, XX(3,ii), 1e-5 );
		}
		for ( core::Size ii = 1; ii <= 4; ++ii ) {
			TS_ASSERT_DELTA( result_yy[ii].x_, YY(1,ii), 1e-5 );
			TS_ASSERT_DELTA( result_yy[ii].y_, YY(2,ii), 1e-5 );
			TS_ASSERT_DELTA( result_yy[ii].z_, YY(3,ii), 1e-5 );
		}

		TS_ASSERT_DELTA( result_uu.xx_, UU(1,1), 1e-5 );
		TS_ASSERT_DELTA( result_uu.xy_, UU(1,2), 1e-5 );
		TS_ASSERT_DELTA( result_uu.xz_, UU(1,3), 1e-5 );
		TS_ASSERT_DELTA( result_uu.yx_, UU(2,1), 1e-5 );
		TS_ASSERT_DELTA( result_uu.yy_, UU(2,2), 1e-5 );
		TS_ASSERT_DELTA( result_uu.yz_, UU(2,3), 1e-5 );
		TS_ASSERT_DELTA( result_uu.zx_, UU(3,1), 1e-5 );
		TS_ASSERT_DELTA( result_uu.zy_, UU(3,2), 1e-5 );
		TS_ASSERT_DELTA( result_uu.zz_, UU(3,3), 1e-5 );

		TS_ASSERT_DELTA( result_sigma3, sigma3, 1e-5 );

#endif
	}


	void test_transform_center_of_mass_to_origin_and_calculate_rmsd_fast_no_translation_to_origin() {
		TS_ASSERT( true ); // for non-gpu builds
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_translate_coords_to_origin = gpu.BuildKernel( "test_translate_coords_to_origin" );
		if ( ! test_translate_coords_to_origin ) {
			std::cerr << "Failed to build test_translate_coords_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel test_calculate_rmsd_fast_no_translation_to_origin = gpu.BuildKernel( "test_calculate_rmsd_fast_no_translation_to_origin" );
		if ( ! test_calculate_rmsd_fast_no_translation_to_origin ) {
			std::cerr << "Failed to build test_calculate_rmsd_fast_no_translation_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		// Need a bunch of coordinates to give to findUU
		utility::vector1< xyzvector > xx( 4 );
		utility::vector1< xyzvector > yy( 4 );
		xx[1].x_ = 1.0; xx[1].y_ = 0.0; xx[1].z_ = 0.0;
		xx[2].x_ = 0.0; xx[2].y_ = 1.0; xx[2].z_ = 0.0;
		xx[3].x_ = 0.0; xx[3].y_ = 0.0; xx[3].z_ = 1.0;
		xx[4].x_ = 2.0; xx[4].y_ = 2.0; xx[4].z_ = 2.0;

		yy[1].x_ = 1.0; yy[1].y_ = 0.0; yy[1].z_ = 0.0;
		yy[2].x_ = 0.0; yy[2].y_ = 0.0; yy[2].z_ = 1.0;
		yy[3].x_ = 0.0; yy[3].y_ =-1.0; yy[3].z_ = 0.0;
		yy[4].x_ = 2.1; yy[4].y_ =-1.9; yy[4].z_ = 2.2;

		utility::vector1< float > ww( 4, 1.0 );

		cl_mem d_xx = gpu.AllocateMemory( sizeof( xyzvector ) * 4, & xx[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_xx ) {
			std::cerr << "Failed to allocate d_xx" << std::endl;
			TS_ASSERT( d_xx );
			return;
		}

		cl_mem d_yy = gpu.AllocateMemory( sizeof( xyzvector ) * 4, & yy[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_yy ) {
			std::cerr << "Failed to allocate d_yy" << std::endl;
			TS_ASSERT( d_yy );
			return;
		}

		cl_mem d_xxmoi = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_xxmoi ) {
			std::cerr << "Failed to allocate d_xxmoi" << std::endl;
			TS_ASSERT( d_xxmoi );
			return;
		}

		cl_mem d_yymoi = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_yymoi ) {
			std::cerr << "Failed to allocate d_yymoi" << std::endl;
			TS_ASSERT( d_yymoi );
			return;
		}

		cl_mem d_rms_computed = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_rms_computed ) {
			std::cerr << "Failed to allocate d_rms_computed" << std::endl;
			TS_ASSERT( d_rms_computed );
			return;
		}

		//cl_mem d_m_moment = gpu.AllocateMemory( sizeof( xyzmatrix ), NULL, CL_MEM_READ_WRITE );
		//if ( ! d_m_moment ) {
		//	std::cerr << "Failed to allocate d_m_moment" << std::endl;
		//	TS_ASSERT( d_m_moment );
		//	return;
		//}

		bool kernel_success;
		kernel_success = gpu.ExecuteKernel( "test_translate_coords_to_origin", 32, 32, 32,
			GPU_DEVMEM, d_xx,
			GPU_DEVMEM, d_xxmoi,
			GPU_INT, 4,
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute test_translate_coords_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		kernel_success = gpu.ExecuteKernel( "test_translate_coords_to_origin", 32, 32, 32,
			GPU_DEVMEM, d_yy,
			GPU_DEVMEM, d_yymoi,
			GPU_INT, 4,
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute test_translate_coords_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		kernel_success = gpu.ExecuteKernel( "test_calculate_rmsd_fast_no_translation_to_origin", 32, 32, 32,
			GPU_DEVMEM, d_xx,
			GPU_DEVMEM, d_yy,
			GPU_INT, 4,
			GPU_DEVMEM, d_xxmoi,
			GPU_DEVMEM, d_yymoi,
			GPU_DEVMEM, d_rms_computed,
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute test_findUU_no_translation_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< xyzvector > result_xx( 4 ), result_yy( 4 );
		gpu.ReadData( & result_xx[1], d_xx, sizeof( xyzvector ) * 4, CL_TRUE );
		gpu.ReadData( & result_yy[1], d_yy, sizeof( xyzvector ) * 4, CL_TRUE );

		float result_rms_computed;
		gpu.ReadData( &result_rms_computed, d_rms_computed, sizeof( float ), CL_TRUE );

		ObjexxFCL::FArray2D< core::Real > XX( 3, 4 ), YY( 3, 4 );
		ObjexxFCL::FArray1D< core::Real > WW( 4, 1.0f );
		ObjexxFCL::FArray2D< core::Real > UU( 3, 3 );
		for ( core::Size ii = 1; ii <= 4; ++ii ) {
			XX(1,ii) = xx[ii].x_;
			XX(2,ii) = xx[ii].y_;
			XX(3,ii) = xx[ii].z_;
			YY(1,ii) = yy[ii].x_;
			YY(2,ii) = yy[ii].y_;
			YY(3,ii) = yy[ii].z_;
		}

		core::Real sigma3;
		numeric::model_quality::findUU( XX, YY, WW, 4, UU, sigma3 );
		float rms_cpu;
		numeric::model_quality::calc_rms_fast( rms_cpu, XX, YY, WW, 4, sigma3 );

		TS_ASSERT_DELTA( result_rms_computed, rms_cpu, 1e-5 );
		//std::cout << "rms: " << result_rms_computed << " vs " << rms_cpu << std::endl;

		for ( core::Size ii = 1; ii <= 4; ++ii ) {
			TS_ASSERT_DELTA( result_xx[ii].x_, XX(1,ii), 1e-5 );
			TS_ASSERT_DELTA( result_xx[ii].y_, XX(2,ii), 1e-5 );
			TS_ASSERT_DELTA( result_xx[ii].z_, XX(3,ii), 1e-5 );
		}
		for ( core::Size ii = 1; ii <= 4; ++ii ) {
			TS_ASSERT_DELTA( result_yy[ii].x_, YY(1,ii), 1e-5 );
			TS_ASSERT_DELTA( result_yy[ii].y_, YY(2,ii), 1e-5 );
			TS_ASSERT_DELTA( result_yy[ii].z_, YY(3,ii), 1e-5 );
		}

#endif
	}


	void test_compute_rmsd_and_clash_scores_() {
		TS_ASSERT( true );
#ifdef USEOPENCL
		utility::vector1< numeric::xyzVector< float > > rmsdhelices = rmsd_coordinates_for_two_helical_pairs(); // helix pair 1 and 10
		utility::vector1< numeric::xyzVector< float > > collisionhelices = collision_coordinates_for_two_helix_pairs(); // helix pair 1 and 10
		run_compute_rmsd_and_clash_scores_test_for_helical_pairs( rmsdhelices, collisionhelices );
#endif
	}


	void run_compute_rmsd_and_clash_scores_test_for_helical_pairs(
		utility::vector1< numeric::xyzVector< float > > & rmsdhelices,
		utility::vector1< numeric::xyzVector< float > > & collisionhelices
	)
	{
#ifdef USEOPENCL
		utility::vector1< int > col_coords_offsets( 2 ); col_coords_offsets[1]=0; col_coords_offsets[2] = 112;
		utility::vector1< int > n_coords_for_col_calc( 2, 112 );
		utility::vector1< int > comparison_coords_ind_list( 4, 0 );
		comparison_coords_ind_list[1] = 29; comparison_coords_ind_list[2] = 85;
		comparison_coords_ind_list[3] = 29; comparison_coords_ind_list[4] = 85;
		utility::vector1< int > comparison_coord_ind_offsets( 2, 0 ); comparison_coord_ind_offsets[2] = 2;
		utility::vector1< int > n_comparison_coords( 2, 2 );
		utility::vector1< unsigned char > calculate_pair_table( 4, 0 ); calculate_pair_table[ 2 ] = 1; // calculate_pair_table[3] = 1;
		utility::vector1< float > null_table( 4, -1234 );
		utility::vector1< int > bundle_inds( 2, 1 ); bundle_inds[2] = 2; // [ 1, 2 ]


		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel translate_rmscoords_to_origin = gpu.BuildKernel( "translate_rmscoords_to_origin" );
		if ( ! translate_rmscoords_to_origin ) {
			std::cerr << "Failed to build translate_rmscoords_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel compute_rmsd_and_clash_scores = gpu.BuildKernel( "compute_rmsd_and_clash_scores" );
		if ( ! compute_rmsd_and_clash_scores ) {
			std::cerr << "Failed to build compute_rmsd_and_clash_scores kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		cl_mem d_bundle_inds = gpu.AllocateMemory( sizeof( int ) * 2, & bundle_inds[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_bundle_inds ) {
			std::cerr << "Failed to allocate d_bundle_inds" << std::endl;
			TS_ASSERT( d_bundle_inds );
			return;
		}
		cl_mem d_rmscoords = gpu.AllocateMemory( sizeof( xyzvector ) * rmsdhelices.size(), & rmsdhelices[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_rmscoords ) {
			std::cerr << "Failed to allocate d_rmscoords" << std::endl;
			TS_ASSERT( d_rmscoords );
			return;
		}
		cl_mem d_rms_moments_of_inertia = gpu.AllocateMemory( sizeof( float ) * 2, NULL, CL_MEM_READ_WRITE );
		if ( ! d_rms_moments_of_inertia ) {
			std::cerr << "Failed to allocate d_rms_moments_of_inertia" << std::endl;
			TS_ASSERT( d_rms_moments_of_inertia );
			return;
		}

		cl_mem d_collision_coords = gpu.AllocateMemory( sizeof( xyzvector ) * collisionhelices.size(), & collisionhelices[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_collision_coords ) {
			std::cerr << "Failed to allocate d_collision_coords" << std::endl;
			TS_ASSERT( d_collision_coords );
			return;
		}

		cl_mem d_col_coords_offsets = gpu.AllocateMemory( sizeof( int ) * 2, & col_coords_offsets[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_col_coords_offsets ) {
			std::cerr << "Failed to allocate d_col_coords_offsets" << std::endl;
			TS_ASSERT( d_col_coords_offsets );
			return;
		}

		cl_mem d_n_coords_for_col_calc = gpu.AllocateMemory( sizeof( int ) * 2, & n_coords_for_col_calc[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_n_coords_for_col_calc ) {
			std::cerr << "Failed to allocate d_n_coords_for_col_calc" << std::endl;
			TS_ASSERT( d_n_coords_for_col_calc );
			return;
		}

		cl_mem d_comparison_coords_ind_list = gpu.AllocateMemory( sizeof( int ) * 4, & comparison_coords_ind_list[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_comparison_coords_ind_list ) {
			std::cerr << "Failed to allocate d_comparison_coords_ind_list" << std::endl;
			TS_ASSERT( d_comparison_coords_ind_list );
			return;
		}

		cl_mem d_comparison_coord_ind_offsets = gpu.AllocateMemory( sizeof( int ) * 2, & comparison_coord_ind_offsets[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_comparison_coord_ind_offsets ) {
			std::cerr << "Failed to allocate d_comparison_coord_ind_offsets" << std::endl;
			TS_ASSERT( d_comparison_coord_ind_offsets );
			return;
		}

		cl_mem d_n_comparison_coords = gpu.AllocateMemory( sizeof( int ) * 2, & n_comparison_coords[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_n_comparison_coords ) {
			std::cerr << "Failed to allocate d_n_comparison_coords" << std::endl;
			TS_ASSERT( d_n_comparison_coords );
			return;
		}

		cl_mem d_rms_table = gpu.AllocateMemory( sizeof( float ) * 4, &null_table[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR);
		if ( ! d_rms_table ) {
			std::cerr << "Failed to allocate d_rms_table" << std::endl;
			TS_ASSERT( d_rms_table );
			return;
		}

		cl_mem d_clash_table = gpu.AllocateMemory( sizeof( float ) * 4, NULL, CL_MEM_READ_WRITE );
		if ( ! d_clash_table ) {
			std::cerr << "Failed to allocate d_clash_table" << std::endl;
			TS_ASSERT( d_clash_table );
			return;
		}

		cl_mem d_calculate_pair_table = gpu.AllocateMemory( sizeof( unsigned char ) * 4, & calculate_pair_table[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_calculate_pair_table ) {
			std::cerr << "Failed to allocate d_calculate_pair_table" << std::endl;
			TS_ASSERT( d_calculate_pair_table );
			return;
		}

		bool kernel_success;
		kernel_success = gpu.ExecuteKernel( "translate_rmscoords_to_origin", 32, 32, 32,
			GPU_INT, 2, // n_nodes
			GPU_INT, 112, // number of coordinates per helical pair
			GPU_DEVMEM, d_rmscoords,
			GPU_DEVMEM, d_rms_moments_of_inertia,
			GPU_DEVMEM, d_collision_coords,
			GPU_DEVMEM, d_col_coords_offsets,
			GPU_DEVMEM, d_n_coords_for_col_calc,
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute translate_rmscoords_to_origin kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< numeric::xyzVector< float > > translated_rms_coords( 224 );
		bool readRMSCoordsSuccess = gpu.ReadData( & translated_rms_coords[1], d_rmscoords, sizeof( xyzvector ) * 224 );
		if ( ! readRMSCoordsSuccess ) { TS_ASSERT( readRMSCoordsSuccess ); return; }

		utility::vector1< numeric::xyzVector< float > > translated_col_coords( 224 );
		bool readColCoordsSuccess = gpu.ReadData( & translated_col_coords[1], d_collision_coords, sizeof( xyzvector ) * 224 );
		if ( ! readColCoordsSuccess ) { TS_ASSERT( readColCoordsSuccess ); return; }

		utility::vector1< float > gpu_mois( 2, 0 );
		bool readMOIsSuccess = gpu.ReadData( & gpu_mois[1], d_rms_moments_of_inertia, sizeof( float ) * 2 );
		if ( ! readMOIsSuccess ) { TS_ASSERT( readMOIsSuccess ); return; }

		numeric::xyzVector< float > com1(0.0f), com2(0.0f);
		for ( core::Size ii = 1; ii <= 112; ++ii ) {
			com1 += rmsdhelices[ ii ];
			com2 += rmsdhelices[ ii + 112 ];
		}
		com1 /= 112; com2 /= 112;
		for ( core::Size ii = 1; ii <= 112; ++ii ) {
			TS_ASSERT( translated_rms_coords[ ii     ].distance( rmsdhelices[ ii     ] - com1 ) < 1e-5 );
			TS_ASSERT( translated_rms_coords[ ii+112 ].distance( rmsdhelices[ ii+112 ] - com2 ) < 1e-5 );
			TS_ASSERT( translated_col_coords[ ii     ].distance( collisionhelices[ ii     ] - com1 ) < 1e-5 );
			TS_ASSERT( translated_col_coords[ ii+112 ].distance( collisionhelices[ ii+112 ] - com2 ) < 1e-5 );
		}
		utility::vector1< float > cpu_mois(2,0);
		for ( core::Size ii = 1; ii <= 112; ++ii ) {
			cpu_mois[1] += translated_rms_coords[ ii ].length_squared();
			cpu_mois[2] += translated_rms_coords[ ii+112 ].length_squared();
		}
		TS_ASSERT_DELTA( gpu_mois[ 1 ], cpu_mois[ 1 ], 1e-2 );
		TS_ASSERT_DELTA( gpu_mois[ 2 ], cpu_mois[ 2 ], 1e-2 );

		// OK: now launch the second kernel as a 2DRangeKernel.
		int block_size( 2 ), n_atoms_in_rms_calc( 112 );
		cl_int errNum;
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 0, sizeof(int), &block_size );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 0, BLOCK_SIZE" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 1, sizeof(int), & block_size  ); // num nodes, block 1
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 1, block1_nnodes" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 2, sizeof(int), & block_size ); // num nodes, block 2
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 2, block2_nnodes" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 3, sizeof(int), & n_atoms_in_rms_calc );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 3, ndat[ii+block1_offset].n_atoms_in_rms_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 4, sizeof(cl_mem), & d_bundle_inds );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 4, block1ptrs.coords_for_rms_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 5, sizeof(cl_mem), & d_bundle_inds );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 5, block2ptrs.coords_for_rms_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 6, sizeof(cl_mem), & d_rmscoords );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 4, block1ptrs.coords_for_rms_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 7, sizeof(cl_mem), & d_rmscoords );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 5, block2ptrs.coords_for_rms_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 8, sizeof(cl_mem), & d_rms_moments_of_inertia );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 6, block1ptrs.rms_coords_moment_of_inertia" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores, 9, sizeof(cl_mem), & d_rms_moments_of_inertia );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 7, block2ptrs.rms_coords_moment_of_inertia" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,10, sizeof(cl_mem), & d_collision_coords );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 8, block1ptrs.coords_for_col_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,11, sizeof(cl_mem), & d_collision_coords );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 9, block2ptrs.coords_for_col_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,12, sizeof(cl_mem), & d_col_coords_offsets );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 10, block1ptrs.coords_for_col_calc_offset" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,13, sizeof(cl_mem), & d_col_coords_offsets );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 11, block2ptrs.coords_for_col_calc_offset" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,14, sizeof(cl_mem), & d_n_coords_for_col_calc );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 12, block1ptrs.n_coords_for_col_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,15, sizeof(cl_mem), & d_n_coords_for_col_calc );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 13, block2ptrs.n_coords_for_col_calc" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,16, sizeof(cl_mem), & d_comparison_coords_ind_list );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 14, block1ptrs.comparison_coord_ind_list" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,17, sizeof(cl_mem), & d_comparison_coords_ind_list );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 15, block2ptrs.comparison_coord_ind_list" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,18, sizeof(cl_mem), & d_comparison_coord_ind_offsets );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 16, block1ptrs.comparison_coord_ind_offset" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,19, sizeof(cl_mem), & d_comparison_coord_ind_offsets );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 17, block2ptrs.comparison_coord_ind_offset" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,20, sizeof(cl_mem), & d_n_comparison_coords );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 18, block1ptrs.n_comparison_coords" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,21, sizeof(cl_mem), & d_n_comparison_coords );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 19, block2ptrs.n_comparison_coords" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,22, sizeof(cl_mem), & d_rms_table );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 20, d_rms_table" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,23, sizeof(cl_mem), & d_clash_table );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 21, d_clash_table" ); }
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,24, sizeof(cl_mem), & d_calculate_pair_table );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 22, d_calculate_pair_table" ); }

		//#define DEBUG_FINDUU_GPU_FUNCTION

#ifdef DEBUG_FINDUU_GPU_FUNCTION
		cl_mem d_UUout = gpu.AllocateMemory( sizeof( float ) * 9, NULL, CL_MEM_READ_WRITE );
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,25, sizeof(cl_mem), & d_UUout );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 23, d_UUout" ); }
		cl_mem d_sigma3out = gpu.AllocateMemory( sizeof( float ), NULL, CL_MEM_READ_WRITE );
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,26, sizeof(cl_mem), & d_sigma3out );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 24, d_sigma3out" ); }
		cl_mem d_m_moment_out = gpu.AllocateMemory( sizeof( float ) * 9, NULL, CL_MEM_READ_WRITE );
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,27, sizeof(cl_mem), & d_m_moment_out );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 25, d_m_moment_out" ); }
#endif

		//#define DEBUG_COLLISION_CALCULATION
#ifdef DEBUG_COLLISION_CALCULATION
		cl_mem d_all_square_distances = gpu.AllocateMemory( sizeof( float ) * 4 * 112, NULL, CL_MEM_READ_WRITE );
		errNum = clSetKernelArg(compute_rmsd_and_clash_scores,23, sizeof(cl_mem), &d_all_square_distances );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to set arg 23, d_all_square_distances" ); }
#endif

		size_t global_work_size[2] = { 8, 8 };
		size_t local_work_size[2] = {8, 8};
		cl_event kernelEvent;

		basic::gpu::GPU_DEV const & gpu_dev = gpu.device();

		errNum = clEnqueueNDRangeKernel( gpu_dev.commandQueue, compute_rmsd_and_clash_scores,
			2, NULL, global_work_size, local_work_size, 0, NULL, &kernelEvent );
		if ( errNum != CL_SUCCESS ) { utility_exit_with_message( "Failed to run calc rms and clash score kernel" ); }

		clWaitForEvents( 1, &kernelEvent );

		ObjexxFCL::FArray2D< float > result_rms(2,2,1234.5f);
		ObjexxFCL::FArray2D< float > result_clash(2,2,1234.5f);

		bool rmsReadSuccess = gpu.ReadData( &result_rms(1,1), d_rms_table, sizeof( float ) * 4, CL_TRUE );
		if ( ! rmsReadSuccess ) { TS_ASSERT( rmsReadSuccess ); return; }
		//for ( int ii = 1; ii <= 2; ++ii ) {
		//	for ( int jj = 1; jj <= 2; ++jj ) {
		//		std::cout << "result_rms: ii " << " jj " << result_rms(jj,ii) << std::endl;
		//	}
		//}
		bool clashReadSuccess = gpu.ReadData( &result_clash(1,1), d_clash_table, sizeof( float ) * 4, CL_TRUE );
		if ( ! clashReadSuccess ) { TS_ASSERT( clashReadSuccess ); return; }

#ifdef DEBUG_FINDUU_GPU_FUNCTION
		ObjexxFCL::FArray2D< float > result_UUout(3,3,0.0f);
		float result_sigma3out;
		bool uuoutReadSuccess = gpu.ReadData( &result_UUout(1,1), d_UUout, sizeof(float)*9, CL_TRUE );
		if ( ! uuoutReadSuccess ) { TS_ASSERT( uuoutReadSuccess ); return; }
		bool sigma3outReadSuccess = gpu.ReadData( &result_sigma3out, d_sigma3out, sizeof(float), CL_TRUE );
		if ( ! sigma3outReadSuccess ) { TS_ASSERT( sigma3outReadSuccess ); return; }
		struct xyzmatrix result_m_moment;
		gpu.ReadData( &result_m_moment, d_m_moment_out, sizeof( xyzmatrix ), CL_TRUE );
		std::cout << "m_moment from GPU\n";
		std::cout << result_m_moment.xx_ << " " << result_m_moment.xy_ << " " << result_m_moment.xz_ << "\n";
		std::cout << result_m_moment.yx_ << " " << result_m_moment.yy_ << " " << result_m_moment.yz_ << "\n";
		std::cout << result_m_moment.zx_ << " " << result_m_moment.zy_ << " " << result_m_moment.zz_ << std::endl;;
#endif

		// calculate RMSD on the CPU
		ObjexxFCL::FArray2D< numeric::Real > p1_coords( 3, 112 );
		ObjexxFCL::FArray2D< numeric::Real > p2_coords( 3, 112 );
		for ( core::Size ii = 1; ii <= 112; ++ii ) {
			for ( core::Size kk = 1; kk <= 3; ++kk ) {
				p1_coords(kk,ii) = rmsdhelices[ ii ](kk);
				p2_coords(kk,ii) = rmsdhelices[ ii+112 ](kk);
			}
		}
		ObjexxFCL::FArray1D< numeric::Real > ww( 112, 1.0 );//weight matrix, all 1 for my purposes
		ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );//transformation matrix
		numeric::Real ctx;
#ifndef DEBUG_FINDUU_GPU_FUNCTION
		numeric::model_quality::findUU( p1_coords, p2_coords, ww, 112, uu, ctx );
#else
		my_findUU( p1_coords, p2_coords, ww, 112, uu, ctx );

		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			for ( core::Size jj = 1; jj <= 3; ++jj ) {
				TS_ASSERT_DELTA( result_UUout(jj,ii), uu(ii,jj), 1e-5 ); // why do I have to reverse the indices?!
			}
		}
		TS_ASSERT_DELTA( result_sigma3out, ctx, 1e-3 );
#endif

		float rms_from_calc_fast_rms;
		numeric::model_quality::calc_rms_fast( rms_from_calc_fast_rms, p1_coords, p2_coords, ww, 112, ctx );

    for ( core::Size ii = 1; ii <= 112; ++ii ) {
			numeric::xyzVector< float > iip1( p1_coords(1,ii), p1_coords(2,ii), p1_coords(3,ii) );
			numeric::xyzVector< float > iip2( p2_coords(1,ii), p2_coords(2,ii), p2_coords(3,ii) );
      TS_ASSERT( translated_rms_coords[ ii     ].distance( iip1 ) < 1e-5 );
      TS_ASSERT( translated_rms_coords[ ii+112 ].distance( iip2 ) < 1e-5 );
    }

		TS_ASSERT_DELTA( result_rms(2,1), rms_from_calc_fast_rms, 1e-2 );

		// not calculated
		TS_ASSERT_DELTA( result_rms(1,1), -1234, 1e-5 );
		TS_ASSERT_DELTA( result_rms(1,2), -1234, 1e-5 );
		TS_ASSERT_DELTA( result_rms(2,2), -1234, 1e-5 );

		/// Pull down node 2's rotated coordinates -- check that they're in the same position as where we woudl rotate them to
		/// on the CPU
		utility::vector1< numeric::xyzVector< float > > clash_coords_after_rotation( 224 );
		gpu.ReadData( &clash_coords_after_rotation[1], d_collision_coords, sizeof(float)*3*224 );

		//for ( core::Size ii = 1; ii <= 112; ++ii ) {
		//	numeric::xyzVector< float > const & iicoord = translated_col_coords[ 112 + ii ];
		//	numeric::xyzVector< float > iicoord_rotated;
		//	iicoord_rotated(1) = ( uu(1,1)*iicoord(1) )+( uu(1,2)*iicoord(2) ) +( uu(1,3)*iicoord(3) );
		//	iicoord_rotated(2) = ( uu(2,1)*iicoord(1) )+( uu(2,2)*iicoord(2) ) +( uu(2,3)*iicoord(3) );
		//	iicoord_rotated(3) = ( uu(3,1)*iicoord(1) )+( uu(3,2)*iicoord(2) ) +( uu(3,3)*iicoord(3) );
		//	if ( iicoord_rotated.distance_squared( clash_coords_after_rotation[ ii + 112 ] ) > 1e-4 ) {
		//		std::cout << "cpu: " << iicoord_rotated.x() << " "<< iicoord_rotated.y() << " "<< iicoord_rotated.z() << " ";
		//		std::cout << " gpu: " << clash_coords_after_rotation[ ii + 112 ].x() << " " << clash_coords_after_rotation[ ii + 112 ].y() << " " << clash_coords_after_rotation[ ii + 112 ].z() << std::endl;
		//	}
		//	TS_ASSERT( iicoord_rotated.distance_squared( clash_coords_after_rotation[ ii + 112 ] ) < 1e-4 );
		//}

#ifdef DEBUG_COLLISION_CALCULATION
		// pull down the collision calculation results
		utility::vector1< float > result_all_square_distances( 4 * 112 );
		bool allD2ReadSuccess = gpu.ReadData( & result_all_square_distances[1], d_all_square_distances, sizeof(float) * 112 * 4 );
		if ( ! allD2ReadSuccess ) { utility_exit_with_message( "failed to read from d_all_square_distances" ); }
		core::Size count_n_square_distances_examined = 0;
#endif


		core::Real closest_contact = 12345;
		{ // scope
			for ( core::Size ii = 1; ii <= 2; ++ii ) {
				numeric::xyzVector< float > const & iicoord = translated_col_coords[ comparison_coords_ind_list[ ii ] + 1 ];
				for ( core::Size jj = 1; jj <= 112; ++jj ) {
					numeric::xyzVector< float > const & jjcoord = translated_col_coords[ 112 + jj ];
					numeric::xyzVector< float > jjcoord_rotated;
          jjcoord_rotated(1) = ( uu(1,1)*jjcoord(1) )+( uu(1,2)*jjcoord(2) ) +( uu(1,3)*jjcoord(3) );
          jjcoord_rotated(2) = ( uu(2,1)*jjcoord(1) )+( uu(2,2)*jjcoord(2) ) +( uu(2,3)*jjcoord(3) );
          jjcoord_rotated(3) = ( uu(3,1)*jjcoord(1) )+( uu(3,2)*jjcoord(2) ) +( uu(3,3)*jjcoord(3) );
					core::Real d2 = iicoord.distance_squared( jjcoord_rotated );
					if ( d2 < closest_contact ) closest_contact = d2;
#ifdef DEBUG_COLLISION_CALCULATION
					++count_n_square_distances_examined;
					TS_ASSERT_DELTA( result_all_square_distances[ count_n_square_distances_examined ], d2, 1e-3 );
#endif
				}
			}
			for ( core::Size ii = 1; ii <= 2; ++ii ) {
				numeric::xyzVector< float > const & iicoord = translated_col_coords[ 112 + comparison_coords_ind_list[ ii+2 ] + 1 ];
				numeric::xyzVector< float > iicoord_rotated;
				iicoord_rotated(1) = ( uu(1,1)*iicoord(1) )+( uu(1,2)*iicoord(2) ) +( uu(1,3)*iicoord(3) );
				iicoord_rotated(2) = ( uu(2,1)*iicoord(1) )+( uu(2,2)*iicoord(2) ) +( uu(2,3)*iicoord(3) );
				iicoord_rotated(3) = ( uu(3,1)*iicoord(1) )+( uu(3,2)*iicoord(2) ) +( uu(3,3)*iicoord(3) );
				for ( core::Size jj = 1; jj <= 112; ++jj ) {
					numeric::xyzVector< float > const & jjcoord = translated_col_coords[ jj ];
					core::Real d2 = iicoord_rotated.distance_squared( jjcoord );
					if ( d2 < closest_contact ) closest_contact = d2;
#ifdef DEBUG_COLLISION_CALCULATION
					++count_n_square_distances_examined;
					TS_ASSERT_DELTA( result_all_square_distances[ count_n_square_distances_examined ], d2, 1e-3 );
#endif
				}
			}
		}
		//std::cout << "closest contact GPU: " << result_clash(2,1) << " CPU: " << closest_contact << std::endl;
		TS_ASSERT_DELTA( result_clash(2,1), closest_contact, 1e-3 );
#endif

	}

	void test_compute_rmsd_and_clash_scores_184_vs_196() {
		TS_ASSERT( true );
#ifdef USEOPENCL
		utility::vector1< numeric::xyzVector< float > > rmsdhelices = rmsd_coords_helices_184_196(); // helix pair 184 and 196
		utility::vector1< numeric::xyzVector< float > > collisionhelices = clash_coords_helices_184_196(); // helix pair 184 and 196
    run_compute_rmsd_and_clash_scores_test_for_helical_pairs( rmsdhelices, collisionhelices );
#endif
	}

	void test_scan_good_rmsd_clash_pairs2() {
		TS_ASSERT( true );
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel scan_good_rmsd_clash_pairs = gpu.BuildKernel( "scan_good_rmsd_clash_pairs" );
		if ( ! scan_good_rmsd_clash_pairs ) {
			std::cerr << "Failed to build scan_good_rmsd_clash_pairs kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		// Need a bunch of rmsd/clash pairs to give to scan_good_rmsd_clash_pairs
		utility::vector1< float > rms_vals( 1024 );
		utility::vector1< float > clash_vals( 1024 );
		utility::vector1< unsigned char > calculate_pair_vals( 1024, 0 );
		int count = 0;
		int n_nodes_1( 30 ), n_nodes_2( 19 );
		for ( int ii = 1; ii <= 32; ++ii ) {
			bool iigood = ii <= n_nodes_1;
			for ( int jj = 1; jj <= 32; ++jj ) {
				bool jjgood = jj <= n_nodes_2;
				++count;
				if ( iigood && jjgood ) calculate_pair_vals[ count ] = 1;
			}
		}
		count = 0;
		for ( core::Size ii = 1; ii <= 1024 / 16; ++ii ) {
			rms_vals[ ++count ] = 1.01; clash_vals[ count ] = 9.01; // good
			rms_vals[ ++count ] = 2.02; clash_vals[ count ] = 9.02; // bad rms
			rms_vals[ ++count ] = 2.03; clash_vals[ count ] = 1.03; // bad both
			rms_vals[ ++count ] = 1.04; clash_vals[ count ] = 9.04; // good

			rms_vals[ ++count ] = 1.05; clash_vals[ count ] = 1.05; // bad clash
			rms_vals[ ++count ] = 1.06; clash_vals[ count ] = 1.06; // bad clash
			rms_vals[ ++count ] = 2.07; clash_vals[ count ] = 1.07; // bad both
			rms_vals[ ++count ] = 1.08; clash_vals[ count ] = 9.08; // good

			rms_vals[ ++count ] = 1.09; clash_vals[ count ] = 9.09; // good
			rms_vals[ ++count ] = 2.00; clash_vals[ count ] = 9.00; // bad rms
			rms_vals[ ++count ] = 2.01; clash_vals[ count ] = 1.01; // bad both
			rms_vals[ ++count ] = 1.02; clash_vals[ count ] = 1.02; // bad clash

			rms_vals[ ++count ] = 1.03; clash_vals[ count ] = 9.03; // good
			rms_vals[ ++count ] = 2.04; clash_vals[ count ] = 9.04; // bad
			rms_vals[ ++count ] = 2.05; clash_vals[ count ] = 1.05; // bad
			rms_vals[ ++count ] = 1.06; clash_vals[ count ] = 9.06; // good
		}


		cl_mem d_calculate_pair = gpu.AllocateMemory( sizeof( float ) * 1024, & calculate_pair_vals[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_calculate_pair ) {
			std::cerr << "Failed to allocate d_calculate_pair" << std::endl;
			TS_ASSERT( d_calculate_pair );
			return;
		}

		cl_mem d_rms = gpu.AllocateMemory( sizeof( float ) * 1024, & rms_vals[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_rms ) {
			std::cerr << "Failed to allocate d_rms" << std::endl;
			TS_ASSERT( d_rms );
			return;
		}

		cl_mem d_clash = gpu.AllocateMemory( sizeof( float ) * 1024, & clash_vals[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_clash ) {
			std::cerr << "Failed to allocate d_clash" << std::endl;
			TS_ASSERT( d_clash );
			return;
		}

		cl_mem d_good_indices = gpu.AllocateMemory( sizeof( int ) * 1024, NULL, CL_MEM_READ_WRITE );
		if ( ! d_good_indices ) {
			std::cerr << "Failed to allocate d_good_indices" << std::endl;
			TS_ASSERT( d_good_indices );
			return;
		}

		cl_mem d_n_good_indices = gpu.AllocateMemory( sizeof( int ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_n_good_indices ) {
			std::cerr << "Failed to allocate d_n_good_indices" << std::endl;
			TS_ASSERT( d_n_good_indices );
			return;
		}

#define DEBUG_SCAN_GOOD_RMS_AND_CLASH
#ifdef DEBUG_SCAN_GOOD_RMS_AND_CLASH
		cl_mem d_scan_results = gpu.AllocateMemory( sizeof( int ) * 1024, NULL, CL_MEM_READ_WRITE );
		if ( ! d_scan_results ) {
			std::cerr << "Failed to allocate d_scan_results" << std::endl;
			TS_ASSERT( d_scan_results );
			return;
		}
		cl_mem d_start_scan_results = gpu.AllocateMemory( sizeof( int ) * 1024, NULL, CL_MEM_READ_WRITE );
		if ( ! d_start_scan_results ) {
			std::cerr << "Failed to allocate d_start_scan_results" << std::endl;
			TS_ASSERT( d_start_scan_results );
			return;
		}
#endif

		bool kernel_success = gpu.ExecuteKernel( "scan_good_rmsd_clash_pairs", 256, 256, 256,
			GPU_INT, 32,
			GPU_INT, n_nodes_1,
			GPU_INT, n_nodes_2,
			GPU_FLOAT, 1.5, // max acceptible rms
			GPU_FLOAT, 4.0, // min tolerated d2
			GPU_DEVMEM, d_calculate_pair,
			GPU_DEVMEM, d_rms,
			GPU_DEVMEM, d_clash,
			GPU_DEVMEM, d_good_indices,
			GPU_DEVMEM, d_n_good_indices,
#ifdef DEBUG_SCAN_GOOD_RMS_AND_CLASH
			GPU_DEVMEM, d_scan_results,
			GPU_DEVMEM, d_start_scan_results,
#endif
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute scan_good_rmsd_clash_pairs kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< float > result_rms( 1024 ), result_clash( 1024 );
		utility::vector1< int > result_good_indices( 1024 );
		int result_n_good_indices;
		gpu.ReadData( & result_rms[1], d_rms, sizeof( float ) * 1024, CL_TRUE );
		gpu.ReadData( & result_clash[1], d_clash, sizeof( float ) * 1024, CL_TRUE );
		gpu.ReadData( & result_good_indices[1], d_good_indices, sizeof( int ) * 1024, CL_TRUE );
		gpu.ReadData( & result_n_good_indices, d_n_good_indices, sizeof( int ), CL_TRUE );

#ifdef DEBUG_SCAN_GOOD_RMS_AND_CLASH
		utility::vector1< int > scan_results( 1024 );
		gpu.ReadData( & scan_results[1], d_scan_results, sizeof(int) * 1024, CL_TRUE );
		utility::vector1< int > start_scan_results( 1024 );
		gpu.ReadData( & start_scan_results[1], d_start_scan_results, sizeof(int) * 1024, CL_TRUE );
#endif

		count = 0;
		int cumulative_sum = 0;
		for ( int ii = 1; ii <= 32; ++ii ) {
			bool ii_in_range = ii <= n_nodes_1;
			for ( int jj = 1; jj <= 32; ++jj ) {
				bool jj_in_range = jj <= n_nodes_2;
				++count;
#ifdef DEBUG_SCAN_GOOD_RMS_AND_CLASH
				TS_ASSERT( cumulative_sum == scan_results[ count ] );
				//std::cout << "count: " << count << " cumsum: " << cumulative_sum << " scan_results: " << scan_results[ count ] << " start_scan_results: " << start_scan_results[ count ] << std::endl;
#endif
				bool iijjgood = ii_in_range && jj_in_range && rms_vals[ count ] <= 1.5 && clash_vals[ count ] >= 4 && calculate_pair_vals[ count ] == 1;
				if ( iijjgood ) {
					//std::cout << "count: " << count << " rms_val: " << rms_vals[ count ] << " result_rms: " << result_rms[ 1+cumulative_sum ];
					//std::cout << " clash_val: " << clash_vals[ count ] << " result_good_indices: " << result_good_indices[ 1+cumulative_sum ] << std::endl;
					TS_ASSERT( result_rms[          1+cumulative_sum ] == rms_vals[   count ] );
					TS_ASSERT( result_clash[        1+cumulative_sum ] == clash_vals[ count ] );
					TS_ASSERT( result_good_indices[ 1+cumulative_sum ] == count-1             );
					++cumulative_sum;
				}
			}
		}
		//std::cout << "final cumulative sum: " << cumulative_sum << " and from GPU: " << result_n_good_indices << std::endl;
		TS_ASSERT( cumulative_sum == result_n_good_indices );
#endif
	}

	void test_scan_good_rmsd_clash_pairs1() {
		TS_ASSERT( true );
#ifdef USEOPENCL
		basic::gpu::GPU gpu;
		std::vector< std::string > programs = gpu_rmsd_test_programs();
		if ( ! gpu.RegisterProgram( programs ) ) {
			std::cerr << "Failed to build program gpu_xyzmatrix.cl" << std::endl;
			TS_ASSERT( false );
			return;
		}
		cl_kernel scan_good_rmsd_clash_pairs = gpu.BuildKernel( "scan_good_rmsd_clash_pairs" );
		if ( ! scan_good_rmsd_clash_pairs ) {
			std::cerr << "Failed to build scan_good_rmsd_clash_pairs kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		// Need a bunch of rmsd/clash pairs to give to scan_good_rmsd_clash_pairs
		utility::vector1< float > rms_vals( 1024 );
		utility::vector1< float > clash_vals( 1024 );
		utility::vector1< unsigned char > calculate_pair_vals( 1024, 0 );
		int count = 0;
		int n_nodes_1( 32 ), n_nodes_2( 32 );
		for ( int ii = 1; ii <= 32; ++ii ) {
			bool iigood = ii <= n_nodes_1;
			for ( int jj = 1; jj <= 32; ++jj ) {
				bool jjgood = jj <= n_nodes_2 && ii < jj;
				++count;
				if ( iigood && jjgood ) calculate_pair_vals[ count ] = 1;
			}
		}
		count = 0;
		for ( core::Size ii = 1; ii <= 1024 / 16; ++ii ) {
			rms_vals[ ++count ] = 1.01; clash_vals[ count ] = 9.01; // good
			rms_vals[ ++count ] = 2.02; clash_vals[ count ] = 9.02; // bad rms
			rms_vals[ ++count ] = 2.03; clash_vals[ count ] = 1.03; // bad both
			rms_vals[ ++count ] = 1.04; clash_vals[ count ] = 9.04; // good

			rms_vals[ ++count ] = 1.05; clash_vals[ count ] = 1.05; // bad clash
			rms_vals[ ++count ] = 1.06; clash_vals[ count ] = 1.06; // bad clash
			rms_vals[ ++count ] = 2.07; clash_vals[ count ] = 1.07; // bad both
			rms_vals[ ++count ] = 1.08; clash_vals[ count ] = 9.08; // good

			rms_vals[ ++count ] = 1.09; clash_vals[ count ] = 9.09; // good
			rms_vals[ ++count ] = 2.00; clash_vals[ count ] = 9.00; // bad rms
			rms_vals[ ++count ] = 2.01; clash_vals[ count ] = 1.01; // bad both
			rms_vals[ ++count ] = 1.02; clash_vals[ count ] = 1.02; // bad clash

			rms_vals[ ++count ] = 1.03; clash_vals[ count ] = 9.03; // good
			rms_vals[ ++count ] = 2.04; clash_vals[ count ] = 9.04; // bad
			rms_vals[ ++count ] = 2.05; clash_vals[ count ] = 1.05; // bad
			rms_vals[ ++count ] = 2.06; clash_vals[ count ] = 9.06; // bad
		}


		cl_mem d_calculate_pair = gpu.AllocateMemory( sizeof( float ) * 1024, & calculate_pair_vals[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_calculate_pair ) {
			std::cerr << "Failed to allocate d_calculate_pair" << std::endl;
			TS_ASSERT( d_calculate_pair );
			return;
		}

		cl_mem d_rms = gpu.AllocateMemory( sizeof( float ) * 1024, & rms_vals[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_rms ) {
			std::cerr << "Failed to allocate d_rms" << std::endl;
			TS_ASSERT( d_rms );
			return;
		}

		cl_mem d_clash = gpu.AllocateMemory( sizeof( float ) * 1024, & clash_vals[1], CL_MEM_READ_WRITE & CL_MEM_COPY_HOST_PTR );
		if ( ! d_clash ) {
			std::cerr << "Failed to allocate d_clash" << std::endl;
			TS_ASSERT( d_clash );
			return;
		}

		cl_mem d_good_indices = gpu.AllocateMemory( sizeof( int ) * 1024, NULL, CL_MEM_READ_WRITE );
		if ( ! d_good_indices ) {
			std::cerr << "Failed to allocate d_good_indices" << std::endl;
			TS_ASSERT( d_good_indices );
			return;
		}

		cl_mem d_n_good_indices = gpu.AllocateMemory( sizeof( int ), NULL, CL_MEM_READ_WRITE );
		if ( ! d_n_good_indices ) {
			std::cerr << "Failed to allocate d_n_good_indices" << std::endl;
			TS_ASSERT( d_n_good_indices );
			return;
		}

		//#define DEBUG_SCAN_GOOD_RMS_AND_CLASH
#ifdef DEBUG_SCAN_GOOD_RMS_AND_CLASH
		cl_mem d_scan_results = gpu.AllocateMemory( sizeof( int ) * 1024, NULL, CL_MEM_READ_WRITE );
		if ( ! d_scan_results ) {
			std::cerr << "Failed to allocate d_scan_results" << std::endl;
			TS_ASSERT( d_scan_results );
			return;
		}
		cl_mem d_start_scan_results = gpu.AllocateMemory( sizeof( int ) * 1024, NULL, CL_MEM_READ_WRITE );
		if ( ! d_start_scan_results ) {
			std::cerr << "Failed to allocate d_start_scan_results" << std::endl;
			TS_ASSERT( d_start_scan_results );
			return;
		}
#endif

		bool kernel_success = gpu.ExecuteKernel( "scan_good_rmsd_clash_pairs", 256, 256, 256,
			GPU_INT, 32,
			GPU_INT, n_nodes_1,
			GPU_INT, n_nodes_2,
			GPU_FLOAT, 1.5, // max acceptible rms
			GPU_FLOAT, 4.0, // min tolerated d2
			GPU_DEVMEM, d_calculate_pair,
			GPU_DEVMEM, d_rms,
			GPU_DEVMEM, d_clash,
			GPU_DEVMEM, d_good_indices,
			GPU_DEVMEM, d_n_good_indices,
#ifdef DEBUG_SCAN_GOOD_RMS_AND_CLASH
			GPU_DEVMEM, d_scan_results,
			GPU_DEVMEM, d_start_scan_results,
#endif
			NULL );
		if ( ! kernel_success ) {
			std::cerr << "Failed to execute scan_good_rmsd_clash_pairs kernel" << std::endl;
			TS_ASSERT( false );
			return;
		}

		utility::vector1< float > result_rms( 1024 ), result_clash( 1024 );
		utility::vector1< int > result_good_indices( 1024 );
		int result_n_good_indices;
		gpu.ReadData( & result_rms[1], d_rms, sizeof( float ) * 1024, CL_TRUE );
		gpu.ReadData( & result_clash[1], d_clash, sizeof( float ) * 1024, CL_TRUE );
		gpu.ReadData( & result_good_indices[1], d_good_indices, sizeof( int ) * 1024, CL_TRUE );
		gpu.ReadData( & result_n_good_indices, d_n_good_indices, sizeof( int ), CL_TRUE );

#ifdef DEBUG_SCAN_GOOD_RMS_AND_CLASH
		utility::vector1< int > scan_results( 1024 );
		gpu.ReadData( & scan_results[1], d_scan_results, sizeof(int) * 1024, CL_TRUE );
		utility::vector1< int > start_scan_results( 1024 );
		gpu.ReadData( & start_scan_results[1], d_start_scan_results, sizeof(int) * 1024, CL_TRUE );
#endif

		count = 0;
		int cumulative_sum = 0;
		for ( int ii = 1; ii <= 32; ++ii ) {
			bool ii_in_range = ii <= n_nodes_1;
			for ( int jj = 1; jj <= 32; ++jj ) {
				bool jj_in_range = jj <= n_nodes_2;
				++count;
#ifdef DEBUG_SCAN_GOOD_RMS_AND_CLASH
				TS_ASSERT( cumulative_sum == scan_results[ count ] );
				std::cout << "count: " << count << " cumsum: " << cumulative_sum << " scan_results: " << scan_results[ count ] << " start_scan_results: " << start_scan_results[ count ] << std::endl;
#endif
				bool iijjgood = ii_in_range && jj_in_range && rms_vals[ count ] <= 1.5 && clash_vals[ count ] >= 4 && calculate_pair_vals[ count ] == 1;
				if ( iijjgood ) {
					std::cout << "count: " << count << " rms_val: " << rms_vals[ count ] << " result_rms: " << result_rms[ 1+cumulative_sum ];
					std::cout << " clash_val: " << clash_vals[ count ] << " result_good_indices: " << result_good_indices[ 1+cumulative_sum ] << std::endl;
					TS_ASSERT( result_rms[          1+cumulative_sum ] == rms_vals[   count ] );
					TS_ASSERT( result_clash[        1+cumulative_sum ] == clash_vals[ count ] );
					TS_ASSERT( result_good_indices[ 1+cumulative_sum ] == count-1             );
					++cumulative_sum;
				}
			}
		}
		//std::cout << "final cumulative sum: " << cumulative_sum << " and from GPU: " << result_n_good_indices << std::endl;
		TS_ASSERT( cumulative_sum == result_n_good_indices );
#endif
	}

	utility::vector1< numeric::xyzVector< float > >
	rmsd_coordinates_for_two_helical_pairs() {
		utility::vector1< numeric::xyzVector< float > > xyzarray( 224 );
		core::Size ii = 0;
		// coordinates taken from Tim's helical assembly test case, nodes 1 and 20.
		xyzarray[ ++ii ].x() = 16.614; xyzarray[ ii ].y() = 6.327; xyzarray[ ii ].z() = 9.328;
		xyzarray[ ++ii ].x() = 17.629; xyzarray[ ii ].y() = 6.362; xyzarray[ ii ].z() = 8.277;
		xyzarray[ ++ii ].x() = 18.32; xyzarray[ ii ].y() = 5.009; xyzarray[ ii ].z() = 8.073;
		xyzarray[ ++ii ].x() = 18.553; xyzarray[ ii ].y() = 4.59; xyzarray[ ii ].z() = 6.929;
		xyzarray[ ++ii ].x() = 18.654; xyzarray[ ii ].y() = 4.321; xyzarray[ ii ].z() = 9.171;
		xyzarray[ ++ii ].x() = 19.157; xyzarray[ ii ].y() = 2.956; xyzarray[ ii ].z() = 9.07;
		xyzarray[ ++ii ].x() = 18.185; xyzarray[ ii ].y() = 2.097; xyzarray[ ii ].z() = 8.257;
		xyzarray[ ++ii ].x() = 18.59; xyzarray[ ii ].y() = 1.291; xyzarray[ ii ].z() = 7.399;
		xyzarray[ ++ii ].x() = 16.894; xyzarray[ ii ].y() = 2.26; xyzarray[ ii ].z() = 8.537;
		xyzarray[ ++ii ].x() = 15.908; xyzarray[ ii ].y() = 1.449; xyzarray[ ii ].z() = 7.858;
		xyzarray[ ++ii ].x() = 15.918; xyzarray[ ii ].y() = 1.707; xyzarray[ ii ].z() = 6.356;
		xyzarray[ ++ii ].x() = 15.878; xyzarray[ ii ].y() = 0.758; xyzarray[ ii ].z() = 5.554;
		xyzarray[ ++ii ].x() = 16.01; xyzarray[ ii ].y() = 2.984; xyzarray[ ii ].z() = 5.977;
		xyzarray[ ++ii ].x() = 16.098; xyzarray[ ii ].y() = 3.349; xyzarray[ ii ].z() = 4.559;
		xyzarray[ ++ii ].x() = 17.351; xyzarray[ ii ].y() = 2.776; xyzarray[ ii ].z() = 3.881;
		xyzarray[ ++ii ].x() = 17.292; xyzarray[ ii ].y() = 2.279; xyzarray[ ii ].z() = 2.745;
		xyzarray[ ++ii ].x() = 18.485; xyzarray[ ii ].y() = 2.873; xyzarray[ ii ].z() = 4.568;
		xyzarray[ ++ii ].x() = 19.766; xyzarray[ ii ].y() = 2.393; xyzarray[ ii ].z() = 4.02;
		xyzarray[ ++ii ].x() = 19.684; xyzarray[ ii ].y() = 0.873; xyzarray[ ii ].z() = 3.871;
		xyzarray[ ++ii ].x() = 20.057; xyzarray[ ii ].y() = 0.328; xyzarray[ ii ].z() = 2.829;
		xyzarray[ ++ii ].x() = 19.121; xyzarray[ ii ].y() = 0.207; xyzarray[ ii ].z() = 4.879;
		xyzarray[ ++ii ].x() = 18.844; xyzarray[ ii ].y() = -1.248; xyzarray[ ii ].z() = 4.795;
		xyzarray[ ++ii ].x() = 18.026; xyzarray[ ii ].y() = -1.658; xyzarray[ ii ].z() = 3.583;
		xyzarray[ ++ii ].x() = 18.348; xyzarray[ ii ].y() = -2.64; xyzarray[ ii ].z() = 2.892;
		xyzarray[ ++ii ].x() = 16.962; xyzarray[ ii ].y() = -0.908; xyzarray[ ii ].z() = 3.33;
		xyzarray[ ++ii ].x() = 16.097; xyzarray[ ii ].y() = -1.125; xyzarray[ ii ].z() = 2.182;
		xyzarray[ ++ii ].x() = 16.869; xyzarray[ ii ].y() = -0.97; xyzarray[ ii ].z() = 0.85;
		xyzarray[ ++ii ].x() = 16.738; xyzarray[ ii ].y() = -1.809; xyzarray[ ii ].z() = -0.05;
		xyzarray[ ++ii ].x() = 17.665; xyzarray[ ii ].y() = 0.096; xyzarray[ ii ].z() = 0.736;
		xyzarray[ ++ii ].x() = 18.475; xyzarray[ ii ].y() = 0.343; xyzarray[ ii ].z() = -0.472;
		xyzarray[ ++ii ].x() = 19.518; xyzarray[ ii ].y() = -0.738; xyzarray[ ii ].z() = -0.716;
		xyzarray[ ++ii ].x() = 19.664; xyzarray[ ii ].y() = -1.207; xyzarray[ ii ].z() = -1.868;
		xyzarray[ ++ii ].x() = 20.228; xyzarray[ ii ].y() = -1.126; xyzarray[ ii ].z() = 0.35;
		xyzarray[ ++ii ].x() = 21.204; xyzarray[ ii ].y() = -2.22; xyzarray[ ii ].z() = 0.257;
		xyzarray[ ++ii ].x() = 20.536; xyzarray[ ii ].y() = -3.519; xyzarray[ ii ].z() = -0.16;
		xyzarray[ ++ii ].x() = 21.082; xyzarray[ ii ].y() = -4.24; xyzarray[ ii ].z() = -0.991;
		xyzarray[ ++ii ].x() = 19.359; xyzarray[ ii ].y() = -3.814; xyzarray[ ii ].z() = 0.397;
		xyzarray[ ++ii ].x() = 18.66; xyzarray[ ii ].y() = -5.056; xyzarray[ ii ].z() = 0.06;
		xyzarray[ ++ii ].x() = 18.223; xyzarray[ ii ].y() = -5.058; xyzarray[ ii ].z() = -1.4;
		xyzarray[ ++ii ].x() = 18.292; xyzarray[ ii ].y() = -6.113; xyzarray[ ii ].z() = -2.078;
		xyzarray[ ++ii ].x() = 17.735; xyzarray[ ii ].y() = -3.904; xyzarray[ ii ].z() = -1.867;
		xyzarray[ ++ii ].x() = 17.333; xyzarray[ ii ].y() = -3.773; xyzarray[ ii ].z() = -3.262;
		xyzarray[ ++ii ].x() = 18.524; xyzarray[ ii ].y() = -3.977; xyzarray[ ii ].z() = -4.22;
		xyzarray[ ++ii ].x() = 18.407; xyzarray[ ii ].y() = -4.699; xyzarray[ ii ].z() = -5.238;
		xyzarray[ ++ii ].x() = 19.653; xyzarray[ ii ].y() = -3.335; xyzarray[ ii ].z() = -3.909;
		xyzarray[ ++ii ].x() = 20.868; xyzarray[ ii ].y() = -3.529; xyzarray[ ii ].z() = -4.709;
		xyzarray[ ++ii ].x() = 21.208; xyzarray[ ii ].y() = -5.012; xyzarray[ ii ].z() = -4.713;
		xyzarray[ ++ii ].x() = 21.425; xyzarray[ ii ].y() = -5.616; xyzarray[ ii ].z() = -5.769;
		xyzarray[ ++ii ].x() = 21.216; xyzarray[ ii ].y() = -5.61; xyzarray[ ii ].z() = -3.534;
		xyzarray[ ++ii ].x() = 21.596; xyzarray[ ii ].y() = -7.016; xyzarray[ ii ].z() = -3.404;
		xyzarray[ ++ii ].x() = 20.683; xyzarray[ ii ].y() = -7.93; xyzarray[ ii ].z() = -4.205;
		xyzarray[ ++ii ].x() = 21.157; xyzarray[ ii ].y() = -8.813; xyzarray[ ii ].z() = -4.937;
		xyzarray[ ++ii ].x() = 19.369; xyzarray[ ii ].y() = -7.728; xyzarray[ ii ].z() = -4.096;
		xyzarray[ ++ii ].x() = 18.435; xyzarray[ ii ].y() = -8.563; xyzarray[ ii ].z() = -4.844;
		xyzarray[ ++ii ].x() = 18.706; xyzarray[ ii ].y() = -8.513; xyzarray[ ii ].z() = -6.37;
		xyzarray[ ++ii ].x() = 18.553; xyzarray[ ii ].y() = -9.516; xyzarray[ ii ].z() = -7.066;
		xyzarray[ ++ii ].x() = 29.401; xyzarray[ ii ].y() = -10.717; xyzarray[ ii ].z() = -8.443;
		xyzarray[ ++ii ].x() = 30.036; xyzarray[ ii ].y() = -10.857; xyzarray[ ii ].z() = -7.131;
		xyzarray[ ++ii ].x() = 30.897; xyzarray[ ii ].y() = -9.691; xyzarray[ ii ].z() = -6.666;
		xyzarray[ ++ii ].x() = 31.085; xyzarray[ ii ].y() = -9.528; xyzarray[ ii ].z() = -5.454;
		xyzarray[ ++ii ].x() = 31.403; xyzarray[ ii ].y() = -8.882; xyzarray[ ii ].z() = -7.598;
		xyzarray[ ++ii ].x() = 32.169; xyzarray[ ii ].y() = -7.697; xyzarray[ ii ].z() = -7.236;
		xyzarray[ ++ii ].x() = 31.258; xyzarray[ ii ].y() = -6.655; xyzarray[ ii ].z() = -6.582;
		xyzarray[ ++ii ].x() = 31.69; xyzarray[ ii ].y() = -5.925; xyzarray[ ii ].z() = -5.685;
		xyzarray[ ++ii ].x() = 29.99; xyzarray[ ii ].y() = -6.618; xyzarray[ ii ].z() = -7.006;
		xyzarray[ ++ii ].x() = 28.972; xyzarray[ ii ].y() = -5.752; xyzarray[ ii ].z() = -6.39;
		xyzarray[ ++ii ].x() = 28.459; xyzarray[ ii ].y() = -6.398; xyzarray[ ii ].z() = -5.106;
		xyzarray[ ++ii ].x() = 28.475; xyzarray[ ii ].y() = -5.787; xyzarray[ ii ].z() = -4.034;
		xyzarray[ ++ii ].x() = 28.028; xyzarray[ ii ].y() = -7.652; xyzarray[ ii ].z() = -5.186;
		xyzarray[ ++ii ].x() = 27.322; xyzarray[ ii ].y() = -8.219; xyzarray[ ii ].z() = -4.023;
		xyzarray[ ++ii ].x() = 28.205; xyzarray[ ii ].y() = -8.388; xyzarray[ ii ].z() = -2.782;
		xyzarray[ ++ii ].x() = 27.74; xyzarray[ ii ].y() = -8.178; xyzarray[ ii ].z() = -1.639;
		xyzarray[ ++ii ].x() = 29.479; xyzarray[ ii ].y() = -8.728; xyzarray[ ii ].z() = -2.988;
		xyzarray[ ++ii ].x() = 30.383; xyzarray[ ii ].y() = -8.922; xyzarray[ ii ].z() = -1.848;
		xyzarray[ ++ii ].x() = 30.641; xyzarray[ ii ].y() = -7.581; xyzarray[ ii ].z() = -1.108;
		xyzarray[ ++ii ].x() = 30.751; xyzarray[ ii ].y() = -7.55; xyzarray[ ii ].z() = 0.13;
		xyzarray[ ++ii ].x() = 30.721; xyzarray[ ii ].y() = -6.499; xyzarray[ ii ].z() = -1.866;
		xyzarray[ ++ii ].x() = 30.879; xyzarray[ ii ].y() = -5.153; xyzarray[ ii ].z() = -1.288;
		xyzarray[ ++ii ].x() = 29.636; xyzarray[ ii ].y() = -4.705; xyzarray[ ii ].z() = -0.541;
		xyzarray[ ++ii ].x() = 29.71; xyzarray[ ii ].y() = -4.153; xyzarray[ ii ].z() = 0.571;
		xyzarray[ ++ii ].x() = 28.488; xyzarray[ ii ].y() = -4.946; xyzarray[ ii ].z() = -1.138;
		xyzarray[ ++ii ].x() = 27.216; xyzarray[ ii ].y() = -4.597; xyzarray[ ii ].z() = -0.494;
		xyzarray[ ++ii ].x() = 27.001; xyzarray[ ii ].y() = -5.415; xyzarray[ ii ].z() = 0.811;
		xyzarray[ ++ii ].x() = 26.562; xyzarray[ ii ].y() = -4.871; xyzarray[ ii ].z() = 1.827;
		xyzarray[ ++ii ].x() = 27.335; xyzarray[ ii ].y() = -6.699; xyzarray[ ii ].z() = 0.782;
		xyzarray[ ++ii ].x() = 27.217; xyzarray[ ii ].y() = -7.556; xyzarray[ ii ].z() = 1.98;
		xyzarray[ ++ii ].x() = 28.107; xyzarray[ ii ].y() = -7.063; xyzarray[ ii ].z() = 3.109;
		xyzarray[ ++ii ].x() = 27.699; xyzarray[ ii ].y() = -7.09; xyzarray[ ii ].z() = 4.28;
		xyzarray[ ++ii ].x() = 29.315; xyzarray[ ii ].y() = -6.618; xyzarray[ ii ].z() = 2.76;
		xyzarray[ ++ii ].x() = 30.234; xyzarray[ ii ].y() = -6.013; xyzarray[ ii ].z() = 3.721;
		xyzarray[ ++ii ].x() = 29.65; xyzarray[ ii ].y() = -4.777; xyzarray[ ii ].z() = 4.365;
		xyzarray[ ++ii ].x() = 29.783; xyzarray[ ii ].y() = -4.609; xyzarray[ ii ].z() = 5.591;
		xyzarray[ ++ii ].x() = 29.014; xyzarray[ ii ].y() = -3.912; xyzarray[ ii ].z() = 3.559;
		xyzarray[ ++ii ].x() = 28.433; xyzarray[ ii ].y() = -2.686; xyzarray[ ii ].z() = 4.098;
		xyzarray[ ++ii ].x() = 27.231; xyzarray[ ii ].y() = -3.045; xyzarray[ ii ].z() = 5.008;
		xyzarray[ ++ii ].x() = 27.09; xyzarray[ ii ].y() = -2.522; xyzarray[ ii ].z() = 6.144;
		xyzarray[ ++ii ].x() = 26.386; xyzarray[ ii ].y() = -3.967; xyzarray[ ii ].z() = 4.537;
		xyzarray[ ++ii ].x() = 25.226; xyzarray[ ii ].y() = -4.368; xyzarray[ ii ].z() = 5.341;
		xyzarray[ ++ii ].x() = 25.651; xyzarray[ ii ].y() = -4.955; xyzarray[ ii ].z() = 6.698;
		xyzarray[ ++ii ].x() = 25.126; xyzarray[ ii ].y() = -4.599; xyzarray[ ii ].z() = 7.756;
		xyzarray[ ++ii ].x() = 26.615; xyzarray[ ii ].y() = -5.86; xyzarray[ ii ].z() = 6.663;
		xyzarray[ ++ii ].x() = 27.081; xyzarray[ ii ].y() = -6.515; xyzarray[ ii ].z() = 7.876;
		xyzarray[ ++ii ].x() = 27.695; xyzarray[ ii ].y() = -5.504; xyzarray[ ii ].z() = 8.82;
		xyzarray[ ++ii ].x() = 27.544; xyzarray[ ii ].y() = -5.615; xyzarray[ ii ].z() = 10.047;
		xyzarray[ ++ii ].x() = 28.411; xyzarray[ ii ].y() = -4.533; xyzarray[ ii ].z() = 8.264;
		xyzarray[ ++ii ].x() = 29.057; xyzarray[ ii ].y() = -3.511; xyzarray[ ii ].z() = 9.11;
		xyzarray[ ++ii ].x() = 28.041; xyzarray[ ii ].y() = -2.668; xyzarray[ ii ].z() = 9.858;
		xyzarray[ ++ii ].x() = 28.115; xyzarray[ ii ].y() = -2.494; xyzarray[ ii ].z() = 11.105;
		xyzarray[ ++ii ].x() = 27.075; xyzarray[ ii ].y() = -2.147; xyzarray[ ii ].z() = 9.105;
		xyzarray[ ++ii ].x() = 26.049; xyzarray[ ii ].y() = -1.302; xyzarray[ ii ].z() = 9.686;
		xyzarray[ ++ii ].x() = 25.188; xyzarray[ ii ].y() = -2.1; xyzarray[ ii ].z() = 10.678;
		xyzarray[ ++ii ].x() = 24.68; xyzarray[ ii ].y() = -1.559; xyzarray[ ii ].z() = 11.644;
		xyzarray[ ++ii ].x() = 29.401; xyzarray[ ii ].y() = -10.717; xyzarray[ ii ].z() = -8.443;
		xyzarray[ ++ii ].x() = 30.036; xyzarray[ ii ].y() = -10.857; xyzarray[ ii ].z() = -7.131;
		xyzarray[ ++ii ].x() = 30.897; xyzarray[ ii ].y() = -9.691; xyzarray[ ii ].z() = -6.666;
		xyzarray[ ++ii ].x() = 31.085; xyzarray[ ii ].y() = -9.528; xyzarray[ ii ].z() = -5.454;
		xyzarray[ ++ii ].x() = 31.403; xyzarray[ ii ].y() = -8.882; xyzarray[ ii ].z() = -7.598;
		xyzarray[ ++ii ].x() = 32.169; xyzarray[ ii ].y() = -7.697; xyzarray[ ii ].z() = -7.236;
		xyzarray[ ++ii ].x() = 31.258; xyzarray[ ii ].y() = -6.655; xyzarray[ ii ].z() = -6.582;
		xyzarray[ ++ii ].x() = 31.69; xyzarray[ ii ].y() = -5.925; xyzarray[ ii ].z() = -5.685;
		xyzarray[ ++ii ].x() = 29.99; xyzarray[ ii ].y() = -6.618; xyzarray[ ii ].z() = -7.006;
		xyzarray[ ++ii ].x() = 28.972; xyzarray[ ii ].y() = -5.752; xyzarray[ ii ].z() = -6.39;
		xyzarray[ ++ii ].x() = 28.459; xyzarray[ ii ].y() = -6.398; xyzarray[ ii ].z() = -5.106;
		xyzarray[ ++ii ].x() = 28.475; xyzarray[ ii ].y() = -5.787; xyzarray[ ii ].z() = -4.034;
		xyzarray[ ++ii ].x() = 28.028; xyzarray[ ii ].y() = -7.652; xyzarray[ ii ].z() = -5.186;
		xyzarray[ ++ii ].x() = 27.322; xyzarray[ ii ].y() = -8.219; xyzarray[ ii ].z() = -4.023;
		xyzarray[ ++ii ].x() = 28.205; xyzarray[ ii ].y() = -8.388; xyzarray[ ii ].z() = -2.782;
		xyzarray[ ++ii ].x() = 27.74; xyzarray[ ii ].y() = -8.178; xyzarray[ ii ].z() = -1.639;
		xyzarray[ ++ii ].x() = 29.479; xyzarray[ ii ].y() = -8.728; xyzarray[ ii ].z() = -2.988;
		xyzarray[ ++ii ].x() = 30.383; xyzarray[ ii ].y() = -8.922; xyzarray[ ii ].z() = -1.848;
		xyzarray[ ++ii ].x() = 30.641; xyzarray[ ii ].y() = -7.581; xyzarray[ ii ].z() = -1.108;
		xyzarray[ ++ii ].x() = 30.751; xyzarray[ ii ].y() = -7.55; xyzarray[ ii ].z() = 0.13;
		xyzarray[ ++ii ].x() = 30.721; xyzarray[ ii ].y() = -6.499; xyzarray[ ii ].z() = -1.866;
		xyzarray[ ++ii ].x() = 30.879; xyzarray[ ii ].y() = -5.153; xyzarray[ ii ].z() = -1.288;
		xyzarray[ ++ii ].x() = 29.636; xyzarray[ ii ].y() = -4.705; xyzarray[ ii ].z() = -0.541;
		xyzarray[ ++ii ].x() = 29.71; xyzarray[ ii ].y() = -4.153; xyzarray[ ii ].z() = 0.571;
		xyzarray[ ++ii ].x() = 28.488; xyzarray[ ii ].y() = -4.946; xyzarray[ ii ].z() = -1.138;
		xyzarray[ ++ii ].x() = 27.216; xyzarray[ ii ].y() = -4.597; xyzarray[ ii ].z() = -0.494;
		xyzarray[ ++ii ].x() = 27.001; xyzarray[ ii ].y() = -5.415; xyzarray[ ii ].z() = 0.811;
		xyzarray[ ++ii ].x() = 26.562; xyzarray[ ii ].y() = -4.871; xyzarray[ ii ].z() = 1.827;
		xyzarray[ ++ii ].x() = 27.335; xyzarray[ ii ].y() = -6.699; xyzarray[ ii ].z() = 0.782;
		xyzarray[ ++ii ].x() = 27.217; xyzarray[ ii ].y() = -7.556; xyzarray[ ii ].z() = 1.98;
		xyzarray[ ++ii ].x() = 28.107; xyzarray[ ii ].y() = -7.063; xyzarray[ ii ].z() = 3.109;
		xyzarray[ ++ii ].x() = 27.699; xyzarray[ ii ].y() = -7.09; xyzarray[ ii ].z() = 4.28;
		xyzarray[ ++ii ].x() = 29.315; xyzarray[ ii ].y() = -6.618; xyzarray[ ii ].z() = 2.76;
		xyzarray[ ++ii ].x() = 30.234; xyzarray[ ii ].y() = -6.013; xyzarray[ ii ].z() = 3.721;
		xyzarray[ ++ii ].x() = 29.65; xyzarray[ ii ].y() = -4.777; xyzarray[ ii ].z() = 4.365;
		xyzarray[ ++ii ].x() = 29.783; xyzarray[ ii ].y() = -4.609; xyzarray[ ii ].z() = 5.591;
		xyzarray[ ++ii ].x() = 29.014; xyzarray[ ii ].y() = -3.912; xyzarray[ ii ].z() = 3.559;
		xyzarray[ ++ii ].x() = 28.433; xyzarray[ ii ].y() = -2.686; xyzarray[ ii ].z() = 4.098;
		xyzarray[ ++ii ].x() = 27.231; xyzarray[ ii ].y() = -3.045; xyzarray[ ii ].z() = 5.008;
		xyzarray[ ++ii ].x() = 27.09; xyzarray[ ii ].y() = -2.522; xyzarray[ ii ].z() = 6.144;
		xyzarray[ ++ii ].x() = 26.386; xyzarray[ ii ].y() = -3.967; xyzarray[ ii ].z() = 4.537;
		xyzarray[ ++ii ].x() = 25.226; xyzarray[ ii ].y() = -4.368; xyzarray[ ii ].z() = 5.341;
		xyzarray[ ++ii ].x() = 25.651; xyzarray[ ii ].y() = -4.955; xyzarray[ ii ].z() = 6.698;
		xyzarray[ ++ii ].x() = 25.126; xyzarray[ ii ].y() = -4.599; xyzarray[ ii ].z() = 7.756;
		xyzarray[ ++ii ].x() = 26.615; xyzarray[ ii ].y() = -5.86; xyzarray[ ii ].z() = 6.663;
		xyzarray[ ++ii ].x() = 27.081; xyzarray[ ii ].y() = -6.515; xyzarray[ ii ].z() = 7.876;
		xyzarray[ ++ii ].x() = 27.695; xyzarray[ ii ].y() = -5.504; xyzarray[ ii ].z() = 8.82;
		xyzarray[ ++ii ].x() = 27.544; xyzarray[ ii ].y() = -5.615; xyzarray[ ii ].z() = 10.047;
		xyzarray[ ++ii ].x() = 28.411; xyzarray[ ii ].y() = -4.533; xyzarray[ ii ].z() = 8.264;
		xyzarray[ ++ii ].x() = 29.057; xyzarray[ ii ].y() = -3.511; xyzarray[ ii ].z() = 9.11;
		xyzarray[ ++ii ].x() = 28.041; xyzarray[ ii ].y() = -2.668; xyzarray[ ii ].z() = 9.858;
		xyzarray[ ++ii ].x() = 28.115; xyzarray[ ii ].y() = -2.494; xyzarray[ ii ].z() = 11.105;
		xyzarray[ ++ii ].x() = 27.075; xyzarray[ ii ].y() = -2.147; xyzarray[ ii ].z() = 9.105;
		xyzarray[ ++ii ].x() = 26.049; xyzarray[ ii ].y() = -1.302; xyzarray[ ii ].z() = 9.686;
		xyzarray[ ++ii ].x() = 25.188; xyzarray[ ii ].y() = -2.1; xyzarray[ ii ].z() = 10.678;
		xyzarray[ ++ii ].x() = 24.68; xyzarray[ ii ].y() = -1.559; xyzarray[ ii ].z() = 11.644;
		xyzarray[ ++ii ].x() = 21.454; xyzarray[ ii ].y() = -0.628; xyzarray[ ii ].z() = -14.376;
		xyzarray[ ++ii ].x() = 22.687; xyzarray[ ii ].y() = 0.067; xyzarray[ ii ].z() = -14.016;
		xyzarray[ ++ii ].x() = 22.354; xyzarray[ ii ].y() = 1.443; xyzarray[ ii ].z() = -13.394;
		xyzarray[ ++ii ].x() = 22.956; xyzarray[ ii ].y() = 1.856; xyzarray[ ii ].z() = -12.376;
		xyzarray[ ++ii ].x() = 21.4; xyzarray[ ii ].y() = 2.132; xyzarray[ ii ].z() = -14.014;
		xyzarray[ ++ii ].x() = 20.965; xyzarray[ ii ].y() = 3.469; xyzarray[ ii ].z() = -13.587;
		xyzarray[ ++ii ].x() = 20.305; xyzarray[ ii ].y() = 3.427; xyzarray[ ii ].z() = -12.241;
		xyzarray[ ++ii ].x() = 20.612; xyzarray[ ii ].y() = 4.236; xyzarray[ ii ].z() = -11.368;
		xyzarray[ ++ii ].x() = 19.371; xyzarray[ ii ].y() = 2.495; xyzarray[ ii ].z() = -12.083;
		xyzarray[ ++ii ].x() = 18.688; xyzarray[ ii ].y() = 2.32; xyzarray[ ii ].z() = -10.809;
		xyzarray[ ++ii ].x() = 19.671; xyzarray[ ii ].y() = 2.061; xyzarray[ ii ].z() = -9.666;
		xyzarray[ ++ii ].x() = 19.638; xyzarray[ ii ].y() = 2.736; xyzarray[ ii ].z() = -8.639;
		xyzarray[ ++ii ].x() = 20.565; xyzarray[ ii ].y() = 1.11; xyzarray[ ii ].z() = -9.867;
		xyzarray[ ++ii ].x() = 21.519; xyzarray[ ii ].y() = 0.753; xyzarray[ ii ].z() = -8.84;
		xyzarray[ ++ii ].x() = 22.493; xyzarray[ ii ].y() = 1.877; xyzarray[ ii ].z() = -8.551;
		xyzarray[ ++ii ].x() = 22.862; xyzarray[ ii ].y() = 2.08; xyzarray[ ii ].z() = -7.397;
		xyzarray[ ++ii ].x() = 22.886; xyzarray[ ii ].y() = 2.615; xyzarray[ ii ].z() = -9.588;
		xyzarray[ ++ii ].x() = 23.818; xyzarray[ ii ].y() = 3.736; xyzarray[ ii ].z() = -9.399;
		xyzarray[ ++ii ].x() = 23.163; xyzarray[ ii ].y() = 4.785; xyzarray[ ii ].z() = -8.479;
		xyzarray[ ++ii ].x() = 23.824; xyzarray[ ii ].y() = 5.345; xyzarray[ ii ].z() = -7.589;
		xyzarray[ ++ii ].x() = 21.871; xyzarray[ ii ].y() = 5.039; xyzarray[ ii ].z() = -8.685;
		xyzarray[ ++ii ].x() = 21.154; xyzarray[ ii ].y() = 5.96; xyzarray[ ii ].z() = -7.83;
		xyzarray[ ++ii ].x() = 20.987; xyzarray[ ii ].y() = 5.458; xyzarray[ ii ].z() = -6.391;
		xyzarray[ ++ii ].x() = 21.087; xyzarray[ ii ].y() = 6.243; xyzarray[ ii ].z() = -5.449;
		xyzarray[ ++ii ].x() = 20.723; xyzarray[ ii ].y() = 4.164; xyzarray[ ii ].z() = -6.223;
		xyzarray[ ++ii ].x() = 20.626; xyzarray[ ii ].y() = 3.607; xyzarray[ ii ].z() = -4.864;
		xyzarray[ ++ii ].x() = 21.973; xyzarray[ ii ].y() = 3.714; xyzarray[ ii ].z() = -4.154;
		xyzarray[ ++ii ].x() = 22.029; xyzarray[ ii ].y() = 4.074; xyzarray[ ii ].z() = -2.977;
		xyzarray[ ++ii ].x() = 23.049; xyzarray[ ii ].y() = 3.395; xyzarray[ ii ].z() = -4.867;
		xyzarray[ ++ii ].x() = 24.4; xyzarray[ ii ].y() = 3.591; xyzarray[ ii ].z() = -4.306;
		xyzarray[ ++ii ].x() = 24.644; xyzarray[ ii ].y() = 5.027; xyzarray[ ii ].z() = -3.883;
		xyzarray[ ++ii ].x() = 25.157; xyzarray[ ii ].y() = 5.281; xyzarray[ ii ].z() = -2.77;
		xyzarray[ ++ii ].x() = 24.251; xyzarray[ ii ].y() = 5.97; xyzarray[ ii ].z() = -4.739;
		xyzarray[ ++ii ].x() = 24.398; xyzarray[ ii ].y() = 7.391; xyzarray[ ii ].z() = -4.404;
		xyzarray[ ++ii ].x() = 23.609; xyzarray[ ii ].y() = 7.745; xyzarray[ ii ].z() = -3.122;
		xyzarray[ ++ii ].x() = 24.102; xyzarray[ ii ].y() = 8.485; xyzarray[ ii ].z() = -2.272;
		xyzarray[ ++ii ].x() = 22.395; xyzarray[ ii ].y() = 7.206; xyzarray[ ii ].z() = -2.998;
		xyzarray[ ++ii ].x() = 21.578; xyzarray[ ii ].y() = 7.423; xyzarray[ ii ].z() = -1.813;
		xyzarray[ ++ii ].x() = 22.308; xyzarray[ ii ].y() = 6.933; xyzarray[ ii ].z() = -0.563;
		xyzarray[ ++ii ].x() = 22.333; xyzarray[ ii ].y() = 7.635; xyzarray[ ii ].z() = 0.464;
		xyzarray[ ++ii ].x() = 22.901; xyzarray[ ii ].y() = 5.738; xyzarray[ ii ].z() = -0.653;
		xyzarray[ ++ii ].x() = 23.595; xyzarray[ ii ].y() = 5.163; xyzarray[ ii ].z() = 0.501;
		xyzarray[ ++ii ].x() = 24.824; xyzarray[ ii ].y() = 6; xyzarray[ ii ].z() = 0.875;
		xyzarray[ ++ii ].x() = 25.02; xyzarray[ ii ].y() = 6.353; xyzarray[ ii ].z() = 2.051;
		xyzarray[ ++ii ].x() = 25.639; xyzarray[ ii ].y() = 6.317; xyzarray[ ii ].z() = -0.127;
		xyzarray[ ++ii ].x() = 26.867; xyzarray[ ii ].y() = 7.07; xyzarray[ ii ].z() = 0.082;
		xyzarray[ ++ii ].x() = 26.595; xyzarray[ ii ].y() = 8.452; xyzarray[ ii ].z() = 0.682;
		xyzarray[ ++ii ].x() = 27.216; xyzarray[ ii ].y() = 8.83; xyzarray[ ii ].z() = 1.684;
		xyzarray[ ++ii ].x() = 25.656; xyzarray[ ii ].y() = 9.185; xyzarray[ ii ].z() = 0.079;
		xyzarray[ ++ii ].x() = 25.231; xyzarray[ ii ].y() = 10.521; xyzarray[ ii ].z() = 0.563;
		xyzarray[ ++ii ].x() = 24.746; xyzarray[ ii ].y() = 10.459; xyzarray[ ii ].z() = 2.037;
		xyzarray[ ++ii ].x() = 25.06; xyzarray[ ii ].y() = 11.345; xyzarray[ ii ].z() = 2.836;
		xyzarray[ ++ii ].x() = 24.003; xyzarray[ ii ].y() = 9.406; xyzarray[ ii ].z() = 2.392;
		xyzarray[ ++ii ].x() = 23.475; xyzarray[ ii ].y() = 9.235; xyzarray[ ii ].z() = 3.757;
		xyzarray[ ++ii ].x() = 24.584; xyzarray[ ii ].y() = 8.937; xyzarray[ ii ].z() = 4.776;
		xyzarray[ ++ii ].x() = 24.633; xyzarray[ ii ].y() = 9.529; xyzarray[ ii ].z() = 5.859;
		return xyzarray;
	}

	utility::vector1< numeric::xyzVector< float > >
	collision_coordinates_for_two_helix_pairs() {

		utility::vector1< numeric::xyzVector< float > > xyzarray( 224 );
		core::Size ii = 0;
		// coordinates taken directly from Tim's helical-assembly test case
		xyzarray[ ++ii ].x() = 18.78; xyzarray[ ii ].y() = -1.096; xyzarray[ ii ].z() = -15.207;
		xyzarray[ ++ii ].x() = 19.242; xyzarray[ ii ].y() = -1.605; xyzarray[ ii ].z() = -13.91;
		xyzarray[ ++ii ].x() = 20.529; xyzarray[ ii ].y() = -0.887; xyzarray[ ii ].z() = -13.448;
		xyzarray[ ++ii ].x() = 20.695; xyzarray[ ii ].y() = -0.604; xyzarray[ ii ].z() = -12.25;
		xyzarray[ ++ii ].x() = 21.454; xyzarray[ ii ].y() = -0.628; xyzarray[ ii ].z() = -14.376;
		xyzarray[ ++ii ].x() = 22.687; xyzarray[ ii ].y() = 0.067; xyzarray[ ii ].z() = -14.016;
		xyzarray[ ++ii ].x() = 22.354; xyzarray[ ii ].y() = 1.443; xyzarray[ ii ].z() = -13.394;
		xyzarray[ ++ii ].x() = 22.956; xyzarray[ ii ].y() = 1.856; xyzarray[ ii ].z() = -12.376;
		xyzarray[ ++ii ].x() = 21.4; xyzarray[ ii ].y() = 2.132; xyzarray[ ii ].z() = -14.014;
		xyzarray[ ++ii ].x() = 20.965; xyzarray[ ii ].y() = 3.469; xyzarray[ ii ].z() = -13.587;
		xyzarray[ ++ii ].x() = 20.305; xyzarray[ ii ].y() = 3.427; xyzarray[ ii ].z() = -12.241;
		xyzarray[ ++ii ].x() = 20.612; xyzarray[ ii ].y() = 4.236; xyzarray[ ii ].z() = -11.368;
		xyzarray[ ++ii ].x() = 19.371; xyzarray[ ii ].y() = 2.495; xyzarray[ ii ].z() = -12.083;
		xyzarray[ ++ii ].x() = 18.688; xyzarray[ ii ].y() = 2.32; xyzarray[ ii ].z() = -10.809;
		xyzarray[ ++ii ].x() = 19.671; xyzarray[ ii ].y() = 2.061; xyzarray[ ii ].z() = -9.666;
		xyzarray[ ++ii ].x() = 19.638; xyzarray[ ii ].y() = 2.736; xyzarray[ ii ].z() = -8.639;
		xyzarray[ ++ii ].x() = 20.565; xyzarray[ ii ].y() = 1.11; xyzarray[ ii ].z() = -9.867;
		xyzarray[ ++ii ].x() = 21.519; xyzarray[ ii ].y() = 0.753; xyzarray[ ii ].z() = -8.84;
		xyzarray[ ++ii ].x() = 22.493; xyzarray[ ii ].y() = 1.877; xyzarray[ ii ].z() = -8.551;
		xyzarray[ ++ii ].x() = 22.862; xyzarray[ ii ].y() = 2.08; xyzarray[ ii ].z() = -7.397;
		xyzarray[ ++ii ].x() = 22.886; xyzarray[ ii ].y() = 2.615; xyzarray[ ii ].z() = -9.588;
		xyzarray[ ++ii ].x() = 23.818; xyzarray[ ii ].y() = 3.736; xyzarray[ ii ].z() = -9.399;
		xyzarray[ ++ii ].x() = 23.163; xyzarray[ ii ].y() = 4.785; xyzarray[ ii ].z() = -8.479;
		xyzarray[ ++ii ].x() = 23.824; xyzarray[ ii ].y() = 5.345; xyzarray[ ii ].z() = -7.589;
		xyzarray[ ++ii ].x() = 21.871; xyzarray[ ii ].y() = 5.039; xyzarray[ ii ].z() = -8.685;
		xyzarray[ ++ii ].x() = 21.154; xyzarray[ ii ].y() = 5.96; xyzarray[ ii ].z() = -7.83;
		xyzarray[ ++ii ].x() = 20.987; xyzarray[ ii ].y() = 5.458; xyzarray[ ii ].z() = -6.391;
		xyzarray[ ++ii ].x() = 21.087; xyzarray[ ii ].y() = 6.243; xyzarray[ ii ].z() = -5.449;
		xyzarray[ ++ii ].x() = 20.723; xyzarray[ ii ].y() = 4.164; xyzarray[ ii ].z() = -6.223;
		xyzarray[ ++ii ].x() = 20.626; xyzarray[ ii ].y() = 3.607; xyzarray[ ii ].z() = -4.864;
		xyzarray[ ++ii ].x() = 21.973; xyzarray[ ii ].y() = 3.714; xyzarray[ ii ].z() = -4.154;
		xyzarray[ ++ii ].x() = 22.029; xyzarray[ ii ].y() = 4.074; xyzarray[ ii ].z() = -2.977;
		xyzarray[ ++ii ].x() = 23.049; xyzarray[ ii ].y() = 3.395; xyzarray[ ii ].z() = -4.867;
		xyzarray[ ++ii ].x() = 24.4; xyzarray[ ii ].y() = 3.591; xyzarray[ ii ].z() = -4.306;
		xyzarray[ ++ii ].x() = 24.644; xyzarray[ ii ].y() = 5.027; xyzarray[ ii ].z() = -3.883;
		xyzarray[ ++ii ].x() = 25.157; xyzarray[ ii ].y() = 5.281; xyzarray[ ii ].z() = -2.77;
		xyzarray[ ++ii ].x() = 24.251; xyzarray[ ii ].y() = 5.97; xyzarray[ ii ].z() = -4.739;
		xyzarray[ ++ii ].x() = 24.398; xyzarray[ ii ].y() = 7.391; xyzarray[ ii ].z() = -4.404;
		xyzarray[ ++ii ].x() = 23.609; xyzarray[ ii ].y() = 7.745; xyzarray[ ii ].z() = -3.122;
		xyzarray[ ++ii ].x() = 24.102; xyzarray[ ii ].y() = 8.485; xyzarray[ ii ].z() = -2.272;
		xyzarray[ ++ii ].x() = 22.395; xyzarray[ ii ].y() = 7.206; xyzarray[ ii ].z() = -2.998;
		xyzarray[ ++ii ].x() = 21.578; xyzarray[ ii ].y() = 7.423; xyzarray[ ii ].z() = -1.813;
		xyzarray[ ++ii ].x() = 22.308; xyzarray[ ii ].y() = 6.933; xyzarray[ ii ].z() = -0.563;
		xyzarray[ ++ii ].x() = 22.333; xyzarray[ ii ].y() = 7.635; xyzarray[ ii ].z() = 0.464;
		xyzarray[ ++ii ].x() = 22.901; xyzarray[ ii ].y() = 5.738; xyzarray[ ii ].z() = -0.653;
		xyzarray[ ++ii ].x() = 23.595; xyzarray[ ii ].y() = 5.163; xyzarray[ ii ].z() = 0.501;
		xyzarray[ ++ii ].x() = 24.824; xyzarray[ ii ].y() = 6; xyzarray[ ii ].z() = 0.875;
		xyzarray[ ++ii ].x() = 25.02; xyzarray[ ii ].y() = 6.353; xyzarray[ ii ].z() = 2.051;
		xyzarray[ ++ii ].x() = 25.639; xyzarray[ ii ].y() = 6.317; xyzarray[ ii ].z() = -0.127;
		xyzarray[ ++ii ].x() = 26.867; xyzarray[ ii ].y() = 7.07; xyzarray[ ii ].z() = 0.082;
		xyzarray[ ++ii ].x() = 26.595; xyzarray[ ii ].y() = 8.452; xyzarray[ ii ].z() = 0.682;
		xyzarray[ ++ii ].x() = 27.216; xyzarray[ ii ].y() = 8.83; xyzarray[ ii ].z() = 1.684;
		xyzarray[ ++ii ].x() = 25.656; xyzarray[ ii ].y() = 9.185; xyzarray[ ii ].z() = 0.079;
		xyzarray[ ++ii ].x() = 25.231; xyzarray[ ii ].y() = 10.521; xyzarray[ ii ].z() = 0.563;
		xyzarray[ ++ii ].x() = 24.746; xyzarray[ ii ].y() = 10.459; xyzarray[ ii ].z() = 2.037;
		xyzarray[ ++ii ].x() = 25.06; xyzarray[ ii ].y() = 11.345; xyzarray[ ii ].z() = 2.836;
		xyzarray[ ++ii ].x() = 34.211; xyzarray[ ii ].y() = 4.736; xyzarray[ ii ].z() = 8.755;
		xyzarray[ ++ii ].x() = 33.305; xyzarray[ ii ].y() = 4.996; xyzarray[ ii ].z() = 7.626;
		xyzarray[ ++ii ].x() = 34.037; xyzarray[ ii ].y() = 4.914; xyzarray[ ii ].z() = 6.279;
		xyzarray[ ++ii ].x() = 33.425; xyzarray[ ii ].y() = 4.569; xyzarray[ ii ].z() = 5.269;
		xyzarray[ ++ii ].x() = 35.346; xyzarray[ ii ].y() = 5.166; xyzarray[ ii ].z() = 6.26;
		xyzarray[ ++ii ].x() = 36.072; xyzarray[ ii ].y() = 5.148; xyzarray[ ii ].z() = 4.98;
		xyzarray[ ++ii ].x() = 36.024; xyzarray[ ii ].y() = 3.765; xyzarray[ ii ].z() = 4.323;
		xyzarray[ ++ii ].x() = 35.895; xyzarray[ ii ].y() = 3.652; xyzarray[ ii ].z() = 3.095;
		xyzarray[ ++ii ].x() = 36.072; xyzarray[ ii ].y() = 2.717; xyzarray[ ii ].z() = 5.143;
		xyzarray[ ++ii ].x() = 35.925; xyzarray[ ii ].y() = 1.365; xyzarray[ ii ].z() = 4.617;
		xyzarray[ ++ii ].x() = 34.641; xyzarray[ ii ].y() = 1.198; xyzarray[ ii ].z() = 3.788;
		xyzarray[ ++ii ].x() = 34.657; xyzarray[ ii ].y() = 0.582; xyzarray[ ii ].z() = 2.691;
		xyzarray[ ++ii ].x() = 33.532; xyzarray[ ii ].y() = 1.721; xyzarray[ ii ].z() = 4.327;
		xyzarray[ ++ii ].x() = 32.238; xyzarray[ ii ].y() = 1.679; xyzarray[ ii ].z() = 3.648;
		xyzarray[ ++ii ].x() = 32.207; xyzarray[ ii ].y() = 2.482; xyzarray[ ii ].z() = 2.355;
		xyzarray[ ++ii ].x() = 31.659; xyzarray[ ii ].y() = 2.003; xyzarray[ ii ].z() = 1.339;
		xyzarray[ ++ii ].x() = 32.79; xyzarray[ ii ].y() = 3.684; xyzarray[ ii ].z() = 2.394;
		xyzarray[ ++ii ].x() = 32.781; xyzarray[ ii ].y() = 4.587; xyzarray[ ii ].z() = 1.247;
		xyzarray[ ++ii ].x() = 33.631; xyzarray[ ii ].y() = 4.007; xyzarray[ ii ].z() = 0.119;
		xyzarray[ ++ii ].x() = 33.276; xyzarray[ ii ].y() = 4.127; xyzarray[ ii ].z() = -1.06;
		xyzarray[ ++ii ].x() = 34.744; xyzarray[ ii ].y() = 3.359; xyzarray[ ii ].z() = 0.487;
		xyzarray[ ++ii ].x() = 35.528; xyzarray[ ii ].y() = 2.641; xyzarray[ ii ].z() = -0.503;
		xyzarray[ ++ii ].x() = 34.775; xyzarray[ ii ].y() = 1.477; xyzarray[ ii ].z() = -1.154;
		xyzarray[ ++ii ].x() = 34.762; xyzarray[ ii ].y() = 1.368; xyzarray[ ii ].z() = -2.384;
		xyzarray[ ++ii ].x() = 34.157; xyzarray[ ii ].y() = 0.617; xyzarray[ ii ].z() = -0.339;
		xyzarray[ ++ii ].x() = 33.395; xyzarray[ ii ].y() = -0.5; xyzarray[ ii ].z() = -0.869;
		xyzarray[ ++ii ].x() = 32.315; xyzarray[ ii ].y() = 0.014; xyzarray[ ii ].z() = -1.829;
		xyzarray[ ++ii ].x() = 32.101; xyzarray[ ii ].y() = -0.587; xyzarray[ ii ].z() = -2.914;
		xyzarray[ ++ii ].x() = 31.68; xyzarray[ ii ].y() = 1.127; xyzarray[ ii ].z() = -1.461;
		xyzarray[ ++ii ].x() = 30.572; xyzarray[ ii ].y() = 1.665; xyzarray[ ii ].z() = -2.285;
		xyzarray[ ++ii ].x() = 31.118; xyzarray[ ii ].y() = 2.266; xyzarray[ ii ].z() = -3.564;
		xyzarray[ ++ii ].x() = 30.494; xyzarray[ ii ].y() = 2.15; xyzarray[ ii ].z() = -4.621;
		xyzarray[ ++ii ].x() = 32.313; xyzarray[ ii ].y() = 2.852; xyzarray[ ii ].z() = -3.496;
		xyzarray[ ++ii ].x() = 32.965; xyzarray[ ii ].y() = 3.336; xyzarray[ ii ].z() = -4.719;
		xyzarray[ ++ii ].x() = 33.302; xyzarray[ ii ].y() = 2.146; xyzarray[ ii ].z() = -5.656;
		xyzarray[ ++ii ].x() = 33.153; xyzarray[ ii ].y() = 2.227; xyzarray[ ii ].z() = -6.893;
		xyzarray[ ++ii ].x() = 33.756; xyzarray[ ii ].y() = 1.035; xyzarray[ ii ].z() = -5.078;
		xyzarray[ ++ii ].x() = 33.989; xyzarray[ ii ].y() = -0.179; xyzarray[ ii ].z() = -5.866;
		xyzarray[ ++ii ].x() = 32.742; xyzarray[ ii ].y() = -0.709; xyzarray[ ii ].z() = -6.553;
		xyzarray[ ++ii ].x() = 32.811; xyzarray[ ii ].y() = -1.185; xyzarray[ ii ].z() = -7.71;
		xyzarray[ ++ii ].x() = 31.598; xyzarray[ ii ].y() = -0.614; xyzarray[ ii ].z() = -5.882;
		xyzarray[ ++ii ].x() = 30.335; xyzarray[ ii ].y() = -0.963; xyzarray[ ii ].z() = -6.545;
		xyzarray[ ++ii ].x() = 30.138; xyzarray[ ii ].y() = -0.109; xyzarray[ ii ].z() = -7.798;
		xyzarray[ ++ii ].x() = 29.796; xyzarray[ ii ].y() = -0.623; xyzarray[ ii ].z() = -8.884;
		xyzarray[ ++ii ].x() = 30.382; xyzarray[ ii ].y() = 1.19; xyzarray[ ii ].z() = -7.67;
		xyzarray[ ++ii ].x() = 30.304; xyzarray[ ii ].y() = 2.088; xyzarray[ ii ].z() = -8.828;
		xyzarray[ ++ii ].x() = 31.296; xyzarray[ ii ].y() = 1.748; xyzarray[ ii ].z() = -9.934;
		xyzarray[ ++ii ].x() = 30.95; xyzarray[ ii ].y() = 1.83; xyzarray[ ii ].z() = -11.133;
		xyzarray[ ++ii ].x() = 32.521; xyzarray[ ii ].y() = 1.368; xyzarray[ ii ].z() = -9.562;
		xyzarray[ ++ii ].x() = 33.52; xyzarray[ ii ].y() = 0.894; xyzarray[ ii ].z() = -10.538;
		xyzarray[ ++ii ].x() = 32.955; xyzarray[ ii ].y() = -0.336; xyzarray[ ii ].z() = -11.275;
		xyzarray[ ++ii ].x() = 33.116; xyzarray[ ii ].y() = -0.479; xyzarray[ ii ].z() = -12.481;
		xyzarray[ ++ii ].x() = 32.299; xyzarray[ ii ].y() = -1.228; xyzarray[ ii ].z() = -10.548;
		xyzarray[ ++ii ].x() = 31.702; xyzarray[ ii ].y() = -2.417; xyzarray[ ii ].z() = -11.196;
		xyzarray[ ++ii ].x() = 30.489; xyzarray[ ii ].y() = -2.092; xyzarray[ ii ].z() = -12.072;
		xyzarray[ ++ii ].x() = 30.35; xyzarray[ ii ].y() = -2.625; xyzarray[ ii ].z() = -13.182;
		xyzarray[ ++ii ].x() = 16.614; xyzarray[ ii ].y() = 6.327; xyzarray[ ii ].z() = 9.328;
		xyzarray[ ++ii ].x() = 17.629; xyzarray[ ii ].y() = 6.362; xyzarray[ ii ].z() = 8.277;
		xyzarray[ ++ii ].x() = 18.32; xyzarray[ ii ].y() = 5.009; xyzarray[ ii ].z() = 8.073;
		xyzarray[ ++ii ].x() = 18.553; xyzarray[ ii ].y() = 4.59; xyzarray[ ii ].z() = 6.929;
		xyzarray[ ++ii ].x() = 18.654; xyzarray[ ii ].y() = 4.321; xyzarray[ ii ].z() = 9.171;
		xyzarray[ ++ii ].x() = 19.157; xyzarray[ ii ].y() = 2.956; xyzarray[ ii ].z() = 9.07;
		xyzarray[ ++ii ].x() = 18.185; xyzarray[ ii ].y() = 2.097; xyzarray[ ii ].z() = 8.257;
		xyzarray[ ++ii ].x() = 18.59; xyzarray[ ii ].y() = 1.291; xyzarray[ ii ].z() = 7.399;
		xyzarray[ ++ii ].x() = 16.894; xyzarray[ ii ].y() = 2.26; xyzarray[ ii ].z() = 8.537;
		xyzarray[ ++ii ].x() = 15.908; xyzarray[ ii ].y() = 1.449; xyzarray[ ii ].z() = 7.858;
		xyzarray[ ++ii ].x() = 15.918; xyzarray[ ii ].y() = 1.707; xyzarray[ ii ].z() = 6.356;
		xyzarray[ ++ii ].x() = 15.878; xyzarray[ ii ].y() = 0.758; xyzarray[ ii ].z() = 5.554;
		xyzarray[ ++ii ].x() = 16.01; xyzarray[ ii ].y() = 2.984; xyzarray[ ii ].z() = 5.977;
		xyzarray[ ++ii ].x() = 16.098; xyzarray[ ii ].y() = 3.349; xyzarray[ ii ].z() = 4.559;
		xyzarray[ ++ii ].x() = 17.351; xyzarray[ ii ].y() = 2.776; xyzarray[ ii ].z() = 3.881;
		xyzarray[ ++ii ].x() = 17.292; xyzarray[ ii ].y() = 2.279; xyzarray[ ii ].z() = 2.745;
		xyzarray[ ++ii ].x() = 18.485; xyzarray[ ii ].y() = 2.873; xyzarray[ ii ].z() = 4.568;
		xyzarray[ ++ii ].x() = 19.766; xyzarray[ ii ].y() = 2.393; xyzarray[ ii ].z() = 4.02;
		xyzarray[ ++ii ].x() = 19.684; xyzarray[ ii ].y() = 0.873; xyzarray[ ii ].z() = 3.871;
		xyzarray[ ++ii ].x() = 20.057; xyzarray[ ii ].y() = 0.328; xyzarray[ ii ].z() = 2.829;
		xyzarray[ ++ii ].x() = 19.121; xyzarray[ ii ].y() = 0.207; xyzarray[ ii ].z() = 4.879;
		xyzarray[ ++ii ].x() = 18.844; xyzarray[ ii ].y() = -1.248; xyzarray[ ii ].z() = 4.795;
		xyzarray[ ++ii ].x() = 18.026; xyzarray[ ii ].y() = -1.658; xyzarray[ ii ].z() = 3.583;
		xyzarray[ ++ii ].x() = 18.348; xyzarray[ ii ].y() = -2.64; xyzarray[ ii ].z() = 2.892;
		xyzarray[ ++ii ].x() = 16.962; xyzarray[ ii ].y() = -0.908; xyzarray[ ii ].z() = 3.33;
		xyzarray[ ++ii ].x() = 16.097; xyzarray[ ii ].y() = -1.125; xyzarray[ ii ].z() = 2.182;
		xyzarray[ ++ii ].x() = 16.869; xyzarray[ ii ].y() = -0.97; xyzarray[ ii ].z() = 0.85;
		xyzarray[ ++ii ].x() = 16.738; xyzarray[ ii ].y() = -1.809; xyzarray[ ii ].z() = -0.05;
		xyzarray[ ++ii ].x() = 17.665; xyzarray[ ii ].y() = 0.096; xyzarray[ ii ].z() = 0.736;
		xyzarray[ ++ii ].x() = 18.475; xyzarray[ ii ].y() = 0.343; xyzarray[ ii ].z() = -0.472;
		xyzarray[ ++ii ].x() = 19.518; xyzarray[ ii ].y() = -0.738; xyzarray[ ii ].z() = -0.716;
		xyzarray[ ++ii ].x() = 19.664; xyzarray[ ii ].y() = -1.207; xyzarray[ ii ].z() = -1.868;
		xyzarray[ ++ii ].x() = 20.228; xyzarray[ ii ].y() = -1.126; xyzarray[ ii ].z() = 0.35;
		xyzarray[ ++ii ].x() = 21.204; xyzarray[ ii ].y() = -2.22; xyzarray[ ii ].z() = 0.257;
		xyzarray[ ++ii ].x() = 20.536; xyzarray[ ii ].y() = -3.519; xyzarray[ ii ].z() = -0.16;
		xyzarray[ ++ii ].x() = 21.082; xyzarray[ ii ].y() = -4.24; xyzarray[ ii ].z() = -0.991;
		xyzarray[ ++ii ].x() = 19.359; xyzarray[ ii ].y() = -3.814; xyzarray[ ii ].z() = 0.397;
		xyzarray[ ++ii ].x() = 18.66; xyzarray[ ii ].y() = -5.056; xyzarray[ ii ].z() = 0.06;
		xyzarray[ ++ii ].x() = 18.223; xyzarray[ ii ].y() = -5.058; xyzarray[ ii ].z() = -1.4;
		xyzarray[ ++ii ].x() = 18.292; xyzarray[ ii ].y() = -6.113; xyzarray[ ii ].z() = -2.078;
		xyzarray[ ++ii ].x() = 17.735; xyzarray[ ii ].y() = -3.904; xyzarray[ ii ].z() = -1.867;
		xyzarray[ ++ii ].x() = 17.333; xyzarray[ ii ].y() = -3.773; xyzarray[ ii ].z() = -3.262;
		xyzarray[ ++ii ].x() = 18.524; xyzarray[ ii ].y() = -3.977; xyzarray[ ii ].z() = -4.22;
		xyzarray[ ++ii ].x() = 18.407; xyzarray[ ii ].y() = -4.699; xyzarray[ ii ].z() = -5.238;
		xyzarray[ ++ii ].x() = 19.653; xyzarray[ ii ].y() = -3.335; xyzarray[ ii ].z() = -3.909;
		xyzarray[ ++ii ].x() = 20.868; xyzarray[ ii ].y() = -3.529; xyzarray[ ii ].z() = -4.709;
		xyzarray[ ++ii ].x() = 21.208; xyzarray[ ii ].y() = -5.012; xyzarray[ ii ].z() = -4.713;
		xyzarray[ ++ii ].x() = 21.425; xyzarray[ ii ].y() = -5.616; xyzarray[ ii ].z() = -5.769;
		xyzarray[ ++ii ].x() = 21.216; xyzarray[ ii ].y() = -5.61; xyzarray[ ii ].z() = -3.534;
		xyzarray[ ++ii ].x() = 21.596; xyzarray[ ii ].y() = -7.016; xyzarray[ ii ].z() = -3.404;
		xyzarray[ ++ii ].x() = 20.683; xyzarray[ ii ].y() = -7.93; xyzarray[ ii ].z() = -4.205;
		xyzarray[ ++ii ].x() = 21.157; xyzarray[ ii ].y() = -8.813; xyzarray[ ii ].z() = -4.937;
		xyzarray[ ++ii ].x() = 19.369; xyzarray[ ii ].y() = -7.728; xyzarray[ ii ].z() = -4.096;
		xyzarray[ ++ii ].x() = 18.435; xyzarray[ ii ].y() = -8.563; xyzarray[ ii ].z() = -4.844;
		xyzarray[ ++ii ].x() = 18.706; xyzarray[ ii ].y() = -8.513; xyzarray[ ii ].z() = -6.37;
		xyzarray[ ++ii ].x() = 18.553; xyzarray[ ii ].y() = -9.516; xyzarray[ ii ].z() = -7.066;
		xyzarray[ ++ii ].x() = 34.211; xyzarray[ ii ].y() = 4.736; xyzarray[ ii ].z() = 8.755;
		xyzarray[ ++ii ].x() = 33.305; xyzarray[ ii ].y() = 4.996; xyzarray[ ii ].z() = 7.626;
		xyzarray[ ++ii ].x() = 34.037; xyzarray[ ii ].y() = 4.914; xyzarray[ ii ].z() = 6.279;
		xyzarray[ ++ii ].x() = 33.425; xyzarray[ ii ].y() = 4.569; xyzarray[ ii ].z() = 5.269;
		xyzarray[ ++ii ].x() = 35.346; xyzarray[ ii ].y() = 5.166; xyzarray[ ii ].z() = 6.26;
		xyzarray[ ++ii ].x() = 36.072; xyzarray[ ii ].y() = 5.148; xyzarray[ ii ].z() = 4.98;
		xyzarray[ ++ii ].x() = 36.024; xyzarray[ ii ].y() = 3.765; xyzarray[ ii ].z() = 4.323;
		xyzarray[ ++ii ].x() = 35.895; xyzarray[ ii ].y() = 3.652; xyzarray[ ii ].z() = 3.095;
		xyzarray[ ++ii ].x() = 36.072; xyzarray[ ii ].y() = 2.717; xyzarray[ ii ].z() = 5.143;
		xyzarray[ ++ii ].x() = 35.925; xyzarray[ ii ].y() = 1.365; xyzarray[ ii ].z() = 4.617;
		xyzarray[ ++ii ].x() = 34.641; xyzarray[ ii ].y() = 1.198; xyzarray[ ii ].z() = 3.788;
		xyzarray[ ++ii ].x() = 34.657; xyzarray[ ii ].y() = 0.582; xyzarray[ ii ].z() = 2.691;
		xyzarray[ ++ii ].x() = 33.532; xyzarray[ ii ].y() = 1.721; xyzarray[ ii ].z() = 4.327;
		xyzarray[ ++ii ].x() = 32.238; xyzarray[ ii ].y() = 1.679; xyzarray[ ii ].z() = 3.648;
		xyzarray[ ++ii ].x() = 32.207; xyzarray[ ii ].y() = 2.482; xyzarray[ ii ].z() = 2.355;
		xyzarray[ ++ii ].x() = 31.659; xyzarray[ ii ].y() = 2.003; xyzarray[ ii ].z() = 1.339;
		xyzarray[ ++ii ].x() = 32.79; xyzarray[ ii ].y() = 3.684; xyzarray[ ii ].z() = 2.394;
		xyzarray[ ++ii ].x() = 32.781; xyzarray[ ii ].y() = 4.587; xyzarray[ ii ].z() = 1.247;
		xyzarray[ ++ii ].x() = 33.631; xyzarray[ ii ].y() = 4.007; xyzarray[ ii ].z() = 0.119;
		xyzarray[ ++ii ].x() = 33.276; xyzarray[ ii ].y() = 4.127; xyzarray[ ii ].z() = -1.06;
		xyzarray[ ++ii ].x() = 34.744; xyzarray[ ii ].y() = 3.359; xyzarray[ ii ].z() = 0.487;
		xyzarray[ ++ii ].x() = 35.528; xyzarray[ ii ].y() = 2.641; xyzarray[ ii ].z() = -0.503;
		xyzarray[ ++ii ].x() = 34.775; xyzarray[ ii ].y() = 1.477; xyzarray[ ii ].z() = -1.154;
		xyzarray[ ++ii ].x() = 34.762; xyzarray[ ii ].y() = 1.368; xyzarray[ ii ].z() = -2.384;
		xyzarray[ ++ii ].x() = 34.157; xyzarray[ ii ].y() = 0.617; xyzarray[ ii ].z() = -0.339;
		xyzarray[ ++ii ].x() = 33.395; xyzarray[ ii ].y() = -0.5; xyzarray[ ii ].z() = -0.869;
		xyzarray[ ++ii ].x() = 32.315; xyzarray[ ii ].y() = 0.014; xyzarray[ ii ].z() = -1.829;
		xyzarray[ ++ii ].x() = 32.101; xyzarray[ ii ].y() = -0.587; xyzarray[ ii ].z() = -2.914;
		xyzarray[ ++ii ].x() = 31.68; xyzarray[ ii ].y() = 1.127; xyzarray[ ii ].z() = -1.461;
		xyzarray[ ++ii ].x() = 30.572; xyzarray[ ii ].y() = 1.665; xyzarray[ ii ].z() = -2.285;
		xyzarray[ ++ii ].x() = 31.118; xyzarray[ ii ].y() = 2.266; xyzarray[ ii ].z() = -3.564;
		xyzarray[ ++ii ].x() = 30.494; xyzarray[ ii ].y() = 2.15; xyzarray[ ii ].z() = -4.621;
		xyzarray[ ++ii ].x() = 32.313; xyzarray[ ii ].y() = 2.852; xyzarray[ ii ].z() = -3.496;
		xyzarray[ ++ii ].x() = 32.965; xyzarray[ ii ].y() = 3.336; xyzarray[ ii ].z() = -4.719;
		xyzarray[ ++ii ].x() = 33.302; xyzarray[ ii ].y() = 2.146; xyzarray[ ii ].z() = -5.656;
		xyzarray[ ++ii ].x() = 33.153; xyzarray[ ii ].y() = 2.227; xyzarray[ ii ].z() = -6.893;
		xyzarray[ ++ii ].x() = 33.756; xyzarray[ ii ].y() = 1.035; xyzarray[ ii ].z() = -5.078;
		xyzarray[ ++ii ].x() = 33.989; xyzarray[ ii ].y() = -0.179; xyzarray[ ii ].z() = -5.866;
		xyzarray[ ++ii ].x() = 32.742; xyzarray[ ii ].y() = -0.709; xyzarray[ ii ].z() = -6.553;
		xyzarray[ ++ii ].x() = 32.811; xyzarray[ ii ].y() = -1.185; xyzarray[ ii ].z() = -7.71;
		xyzarray[ ++ii ].x() = 31.598; xyzarray[ ii ].y() = -0.614; xyzarray[ ii ].z() = -5.882;
		xyzarray[ ++ii ].x() = 30.335; xyzarray[ ii ].y() = -0.963; xyzarray[ ii ].z() = -6.545;
		xyzarray[ ++ii ].x() = 30.138; xyzarray[ ii ].y() = -0.109; xyzarray[ ii ].z() = -7.798;
		xyzarray[ ++ii ].x() = 29.796; xyzarray[ ii ].y() = -0.623; xyzarray[ ii ].z() = -8.884;
		xyzarray[ ++ii ].x() = 30.382; xyzarray[ ii ].y() = 1.19; xyzarray[ ii ].z() = -7.67;
		xyzarray[ ++ii ].x() = 30.304; xyzarray[ ii ].y() = 2.088; xyzarray[ ii ].z() = -8.828;
		xyzarray[ ++ii ].x() = 31.296; xyzarray[ ii ].y() = 1.748; xyzarray[ ii ].z() = -9.934;
		xyzarray[ ++ii ].x() = 30.95; xyzarray[ ii ].y() = 1.83; xyzarray[ ii ].z() = -11.133;
		xyzarray[ ++ii ].x() = 32.521; xyzarray[ ii ].y() = 1.368; xyzarray[ ii ].z() = -9.562;
		xyzarray[ ++ii ].x() = 33.52; xyzarray[ ii ].y() = 0.894; xyzarray[ ii ].z() = -10.538;
		xyzarray[ ++ii ].x() = 32.955; xyzarray[ ii ].y() = -0.336; xyzarray[ ii ].z() = -11.275;
		xyzarray[ ++ii ].x() = 33.116; xyzarray[ ii ].y() = -0.479; xyzarray[ ii ].z() = -12.481;
		xyzarray[ ++ii ].x() = 32.299; xyzarray[ ii ].y() = -1.228; xyzarray[ ii ].z() = -10.548;
		xyzarray[ ++ii ].x() = 31.702; xyzarray[ ii ].y() = -2.417; xyzarray[ ii ].z() = -11.196;
		xyzarray[ ++ii ].x() = 30.489; xyzarray[ ii ].y() = -2.092; xyzarray[ ii ].z() = -12.072;
		xyzarray[ ++ii ].x() = 30.35; xyzarray[ ii ].y() = -2.625; xyzarray[ ii ].z() = -13.182;
		return xyzarray;
	}

	utility::vector1< numeric::xyzVector< float > >
	rmsd_coords_helices_184_196() {
		utility::vector1< numeric::xyzVector< float > > xyzarray( 224 );
		core::Size ii = 0;
		xyzarray[ ++ii ].x() = 18.654; xyzarray[ ii ].y() = 4.321; xyzarray[ ii ].z() = 9.171;
		xyzarray[ ++ii ].x() = 19.157; xyzarray[ ii ].y() = 2.956; xyzarray[ ii ].z() = 9.07;
		xyzarray[ ++ii ].x() = 18.185; xyzarray[ ii ].y() = 2.097; xyzarray[ ii ].z() = 8.257;
		xyzarray[ ++ii ].x() = 18.59; xyzarray[ ii ].y() = 1.291; xyzarray[ ii ].z() = 7.399;
		xyzarray[ ++ii ].x() = 16.894; xyzarray[ ii ].y() = 2.26; xyzarray[ ii ].z() = 8.537;
		xyzarray[ ++ii ].x() = 15.908; xyzarray[ ii ].y() = 1.449; xyzarray[ ii ].z() = 7.858;
		xyzarray[ ++ii ].x() = 15.918; xyzarray[ ii ].y() = 1.707; xyzarray[ ii ].z() = 6.356;
		xyzarray[ ++ii ].x() = 15.878; xyzarray[ ii ].y() = 0.758; xyzarray[ ii ].z() = 5.554;
		xyzarray[ ++ii ].x() = 16.01; xyzarray[ ii ].y() = 2.984; xyzarray[ ii ].z() = 5.977;
		xyzarray[ ++ii ].x() = 16.098; xyzarray[ ii ].y() = 3.349; xyzarray[ ii ].z() = 4.559;
		xyzarray[ ++ii ].x() = 17.351; xyzarray[ ii ].y() = 2.776; xyzarray[ ii ].z() = 3.881;
		xyzarray[ ++ii ].x() = 17.292; xyzarray[ ii ].y() = 2.279; xyzarray[ ii ].z() = 2.745;
		xyzarray[ ++ii ].x() = 18.485; xyzarray[ ii ].y() = 2.873; xyzarray[ ii ].z() = 4.568;
		xyzarray[ ++ii ].x() = 19.766; xyzarray[ ii ].y() = 2.393; xyzarray[ ii ].z() = 4.02;
		xyzarray[ ++ii ].x() = 19.684; xyzarray[ ii ].y() = 0.873; xyzarray[ ii ].z() = 3.871;
		xyzarray[ ++ii ].x() = 20.057; xyzarray[ ii ].y() = 0.328; xyzarray[ ii ].z() = 2.829;
		xyzarray[ ++ii ].x() = 19.121; xyzarray[ ii ].y() = 0.207; xyzarray[ ii ].z() = 4.879;
		xyzarray[ ++ii ].x() = 18.844; xyzarray[ ii ].y() = -1.248; xyzarray[ ii ].z() = 4.795;
		xyzarray[ ++ii ].x() = 18.026; xyzarray[ ii ].y() = -1.658; xyzarray[ ii ].z() = 3.583;
		xyzarray[ ++ii ].x() = 18.348; xyzarray[ ii ].y() = -2.64; xyzarray[ ii ].z() = 2.892;
		xyzarray[ ++ii ].x() = 16.962; xyzarray[ ii ].y() = -0.908; xyzarray[ ii ].z() = 3.33;
		xyzarray[ ++ii ].x() = 16.097; xyzarray[ ii ].y() = -1.125; xyzarray[ ii ].z() = 2.182;
		xyzarray[ ++ii ].x() = 16.869; xyzarray[ ii ].y() = -0.97; xyzarray[ ii ].z() = 0.85;
		xyzarray[ ++ii ].x() = 16.738; xyzarray[ ii ].y() = -1.809; xyzarray[ ii ].z() = -0.05;
		xyzarray[ ++ii ].x() = 17.665; xyzarray[ ii ].y() = 0.096; xyzarray[ ii ].z() = 0.736;
		xyzarray[ ++ii ].x() = 18.475; xyzarray[ ii ].y() = 0.343; xyzarray[ ii ].z() = -0.472;
		xyzarray[ ++ii ].x() = 19.518; xyzarray[ ii ].y() = -0.738; xyzarray[ ii ].z() = -0.716;
		xyzarray[ ++ii ].x() = 19.664; xyzarray[ ii ].y() = -1.207; xyzarray[ ii ].z() = -1.868;
		xyzarray[ ++ii ].x() = 20.228; xyzarray[ ii ].y() = -1.126; xyzarray[ ii ].z() = 0.35;
		xyzarray[ ++ii ].x() = 21.204; xyzarray[ ii ].y() = -2.22; xyzarray[ ii ].z() = 0.257;
		xyzarray[ ++ii ].x() = 20.536; xyzarray[ ii ].y() = -3.519; xyzarray[ ii ].z() = -0.16;
		xyzarray[ ++ii ].x() = 21.082; xyzarray[ ii ].y() = -4.24; xyzarray[ ii ].z() = -0.991;
		xyzarray[ ++ii ].x() = 19.359; xyzarray[ ii ].y() = -3.814; xyzarray[ ii ].z() = 0.397;
		xyzarray[ ++ii ].x() = 18.66; xyzarray[ ii ].y() = -5.056; xyzarray[ ii ].z() = 0.06;
		xyzarray[ ++ii ].x() = 18.223; xyzarray[ ii ].y() = -5.058; xyzarray[ ii ].z() = -1.4;
		xyzarray[ ++ii ].x() = 18.292; xyzarray[ ii ].y() = -6.113; xyzarray[ ii ].z() = -2.078;
		xyzarray[ ++ii ].x() = 17.735; xyzarray[ ii ].y() = -3.904; xyzarray[ ii ].z() = -1.867;
		xyzarray[ ++ii ].x() = 17.333; xyzarray[ ii ].y() = -3.773; xyzarray[ ii ].z() = -3.262;
		xyzarray[ ++ii ].x() = 18.524; xyzarray[ ii ].y() = -3.977; xyzarray[ ii ].z() = -4.22;
		xyzarray[ ++ii ].x() = 18.407; xyzarray[ ii ].y() = -4.699; xyzarray[ ii ].z() = -5.238;
		xyzarray[ ++ii ].x() = 19.653; xyzarray[ ii ].y() = -3.335; xyzarray[ ii ].z() = -3.909;
		xyzarray[ ++ii ].x() = 20.868; xyzarray[ ii ].y() = -3.529; xyzarray[ ii ].z() = -4.709;
		xyzarray[ ++ii ].x() = 21.208; xyzarray[ ii ].y() = -5.012; xyzarray[ ii ].z() = -4.713;
		xyzarray[ ++ii ].x() = 21.425; xyzarray[ ii ].y() = -5.616; xyzarray[ ii ].z() = -5.769;
		xyzarray[ ++ii ].x() = 21.216; xyzarray[ ii ].y() = -5.61; xyzarray[ ii ].z() = -3.534;
		xyzarray[ ++ii ].x() = 21.596; xyzarray[ ii ].y() = -7.016; xyzarray[ ii ].z() = -3.404;
		xyzarray[ ++ii ].x() = 20.683; xyzarray[ ii ].y() = -7.93; xyzarray[ ii ].z() = -4.205;
		xyzarray[ ++ii ].x() = 21.157; xyzarray[ ii ].y() = -8.813; xyzarray[ ii ].z() = -4.937;
		xyzarray[ ++ii ].x() = 19.369; xyzarray[ ii ].y() = -7.728; xyzarray[ ii ].z() = -4.096;
		xyzarray[ ++ii ].x() = 18.435; xyzarray[ ii ].y() = -8.563; xyzarray[ ii ].z() = -4.844;
		xyzarray[ ++ii ].x() = 18.706; xyzarray[ ii ].y() = -8.513; xyzarray[ ii ].z() = -6.37;
		xyzarray[ ++ii ].x() = 18.553; xyzarray[ ii ].y() = -9.516; xyzarray[ ii ].z() = -7.066;
		xyzarray[ ++ii ].x() = 19.134; xyzarray[ ii ].y() = -7.353; xyzarray[ ii ].z() = -6.871;
		xyzarray[ ++ii ].x() = 19.336; xyzarray[ ii ].y() = -7.15; xyzarray[ ii ].z() = -8.301;
		xyzarray[ ++ii ].x() = 20.728; xyzarray[ ii ].y() = -7.594; xyzarray[ ii ].z() = -8.797;
		xyzarray[ ++ii ].x() = 20.957; xyzarray[ ii ].y() = -7.666; xyzarray[ ii ].z() = -10.005;
		xyzarray[ ++ii ].x() = 28.028; xyzarray[ ii ].y() = -7.652; xyzarray[ ii ].z() = -5.186;
		xyzarray[ ++ii ].x() = 27.322; xyzarray[ ii ].y() = -8.219; xyzarray[ ii ].z() = -4.023;
		xyzarray[ ++ii ].x() = 28.205; xyzarray[ ii ].y() = -8.388; xyzarray[ ii ].z() = -2.782;
		xyzarray[ ++ii ].x() = 27.74; xyzarray[ ii ].y() = -8.178; xyzarray[ ii ].z() = -1.639;
		xyzarray[ ++ii ].x() = 29.479; xyzarray[ ii ].y() = -8.728; xyzarray[ ii ].z() = -2.988;
		xyzarray[ ++ii ].x() = 30.383; xyzarray[ ii ].y() = -8.922; xyzarray[ ii ].z() = -1.848;
		xyzarray[ ++ii ].x() = 30.641; xyzarray[ ii ].y() = -7.581; xyzarray[ ii ].z() = -1.108;
		xyzarray[ ++ii ].x() = 30.751; xyzarray[ ii ].y() = -7.55; xyzarray[ ii ].z() = 0.13;
		xyzarray[ ++ii ].x() = 30.721; xyzarray[ ii ].y() = -6.499; xyzarray[ ii ].z() = -1.866;
		xyzarray[ ++ii ].x() = 30.879; xyzarray[ ii ].y() = -5.153; xyzarray[ ii ].z() = -1.288;
		xyzarray[ ++ii ].x() = 29.636; xyzarray[ ii ].y() = -4.705; xyzarray[ ii ].z() = -0.541;
		xyzarray[ ++ii ].x() = 29.71; xyzarray[ ii ].y() = -4.153; xyzarray[ ii ].z() = 0.571;
		xyzarray[ ++ii ].x() = 28.488; xyzarray[ ii ].y() = -4.946; xyzarray[ ii ].z() = -1.138;
		xyzarray[ ++ii ].x() = 27.216; xyzarray[ ii ].y() = -4.597; xyzarray[ ii ].z() = -0.494;
		xyzarray[ ++ii ].x() = 27.001; xyzarray[ ii ].y() = -5.415; xyzarray[ ii ].z() = 0.811;
		xyzarray[ ++ii ].x() = 26.562; xyzarray[ ii ].y() = -4.871; xyzarray[ ii ].z() = 1.827;
		xyzarray[ ++ii ].x() = 27.335; xyzarray[ ii ].y() = -6.699; xyzarray[ ii ].z() = 0.782;
		xyzarray[ ++ii ].x() = 27.217; xyzarray[ ii ].y() = -7.556; xyzarray[ ii ].z() = 1.98;
		xyzarray[ ++ii ].x() = 28.107; xyzarray[ ii ].y() = -7.063; xyzarray[ ii ].z() = 3.109;
		xyzarray[ ++ii ].x() = 27.699; xyzarray[ ii ].y() = -7.09; xyzarray[ ii ].z() = 4.28;
		xyzarray[ ++ii ].x() = 29.315; xyzarray[ ii ].y() = -6.618; xyzarray[ ii ].z() = 2.76;
		xyzarray[ ++ii ].x() = 30.234; xyzarray[ ii ].y() = -6.013; xyzarray[ ii ].z() = 3.721;
		xyzarray[ ++ii ].x() = 29.65; xyzarray[ ii ].y() = -4.777; xyzarray[ ii ].z() = 4.365;
		xyzarray[ ++ii ].x() = 29.783; xyzarray[ ii ].y() = -4.609; xyzarray[ ii ].z() = 5.591;
		xyzarray[ ++ii ].x() = 29.014; xyzarray[ ii ].y() = -3.912; xyzarray[ ii ].z() = 3.559;
		xyzarray[ ++ii ].x() = 28.433; xyzarray[ ii ].y() = -2.686; xyzarray[ ii ].z() = 4.098;
		xyzarray[ ++ii ].x() = 27.231; xyzarray[ ii ].y() = -3.045; xyzarray[ ii ].z() = 5.008;
		xyzarray[ ++ii ].x() = 27.09; xyzarray[ ii ].y() = -2.522; xyzarray[ ii ].z() = 6.144;
		xyzarray[ ++ii ].x() = 26.386; xyzarray[ ii ].y() = -3.967; xyzarray[ ii ].z() = 4.537;
		xyzarray[ ++ii ].x() = 25.226; xyzarray[ ii ].y() = -4.368; xyzarray[ ii ].z() = 5.341;
		xyzarray[ ++ii ].x() = 25.651; xyzarray[ ii ].y() = -4.955; xyzarray[ ii ].z() = 6.698;
		xyzarray[ ++ii ].x() = 25.126; xyzarray[ ii ].y() = -4.599; xyzarray[ ii ].z() = 7.756;
		xyzarray[ ++ii ].x() = 26.615; xyzarray[ ii ].y() = -5.86; xyzarray[ ii ].z() = 6.663;
		xyzarray[ ++ii ].x() = 27.081; xyzarray[ ii ].y() = -6.515; xyzarray[ ii ].z() = 7.876;
		xyzarray[ ++ii ].x() = 27.695; xyzarray[ ii ].y() = -5.504; xyzarray[ ii ].z() = 8.82;
		xyzarray[ ++ii ].x() = 27.544; xyzarray[ ii ].y() = -5.615; xyzarray[ ii ].z() = 10.047;
		xyzarray[ ++ii ].x() = 28.411; xyzarray[ ii ].y() = -4.533; xyzarray[ ii ].z() = 8.264;
		xyzarray[ ++ii ].x() = 29.057; xyzarray[ ii ].y() = -3.511; xyzarray[ ii ].z() = 9.11;
		xyzarray[ ++ii ].x() = 28.041; xyzarray[ ii ].y() = -2.668; xyzarray[ ii ].z() = 9.858;
		xyzarray[ ++ii ].x() = 28.115; xyzarray[ ii ].y() = -2.494; xyzarray[ ii ].z() = 11.105;
		xyzarray[ ++ii ].x() = 27.075; xyzarray[ ii ].y() = -2.147; xyzarray[ ii ].z() = 9.105;
		xyzarray[ ++ii ].x() = 26.049; xyzarray[ ii ].y() = -1.302; xyzarray[ ii ].z() = 9.686;
		xyzarray[ ++ii ].x() = 25.188; xyzarray[ ii ].y() = -2.1; xyzarray[ ii ].z() = 10.678;
		xyzarray[ ++ii ].x() = 24.68; xyzarray[ ii ].y() = -1.559; xyzarray[ ii ].z() = 11.644;
		xyzarray[ ++ii ].x() = 24.997; xyzarray[ ii ].y() = -3.379; xyzarray[ ii ].z() = 10.413;
		xyzarray[ ++ii ].x() = 24.176; xyzarray[ ii ].y() = -4.222; xyzarray[ ii ].z() = 11.273;
		xyzarray[ ++ii ].x() = 24.804; xyzarray[ ii ].y() = -4.288; xyzarray[ ii ].z() = 12.666;
		xyzarray[ ++ii ].x() = 24.111; xyzarray[ ii ].y() = -4.251; xyzarray[ ii ].z() = 13.682;
		xyzarray[ ++ii ].x() = 26.124; xyzarray[ ii ].y() = -4.376; xyzarray[ ii ].z() = 12.722;
		xyzarray[ ++ii ].x() = 26.851; xyzarray[ ii ].y() = -4.358; xyzarray[ ii ].z() = 13.994;
		xyzarray[ ++ii ].x() = 26.77; xyzarray[ ii ].y() = -2.98; xyzarray[ ii ].z() = 14.662;
		xyzarray[ ++ii ].x() = 26.566; xyzarray[ ii ].y() = -2.885; xyzarray[ ii ].z() = 15.89;
		xyzarray[ ++ii ].x() = 26.95; xyzarray[ ii ].y() = -1.916; xyzarray[ ii ].z() = 13.88;
		xyzarray[ ++ii ].x() = 26.964; xyzarray[ ii ].y() = -0.568; xyzarray[ ii ].z() = 14.482;
		xyzarray[ ++ii ].x() = 25.616; xyzarray[ ii ].y() = -0.173; xyzarray[ ii ].z() = 15.086;
		xyzarray[ ++ii ].x() = 25.556; xyzarray[ ii ].y() = 0.366; xyzarray[ ii ].z() = 16.204;
		xyzarray[ ++ii ].x() = 18.654; xyzarray[ ii ].y() = 4.321; xyzarray[ ii ].z() = 9.171;
		xyzarray[ ++ii ].x() = 19.157; xyzarray[ ii ].y() = 2.956; xyzarray[ ii ].z() = 9.07;
		xyzarray[ ++ii ].x() = 18.185; xyzarray[ ii ].y() = 2.097; xyzarray[ ii ].z() = 8.257;
		xyzarray[ ++ii ].x() = 18.59; xyzarray[ ii ].y() = 1.291; xyzarray[ ii ].z() = 7.399;
		xyzarray[ ++ii ].x() = 16.894; xyzarray[ ii ].y() = 2.26; xyzarray[ ii ].z() = 8.537;
		xyzarray[ ++ii ].x() = 15.908; xyzarray[ ii ].y() = 1.449; xyzarray[ ii ].z() = 7.858;
		xyzarray[ ++ii ].x() = 15.918; xyzarray[ ii ].y() = 1.707; xyzarray[ ii ].z() = 6.356;
		xyzarray[ ++ii ].x() = 15.878; xyzarray[ ii ].y() = 0.758; xyzarray[ ii ].z() = 5.554;
		xyzarray[ ++ii ].x() = 16.01; xyzarray[ ii ].y() = 2.984; xyzarray[ ii ].z() = 5.977;
		xyzarray[ ++ii ].x() = 16.098; xyzarray[ ii ].y() = 3.349; xyzarray[ ii ].z() = 4.559;
		xyzarray[ ++ii ].x() = 17.351; xyzarray[ ii ].y() = 2.776; xyzarray[ ii ].z() = 3.881;
		xyzarray[ ++ii ].x() = 17.292; xyzarray[ ii ].y() = 2.279; xyzarray[ ii ].z() = 2.745;
		xyzarray[ ++ii ].x() = 18.485; xyzarray[ ii ].y() = 2.873; xyzarray[ ii ].z() = 4.568;
		xyzarray[ ++ii ].x() = 19.766; xyzarray[ ii ].y() = 2.393; xyzarray[ ii ].z() = 4.02;
		xyzarray[ ++ii ].x() = 19.684; xyzarray[ ii ].y() = 0.873; xyzarray[ ii ].z() = 3.871;
		xyzarray[ ++ii ].x() = 20.057; xyzarray[ ii ].y() = 0.328; xyzarray[ ii ].z() = 2.829;
		xyzarray[ ++ii ].x() = 19.121; xyzarray[ ii ].y() = 0.207; xyzarray[ ii ].z() = 4.879;
		xyzarray[ ++ii ].x() = 18.844; xyzarray[ ii ].y() = -1.248; xyzarray[ ii ].z() = 4.795;
		xyzarray[ ++ii ].x() = 18.026; xyzarray[ ii ].y() = -1.658; xyzarray[ ii ].z() = 3.583;
		xyzarray[ ++ii ].x() = 18.348; xyzarray[ ii ].y() = -2.64; xyzarray[ ii ].z() = 2.892;
		xyzarray[ ++ii ].x() = 16.962; xyzarray[ ii ].y() = -0.908; xyzarray[ ii ].z() = 3.33;
		xyzarray[ ++ii ].x() = 16.097; xyzarray[ ii ].y() = -1.125; xyzarray[ ii ].z() = 2.182;
		xyzarray[ ++ii ].x() = 16.869; xyzarray[ ii ].y() = -0.97; xyzarray[ ii ].z() = 0.85;
		xyzarray[ ++ii ].x() = 16.738; xyzarray[ ii ].y() = -1.809; xyzarray[ ii ].z() = -0.05;
		xyzarray[ ++ii ].x() = 17.665; xyzarray[ ii ].y() = 0.096; xyzarray[ ii ].z() = 0.736;
		xyzarray[ ++ii ].x() = 18.475; xyzarray[ ii ].y() = 0.343; xyzarray[ ii ].z() = -0.472;
		xyzarray[ ++ii ].x() = 19.518; xyzarray[ ii ].y() = -0.738; xyzarray[ ii ].z() = -0.716;
		xyzarray[ ++ii ].x() = 19.664; xyzarray[ ii ].y() = -1.207; xyzarray[ ii ].z() = -1.868;
		xyzarray[ ++ii ].x() = 20.228; xyzarray[ ii ].y() = -1.126; xyzarray[ ii ].z() = 0.35;
		xyzarray[ ++ii ].x() = 21.204; xyzarray[ ii ].y() = -2.22; xyzarray[ ii ].z() = 0.257;
		xyzarray[ ++ii ].x() = 20.536; xyzarray[ ii ].y() = -3.519; xyzarray[ ii ].z() = -0.16;
		xyzarray[ ++ii ].x() = 21.082; xyzarray[ ii ].y() = -4.24; xyzarray[ ii ].z() = -0.991;
		xyzarray[ ++ii ].x() = 19.359; xyzarray[ ii ].y() = -3.814; xyzarray[ ii ].z() = 0.397;
		xyzarray[ ++ii ].x() = 18.66; xyzarray[ ii ].y() = -5.056; xyzarray[ ii ].z() = 0.06;
		xyzarray[ ++ii ].x() = 18.223; xyzarray[ ii ].y() = -5.058; xyzarray[ ii ].z() = -1.4;
		xyzarray[ ++ii ].x() = 18.292; xyzarray[ ii ].y() = -6.113; xyzarray[ ii ].z() = -2.078;
		xyzarray[ ++ii ].x() = 17.735; xyzarray[ ii ].y() = -3.904; xyzarray[ ii ].z() = -1.867;
		xyzarray[ ++ii ].x() = 17.333; xyzarray[ ii ].y() = -3.773; xyzarray[ ii ].z() = -3.262;
		xyzarray[ ++ii ].x() = 18.524; xyzarray[ ii ].y() = -3.977; xyzarray[ ii ].z() = -4.22;
		xyzarray[ ++ii ].x() = 18.407; xyzarray[ ii ].y() = -4.699; xyzarray[ ii ].z() = -5.238;
		xyzarray[ ++ii ].x() = 19.653; xyzarray[ ii ].y() = -3.335; xyzarray[ ii ].z() = -3.909;
		xyzarray[ ++ii ].x() = 20.868; xyzarray[ ii ].y() = -3.529; xyzarray[ ii ].z() = -4.709;
		xyzarray[ ++ii ].x() = 21.208; xyzarray[ ii ].y() = -5.012; xyzarray[ ii ].z() = -4.713;
		xyzarray[ ++ii ].x() = 21.425; xyzarray[ ii ].y() = -5.616; xyzarray[ ii ].z() = -5.769;
		xyzarray[ ++ii ].x() = 21.216; xyzarray[ ii ].y() = -5.61; xyzarray[ ii ].z() = -3.534;
		xyzarray[ ++ii ].x() = 21.596; xyzarray[ ii ].y() = -7.016; xyzarray[ ii ].z() = -3.404;
		xyzarray[ ++ii ].x() = 20.683; xyzarray[ ii ].y() = -7.93; xyzarray[ ii ].z() = -4.205;
		xyzarray[ ++ii ].x() = 21.157; xyzarray[ ii ].y() = -8.813; xyzarray[ ii ].z() = -4.937;
		xyzarray[ ++ii ].x() = 19.369; xyzarray[ ii ].y() = -7.728; xyzarray[ ii ].z() = -4.096;
		xyzarray[ ++ii ].x() = 18.435; xyzarray[ ii ].y() = -8.563; xyzarray[ ii ].z() = -4.844;
		xyzarray[ ++ii ].x() = 18.706; xyzarray[ ii ].y() = -8.513; xyzarray[ ii ].z() = -6.37;
		xyzarray[ ++ii ].x() = 18.553; xyzarray[ ii ].y() = -9.516; xyzarray[ ii ].z() = -7.066;
		xyzarray[ ++ii ].x() = 19.134; xyzarray[ ii ].y() = -7.353; xyzarray[ ii ].z() = -6.871;
		xyzarray[ ++ii ].x() = 19.336; xyzarray[ ii ].y() = -7.15; xyzarray[ ii ].z() = -8.301;
		xyzarray[ ++ii ].x() = 20.728; xyzarray[ ii ].y() = -7.594; xyzarray[ ii ].z() = -8.797;
		xyzarray[ ++ii ].x() = 20.957; xyzarray[ ii ].y() = -7.666; xyzarray[ ii ].z() = -10.005;
		xyzarray[ ++ii ].x() = 28.028; xyzarray[ ii ].y() = -7.652; xyzarray[ ii ].z() = -5.186;
		xyzarray[ ++ii ].x() = 27.322; xyzarray[ ii ].y() = -8.219; xyzarray[ ii ].z() = -4.023;
		xyzarray[ ++ii ].x() = 28.205; xyzarray[ ii ].y() = -8.388; xyzarray[ ii ].z() = -2.782;
		xyzarray[ ++ii ].x() = 27.74; xyzarray[ ii ].y() = -8.178; xyzarray[ ii ].z() = -1.639;
		xyzarray[ ++ii ].x() = 29.479; xyzarray[ ii ].y() = -8.728; xyzarray[ ii ].z() = -2.988;
		xyzarray[ ++ii ].x() = 30.383; xyzarray[ ii ].y() = -8.922; xyzarray[ ii ].z() = -1.848;
		xyzarray[ ++ii ].x() = 30.641; xyzarray[ ii ].y() = -7.581; xyzarray[ ii ].z() = -1.108;
		xyzarray[ ++ii ].x() = 30.751; xyzarray[ ii ].y() = -7.55; xyzarray[ ii ].z() = 0.13;
		xyzarray[ ++ii ].x() = 30.721; xyzarray[ ii ].y() = -6.499; xyzarray[ ii ].z() = -1.866;
		xyzarray[ ++ii ].x() = 30.879; xyzarray[ ii ].y() = -5.153; xyzarray[ ii ].z() = -1.288;
		xyzarray[ ++ii ].x() = 29.636; xyzarray[ ii ].y() = -4.705; xyzarray[ ii ].z() = -0.541;
		xyzarray[ ++ii ].x() = 29.71; xyzarray[ ii ].y() = -4.153; xyzarray[ ii ].z() = 0.571;
		xyzarray[ ++ii ].x() = 28.488; xyzarray[ ii ].y() = -4.946; xyzarray[ ii ].z() = -1.138;
		xyzarray[ ++ii ].x() = 27.216; xyzarray[ ii ].y() = -4.597; xyzarray[ ii ].z() = -0.494;
		xyzarray[ ++ii ].x() = 27.001; xyzarray[ ii ].y() = -5.415; xyzarray[ ii ].z() = 0.811;
		xyzarray[ ++ii ].x() = 26.562; xyzarray[ ii ].y() = -4.871; xyzarray[ ii ].z() = 1.827;
		xyzarray[ ++ii ].x() = 27.335; xyzarray[ ii ].y() = -6.699; xyzarray[ ii ].z() = 0.782;
		xyzarray[ ++ii ].x() = 27.217; xyzarray[ ii ].y() = -7.556; xyzarray[ ii ].z() = 1.98;
		xyzarray[ ++ii ].x() = 28.107; xyzarray[ ii ].y() = -7.063; xyzarray[ ii ].z() = 3.109;
		xyzarray[ ++ii ].x() = 27.699; xyzarray[ ii ].y() = -7.09; xyzarray[ ii ].z() = 4.28;
		xyzarray[ ++ii ].x() = 29.315; xyzarray[ ii ].y() = -6.618; xyzarray[ ii ].z() = 2.76;
		xyzarray[ ++ii ].x() = 30.234; xyzarray[ ii ].y() = -6.013; xyzarray[ ii ].z() = 3.721;
		xyzarray[ ++ii ].x() = 29.65; xyzarray[ ii ].y() = -4.777; xyzarray[ ii ].z() = 4.365;
		xyzarray[ ++ii ].x() = 29.783; xyzarray[ ii ].y() = -4.609; xyzarray[ ii ].z() = 5.591;
		xyzarray[ ++ii ].x() = 29.014; xyzarray[ ii ].y() = -3.912; xyzarray[ ii ].z() = 3.559;
		xyzarray[ ++ii ].x() = 28.433; xyzarray[ ii ].y() = -2.686; xyzarray[ ii ].z() = 4.098;
		xyzarray[ ++ii ].x() = 27.231; xyzarray[ ii ].y() = -3.045; xyzarray[ ii ].z() = 5.008;
		xyzarray[ ++ii ].x() = 27.09; xyzarray[ ii ].y() = -2.522; xyzarray[ ii ].z() = 6.144;
		xyzarray[ ++ii ].x() = 26.386; xyzarray[ ii ].y() = -3.967; xyzarray[ ii ].z() = 4.537;
		xyzarray[ ++ii ].x() = 25.226; xyzarray[ ii ].y() = -4.368; xyzarray[ ii ].z() = 5.341;
		xyzarray[ ++ii ].x() = 25.651; xyzarray[ ii ].y() = -4.955; xyzarray[ ii ].z() = 6.698;
		xyzarray[ ++ii ].x() = 25.126; xyzarray[ ii ].y() = -4.599; xyzarray[ ii ].z() = 7.756;
		xyzarray[ ++ii ].x() = 26.615; xyzarray[ ii ].y() = -5.86; xyzarray[ ii ].z() = 6.663;
		xyzarray[ ++ii ].x() = 27.081; xyzarray[ ii ].y() = -6.515; xyzarray[ ii ].z() = 7.876;
		xyzarray[ ++ii ].x() = 27.695; xyzarray[ ii ].y() = -5.504; xyzarray[ ii ].z() = 8.82;
		xyzarray[ ++ii ].x() = 27.544; xyzarray[ ii ].y() = -5.615; xyzarray[ ii ].z() = 10.047;
		xyzarray[ ++ii ].x() = 28.411; xyzarray[ ii ].y() = -4.533; xyzarray[ ii ].z() = 8.264;
		xyzarray[ ++ii ].x() = 29.057; xyzarray[ ii ].y() = -3.511; xyzarray[ ii ].z() = 9.11;
		xyzarray[ ++ii ].x() = 28.041; xyzarray[ ii ].y() = -2.668; xyzarray[ ii ].z() = 9.858;
		xyzarray[ ++ii ].x() = 28.115; xyzarray[ ii ].y() = -2.494; xyzarray[ ii ].z() = 11.105;
		xyzarray[ ++ii ].x() = 27.075; xyzarray[ ii ].y() = -2.147; xyzarray[ ii ].z() = 9.105;
		xyzarray[ ++ii ].x() = 26.049; xyzarray[ ii ].y() = -1.302; xyzarray[ ii ].z() = 9.686;
		xyzarray[ ++ii ].x() = 25.188; xyzarray[ ii ].y() = -2.1; xyzarray[ ii ].z() = 10.678;
		xyzarray[ ++ii ].x() = 24.68; xyzarray[ ii ].y() = -1.559; xyzarray[ ii ].z() = 11.644;
		xyzarray[ ++ii ].x() = 24.997; xyzarray[ ii ].y() = -3.379; xyzarray[ ii ].z() = 10.413;
		xyzarray[ ++ii ].x() = 24.176; xyzarray[ ii ].y() = -4.222; xyzarray[ ii ].z() = 11.273;
		xyzarray[ ++ii ].x() = 24.804; xyzarray[ ii ].y() = -4.288; xyzarray[ ii ].z() = 12.666;
		xyzarray[ ++ii ].x() = 24.111; xyzarray[ ii ].y() = -4.251; xyzarray[ ii ].z() = 13.682;
		xyzarray[ ++ii ].x() = 26.124; xyzarray[ ii ].y() = -4.376; xyzarray[ ii ].z() = 12.722;
		xyzarray[ ++ii ].x() = 26.851; xyzarray[ ii ].y() = -4.358; xyzarray[ ii ].z() = 13.994;
		xyzarray[ ++ii ].x() = 26.77; xyzarray[ ii ].y() = -2.98; xyzarray[ ii ].z() = 14.662;
		xyzarray[ ++ii ].x() = 26.566; xyzarray[ ii ].y() = -2.885; xyzarray[ ii ].z() = 15.89;
		xyzarray[ ++ii ].x() = 26.95; xyzarray[ ii ].y() = -1.916; xyzarray[ ii ].z() = 13.88;
		xyzarray[ ++ii ].x() = 26.964; xyzarray[ ii ].y() = -0.568; xyzarray[ ii ].z() = 14.482;
		xyzarray[ ++ii ].x() = 25.616; xyzarray[ ii ].y() = -0.173; xyzarray[ ii ].z() = 15.086;
		xyzarray[ ++ii ].x() = 25.556; xyzarray[ ii ].y() = 0.366; xyzarray[ ii ].z() = 16.204;
		return xyzarray;
	}

	utility::vector1< numeric::xyzVector< float > >
	clash_coords_helices_184_196() {
		utility::vector1< numeric::xyzVector< float > > xyzarray( 224 );
		core::Size ii = 0;
		xyzarray[ ++ii ].x() = 19.371; xyzarray[ ii ].y() = 2.495; xyzarray[ ii ].z() = -12.083;
		xyzarray[ ++ii ].x() = 18.688; xyzarray[ ii ].y() = 2.32; xyzarray[ ii ].z() = -10.809;
		xyzarray[ ++ii ].x() = 19.671; xyzarray[ ii ].y() = 2.061; xyzarray[ ii ].z() = -9.666;
		xyzarray[ ++ii ].x() = 19.638; xyzarray[ ii ].y() = 2.736; xyzarray[ ii ].z() = -8.639;
		xyzarray[ ++ii ].x() = 20.565; xyzarray[ ii ].y() = 1.11; xyzarray[ ii ].z() = -9.867;
		xyzarray[ ++ii ].x() = 21.519; xyzarray[ ii ].y() = 0.753; xyzarray[ ii ].z() = -8.84;
		xyzarray[ ++ii ].x() = 22.493; xyzarray[ ii ].y() = 1.877; xyzarray[ ii ].z() = -8.551;
		xyzarray[ ++ii ].x() = 22.862; xyzarray[ ii ].y() = 2.08; xyzarray[ ii ].z() = -7.397;
		xyzarray[ ++ii ].x() = 22.886; xyzarray[ ii ].y() = 2.615; xyzarray[ ii ].z() = -9.588;
		xyzarray[ ++ii ].x() = 23.818; xyzarray[ ii ].y() = 3.736; xyzarray[ ii ].z() = -9.399;
		xyzarray[ ++ii ].x() = 23.163; xyzarray[ ii ].y() = 4.785; xyzarray[ ii ].z() = -8.479;
		xyzarray[ ++ii ].x() = 23.824; xyzarray[ ii ].y() = 5.345; xyzarray[ ii ].z() = -7.589;
		xyzarray[ ++ii ].x() = 21.871; xyzarray[ ii ].y() = 5.039; xyzarray[ ii ].z() = -8.685;
		xyzarray[ ++ii ].x() = 21.154; xyzarray[ ii ].y() = 5.96; xyzarray[ ii ].z() = -7.83;
		xyzarray[ ++ii ].x() = 20.987; xyzarray[ ii ].y() = 5.458; xyzarray[ ii ].z() = -6.391;
		xyzarray[ ++ii ].x() = 21.087; xyzarray[ ii ].y() = 6.243; xyzarray[ ii ].z() = -5.449;
		xyzarray[ ++ii ].x() = 20.723; xyzarray[ ii ].y() = 4.164; xyzarray[ ii ].z() = -6.223;
		xyzarray[ ++ii ].x() = 20.626; xyzarray[ ii ].y() = 3.607; xyzarray[ ii ].z() = -4.864;
		xyzarray[ ++ii ].x() = 21.973; xyzarray[ ii ].y() = 3.714; xyzarray[ ii ].z() = -4.154;
		xyzarray[ ++ii ].x() = 22.029; xyzarray[ ii ].y() = 4.074; xyzarray[ ii ].z() = -2.977;
		xyzarray[ ++ii ].x() = 23.049; xyzarray[ ii ].y() = 3.395; xyzarray[ ii ].z() = -4.867;
		xyzarray[ ++ii ].x() = 24.4; xyzarray[ ii ].y() = 3.591; xyzarray[ ii ].z() = -4.306;
		xyzarray[ ++ii ].x() = 24.644; xyzarray[ ii ].y() = 5.027; xyzarray[ ii ].z() = -3.883;
		xyzarray[ ++ii ].x() = 25.157; xyzarray[ ii ].y() = 5.281; xyzarray[ ii ].z() = -2.77;
		xyzarray[ ++ii ].x() = 24.251; xyzarray[ ii ].y() = 5.97; xyzarray[ ii ].z() = -4.739;
		xyzarray[ ++ii ].x() = 24.398; xyzarray[ ii ].y() = 7.391; xyzarray[ ii ].z() = -4.404;
		xyzarray[ ++ii ].x() = 23.609; xyzarray[ ii ].y() = 7.745; xyzarray[ ii ].z() = -3.122;
		xyzarray[ ++ii ].x() = 24.102; xyzarray[ ii ].y() = 8.485; xyzarray[ ii ].z() = -2.272;
		xyzarray[ ++ii ].x() = 22.395; xyzarray[ ii ].y() = 7.206; xyzarray[ ii ].z() = -2.998;
		xyzarray[ ++ii ].x() = 21.578; xyzarray[ ii ].y() = 7.423; xyzarray[ ii ].z() = -1.813;
		xyzarray[ ++ii ].x() = 22.308; xyzarray[ ii ].y() = 6.933; xyzarray[ ii ].z() = -0.563;
		xyzarray[ ++ii ].x() = 22.333; xyzarray[ ii ].y() = 7.635; xyzarray[ ii ].z() = 0.464;
		xyzarray[ ++ii ].x() = 22.901; xyzarray[ ii ].y() = 5.738; xyzarray[ ii ].z() = -0.653;
		xyzarray[ ++ii ].x() = 23.595; xyzarray[ ii ].y() = 5.163; xyzarray[ ii ].z() = 0.501;
		xyzarray[ ++ii ].x() = 24.824; xyzarray[ ii ].y() = 6; xyzarray[ ii ].z() = 0.875;
		xyzarray[ ++ii ].x() = 25.02; xyzarray[ ii ].y() = 6.353; xyzarray[ ii ].z() = 2.051;
		xyzarray[ ++ii ].x() = 25.639; xyzarray[ ii ].y() = 6.317; xyzarray[ ii ].z() = -0.127;
		xyzarray[ ++ii ].x() = 26.867; xyzarray[ ii ].y() = 7.07; xyzarray[ ii ].z() = 0.082;
		xyzarray[ ++ii ].x() = 26.595; xyzarray[ ii ].y() = 8.452; xyzarray[ ii ].z() = 0.682;
		xyzarray[ ++ii ].x() = 27.216; xyzarray[ ii ].y() = 8.83; xyzarray[ ii ].z() = 1.684;
		xyzarray[ ++ii ].x() = 25.656; xyzarray[ ii ].y() = 9.185; xyzarray[ ii ].z() = 0.079;
		xyzarray[ ++ii ].x() = 25.231; xyzarray[ ii ].y() = 10.521; xyzarray[ ii ].z() = 0.563;
		xyzarray[ ++ii ].x() = 24.746; xyzarray[ ii ].y() = 10.459; xyzarray[ ii ].z() = 2.037;
		xyzarray[ ++ii ].x() = 25.06; xyzarray[ ii ].y() = 11.345; xyzarray[ ii ].z() = 2.836;
		xyzarray[ ++ii ].x() = 24.003; xyzarray[ ii ].y() = 9.406; xyzarray[ ii ].z() = 2.392;
		xyzarray[ ++ii ].x() = 23.475; xyzarray[ ii ].y() = 9.235; xyzarray[ ii ].z() = 3.757;
		xyzarray[ ++ii ].x() = 24.584; xyzarray[ ii ].y() = 8.937; xyzarray[ ii ].z() = 4.776;
		xyzarray[ ++ii ].x() = 24.633; xyzarray[ ii ].y() = 9.529; xyzarray[ ii ].z() = 5.859;
		xyzarray[ ++ii ].x() = 25.513; xyzarray[ ii ].y() = 8.057; xyzarray[ ii ].z() = 4.415;
		xyzarray[ ++ii ].x() = 26.642; xyzarray[ ii ].y() = 7.807; xyzarray[ ii ].z() = 5.294;
		xyzarray[ ++ii ].x() = 27.463; xyzarray[ ii ].y() = 9.08; xyzarray[ ii ].z() = 5.629;
		xyzarray[ ++ii ].x() = 27.8; xyzarray[ ii ].y() = 9.327; xyzarray[ ii ].z() = 6.816;
		xyzarray[ ++ii ].x() = 27.764; xyzarray[ ii ].y() = 9.899; xyzarray[ ii ].z() = 4.627;
		xyzarray[ ++ii ].x() = 28.581; xyzarray[ ii ].y() = 11.104; xyzarray[ ii ].z() = 4.895;
		xyzarray[ ++ii ].x() = 27.827; xyzarray[ ii ].y() = 12.143; xyzarray[ ii ].z() = 5.711;
		xyzarray[ ++ii ].x() = 28.441; xyzarray[ ii ].y() = 12.951; xyzarray[ ii ].z() = 6.413;
		xyzarray[ ++ii ].x() = 34.211; xyzarray[ ii ].y() = 4.736; xyzarray[ ii ].z() = 8.755;
		xyzarray[ ++ii ].x() = 33.305; xyzarray[ ii ].y() = 4.996; xyzarray[ ii ].z() = 7.626;
		xyzarray[ ++ii ].x() = 34.037; xyzarray[ ii ].y() = 4.914; xyzarray[ ii ].z() = 6.279;
		xyzarray[ ++ii ].x() = 33.425; xyzarray[ ii ].y() = 4.569; xyzarray[ ii ].z() = 5.269;
		xyzarray[ ++ii ].x() = 35.346; xyzarray[ ii ].y() = 5.166; xyzarray[ ii ].z() = 6.26;
		xyzarray[ ++ii ].x() = 36.072; xyzarray[ ii ].y() = 5.148; xyzarray[ ii ].z() = 4.98;
		xyzarray[ ++ii ].x() = 36.024; xyzarray[ ii ].y() = 3.765; xyzarray[ ii ].z() = 4.323;
		xyzarray[ ++ii ].x() = 35.895; xyzarray[ ii ].y() = 3.652; xyzarray[ ii ].z() = 3.095;
		xyzarray[ ++ii ].x() = 36.072; xyzarray[ ii ].y() = 2.717; xyzarray[ ii ].z() = 5.143;
		xyzarray[ ++ii ].x() = 35.925; xyzarray[ ii ].y() = 1.365; xyzarray[ ii ].z() = 4.617;
		xyzarray[ ++ii ].x() = 34.641; xyzarray[ ii ].y() = 1.198; xyzarray[ ii ].z() = 3.788;
		xyzarray[ ++ii ].x() = 34.657; xyzarray[ ii ].y() = 0.582; xyzarray[ ii ].z() = 2.691;
		xyzarray[ ++ii ].x() = 33.532; xyzarray[ ii ].y() = 1.721; xyzarray[ ii ].z() = 4.327;
		xyzarray[ ++ii ].x() = 32.238; xyzarray[ ii ].y() = 1.679; xyzarray[ ii ].z() = 3.648;
		xyzarray[ ++ii ].x() = 32.207; xyzarray[ ii ].y() = 2.482; xyzarray[ ii ].z() = 2.355;
		xyzarray[ ++ii ].x() = 31.659; xyzarray[ ii ].y() = 2.003; xyzarray[ ii ].z() = 1.339;
		xyzarray[ ++ii ].x() = 32.79; xyzarray[ ii ].y() = 3.684; xyzarray[ ii ].z() = 2.394;
		xyzarray[ ++ii ].x() = 32.781; xyzarray[ ii ].y() = 4.587; xyzarray[ ii ].z() = 1.247;
		xyzarray[ ++ii ].x() = 33.631; xyzarray[ ii ].y() = 4.007; xyzarray[ ii ].z() = 0.119;
		xyzarray[ ++ii ].x() = 33.276; xyzarray[ ii ].y() = 4.127; xyzarray[ ii ].z() = -1.06;
		xyzarray[ ++ii ].x() = 34.744; xyzarray[ ii ].y() = 3.359; xyzarray[ ii ].z() = 0.487;
		xyzarray[ ++ii ].x() = 35.528; xyzarray[ ii ].y() = 2.641; xyzarray[ ii ].z() = -0.503;
		xyzarray[ ++ii ].x() = 34.775; xyzarray[ ii ].y() = 1.477; xyzarray[ ii ].z() = -1.154;
		xyzarray[ ++ii ].x() = 34.762; xyzarray[ ii ].y() = 1.368; xyzarray[ ii ].z() = -2.384;
		xyzarray[ ++ii ].x() = 34.157; xyzarray[ ii ].y() = 0.617; xyzarray[ ii ].z() = -0.339;
		xyzarray[ ++ii ].x() = 33.395; xyzarray[ ii ].y() = -0.5; xyzarray[ ii ].z() = -0.869;
		xyzarray[ ++ii ].x() = 32.315; xyzarray[ ii ].y() = 0.014; xyzarray[ ii ].z() = -1.829;
		xyzarray[ ++ii ].x() = 32.101; xyzarray[ ii ].y() = -0.587; xyzarray[ ii ].z() = -2.914;
		xyzarray[ ++ii ].x() = 31.68; xyzarray[ ii ].y() = 1.127; xyzarray[ ii ].z() = -1.461;
		xyzarray[ ++ii ].x() = 30.572; xyzarray[ ii ].y() = 1.665; xyzarray[ ii ].z() = -2.285;
		xyzarray[ ++ii ].x() = 31.118; xyzarray[ ii ].y() = 2.266; xyzarray[ ii ].z() = -3.564;
		xyzarray[ ++ii ].x() = 30.494; xyzarray[ ii ].y() = 2.15; xyzarray[ ii ].z() = -4.621;
		xyzarray[ ++ii ].x() = 32.313; xyzarray[ ii ].y() = 2.852; xyzarray[ ii ].z() = -3.496;
		xyzarray[ ++ii ].x() = 32.965; xyzarray[ ii ].y() = 3.336; xyzarray[ ii ].z() = -4.719;
		xyzarray[ ++ii ].x() = 33.302; xyzarray[ ii ].y() = 2.146; xyzarray[ ii ].z() = -5.656;
		xyzarray[ ++ii ].x() = 33.153; xyzarray[ ii ].y() = 2.227; xyzarray[ ii ].z() = -6.893;
		xyzarray[ ++ii ].x() = 33.756; xyzarray[ ii ].y() = 1.035; xyzarray[ ii ].z() = -5.078;
		xyzarray[ ++ii ].x() = 33.989; xyzarray[ ii ].y() = -0.179; xyzarray[ ii ].z() = -5.866;
		xyzarray[ ++ii ].x() = 32.742; xyzarray[ ii ].y() = -0.709; xyzarray[ ii ].z() = -6.553;
		xyzarray[ ++ii ].x() = 32.811; xyzarray[ ii ].y() = -1.185; xyzarray[ ii ].z() = -7.71;
		xyzarray[ ++ii ].x() = 31.598; xyzarray[ ii ].y() = -0.614; xyzarray[ ii ].z() = -5.882;
		xyzarray[ ++ii ].x() = 30.335; xyzarray[ ii ].y() = -0.963; xyzarray[ ii ].z() = -6.545;
		xyzarray[ ++ii ].x() = 30.138; xyzarray[ ii ].y() = -0.109; xyzarray[ ii ].z() = -7.798;
		xyzarray[ ++ii ].x() = 29.796; xyzarray[ ii ].y() = -0.623; xyzarray[ ii ].z() = -8.884;
		xyzarray[ ++ii ].x() = 30.382; xyzarray[ ii ].y() = 1.19; xyzarray[ ii ].z() = -7.67;
		xyzarray[ ++ii ].x() = 30.304; xyzarray[ ii ].y() = 2.088; xyzarray[ ii ].z() = -8.828;
		xyzarray[ ++ii ].x() = 31.296; xyzarray[ ii ].y() = 1.748; xyzarray[ ii ].z() = -9.934;
		xyzarray[ ++ii ].x() = 30.95; xyzarray[ ii ].y() = 1.83; xyzarray[ ii ].z() = -11.133;
		xyzarray[ ++ii ].x() = 32.521; xyzarray[ ii ].y() = 1.368; xyzarray[ ii ].z() = -9.562;
		xyzarray[ ++ii ].x() = 33.52; xyzarray[ ii ].y() = 0.894; xyzarray[ ii ].z() = -10.538;
		xyzarray[ ++ii ].x() = 32.955; xyzarray[ ii ].y() = -0.336; xyzarray[ ii ].z() = -11.275;
		xyzarray[ ++ii ].x() = 33.116; xyzarray[ ii ].y() = -0.479; xyzarray[ ii ].z() = -12.481;
		xyzarray[ ++ii ].x() = 32.299; xyzarray[ ii ].y() = -1.228; xyzarray[ ii ].z() = -10.548;
		xyzarray[ ++ii ].x() = 31.702; xyzarray[ ii ].y() = -2.417; xyzarray[ ii ].z() = -11.196;
		xyzarray[ ++ii ].x() = 30.489; xyzarray[ ii ].y() = -2.092; xyzarray[ ii ].z() = -12.072;
		xyzarray[ ++ii ].x() = 30.35; xyzarray[ ii ].y() = -2.625; xyzarray[ ii ].z() = -13.182;
		xyzarray[ ++ii ].x() = 19.371; xyzarray[ ii ].y() = 2.495; xyzarray[ ii ].z() = -12.083;
		xyzarray[ ++ii ].x() = 18.688; xyzarray[ ii ].y() = 2.32; xyzarray[ ii ].z() = -10.809;
		xyzarray[ ++ii ].x() = 19.671; xyzarray[ ii ].y() = 2.061; xyzarray[ ii ].z() = -9.666;
		xyzarray[ ++ii ].x() = 19.638; xyzarray[ ii ].y() = 2.736; xyzarray[ ii ].z() = -8.639;
		xyzarray[ ++ii ].x() = 20.565; xyzarray[ ii ].y() = 1.11; xyzarray[ ii ].z() = -9.867;
		xyzarray[ ++ii ].x() = 21.519; xyzarray[ ii ].y() = 0.753; xyzarray[ ii ].z() = -8.84;
		xyzarray[ ++ii ].x() = 22.493; xyzarray[ ii ].y() = 1.877; xyzarray[ ii ].z() = -8.551;
		xyzarray[ ++ii ].x() = 22.862; xyzarray[ ii ].y() = 2.08; xyzarray[ ii ].z() = -7.397;
		xyzarray[ ++ii ].x() = 22.886; xyzarray[ ii ].y() = 2.615; xyzarray[ ii ].z() = -9.588;
		xyzarray[ ++ii ].x() = 23.818; xyzarray[ ii ].y() = 3.736; xyzarray[ ii ].z() = -9.399;
		xyzarray[ ++ii ].x() = 23.163; xyzarray[ ii ].y() = 4.785; xyzarray[ ii ].z() = -8.479;
		xyzarray[ ++ii ].x() = 23.824; xyzarray[ ii ].y() = 5.345; xyzarray[ ii ].z() = -7.589;
		xyzarray[ ++ii ].x() = 21.871; xyzarray[ ii ].y() = 5.039; xyzarray[ ii ].z() = -8.685;
		xyzarray[ ++ii ].x() = 21.154; xyzarray[ ii ].y() = 5.96; xyzarray[ ii ].z() = -7.83;
		xyzarray[ ++ii ].x() = 20.987; xyzarray[ ii ].y() = 5.458; xyzarray[ ii ].z() = -6.391;
		xyzarray[ ++ii ].x() = 21.087; xyzarray[ ii ].y() = 6.243; xyzarray[ ii ].z() = -5.449;
		xyzarray[ ++ii ].x() = 20.723; xyzarray[ ii ].y() = 4.164; xyzarray[ ii ].z() = -6.223;
		xyzarray[ ++ii ].x() = 20.626; xyzarray[ ii ].y() = 3.607; xyzarray[ ii ].z() = -4.864;
		xyzarray[ ++ii ].x() = 21.973; xyzarray[ ii ].y() = 3.714; xyzarray[ ii ].z() = -4.154;
		xyzarray[ ++ii ].x() = 22.029; xyzarray[ ii ].y() = 4.074; xyzarray[ ii ].z() = -2.977;
		xyzarray[ ++ii ].x() = 23.049; xyzarray[ ii ].y() = 3.395; xyzarray[ ii ].z() = -4.867;
		xyzarray[ ++ii ].x() = 24.4; xyzarray[ ii ].y() = 3.591; xyzarray[ ii ].z() = -4.306;
		xyzarray[ ++ii ].x() = 24.644; xyzarray[ ii ].y() = 5.027; xyzarray[ ii ].z() = -3.883;
		xyzarray[ ++ii ].x() = 25.157; xyzarray[ ii ].y() = 5.281; xyzarray[ ii ].z() = -2.77;
		xyzarray[ ++ii ].x() = 24.251; xyzarray[ ii ].y() = 5.97; xyzarray[ ii ].z() = -4.739;
		xyzarray[ ++ii ].x() = 24.398; xyzarray[ ii ].y() = 7.391; xyzarray[ ii ].z() = -4.404;
		xyzarray[ ++ii ].x() = 23.609; xyzarray[ ii ].y() = 7.745; xyzarray[ ii ].z() = -3.122;
		xyzarray[ ++ii ].x() = 24.102; xyzarray[ ii ].y() = 8.485; xyzarray[ ii ].z() = -2.272;
		xyzarray[ ++ii ].x() = 22.395; xyzarray[ ii ].y() = 7.206; xyzarray[ ii ].z() = -2.998;
		xyzarray[ ++ii ].x() = 21.578; xyzarray[ ii ].y() = 7.423; xyzarray[ ii ].z() = -1.813;
		xyzarray[ ++ii ].x() = 22.308; xyzarray[ ii ].y() = 6.933; xyzarray[ ii ].z() = -0.563;
		xyzarray[ ++ii ].x() = 22.333; xyzarray[ ii ].y() = 7.635; xyzarray[ ii ].z() = 0.464;
		xyzarray[ ++ii ].x() = 22.901; xyzarray[ ii ].y() = 5.738; xyzarray[ ii ].z() = -0.653;
		xyzarray[ ++ii ].x() = 23.595; xyzarray[ ii ].y() = 5.163; xyzarray[ ii ].z() = 0.501;
		xyzarray[ ++ii ].x() = 24.824; xyzarray[ ii ].y() = 6; xyzarray[ ii ].z() = 0.875;
		xyzarray[ ++ii ].x() = 25.02; xyzarray[ ii ].y() = 6.353; xyzarray[ ii ].z() = 2.051;
		xyzarray[ ++ii ].x() = 25.639; xyzarray[ ii ].y() = 6.317; xyzarray[ ii ].z() = -0.127;
		xyzarray[ ++ii ].x() = 26.867; xyzarray[ ii ].y() = 7.07; xyzarray[ ii ].z() = 0.082;
		xyzarray[ ++ii ].x() = 26.595; xyzarray[ ii ].y() = 8.452; xyzarray[ ii ].z() = 0.682;
		xyzarray[ ++ii ].x() = 27.216; xyzarray[ ii ].y() = 8.83; xyzarray[ ii ].z() = 1.684;
		xyzarray[ ++ii ].x() = 25.656; xyzarray[ ii ].y() = 9.185; xyzarray[ ii ].z() = 0.079;
		xyzarray[ ++ii ].x() = 25.231; xyzarray[ ii ].y() = 10.521; xyzarray[ ii ].z() = 0.563;
		xyzarray[ ++ii ].x() = 24.746; xyzarray[ ii ].y() = 10.459; xyzarray[ ii ].z() = 2.037;
		xyzarray[ ++ii ].x() = 25.06; xyzarray[ ii ].y() = 11.345; xyzarray[ ii ].z() = 2.836;
		xyzarray[ ++ii ].x() = 24.003; xyzarray[ ii ].y() = 9.406; xyzarray[ ii ].z() = 2.392;
		xyzarray[ ++ii ].x() = 23.475; xyzarray[ ii ].y() = 9.235; xyzarray[ ii ].z() = 3.757;
		xyzarray[ ++ii ].x() = 24.584; xyzarray[ ii ].y() = 8.937; xyzarray[ ii ].z() = 4.776;
		xyzarray[ ++ii ].x() = 24.633; xyzarray[ ii ].y() = 9.529; xyzarray[ ii ].z() = 5.859;
		xyzarray[ ++ii ].x() = 25.513; xyzarray[ ii ].y() = 8.057; xyzarray[ ii ].z() = 4.415;
		xyzarray[ ++ii ].x() = 26.642; xyzarray[ ii ].y() = 7.807; xyzarray[ ii ].z() = 5.294;
		xyzarray[ ++ii ].x() = 27.463; xyzarray[ ii ].y() = 9.08; xyzarray[ ii ].z() = 5.629;
		xyzarray[ ++ii ].x() = 27.8; xyzarray[ ii ].y() = 9.327; xyzarray[ ii ].z() = 6.816;
		xyzarray[ ++ii ].x() = 27.764; xyzarray[ ii ].y() = 9.899; xyzarray[ ii ].z() = 4.627;
		xyzarray[ ++ii ].x() = 28.581; xyzarray[ ii ].y() = 11.104; xyzarray[ ii ].z() = 4.895;
		xyzarray[ ++ii ].x() = 27.827; xyzarray[ ii ].y() = 12.143; xyzarray[ ii ].z() = 5.711;
		xyzarray[ ++ii ].x() = 28.441; xyzarray[ ii ].y() = 12.951; xyzarray[ ii ].z() = 6.413;
		xyzarray[ ++ii ].x() = 32.79; xyzarray[ ii ].y() = 3.684; xyzarray[ ii ].z() = 2.394;
		xyzarray[ ++ii ].x() = 32.781; xyzarray[ ii ].y() = 4.587; xyzarray[ ii ].z() = 1.247;
		xyzarray[ ++ii ].x() = 33.631; xyzarray[ ii ].y() = 4.007; xyzarray[ ii ].z() = 0.119;
		xyzarray[ ++ii ].x() = 33.276; xyzarray[ ii ].y() = 4.127; xyzarray[ ii ].z() = -1.06;
		xyzarray[ ++ii ].x() = 34.744; xyzarray[ ii ].y() = 3.359; xyzarray[ ii ].z() = 0.487;
		xyzarray[ ++ii ].x() = 35.528; xyzarray[ ii ].y() = 2.641; xyzarray[ ii ].z() = -0.503;
		xyzarray[ ++ii ].x() = 34.775; xyzarray[ ii ].y() = 1.477; xyzarray[ ii ].z() = -1.154;
		xyzarray[ ++ii ].x() = 34.762; xyzarray[ ii ].y() = 1.368; xyzarray[ ii ].z() = -2.384;
		xyzarray[ ++ii ].x() = 34.157; xyzarray[ ii ].y() = 0.617; xyzarray[ ii ].z() = -0.339;
		xyzarray[ ++ii ].x() = 33.395; xyzarray[ ii ].y() = -0.5; xyzarray[ ii ].z() = -0.869;
		xyzarray[ ++ii ].x() = 32.315; xyzarray[ ii ].y() = 0.014; xyzarray[ ii ].z() = -1.829;
		xyzarray[ ++ii ].x() = 32.101; xyzarray[ ii ].y() = -0.587; xyzarray[ ii ].z() = -2.914;
		xyzarray[ ++ii ].x() = 31.68; xyzarray[ ii ].y() = 1.127; xyzarray[ ii ].z() = -1.461;
		xyzarray[ ++ii ].x() = 30.572; xyzarray[ ii ].y() = 1.665; xyzarray[ ii ].z() = -2.285;
		xyzarray[ ++ii ].x() = 31.118; xyzarray[ ii ].y() = 2.266; xyzarray[ ii ].z() = -3.564;
		xyzarray[ ++ii ].x() = 30.494; xyzarray[ ii ].y() = 2.15; xyzarray[ ii ].z() = -4.621;
		xyzarray[ ++ii ].x() = 32.313; xyzarray[ ii ].y() = 2.852; xyzarray[ ii ].z() = -3.496;
		xyzarray[ ++ii ].x() = 32.965; xyzarray[ ii ].y() = 3.336; xyzarray[ ii ].z() = -4.719;
		xyzarray[ ++ii ].x() = 33.302; xyzarray[ ii ].y() = 2.146; xyzarray[ ii ].z() = -5.656;
		xyzarray[ ++ii ].x() = 33.153; xyzarray[ ii ].y() = 2.227; xyzarray[ ii ].z() = -6.893;
		xyzarray[ ++ii ].x() = 33.756; xyzarray[ ii ].y() = 1.035; xyzarray[ ii ].z() = -5.078;
		xyzarray[ ++ii ].x() = 33.989; xyzarray[ ii ].y() = -0.179; xyzarray[ ii ].z() = -5.866;
		xyzarray[ ++ii ].x() = 32.742; xyzarray[ ii ].y() = -0.709; xyzarray[ ii ].z() = -6.553;
		xyzarray[ ++ii ].x() = 32.811; xyzarray[ ii ].y() = -1.185; xyzarray[ ii ].z() = -7.71;
		xyzarray[ ++ii ].x() = 31.598; xyzarray[ ii ].y() = -0.614; xyzarray[ ii ].z() = -5.882;
		xyzarray[ ++ii ].x() = 30.335; xyzarray[ ii ].y() = -0.963; xyzarray[ ii ].z() = -6.545;
		xyzarray[ ++ii ].x() = 30.138; xyzarray[ ii ].y() = -0.109; xyzarray[ ii ].z() = -7.798;
		xyzarray[ ++ii ].x() = 29.796; xyzarray[ ii ].y() = -0.623; xyzarray[ ii ].z() = -8.884;
		xyzarray[ ++ii ].x() = 30.382; xyzarray[ ii ].y() = 1.19; xyzarray[ ii ].z() = -7.67;
		xyzarray[ ++ii ].x() = 30.304; xyzarray[ ii ].y() = 2.088; xyzarray[ ii ].z() = -8.828;
		xyzarray[ ++ii ].x() = 31.296; xyzarray[ ii ].y() = 1.748; xyzarray[ ii ].z() = -9.934;
		xyzarray[ ++ii ].x() = 30.95; xyzarray[ ii ].y() = 1.83; xyzarray[ ii ].z() = -11.133;
		xyzarray[ ++ii ].x() = 32.521; xyzarray[ ii ].y() = 1.368; xyzarray[ ii ].z() = -9.562;
		xyzarray[ ++ii ].x() = 33.52; xyzarray[ ii ].y() = 0.894; xyzarray[ ii ].z() = -10.538;
		xyzarray[ ++ii ].x() = 32.955; xyzarray[ ii ].y() = -0.336; xyzarray[ ii ].z() = -11.275;
		xyzarray[ ++ii ].x() = 33.116; xyzarray[ ii ].y() = -0.479; xyzarray[ ii ].z() = -12.481;
		xyzarray[ ++ii ].x() = 32.299; xyzarray[ ii ].y() = -1.228; xyzarray[ ii ].z() = -10.548;
		xyzarray[ ++ii ].x() = 31.702; xyzarray[ ii ].y() = -2.417; xyzarray[ ii ].z() = -11.196;
		xyzarray[ ++ii ].x() = 30.489; xyzarray[ ii ].y() = -2.092; xyzarray[ ii ].z() = -12.072;
		xyzarray[ ++ii ].x() = 30.35; xyzarray[ ii ].y() = -2.625; xyzarray[ ii ].z() = -13.182;
		xyzarray[ ++ii ].x() = 29.621; xyzarray[ ii ].y() = -1.206; xyzarray[ ii ].z() = -11.59;
		xyzarray[ ++ii ].x() = 28.493; xyzarray[ ii ].y() = -0.755; xyzarray[ ii ].z() = -12.419;
		xyzarray[ ++ii ].x() = 28.953; xyzarray[ ii ].y() = -0.164; xyzarray[ ii ].z() = -13.731;
		xyzarray[ ++ii ].x() = 28.315; xyzarray[ ii ].y() = -0.405; xyzarray[ ii ].z() = -14.779;
		xyzarray[ ++ii ].x() = 30.031; xyzarray[ ii ].y() = 0.611; xyzarray[ ii ].z() = -13.672;
		xyzarray[ ++ii ].x() = 30.602; xyzarray[ ii ].y() = 1.252; xyzarray[ ii ].z() = -14.859;
		xyzarray[ ++ii ].x() = 31.11; xyzarray[ ii ].y() = 0.215; xyzarray[ ii ].z() = -15.869;
		xyzarray[ ++ii ].x() = 30.888; xyzarray[ ii ].y() = 0.391; xyzarray[ ii ].z() = -17.072;
		xyzarray[ ++ii ].x() = 31.812; xyzarray[ ii ].y() = -0.818; xyzarray[ ii ].z() = -15.397;
		xyzarray[ ++ii ].x() = 32.33; xyzarray[ ii ].y() = -1.859; xyzarray[ ii ].z() = -16.269;
		xyzarray[ ++ii ].x() = 31.177; xyzarray[ ii ].y() = -2.626; xyzarray[ ii ].z() = -16.93;
		xyzarray[ ++ii ].x() = 31.241; xyzarray[ ii ].y() = -3.001; xyzarray[ ii ].z() = -18.12;
		xyzarray[ ++ii ].x() = 30.129; xyzarray[ ii ].y() = -2.872; xyzarray[ ii ].z() = -16.163;
		xyzarray[ ++ii ].x() = 28.997; xyzarray[ ii ].y() = -3.685; xyzarray[ ii ].z() = -16.643;
		xyzarray[ ++ii ].x() = 28.262; xyzarray[ ii ].y() = -2.839; xyzarray[ ii ].z() = -17.673;
		xyzarray[ ++ii ].x() = 27.929; xyzarray[ ii ].y() = -3.322; xyzarray[ ii ].z() = -18.745;
		return xyzarray;
	}


};


