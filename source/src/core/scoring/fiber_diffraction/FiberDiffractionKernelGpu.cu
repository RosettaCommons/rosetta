// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffractiobn/FiberDiffractionKernel.cu
/// @brief  FiberDiffraction GPU support
/// @author Wojtek Potrzebowski and Ingemar Andre

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>

#include <stdlib.h>
#include <iostream>
#include <cuda.h>
#include <math.h>
#include <core/types.hh>
#include <sys/time.h>

#include "cutil_math.h"

namespace core {
namespace scoring {
namespace fiber_diffraction {

static basic::Tracer TR("core.scoring.fiber_diffraction.FiberDiffractionKernelGpu");

__global__
void calculate_bessels_kernel(
	int l, 
	int n, 
	int abs_n, 
	int const natoms,
	int const legal_R_values,
	float * d_layer_lines_R,
	Real * d_r, float * d_bessel)
{
	unsigned int atom1 = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int R = blockDim.y * blockIdx.y + threadIdx.y;
	if (atom1 < natoms && R < legal_R_values) {
		Real Rinv(d_layer_lines_R[R]);
		Real x_factor( 2*M_PI*Rinv );
		Real X1 (x_factor*d_r[atom1]);
		if ( abs_n <= X1 +2 ) {
			if (n==0) d_bessel[natoms*R+atom1]=j0f(X1);
			if (n==1) d_bessel[natoms*R+atom1]=j1f(X1);
			if (n==-1) d_bessel[natoms*R+atom1]=-j1f(X1);
			if (n>1) d_bessel[natoms*R+atom1]=jnf(n,X1);
			if (n<-1) d_bessel[natoms*R+atom1]=powf(-1.0,(Real)n)*jnf(-n,X1);
		}//abs_n
	}//atom1 R	
}


__global__
void calculate_bessels_derivatives_kernel(
	int l, 
	int n, 
	int abs_n, 
	int const natoms,
	int const start_R_index, 
	int const legal_R_values, 
	float * d_layer_lines_R, Real * d_r, 
	float * d_bessel, 
	float * d_bessel_plus_1)
{
	unsigned int atom1 = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int R = blockDim.y * blockIdx.y + threadIdx.y;
	if (atom1 < natoms && R < legal_R_values) {
		Real Rinv(d_layer_lines_R[start_R_index+R]);
		Real x_factor( 2*M_PI*Rinv );
		Real X1 (x_factor*d_r[atom1]);
		if ( abs_n <= X1 +2 ) {
			if (n==0) d_bessel[natoms*R+atom1]=j0f(X1);
			if (n==1) d_bessel[natoms*R+atom1]=j1f(X1);
			if (n==-1) d_bessel[natoms*R+atom1]=-j1f(X1);
			if (n>1) d_bessel[natoms*R+atom1]=jnf(n,X1);
			if (n<-1) d_bessel[natoms*R+atom1]=powf(-1.0,(float)n)*jnf(-n,X1);
			//Calculating n+1 bessels
			int l = n+1;
			if (l==0) d_bessel_plus_1[natoms*R+atom1]=j0f(X1);
			if (l==1) d_bessel_plus_1[natoms*R+atom1]=j1f(X1);
			if (l==-1) d_bessel_plus_1[natoms*R+atom1]=-j1f(X1);
			if (l>1) d_bessel_plus_1[natoms*R+atom1]=jnf(l,X1);
			if (l<-1) d_bessel_plus_1[natoms*R+atom1]=powf(-1.0,(float)l)*jnf(-l,X1);
		}//abs_n
	}//atom1 R      
}


__global__
void calculate_phase_kernel(
	int l, 
	int n, 
	int abs_n, 
	int const natoms,
  float const c_, 
	Real * d_phi, 
	Real * d_z, 
	Real * d_r, 
	float * d_phase)
{
	//Main thread indexes for R, atom1, atom2
	unsigned int atom1 = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int atom2 = blockDim.y * blockIdx.y + threadIdx.y;

	if (atom1 < natoms && atom2 <= atom1) {
		d_phase[natoms*atom2+atom1] = (cosf(n*(d_phi[atom2]-d_phi[atom1])+2*M_PI*l/c_*(d_z[atom1]-d_z[atom2])));
	}
}


__global__
void calculate_phase_derivatives(
	int l, 
	int n, 
	int abs_n, 
	int const natoms,
 	float const c_,
	Real * d_phi,
	Real * d_z, 	
	Real * d_r,
	float * d_phase, 
	float * d_phase_prime)
{
	//Main thread indexes for R, atom1, atom2
	unsigned int atom1 = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int atom2 = blockDim.y * blockIdx.y + threadIdx.y;
	if (atom1 < natoms && atom2 <= atom1) {
		float phase_arg = n*(d_phi[atom2]-d_phi[atom1] )+2*M_PI*l/c_*( d_z[atom1]-d_z[atom2]);
		d_phase[natoms*atom2+atom1] = cosf(phase_arg);
		d_phase[natoms*atom1+atom2] = cosf(phase_arg);
		d_phase_prime[natoms*atom2+atom1] = -sinf(phase_arg);
		d_phase_prime[natoms*atom1+atom2] = sinf(phase_arg);
	}
}


__global__
void calculate_intensity_kernel(
	int l, 
	int n, 
	int abs_n, 
	int const natoms,
	int const legal_R_values,
	float const c_, 
	Size * d_atom_type_number,  
	float * d_layer_lines_R, 
	float * d_form_factors,
	Real * d_phi, 
	Real * d_z,
	Real * d_r, 
	float * d_bessel, 
	float * d_phase, 
	float * d_I) 
{
	//Main thread indexes for R, atom1, atom2
	unsigned int atom1 = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int atom2 = blockDim.y * blockIdx.y + threadIdx.y;
	unsigned int R = blockDim.z * blockIdx.z + threadIdx.z;
	
	if (atom1 < natoms && atom2 <= atom1 && R < legal_R_values) {
		//float d_phase(0);
		float Rinv(d_layer_lines_R[R]);
		float x_factor( 2*M_PI*Rinv );
		float X1 (x_factor*d_r[atom1]);
		float X2 (x_factor*d_r[atom2]);
		if ( abs_n <= X1 +2 && abs_n <= X2 +2 ) {
			int atom_type_index1 ( d_atom_type_number[atom1] );
			int atom_type_index2 ( d_atom_type_number[atom2] );
	       
			float f_atom1( d_form_factors[legal_R_values*(atom_type_index1-1)+R] );
			float f_atom2( d_form_factors[legal_R_values*(atom_type_index2-1)+R] );
			//__syncthreads();
			float dummyI( f_atom1*f_atom2*d_bessel[natoms*R+atom1]*d_bessel[natoms*R+atom2]*d_phase[natoms*atom2+atom1] );
			d_I[natoms*(natoms*R+atom1)+atom2]+= (atom1==atom2) ? dummyI : 2*dummyI;
		} //end abs_n condition
	} //end main if statement
}

//This function works as well but it is bit slower
__global__
void calculate_intensity_3d_kernel(int l, int n, int abs_n, int const natoms,
	int const legal_R_values, 
	float const c_, Size * d_atom_type_number,  
	float * d_layer_lines_R, float * d_form_factors,
	Real * d_phi, Real * d_z,
	Real * d_r, float * d_bessel, float * d_phase, float * d_I) 
{
	for (int atom1=blockIdx.x; atom1 < natoms; atom1 += gridDim.x){
		for (int R=threadIdx.y; R<legal_R_values; R+= blockDim.y){
			for (int atom2=threadIdx.x; atom2<=atom1; atom2 += blockDim.x){
      	float Rinv(d_layer_lines_R[R]);
      	float x_factor( 2*M_PI*Rinv );
      	float X1 (x_factor*d_r[atom1]);
      	float X2 (x_factor*d_r[atom2]);
      	if ( abs_n <= X1 +2 && abs_n <= X2 +2 ) {
					int atom_type_index1 ( d_atom_type_number[atom1] );
					int atom_type_index2 ( d_atom_type_number[atom2] );

					float f_atom1( d_form_factors[legal_R_values*(atom_type_index1-1)+R] );
					float f_atom2( d_form_factors[legal_R_values*(atom_type_index2-1)+R] );
					//__syncthreads();
					float dummyI( f_atom1*f_atom2*d_bessel[natoms*R+atom1]*d_bessel[natoms*R+atom2]*d_phase[natoms*atom2+atom1] );
					d_I[natoms*(natoms*R+atom1)+atom2]+= (atom1==atom2) ? dummyI : 2*dummyI;
        }
			}
		}
	}
}
//template< class T >
__global__
void calculate_derivatives_kernel(int l, int n, int abs_n, 
	int const natoms, int const start_R_index, 
	int const legal_R_values, int const max_R_values,
	float const scale_factor_, float const square_obs_,
	float const c_, Real * d_phi, Real * d_z, Real * d_r,
	float * d_layer_lines_R, float * d_layer_lines_I,
	float * d_form_factors, Size * d_atom_type_number,
	float * d_bessel, float * d_bessel_plus_1, float * d_I,
	float * d_phases, float * d_phases_prime, 
	float3 * D, float3 * D_cross_R,
	bool rfactor_refinement) 
{
        //Main thread indexes for R, atom1, atom2
	unsigned int atom1 = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int atom2 = blockDim.y * blockIdx.y + threadIdx.y;
	unsigned int R = blockDim.z * blockIdx.z + threadIdx.z;

	if (atom1 < natoms && atom2 < natoms && R < legal_R_values) {
		float Rinv(d_layer_lines_R[start_R_index+R]);
		float x_factor( 2*M_PI*Rinv );
		float X1 (x_factor*d_r[atom1]);
		float X2 (x_factor*d_r[atom2]);
		if ( abs_n <= X1 +2 && abs_n <= X2 +2 ) {	
			//Bessel values retreived
			float jn1( d_bessel[natoms*R+atom1] );
			float jn1_plus_1( d_bessel_plus_1[natoms*R+atom1] );
			float jn2( d_bessel[natoms*R+atom2] );
		
			float d_phase( d_phases[natoms*atom2+atom1] );
			float d_phase_prime( d_phases_prime[natoms*atom2+atom1] );
                	//__syncthreads();
			int atom_type_index1 ( d_atom_type_number[atom1] );
			float Jn_deriv_atom1( -x_factor*jn1_plus_1+n*jn1/d_r[atom1] );
			if ( fabs(d_r[atom1] ) < 1e-2 ) Jn_deriv_atom1=0.0f;
                                        
			float3 cartesian_coord_atom1 = {d_r[atom1]*cosf(d_phi[atom1]), d_r[atom1]*sinf(d_phi[atom1]), d_z[atom1] };
			float3 unit_r = { cosf(d_phi[atom1]), sinf(d_phi[atom1]), 0.0f };
			float3 unit_x = { 1.0f, 0.0f, 0.0f };
			float3 unit_z = { 0.0f, 0.0f, 1.0f };
			
			//Temporary kernel definition
			float3 D_ker= { 0.0f, 0.0f, 0.0f };
			float3 D_cross_R_ker = { 0.0f, 0.0f, 0.0f };
			
			if (atom1==atom2) {
				float tmp( 2*d_form_factors[max_R_values*(atom_type_index1-1)+start_R_index+R]\
						*d_form_factors[max_R_values*(atom_type_index1-1)+start_R_index+R]\
						*jn1*Jn_deriv_atom1 );
				D_ker = tmp*unit_r;
				D_cross_R_ker = cross(D_ker,cartesian_coord_atom1);
				//__syncthreads();
			}
			
			if (atom1!=atom2) {
				int atom_type_index2 ( d_atom_type_number[atom2] );
				float fact( d_form_factors[max_R_values*(atom_type_index1-1)+start_R_index+R]\
							*d_form_factors[max_R_values*(atom_type_index2-1)+start_R_index+R]*jn2 );

				float dr( Jn_deriv_atom1*fact*d_phase );
				float dphi( n*jn1*fact*d_phase_prime );                
				float dz ( 2*M_PI*l/c_*jn1*fact*d_phase_prime );
                		
				float3 D_tmp = { 0.0f,0.0f,0.0f };
				float3 dphi_vec = { 0.0f,0.0f,0.0f };
				if ( d_r[atom1] >= 1e-2 ) {
					if (fabs(sinf(d_phi[atom1])) > 1e-3)
 						dphi_vec = (unit_x-cosf(d_phi[atom1])*unit_r)/(sinf(d_phi[atom1])*d_r[atom1]);
					else 	{
						float3 dummy_dphi = {0.0f,(-1.0f/d_r[atom1]),0.0f};	
						dphi_vec = dummy_dphi;
					}
				}
				D_tmp = 2.0f*(dr*unit_r + dz*unit_z + dphi*dphi_vec);
				D_ker = D_ker+D_tmp;
				//TODO: Check if it is ok with math formula
				D_cross_R_ker=D_cross_R_ker+cross(D_tmp,cartesian_coord_atom1);
				//__syncthreads();
			}

			float dummy_factor(0); 
			if (rfactor_refinement) {
				if (d_I[start_R_index+R] > 0) {
					float F_diff (scale_factor_*sqrt(d_I[start_R_index+R]) - fabs(d_layer_lines_I[start_R_index+R]));
					dummy_factor = scale_factor_*F_diff/(2*square_obs_*fabs(F_diff)*sqrt(d_I[start_R_index+R])); 
				}
			} 
			else { 
				float dummy3( d_layer_lines_I[start_R_index+R]*d_layer_lines_I[start_R_index+R] );
				float I_diff( scale_factor_*d_I[start_R_index+R] - dummy3 );
				dummy_factor = 2*scale_factor_*I_diff/square_obs_;
			}
			D[natoms*(natoms*R+atom1)+atom2] += dummy_factor*D_ker;
			D_cross_R[natoms*(natoms*R+atom1)+atom2] += dummy_factor*D_cross_R_ker;
		} //end abs_n condition
	} //end main if statement
}


__global__
void sum_intensity_kernel(int const natoms, int const legal_R_values,
	float * d_I, float * d_I_R)
{
	//Main thread indexes for R, atom1, atom2
	unsigned int R = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int atom1 = blockDim.y * blockIdx.y + threadIdx.y;
	unsigned int atom2 = blockDim.z * blockIdx.z + threadIdx.z;
	if (atom1 < natoms && atom2 < natoms && atom2 <= atom1 && R < legal_R_values)
		atomicAdd(&d_I_R[R],d_I[natoms*(natoms*R+atom1)+atom2]);
}

__global__
void sum_intensity_3d_kernel(int const natoms, int const legal_R_values,
														float * d_I_l, float * d_I_R)
{  
	for (int R=blockIdx.x; R< legal_R_values; R += gridDim.x){
		float tmp_I_R(0);
		for (int atom1=threadIdx.y; atom1<natoms; atom1+= blockDim.y){
			for (int atom2=threadIdx.x; atom2<=atom1; atom2 += blockDim.x){
				tmp_I_R += d_I_l[natoms*(natoms*R+atom1)+atom2];
			}
		}
		atomicAdd(&d_I_R[R], tmp_I_R);
	}
}


__global__
void sum_derivatives_chi_3d_kernel(int const natoms,
	int const legal_R_values,
	float3 * d_dchi2_d, float3 * d_dchi2_d_cross_R, 
	float3 * D, float3 * D_cross_R)
{
    
    
	for (int atom1=blockIdx.x; atom1< natoms; atom1 += gridDim.x){
    float3 tmp_dchi2_d = { 0.0f, 0.0f, 0.0f };
    float3 tmp_dchi2_d_cross_R = { 0.0f, 0.0f, 0.0f };
    for (int R=threadIdx.y; R<legal_R_values; R+= blockDim.y){
			for (int atom2=threadIdx.x; atom2<natoms; atom2 += blockDim.x){
				tmp_dchi2_d += D[natoms*(natoms*R+atom1)+atom2];
				tmp_dchi2_d_cross_R+=D_cross_R[natoms*(natoms*R+atom1)+atom2];
			}
		}
		atomicAdd( &d_dchi2_d[atom1].x,tmp_dchi2_d.x );
		atomicAdd( &d_dchi2_d[atom1].y,tmp_dchi2_d.y );
		atomicAdd( &d_dchi2_d[atom1].z,tmp_dchi2_d.z );
		atomicAdd( &d_dchi2_d_cross_R[atom1].x,tmp_dchi2_d_cross_R.x );
		atomicAdd( &d_dchi2_d_cross_R[atom1].y,tmp_dchi2_d_cross_R.y );
		atomicAdd( &d_dchi2_d_cross_R[atom1].z,tmp_dchi2_d_cross_R.z );
	}
}

__global__
void sum_derivatives_chi_kernel(int const natoms, int const legal_R_values,
	float3 * d_dchi2_d, float3 * d_dchi2_d_cross_R, 
	float3 * D, float3 * D_cross_R)
{
	//Main thread indexes for R, atom1, atom2
	unsigned int atom1 = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int atom2 = blockDim.y * blockIdx.y + threadIdx.y;
	unsigned int R = blockDim.z * blockIdx.z + threadIdx.z;

	if (atom1 < natoms && atom2 < natoms && R < legal_R_values) {
		atomicAdd( &d_dchi2_d[atom1].x,D[natoms*(natoms*R+atom1)+atom2].x );
		atomicAdd( &d_dchi2_d[atom1].y,D[natoms*(natoms*R+atom1)+atom2].y );
		atomicAdd( &d_dchi2_d[atom1].z,D[natoms*(natoms*R+atom1)+atom2].z );
		atomicAdd( &d_dchi2_d_cross_R[atom1].x,D_cross_R[natoms*(natoms*R+atom1)+atom2].x );
		atomicAdd( &d_dchi2_d_cross_R[atom1].y,D_cross_R[natoms*(natoms*R+atom1)+atom2].y );
		atomicAdd( &d_dchi2_d_cross_R[atom1].z,D_cross_R[natoms*(natoms*R+atom1)+atom2].z );
	}
}


// Utility function for checking CUDA runtime API results
// can be wrapped around any runtime API call.
void checkCuda(cudaError_t result)
{
	if (result != cudaSuccess) {
		TR.Error<<"CUDA Runtime Error: "<<cudaGetErrorString(result);
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
  }
}

void  calculate_intensity_gpu(
	Size const l_max,
	Size const natoms,
	utility::vector0 < utility::vector0 < int > >::iterator & nvals,
	utility::vector0< utility::vector1< Real > >::iterator & layer_lines_R,	
	utility::vector0< utility::vector1< Real > > & I,
	utility::vector0< utility::vector1< utility::vector1< Real > > >::iterator & form_factors,
	utility::vector1< Real > & phi,
	utility::vector1< Real > & z,
	utility::vector1< Real > & r,
	utility::vector1< Size > & atom_type_number,
	Real const c_, 
	Real const res_cutoff_low_, 
	Real const res_cutoff_high_,
	int const gpu_processor_)
{
	//Cuda error retruned by checkCuda function
	cudaError_t error;
	
	checkCuda ( cudaSetDevice( gpu_processor_ ) );
	//TODO: observe!
	//checkCuda ( cudaDeviceReset() );
	
	//GPU device variables 
	Real * d_phi;
	Real * d_z;
	Real * d_r;
	Size * d_atom_type_number;


	//Input variables initialized by H2D memcpy. 
	checkCuda ( cudaMalloc((void **)&d_phi, sizeof( Real ) * natoms) );
	checkCuda ( cudaMalloc((void **)&d_z, sizeof( Real ) * natoms) );
	checkCuda ( cudaMalloc((void **)&d_atom_type_number, sizeof( Size ) *  natoms) );
	checkCuda ( cudaMalloc((void **)&d_r, sizeof( Real ) * natoms) );

	//Copying values from CPU vectors to previously initialized variables on the GPU
	checkCuda ( cudaMemcpy(d_phi, &phi[ 1 ], sizeof( Real ) * natoms, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy(d_z, &z[ 1 ], sizeof( Real ) * natoms, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy(d_r, &r[ 1 ], sizeof( Real ) * natoms, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy(d_atom_type_number, &atom_type_number[ 1 ], sizeof( Size ) * natoms, cudaMemcpyHostToDevice) );

	float totalElapsedTime_bessel(0);
	float totalElapsedTime_main(0);
	float totalElapsedTime_reduct(0);
	float totalElapsedTime_phase(0);
	float Cycle;
	cudaEvent_t t1_bessel,t2_bessel, t1_main, t2_main, t1_reduct, t2_reduct, t1_phase, t2_phase;
	cudaEventCreate(&t1_bessel);
	cudaEventCreate(&t2_bessel);
	cudaEventCreate(&t1_main);
	cudaEventCreate(&t2_main);
	cudaEventCreate(&t1_reduct);
	cudaEventCreate(&t2_reduct);
	cudaEventCreate(&t1_phase);
	cudaEventCreate(&t2_phase);
        
        
	for ( Size l=0; l <= l_max; ++l ) {
		
		Size max_b_order( nvals[l].size() );
		Size max_R_values( layer_lines_R[l].size());
                
		//Intensity per each layer line
		float * h_I_R;
		float * h_layer_lines_R_l;
		float * h_form_factors_l;
             
		if (max_R_values==0) continue;
                
		h_I_R =  (float * ) malloc( max_R_values * sizeof( float ));
		h_layer_lines_R_l =  ( float * ) malloc(max_R_values *  sizeof( float ));
		h_form_factors_l =  ( float * ) malloc(5 * max_R_values *  sizeof( float ));	
		
		Size t_count(0);		
		for ( Size atom=0; atom<5; ++atom ) {
			t_count = 0;
			for ( Size R=0; R<max_R_values; ++R ) {
				h_form_factors_l[max_R_values*atom+t_count]=form_factors[l][atom+1][R+1];
				t_count++;
			}
		}
		t_count =0;	
		for ( Size R=0; R< max_R_values; ++R ) {
			h_layer_lines_R_l[t_count]=layer_lines_R[l][R+1];
			t_count++;
		}
		
		//GPU device variables. It's bit faster when using float instead of Real	
		float * d_I_l;
		float * d_I_R;
		float * d_layer_lines_R_l;
		float * d_form_factors_l;
		float * d_bessels;
		float * d_phases;
	
		checkCuda ( cudaMalloc((void **)&d_layer_lines_R_l, sizeof( float ) * max_R_values) );
		checkCuda ( cudaMalloc((void **)&d_form_factors_l, sizeof( float ) * 5 * max_R_values) );
		checkCuda ( cudaMemcpy(d_layer_lines_R_l, h_layer_lines_R_l, sizeof( float ) * max_R_values, cudaMemcpyHostToDevice) );
		checkCuda ( cudaMemcpy(d_form_factors_l, h_form_factors_l, sizeof( float ) * 5 * max_R_values, cudaMemcpyHostToDevice) );	
		
		checkCuda ( cudaMalloc((void **)&d_bessels, sizeof( float ) * natoms * max_R_values) );
		checkCuda ( cudaMalloc((void **)&d_phases, sizeof( float ) * natoms * natoms) );
		
		//Output intensity
		checkCuda ( cudaMalloc((void **)&d_I_l, sizeof( float ) * max_R_values * natoms * natoms ) );
		checkCuda ( cudaMemset(d_I_l, 0.0, max_R_values * natoms * natoms * sizeof( float )));    
		checkCuda ( cudaMalloc((void **)&d_I_R, sizeof( float ) * max_R_values) );
		checkCuda ( cudaMemset(d_I_R, 0.0, max_R_values * sizeof( float )));
		//cudaDeviceSynchronize();
		
               
		for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {

			int n( nvals[l][b_order-1] );
			int abs_n( abs(n) );

			cudaEventRecord(t1_phase,0);
			//Bessel function calculation (these values are precalcualated and store in GPU mem)
			dim3 threadsPerBlockPhase(8,16);
			dim3 numBlocksPhase((natoms+threadsPerBlockPhase.x-1)/threadsPerBlockPhase.x,
                     (natoms+threadsPerBlockPhase.y-1)/threadsPerBlockPhase.y);
			calculate_phase_kernel<<<numBlocksPhase, threadsPerBlockPhase >>>(l, n, abs_n, natoms, c_,
													d_phi, d_z, d_r, d_phases);
     error = cudaGetLastError();
     if ( error != 0 ) {
				TR.Error<<"Problem running phase kernel! "<<cudaGetErrorString(error)<<std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}

			cudaEventRecord(t2_phase,0);
			cudaEventSynchronize(t2_phase);
			cudaEventElapsedTime(&Cycle,t1_phase,t2_phase);
			totalElapsedTime_phase +=Cycle;

			cudaEventRecord(t1_bessel,0);

			//Bessel function calculation (these values are precalcualated and store in GPU mem)
			dim3 threadsPerBlockBessel(8,16);
			dim3 numBlocksBessel((natoms+threadsPerBlockBessel.x-1)/threadsPerBlockBessel.x,
			(max_R_values+threadsPerBlockBessel.y-1)/threadsPerBlockBessel.y);
			calculate_bessels_kernel<<<numBlocksBessel, threadsPerBlockBessel>>>(l, n, abs_n, natoms,
                                                                       	max_R_values, d_layer_lines_R_l,
                                                                        d_r, d_bessels);
			
			error = cudaGetLastError();
			if ( error != 0 ) {
				TR.Error<<"Problem running bessel kernel! "<<cudaGetErrorString(error)<<std::endl; 
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}
			//////////////////
			cudaEventRecord(t2_bessel,0);
			cudaEventSynchronize(t2_bessel);
			cudaEventElapsedTime(&Cycle,t1_bessel,t2_bessel);
			totalElapsedTime_bessel +=Cycle;
			/////////////////
			cudaEventRecord(t1_main,0);

			//Main GPU function - calculating intensity			
			dim3 threadsPerBlock(2, 8, 8);
			dim3 numBlocks((natoms+threadsPerBlock.x-1)/threadsPerBlock.x,
				(natoms+threadsPerBlock.y-1)/threadsPerBlock.y,
				(max_R_values+threadsPerBlock.z-1)/threadsPerBlock.z);

			calculate_intensity_kernel<<<numBlocks, threadsPerBlock>>>(l, n, abs_n, natoms,
									max_R_values,
									c_, d_atom_type_number, d_layer_lines_R_l, 
									 d_form_factors_l,
									d_phi, d_z, d_r, d_bessels, d_phases, d_I_l);
                        
			error = cudaGetLastError();
			if ( error != 0 ) {
				TR.Error<<"Problem running intensity kernel! "<<cudaGetErrorString(error)<<std::endl; 
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}

			cudaEventRecord(t2_main,0);
			cudaEventSynchronize(t2_main);
			cudaEventElapsedTime(&Cycle,t1_main,t2_main);
			totalElapsedTime_main +=Cycle;
		}//end b_order
                
		cudaEventRecord(t1_reduct,0);
                
		//Horizontal and vertical summing of intensity
		/*dim3 threadsPerBlockSum(32, 4, 4);
		dim3 numBlocksSum((legal_R_values+threadsPerBlockSum.x-1)/threadsPerBlockSum.x,
		(natoms+threadsPerBlockSum.y-1)/threadsPerBlockSum.y,
		(natoms+threadsPerBlockSum.z-1)/threadsPerBlockSum.z);*/

		dim3 threadsPerBlockSum(8, 64);
		dim3 numBlocksSum(max_R_values);
		sum_intensity_3d_kernel<<<numBlocksSum, threadsPerBlockSum>>>(natoms, max_R_values, d_I_l, d_I_R);
		error = cudaGetLastError();
		if ( error != 0 ) {
			TR.Error<<"Problem running kernel! "<<cudaGetErrorString(error)<<std::endl; 
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
    
    cudaEventRecord(t2_reduct,0);
    cudaEventSynchronize(t2_reduct);
    cudaEventElapsedTime(&Cycle,t1_reduct,t2_reduct);
    totalElapsedTime_reduct +=Cycle;
                
		checkCuda( cudaMemcpy(h_I_R, d_I_R ,max_R_values*sizeof(float),cudaMemcpyDeviceToHost) );
		cudaDeviceSynchronize();
		t_count =0;
		for ( Size R=0; R< max_R_values; ++R ) {
			I[l][R+1] = h_I_R[t_count];
			t_count++;
		}		

		cudaFree( d_layer_lines_R_l );
		cudaFree( d_form_factors_l );
		cudaFree( d_I_l );
		cudaFree( d_I_R );
		cudaFree( d_bessels );
		cudaFree( d_phases );
		free(h_I_R);
		free(h_layer_lines_R_l);
		free(h_form_factors_l);
	} //end l
	
	TR << " Phase, bessel,  main, reduction times " <<totalElapsedTime_phase \
                <<", "<<totalElapsedTime_bessel<<", "<< totalElapsedTime_main \
                <<", "<< totalElapsedTime_reduct<< " ms.\n";
	//Freeing memory for the device variables
	cudaFree( d_phi );
	cudaFree( d_z );
	cudaFree( d_r);
	cudaFree( d_atom_type_number );
	//free(streams);
}

void  calculate_derivatives_gpu(
	Size const l_max,
	Size const natoms,
	utility::vector0< utility::vector0 < int > >::iterator & nvals,
	utility::vector0< utility::vector1< Real > >::iterator & layer_lines_R,
	utility::vector0< utility::vector1< Real > >::iterator & layer_lines_I,
	utility::vector0< utility::vector1< Real > > & I,
	utility::vector0< utility::vector1< utility::vector1< Real > > >::iterator & form_factors,
	utility::vector1< Real > & phi,
	utility::vector1< Real > & z,
	utility::vector1< Real > & r,
	utility::vector1< Size > & atom_type_number,
	utility::vector1< numeric::xyzVector< core::Real > > & dchi2_d, 
	utility::vector1< numeric::xyzVector< core::Real > > & dchi2_d_cross_R,
	Real const c_,
	Real const res_cutoff_low_,
	Real const res_cutoff_high_,
	Real const scale_factor_, 
	Real const square_obs_,
	int const gpu_processor_,
	bool rfactor_refinement)
{
	//Cuda error retruned by checkCuda function
	cudaError_t error;
	checkCuda ( cudaSetDevice( gpu_processor_ ) );        

	//It is not necessary to reset but if there is sth wrong the it will raise up an error
	checkCuda ( cudaDeviceReset() );

	//GPU device variables 
	Real * d_phi;
	Real * d_z;
	Real * d_r;
	Size * d_atom_type_number;

	//Input variables initialized by H2D memcpy. 
	checkCuda ( cudaMalloc((void **)&d_phi, sizeof( Real ) * natoms) );
	checkCuda ( cudaMalloc((void **)&d_z, sizeof( Real ) * natoms) );
	checkCuda ( cudaMalloc((void **)&d_atom_type_number, sizeof( Size ) *  natoms) );
	checkCuda ( cudaMalloc((void **)&d_r, sizeof( Real ) * natoms) );

	//Copying values from CPU vectors to previously initialized variables on the GPU
	checkCuda ( cudaMemcpy(d_phi, &phi[ 1 ], sizeof( Real ) * natoms, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy(d_z, &z[ 1 ], sizeof( Real ) * natoms, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy(d_r, &r[ 1 ], sizeof( Real ) * natoms, cudaMemcpyHostToDevice) );
	checkCuda ( cudaMemcpy(d_atom_type_number, &atom_type_number[ 1 ], sizeof( Size ) * natoms, cudaMemcpyHostToDevice) );

	float totalElapsedTime_reduct(0);
	float Cycle;
	cudaEvent_t t1_reduct, t2_reduct;
	cudaEventCreate(&t1_reduct);
	cudaEventCreate(&t2_reduct);

	for ( Size l=0; l <= l_max; ++l ) {
		Size max_b_order( nvals[l].size() );
		Size max_R_values( layer_lines_R[l].size());
		
		//d_dchi and d_dchi_cross_R for each l
		float3 * h_dchi_l;
		float3 * h_dchi_cross_R_l;
		h_dchi_l = (float3 * ) malloc( natoms * sizeof( float3 ));
		h_dchi_cross_R_l = (float3 * ) malloc( natoms * sizeof( float3 ));

		//Intensity per each layer line
		float * h_I_l;
		float * h_layer_lines_R_l;
		float * h_layer_lines_I_l;
		float * h_form_factors_l;
                
		if (max_R_values==0) continue;
                
		h_I_l =  (float * ) malloc( max_R_values * sizeof( float ));
		h_layer_lines_R_l =  ( float * ) malloc(max_R_values *  sizeof( float ));
		h_layer_lines_I_l =  ( float * ) malloc(max_R_values *  sizeof( float ));
		h_form_factors_l =  ( float * ) malloc(5 * max_R_values *  sizeof( float ));

		Size t_count(0);
		for ( Size atom=0; atom<5; ++atom ) {
        t_count = 0;
        for ( Size R=0; R<max_R_values; ++R ) {
					h_form_factors_l[max_R_values*atom+t_count]=form_factors[l][atom+1][R+1];
					t_count++;
        }
		}
		t_count =0;
		for ( Size R=0; R< max_R_values; ++R ) {
			h_layer_lines_R_l[t_count]=layer_lines_R[l][R+1];
			h_layer_lines_I_l[t_count]=layer_lines_I[l][R+1];
			h_I_l[t_count] = I[l][R+1];
			t_count++;
		}
		
		//GPU device variables. It's bit faster when using float instead of Real        
		float3 * d_dchi_l;
		float3 * d_dchi_cross_R_l;
		float * d_I_l;
		float * d_layer_lines_R_l;
		float * d_layer_lines_I_l;
 		float * d_form_factors_l;
		float * d_phases;
		float * d_phases_prime;

		checkCuda ( cudaMalloc((void **)&d_layer_lines_R_l, sizeof( float ) * max_R_values) );
		checkCuda ( cudaMalloc((void **)&d_layer_lines_I_l, sizeof( float ) * max_R_values) );
		checkCuda ( cudaMalloc((void **)&d_I_l, sizeof( float ) * max_R_values) );
		checkCuda ( cudaMalloc((void **)&d_form_factors_l, sizeof( float ) * 5 * max_R_values) );
		checkCuda ( cudaMemcpy(d_layer_lines_R_l, h_layer_lines_R_l, sizeof( float ) * max_R_values, cudaMemcpyHostToDevice) );
		checkCuda ( cudaMemcpy(d_layer_lines_I_l, h_layer_lines_I_l, sizeof( float ) * max_R_values, cudaMemcpyHostToDevice) );
		checkCuda ( cudaMemcpy(d_I_l, h_I_l, sizeof( float ) * max_R_values, cudaMemcpyHostToDevice) );
		checkCuda ( cudaMemcpy(d_form_factors_l, h_form_factors_l, sizeof( float ) * 5 * max_R_values, cudaMemcpyHostToDevice) );

		checkCuda ( cudaMalloc((void **)&d_phases, sizeof( float ) * natoms * natoms) );
		checkCuda ( cudaMalloc((void **)&d_phases_prime, sizeof( float ) * natoms * natoms) );

		checkCuda ( cudaMalloc((void **)&d_dchi_l, sizeof( float3 ) * natoms ) );
		checkCuda ( cudaMemset(d_dchi_l, 0.0, natoms * sizeof( float3 )));
		checkCuda ( cudaMalloc((void **)&d_dchi_cross_R_l, sizeof( float3 ) * natoms) );
		checkCuda ( cudaMemset(d_dchi_cross_R_l, 0.0, natoms * sizeof( float3 )));

		cudaDeviceSynchronize();
	
		size_t gpu_mem_tot = 0;
		size_t gpu_mem_free = 0;
		cudaMemGetInfo(&gpu_mem_free, &gpu_mem_tot) ;
		error = cudaGetLastError();
    if ( error != 0 ) {
				TR.Error<<"Cannot check memory! "<<cudaGetErrorString(error)<<std::endl; 
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		//Chunk size is set based on R values, 500MB and memory for bessels are left as a buffer and 2 variables has to be allocated 
		Size chunk_size =float(gpu_mem_free - 2*sizeof( float ) * natoms * max_R_values - 524288000)/float((sizeof( float3 ) * natoms * natoms ))/2 ;
		Size number_of_chunks(ceil(float(max_R_values)/float(chunk_size)));
		Size start_R_index;
		Size legal_R_values(floor(float(max_R_values)/float(number_of_chunks)));
		for (Size chunk=0; chunk<number_of_chunks; ++chunk) {
			start_R_index = chunk*legal_R_values;
			//Last chunk takes what is left
			if (chunk==number_of_chunks-1) legal_R_values = max_R_values - chunk*legal_R_values; 
			//D and D_cross_R init are splitted into chunks depending on memory limits
			float * d_bessels;
			float * d_bessels_plus_1;	
			float3 * d_D_l;
			float3 * d_D_cross_R_l;
			
			checkCuda ( cudaMalloc((void **)&d_bessels, sizeof( float ) * natoms * legal_R_values) );
			checkCuda ( cudaMalloc((void **)&d_bessels_plus_1, sizeof( float ) * natoms * legal_R_values) );
			checkCuda ( cudaMalloc((void **)&d_D_l, sizeof( float3 ) * legal_R_values * natoms * natoms ) );
			checkCuda ( cudaMemset(d_D_l, 0.0, legal_R_values * natoms * natoms * sizeof( float3 )));
			checkCuda ( cudaMalloc((void **)&d_D_cross_R_l, sizeof( float3 ) * natoms * natoms * legal_R_values) );
			checkCuda ( cudaMemset(d_D_cross_R_l, 0.0, legal_R_values * natoms * natoms * sizeof( float3 )));

			for ( Size b_order=1; b_order <= max_b_order; ++b_order ) {
				int n( nvals[l][b_order-1] );
				int abs_n( abs(n) );

				//TODO: Phase and bessels might be done outside the chunk loop. Time benfit might be slight but memory required
				dim3 threadsPerBlockPhase(8,16);
				dim3 numBlocksPhase((natoms+threadsPerBlockPhase.x-1)/threadsPerBlockPhase.x,
        (natoms+threadsPerBlockPhase.y-1)/threadsPerBlockPhase.y);
				calculate_phase_derivatives<<<numBlocksPhase, threadsPerBlockPhase>>>(l, n, abs_n, natoms, c_,
                                                        d_phi, d_z, d_r,
                                                        d_phases, d_phases_prime);
				cudaDeviceSynchronize();
				error = cudaGetLastError();
				if ( error != 0 ) {
					TR.Error<<"Problem running kernel! "<<cudaGetErrorString(error)<<std::endl; 
					utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
				}

				//Bessel function calculation (these values are precalcualated and store in GPU mem)
				dim3 threadsPerBlockBessel(8,16);
				dim3 numBlocksBessel((natoms+threadsPerBlockBessel.x-1)/threadsPerBlockBessel.x,
        (legal_R_values+threadsPerBlockBessel.y-1)/threadsPerBlockBessel.y);
				calculate_bessels_derivatives_kernel<<<numBlocksBessel, threadsPerBlockBessel>>>(l, n, abs_n, natoms,
                                        start_R_index, legal_R_values, d_layer_lines_R_l,
                                        d_r, d_bessels, d_bessels_plus_1);
				cudaDeviceSynchronize();
				error = cudaGetLastError();
				if ( error != 0 ) {
					TR.Error<<"Problem running kernel! "<<cudaGetErrorString(error)<<std::endl; 
					utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
				}

				//Main GPU function - calculating intensity                     
				dim3 threadsPerBlock(2, 8, 8);
				dim3 numBlocks((natoms+threadsPerBlock.x-1)/threadsPerBlock.x,
      			(natoms+threadsPerBlock.y-1)/threadsPerBlock.y,
      			(legal_R_values+threadsPerBlock.z-1)/threadsPerBlock.z);

				calculate_derivatives_kernel<<<numBlocks, threadsPerBlock>>>(l, n, abs_n, natoms, 
									start_R_index, legal_R_values, max_R_values,
									scale_factor_, square_obs_,
									c_, d_phi, d_z, d_r,
									d_layer_lines_R_l, d_layer_lines_I_l,
									d_form_factors_l, d_atom_type_number,
									d_bessels, d_bessels_plus_1, d_I_l,
									d_phases, d_phases_prime,
									d_D_l, d_D_cross_R_l,
									rfactor_refinement);
									//d_dchi_all, d_dchi_cross_R_all);
				cudaDeviceSynchronize();
				error = cudaGetLastError();
				if ( error != 0 ) {
					TR.Error<<"Problem running main kernel! "<<cudaGetErrorString(error)<<std::endl; 
					utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
				}

			}//end b_order
			//Horizontal and vertical summing of intensity
			cudaEventRecord(t1_reduct,0);

			dim3 threadsPerBlockSum(16, 8);
			dim3 numBlocksSum((natoms+threadsPerBlockSum.x-1)/threadsPerBlockSum.x);
			//dim3 numBlocksSum(natoms);
			sum_derivatives_chi_3d_kernel<<<numBlocksSum,threadsPerBlockSum>>>(natoms, legal_R_values, 
					d_dchi_l, d_dchi_cross_R_l, d_D_l, d_D_cross_R_l);
			cudaDeviceSynchronize();
			error = cudaGetLastError();
			if ( error != 0 ) {
				TR.Error<<"Problem running summation kernel! "<<cudaGetErrorString(error)<<std::endl; 
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}

			cudaEventRecord(t2_reduct,0);
			cudaEventSynchronize(t2_reduct);
			cudaEventElapsedTime(&Cycle,t1_reduct,t2_reduct);
			totalElapsedTime_reduct +=Cycle;
			cudaFree( d_D_l );
			cudaFree( d_D_cross_R_l );
			cudaFree( d_bessels );
			cudaFree( d_bessels_plus_1 );
		}//end chunks
		checkCuda( cudaMemcpy(h_dchi_l, d_dchi_l ,natoms*sizeof(float3),cudaMemcpyDeviceToHost) );
		checkCuda( cudaMemcpy(h_dchi_cross_R_l, d_dchi_cross_R_l, natoms*sizeof(float3),cudaMemcpyDeviceToHost) );
             
		for ( Size atom1=0; atom1< natoms; ++atom1 ) {
			numeric::xyzVector< core::Real> dummy_dchi( h_dchi_l[atom1].x, h_dchi_l[atom1].y, h_dchi_l[atom1].z);
			numeric::xyzVector< core::Real> dummy_dchi_cross_R( h_dchi_cross_R_l[atom1].x, h_dchi_cross_R_l[atom1].y, h_dchi_cross_R_l[atom1].z );
			dchi2_d[atom1+1] += dummy_dchi;
			dchi2_d_cross_R[atom1+1] += dummy_dchi_cross_R; 
		}
	
		cudaFree( d_layer_lines_R_l );
		cudaFree( d_layer_lines_I_l );
		cudaFree( d_form_factors_l );
		cudaFree( d_I_l );
		cudaFree( d_dchi_l );
		cudaFree( d_dchi_cross_R_l );
		cudaFree( d_phases );
		cudaFree( d_phases_prime );
		free(h_dchi_l);
		free(h_dchi_cross_R_l);
		free(h_I_l);
		free(h_layer_lines_R_l);
		free(h_form_factors_l);
	}//end l
	TR << "reduction times " << totalElapsedTime_reduct<< " ms.\n";
	cudaFree( d_phi);
	cudaFree( d_r);
	cudaFree( d_z);
	cudaFree( d_atom_type_number );
}

} // namespace fiber_diffraction
} // namespace scoring
} // namespace core

