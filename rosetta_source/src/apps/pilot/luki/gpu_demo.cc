// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/luki/gpu_demo.cc
/// @brief  A simple GPU demo for RosettaCon tutorial session
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

/*
 * quick build:
g++ -o build/src/release/linux/2.6/64/x86/gcc/4.4/opencl/apps/pilot/luki/gpu_demo.o -c -std=c++98 -isystem external/boost_1_46_1/ -pipe -ffor-scope -Wall -Wextra -pedantic -Wno-long-long -O3 -ffast-math -funroll-loops -finline-functions -s -Wno-unused-variable -finline-limit=487 -DNDEBUG -DUSEOPENCL -Isrc -Iexternal/include -Isrc/platform/linux/64/gcc/4.4 -Isrc/platform/linux/64/gcc -Isrc/platform/linux/64 -Isrc/platform/linux -Iexternal/boost_1_46_1 -Iexternal/dbio -I/usr/include -I/usr/local/include -I/usr/local/cuda/include src/apps/pilot/luki/gpu_demo.cc
g++ -o build/src/release/linux/2.6/64/x86/gcc/4.4/opencl/gpu_demo.opencl.linuxgccrelease -Wl,-rpath=/home/luki/rosetta/mini-darc/build/src/release/linux/2.6/64/x86/gcc/4.4/opencl -Wl,-rpath=/home/luki/rosetta/mini-darc/build/external/release/linux/2.6/64/x86/gcc/4.4/opencl build/src/release/linux/2.6/64/x86/gcc/4.4/opencl/apps/pilot/luki/gpu_demo.o -Lexternal/lib -Lbuild/src/release/linux/2.6/64/x86/gcc/4.4/opencl -Lsrc -Lbuild/external/release/linux/2.6/64/x86/gcc/4.4/opencl -Lexternal -L/usr/lib -L/usr/local/lib -L/usr/local/cuda/lib64 -L/usr/local/cuda/lib -ldevel -lprotocols.7 -lprotocols.6 -lprotocols_f.5 -lprotocols_e.5 -lprotocols_d.5 -lprotocols_c.5 -lprotocols_b.5 -lprotocols_a.5 -lprotocols_h.4 -lprotocols_g.4 -lprotocols_f.4 -lprotocols_e.4 -lprotocols_d.4 -lprotocols_c.4 -lprotocols_b.4 -lprotocols_a.4 -lprotocols.3 -lprotocols_b.2 -lprotocols_a.2 -lprotocols.1 -lcore.5 -lcore.4 -lcore.3 -lcore.2 -lcore.1 -lbasic -lnumeric -lutility -lObjexxFCL -lz -lcppdb -lsqlite3 -lOpenCL
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <devel/init.hh>

#include <core/scoring/sasa.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <numeric/constants.hh>

#ifdef USEOPENCL
#include <utility/GPU.hh>

using namespace std;
using namespace utility;

///////////////////////////////////////////////////////////////////////////////

class Timer {
	timespec start, end;
	const char *tag_;
	public:
	Timer(const char *tag =NULL) {
		tag_ = tag;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
	}
	~Timer() {
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
		double time = (end.tv_sec - start.tv_sec)*1000000 + (end.tv_nsec - start.tv_nsec)/1000;
		if(tag_)
			cout << "Time [" << tag_ << "]: ";
		else
			cout << "Time: ";
		cout << (time/1000) << " ms" << endl;
	}
};

///////////////////////////////////////////////////////////////////////////////

#define NARRAY(x) (sizeof(x)/sizeof(*x))

void simple_add() {
	cout << "========== SIMPLE ADD =========" << endl;

	GPU g1;
	g1.RegisterProgram("demo.cl");

	float a[] = { 5, 7, 13, 21, 66 };
	float b[] = { 0, 2, 5, 17, 77 };
	float c[NARRAY(a)];
	int n = NARRAY(a);

	g1.ExecuteKernel("add", n, 512, 32,
		GPU_IN, sizeof(a), a,
		GPU_IN, sizeof(b), b,
		GPU_OUT, sizeof(c), c,
		GPU_FLOAT, 0.7,
		NULL);

	cout << "Kernel run time: " << g1.lastKernelRuntime() << " ms " << endl;

	for(int i=0; i < n; i++)
		printf("%7.3f + %7.3f + 0.7 = %7.3f\n", a[i], b[i], c[i]);
}

void angle_to_xy(int n) {
	cout << "========== ANGLE to X, Y =========" << endl;

	float *angle = new float[n];
	float *x = new float[n];
	float *y = new float[n];

	for(int i=0; i < n; i++)
		angle[i] = numeric::constants::f::pi/180 * i;

	{
		Timer t1("On GPU");
		GPU g1;
		g1.RegisterProgram("demo.cl");

		{
		Timer t2("Exec Kernel");
		if(!g1.ExecuteKernel("angle_xy", n, 512, 32,
			GPU_IN, sizeof(*angle)*n, angle,
			GPU_OUT, sizeof(*x)*n, x,
			GPU_OUT, sizeof(*y)*n, y,
			NULL))
			  return;
		}

		cout << "Kernel run time: " << g1.lastKernelRuntime() << " ms " << endl;
	}

	for(int i=0; i < n; i++) {
		printf("%7.3f => %7.3f, %7.3f\n", angle[i], x[i], y[i]);
		if(i > 10) { cout << "..." << endl; break; }
	}


	{
		Timer t3("On CPU");
		for(int i=0; i < n; i++) {
			x[i] = sinf(angle[i]);
			y[i] = cosf(angle[i]);
		}
	}

	delete angle;
	delete x;
	delete y;
}

void two_gpu_demo() {

	cout << "========== TWO GPU DEMO =========" << endl;

	GPU g1(1), g2(2);
	g1.RegisterProgram("demo.cl");

	float a[] = { 5, 7, 13, 21, 66 };
	float b[] = { 0, 2, 5, 17, 77 };
	float c[NARRAY(a)];
	int n = NARRAY(a);

	// Allocate GPU memory and copy a, b to GPU memory
	cl_mem d_a = g1.AllocateMemory(sizeof(a), a);
	cl_mem d_b = g1.AllocateMemory(sizeof(b), b);
	cl_mem d_c = g1.AllocateMemory(sizeof(c));

	// Execute kernel (no data transfer)
	g1.ExecuteKernel("add", n, 512, 32,
		GPU_DEVMEM, d_a,
		GPU_DEVMEM, d_b,
		GPU_DEVMEM, d_c,
		GPU_FLOAT, 0.0,
		NULL);

	cout << "First kernel run time: " << g1.lastKernelRuntime() << " ms " << endl;

	// Execute second kernel (no data transfer)
	g2.ExecuteKernel("compute_sin", n, 512, 32,
		GPU_DEVMEM, d_c,
		NULL);

	cout << "Second kernel run time: " << g1.lastKernelRuntime() << " ms " << endl;

	// Read back results
	g2.ReadData(c, d_c, sizeof(c));

	for(int i=0; i < n; i++)
		printf("sin(%7.3f + %7.3f) = %7.3f == %7.3f\n", a[i], b[i], c[i], sin(a[i]+b[i]));
}

void Check_for_intersection(
	const float4 *rays,
	const float4 *atoms,
	float *distances,
	const unsigned int NUM_ATOMS,
	const unsigned int NUM_RAYS
)
{
	const int large_dist = 9999.0;
	unsigned int rayID;
	int unsigned atomID, i, j, l;
	float dirX, dirY, dirZ;
	float distance = 9999.0;
	float4 ray;

	for(rayID = 0; rayID < NUM_RAYS; rayID++) {

	// CUT AND PASTE GPU CODE BELOW

	distance = 9999.0;

	ray = rays[rayID];
	dirX = large_dist*sin(ray.x)*cos(ray.y);
	dirY = large_dist*sin(ray.x)*sin(ray.y);
	dirZ = large_dist*cos(ray.x);

	// setup our quadratic equation
	float a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);

	for(atomID = 0; atomID < NUM_ATOMS; atomID ++) {

			float4 atom = atoms[atomID];
			float b = 2.0 * ( (dirX*(-atom.x)) + (dirY*(-atom.y)) + (dirZ*(-atom.z)) );
			float c = atom.x*atom.x + atom.y*atom.y + atom.z*atom.z - (atom.w * atom.w);

			// test for intersection
			float inside_sq = ( b * b ) - ( 4.0 * a * c );
			if (inside_sq > 0) {
				float inside = sqrt(inside_sq);

				float mu1 = -(b-inside) / ( 2.0 * a);
				float x1 =  mu1 * dirX;
				float y1 =  mu1 * dirY;
				float z1 =  mu1 * dirZ;
				float dist1_sq = x1*x1 + y1*y1 + z1*z1;

				float mu2 = -(b+inside) / ( 2.0 * a);
				float x2 = mu2 * dirX;
				float y2 = mu2 * dirY;
				float z2 = mu2 * dirZ;
				float dist2_sq = x2*x2 + y2*y2 + z2*z2;

				if(dist1_sq < distance)
					distance = dist1_sq;
				if(dist2_sq < distance)
					distance = dist2_sq;
			}
	}

	if(distance < 9999.0) {
		distance = sqrt(distance) - ray.z;
		if(distance < 0)
			distance = -distance;
		else
			distance *= 3.12; // check this
	} else
		distance = 0;

	distances[rayID] = distance;

	}
}


// Limits
#define MAX_NUM_RAYS 10000
#define MAX_NUM_ATOMS 100

void darc_test() {

	cout << "========== DARC TEST =========" << endl;

	ifstream inputFile;
	float4 atoms[MAX_NUM_ATOMS];
	float4 rays[MAX_NUM_RAYS];
	float distances[MAX_NUM_RAYS];
	int NUM_ATOMS =0, NUM_RAYS =0;

	/// DEMO CODE from KAREN ///

	// Load rays
	inputFile.open("rays.txt");
	while(inputFile.good()) {
		if ( NUM_RAYS >= MAX_NUM_RAYS ) { cout << "Error, too many rays." << endl; return; }
		float4 ray;
		// phi, psi, rho
		inputFile >> ray.x >> ray.y >> ray.z;
		if(!inputFile.fail() && ray.z != 9999)
			rays[NUM_RAYS++] = ray;
	}
	inputFile.close();
	cout << "Done reading " << NUM_RAYS << " rays." << endl;

	// Load atoms
	inputFile.open("atoms.txt");
	while(inputFile.good()) {
		if ( NUM_ATOMS >= MAX_NUM_ATOMS ) { cout << "Error, too many atoms." << endl; return; }
		float4 atom;
		inputFile >> atom.x >> atom.y >> atom.z >> atom.w;
		if(!inputFile.fail() && atom.w > 0)
			atoms[NUM_ATOMS++] = atom;
	}
	inputFile.close();
	cout << "Done reading " << NUM_ATOMS << " atoms." << endl;

	/// CPU CODE  ///
	{
		Timer t1("DARC CPU Total");
		float darc_score =0, gpu_time =0;

		for(int cycle = 0; cycle < 10000; cycle++) {

			Check_for_intersection(rays, atoms, distances, NUM_ATOMS, NUM_RAYS);

			// Compute total score for ligand
			double temp_gpu_total_dist = 0;
			for(int k = 0; k < NUM_RAYS; ++k) {
				//if(k < 20) cout << k << ", " << distances[k] << endl;
				temp_gpu_total_dist += distances[k];
			}
			temp_gpu_total_dist /= NUM_RAYS;

			if(cycle == 0)
				darc_score = temp_gpu_total_dist;

			// Trivial particle mover for demo purposes
	        for(int j=0;j<NUM_ATOMS;j++){
            	atoms[j].x += 0.05;
            	atoms[j].y += 0.05;
            	atoms[j].z += 0.05;
            }
		}

		// Calc run time for file
		cout << "DARC score for initial conformation: " << darc_score << endl;
	}
}

void darc_test_gpu() {

	cout << "========== DARC GPU TEST =========" << endl;

	// Limits
	#define MAX_NUM_RAYS 10000
	#define MAX_NUM_ATOMS 100

	ifstream inputFile;
	float4 atoms[MAX_NUM_ATOMS];
	float4 rays[MAX_NUM_RAYS];
	float distances[MAX_NUM_RAYS];
	int NUM_ATOMS =0, NUM_RAYS =0;

	/// DEMO CODE from KAREN ///

	// Load rays
	inputFile.open("rays.txt");
	while(inputFile.good()) {
		if ( NUM_RAYS >= MAX_NUM_RAYS ) { cout << "Error, too many rays." << endl; return; }
		float4 ray;
		// phi, psi, rho
		inputFile >> ray.x >> ray.y >> ray.z;
		if(!inputFile.fail() && ray.z != 9999)
			rays[NUM_RAYS++] = ray;
	}
	inputFile.close();
	cout << "Done reading " << NUM_RAYS << " rays." << endl;

	// Load atoms
	inputFile.open("atoms.txt");
	while(inputFile.good()) {
		if ( NUM_ATOMS >= MAX_NUM_ATOMS ) { cout << "Error, too many atoms." << endl; return; }
		float4 atom;
		inputFile >> atom.x >> atom.y >> atom.z >> atom.w;
		if(!inputFile.fail() && atom.w > 0)
			atoms[NUM_ATOMS++] = atom;
	}
	inputFile.close();
	cout << "Done reading " << NUM_ATOMS << " atoms." << endl;

	// GPU stuff
	GPU gpu;
	gpu.profiling(0);
	gpu.Init();
	gpu.RegisterProgram("darc.cl");

	// Allocate persistent memory for rays
	cl_mem d_rays = gpu.AllocateMemory(sizeof(*rays) * NUM_RAYS);
	cl_mem d_atoms = gpu.AllocateMemory(sizeof(*atoms) * MAX_NUM_ATOMS);
	cl_mem d_distances = gpu.AllocateMemory(sizeof(*distances) * NUM_RAYS);

	// Copy ray data to GPU memory (done above)
	gpu.WriteData(d_rays, rays, sizeof(*rays) * NUM_RAYS);

	/// GPU CODE  ///
	{
		Timer t1("DARC GPU Total");
		float darc_score =0, gpu_time =0;

		for(int cycle = 0; cycle < 10000; cycle++) {

			// Load atoms coords to GPU memory
			gpu.WriteData(d_atoms, atoms, sizeof(*atoms) * NUM_ATOMS);

			// Run GPU kernel
			if(!gpu.ExecuteKernel("Check_for_intersection", NUM_RAYS, 256, 64,
				GPU_DEVMEM, d_rays,
				GPU_DEVMEM, d_atoms,
				GPU_DEVMEM, d_distances,
				GPU_INT, NUM_ATOMS,
				NULL)) {
					cerr << "Failed to launch kernel: " << gpu.lastErrorStr() << endl;
					return;
			}

			gpu_time += gpu.lastKernelRuntime();

			// Read computed distances back from GPU
			gpu.ReadData(distances, d_distances, sizeof(*distances) * NUM_RAYS);

			// Compute total score for ligand
			double temp_gpu_total_dist = 0;
			for(int k = 0; k < NUM_RAYS; ++k) {
				//if(k < 20) cout << k << ", " << distances[k] << endl;
				temp_gpu_total_dist += distances[k];
			}
			temp_gpu_total_dist /= NUM_RAYS;

			if(cycle == 0)
				darc_score = temp_gpu_total_dist;

			// Trivial particle mover for demo purposes
	        for(int j=0;j<NUM_ATOMS;j++){
            	atoms[j].x += 0.05;
            	atoms[j].y += 0.05;
            	atoms[j].z += 0.05;
            }
		}

		// Calc run time for file
		cout << "DARC score for initial conformation: " << darc_score << endl;
		cout << "Total GPU time: " << gpu_time << " ms" << endl;
	}
}

void calc_sasa() {

	cout << "========== SASA CALC =========" << endl;

	core::pose::Pose p;
	core::import_pose::pose_from_pdb(p, "demo.pdb");

	{
	Timer t("SASA calc #1");
	core::Real sasa = core::scoring::calc_total_sasa(p, 1.0);
	cout << "SASA: " << sasa << endl;
	}

	{
	Timer t("SASA calc #2");
	core::Real sasa = core::scoring::calc_total_sasa(p, 1.0);
	cout << "SASA: " << sasa << endl;
	}

	{
	Timer t("SASA calc #3");
	core::Real sasa = core::scoring::calc_total_sasa(p, 1.0);
	cout << "SASA: " << sasa << endl;
	}
}


///////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
	devel::init(argc, argv);

	simple_add();
	two_gpu_demo();
	angle_to_xy(1000000);
	//calc_sasa();
	darc_test();
	darc_test_gpu();

	return 0;
}

#else

int main(int argc, char *argv[]) {
	devel::init(argc, argv);
	std::cout << "This build does not contain GPU support." << std::endl;
	return 0;
}

#endif
