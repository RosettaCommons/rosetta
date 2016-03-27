// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/gpu/GPU.hh
/// @brief  OpenCL-based GPU scheduler class
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

#ifndef INCLUDED_basic_gpu_GPU_hh
#define INCLUDED_basic_gpu_GPU_hh

#ifdef USEOPENCL

#include <string>
#include <map>
#include <vector>

#ifdef MACOPENCL
#include <OpenCL/cl_platform.h>
#include <OpenCL/opencl.h>
#include <OpenCL/cl.h>
#else
#include <CL/cl_platform.h>
#include <CL/opencl.h>
#include <CL/cl.h>
#endif

#endif // USEOPENCL


namespace basic {
namespace gpu {


#ifndef USEOPENCL

// provide a non-GPU definition of basic::gpu::float4 for function definitions that can be GPU or non-GPU
struct float4
{
	float x, y, z, w;
};

#else // USEOPENCL

struct float4
{
	cl_float x, y, z, w;
};

#define GPU_MEM_RO (CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR)
#define GPU_MEM_RW (CL_MEM_READ_WRITE)
#define GPU_IN 0x01
#define GPU_OUT 0x02
#define GPU_IN_OUT (GPU_IN|GPU_OUT)

#define GPU_ARG_TYPE 0xFF00
#define GPU_MEM 0x0000
#define GPU_DEVMEM 0x0100
#define GPU_INT 0x0200
#define GPU_FLOAT 0x0300
#define GPU_DOUBLE 0x0400

struct GPU_DEV
{
	int ndevice;
	char name[128];
	size_t threads;
	cl_uint clockRate;
	cl_uint multiProcessorCount;

	int initialized;
	int usage;

	cl_device_id device;
	cl_command_queue commandQueue;

};

struct GPU_KERNEL_ARG
{
	int type;		// argument type GPU_xxx
	int size;		// argument size (memory size)
	union {
		float f;	// direct float
		int i;		// direct int
		void *ptr;
	};
	cl_mem mem;		// device memory
	size_t k_size;	// kernel argument size
	void *k_p;		// kernel argument pointer
};

class GPU {

private:
	// Device resources
	static GPU_DEV devices_[8];	// max GPUs per system
	static cl_context context_;
	static std::map<std::string,cl_kernel> kernels_;
	static std::map<std::string,cl_program> programs_;

	// Private memory objects are discared in destructor
	static std::vector<cl_mem> privateMemoryObjects_;
	// Shared memory objects are shared between instances and retained
	static std::map<std::string,cl_mem> sharedMemoryObjects_;

	int ndevice_;
	int initialized_;
	cl_int errNum_;
	float kernelRuntime_;
	int profiling_;

public:
	GPU(int ndevice =0);
	~GPU();

	int use();

	static const char *errstr(int errorCode);

	int Init();
	int Release();
	static int Free();

	int RegisterProgram(const char *filename);
	int RegisterProgram(std::string& filename);
	int RegisterProgram(std::vector<std::string> & files);

	cl_kernel BuildKernel(const char *kernel_name);

	cl_mem AllocateMemory(unsigned int size, void *data =NULL, int flags =0, const char *name =NULL);
	cl_mem AllocateMemoryReuse(cl_mem &old_mem, unsigned int &old_size, unsigned int new_size, int flags =0);
	cl_mem GetSharedMemory(const char *name);
	void Free(cl_mem h);

	int ExecuteKernel(const char *kernel_name, int total_threads, int max_conc_threads_high, int max_conc_threads_low, ...);

	int ReadData(void *dst, cl_mem src, unsigned int size, int blocking =CL_TRUE);
	int WriteData(cl_mem dst, void *src, unsigned int size, int blocking =CL_TRUE);
	int setKernelArg(cl_kernel kernel, 	cl_uint i, std::size_t size, const void *p);

	// Getters & setters
	int initialized() { return initialized_; }
	const GPU_DEV &device() { return devices_[ndevice_]; }
	int lastError() { return errNum_; }
	const char *lastErrorStr() { return errstr(errNum_); }
	float lastKernelRuntime() { return kernelRuntime_; }
	int profiling() { return profiling_; }
	void profiling(int i) { profiling_ = i; }

protected:
	GPU_DEV &this_device() { return devices_[ndevice_]; }
	int _ExecuteKernel(const char *kernel_name, int total_threads, int max_conc_threads_high, int max_conc_threads_low, GPU_KERNEL_ARG *args, int async =0);
};

#endif // USEOPENCL

} // gpu
} // basic

#endif // INCLUDED_basic_gpu_GPU_hh
