// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/gpu/GPU.cc
/// @brief  OpenCL-based GPU scheduler class
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

#ifdef USEOPENCL
#ifndef INCLUDED_basic_gpu_GPU_cc
#define INCLUDED_basic_gpu_GPU_cc

#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>

#include <basic/gpu/GPU.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/gpu.OptionKeys.gen.hh>

#define UPPER_MULTIPLE(n,d) (((n)%(d)) ? (((n)/(d)+1)*(d)) : (n))

static thread_local basic::Tracer TR( "basic.gpu" );

using namespace std;

namespace basic {
namespace gpu {

GPU_DEV GPU::devices_[8];
cl_context GPU::context_ = NULL;
map<string,cl_kernel> GPU::kernels_;
map<string,cl_program> GPU::programs_;
map<string,cl_mem> GPU::sharedMemoryObjects_;
std::vector<cl_mem> GPU::privateMemoryObjects_;

const char *GPU::errstr(int errorCode)
{
	switch(errorCode) {
		case CL_DEVICE_NOT_FOUND:
			return "CL_DEVICE_NOT_FOUND";
		case CL_DEVICE_NOT_AVAILABLE:
			return "CL_DEVICE_NOT_AVAILABLE";
		case CL_COMPILER_NOT_AVAILABLE:
			return "CL_COMPILER_NOT_AVAILABLE";
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		case CL_OUT_OF_RESOURCES:
			return "CL_OUT_OF_RESOURCES";
		case CL_OUT_OF_HOST_MEMORY:
			return "CL_OUT_OF_HOST_MEMORY";
		case CL_PROFILING_INFO_NOT_AVAILABLE:
			return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case CL_MEM_COPY_OVERLAP:
			return "CL_MEM_COPY_OVERLAP";
		case CL_IMAGE_FORMAT_MISMATCH:
			return "CL_IMAGE_FORMAT_MISMATCH";
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:
			return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case CL_BUILD_PROGRAM_FAILURE:
			return "CL_BUILD_PROGRAM_FAILURE";
		case CL_MAP_FAILURE:
			return "CL_MAP_FAILURE";
		case CL_MISALIGNED_SUB_BUFFER_OFFSET:
			return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
		case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
			return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
		case CL_INVALID_VALUE:
			return "CL_INVALID_VALUE";
		case CL_INVALID_DEVICE_TYPE:
			return "CL_INVALID_DEVICE_TYPE";
		case CL_INVALID_PLATFORM:
			return "CL_INVALID_PLATFORM";
		case CL_INVALID_DEVICE:
			return "CL_INVALID_DEVICE";
		case CL_INVALID_CONTEXT:
			return "CL_INVALID_CONTEXT";
		case CL_INVALID_QUEUE_PROPERTIES:
			return "CL_INVALID_QUEUE_PROPERTIES";
		case CL_INVALID_COMMAND_QUEUE:
			return "CL_INVALID_COMMAND_QUEUE";
		case CL_INVALID_HOST_PTR:
			return "CL_INVALID_HOST_PTR";
		case CL_INVALID_MEM_OBJECT:
			return "CL_INVALID_MEM_OBJECT";
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
			return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case CL_INVALID_IMAGE_SIZE:
			return "CL_INVALID_IMAGE_SIZE";
		case CL_INVALID_SAMPLER:
			return "CL_INVALID_SAMPLER";
		case CL_INVALID_BINARY:
			return "CL_INVALID_BINARY";
		case CL_INVALID_BUILD_OPTIONS:
			return "CL_INVALID_BUILD_OPTIONS";
		case CL_INVALID_PROGRAM:
			return "CL_INVALID_PROGRAM";
		case CL_INVALID_PROGRAM_EXECUTABLE:
			return "CL_INVALID_PROGRAM_EXECUTABLE";
		case CL_INVALID_KERNEL_NAME:
			return "CL_INVALID_KERNEL_NAME";
		case CL_INVALID_KERNEL_DEFINITION:
			return "CL_INVALID_KERNEL_DEFINITION";
		case CL_INVALID_KERNEL:
			return "CL_INVALID_KERNEL";
		case CL_INVALID_ARG_INDEX:
			return "CL_INVALID_ARG_INDEX";
		case CL_INVALID_ARG_VALUE:
			return "CL_INVALID_ARG_VALUE";
		case CL_INVALID_ARG_SIZE:
			return "CL_INVALID_ARG_SIZE";
		case CL_INVALID_KERNEL_ARGS:
			return "CL_INVALID_KERNEL_ARGS";
		case CL_INVALID_WORK_DIMENSION:
			return "CL_INVALID_WORK_DIMENSION";
		case CL_INVALID_WORK_GROUP_SIZE:
			return "CL_INVALID_WORK_GROUP_SIZE";
		case CL_INVALID_WORK_ITEM_SIZE:
			return "CL_INVALID_WORK_ITEM_SIZE";
		case CL_INVALID_GLOBAL_OFFSET:
			return "CL_INVALID_GLOBAL_OFFSET";
		case CL_INVALID_EVENT_WAIT_LIST:
			return "CL_INVALID_EVENT_WAIT_LIST";
		case CL_INVALID_EVENT:
			return "CL_INVALID_EVENT";
		case CL_INVALID_OPERATION:
			return "CL_INVALID_OPERATION";
		case CL_INVALID_GL_OBJECT:
			return "CL_INVALID_GL_OBJECT";
		case CL_INVALID_BUFFER_SIZE:
			return "CL_INVALID_BUFFER_SIZE";
		case CL_INVALID_MIP_LEVEL:
			return "CL_INVALID_MIP_LEVEL";
		case CL_INVALID_GLOBAL_WORK_SIZE:
			return "CL_INVALID_GLOBAL_WORK_SIZE";
		case CL_SUCCESS:
			return "(No error)";
			//		case CL_PLATFORM_NOT_FOUND_KHR:
			//			return "CL_PLATFORM_NOT_FOUND_KHR";
		default:
			return "Unknown OpenCL error";
	}
}

GPU::GPU(int ndevice)
{
	initialized_ = 0;
	errNum_ = 0;
	kernelRuntime_ = 0.0;
	profiling_ = 0;

	if(!ndevice)
		ndevice = use();
	ndevice_ = ndevice -1;
}

GPU::~GPU()
{
	if(!initialized_)
		return;

	for(vector<cl_mem>::iterator i = privateMemoryObjects_.begin();
		i != privateMemoryObjects_.end();
		i++)
		clReleaseMemObject(*i);

	privateMemoryObjects_.clear();
	this_device().usage--;
}

// Returns true if user wants to use the GPU
// Options -gpu <device number> or 0 to turn off
int GPU::use()
{
	using namespace basic::options;
	if(!option [ OptionKeys::gpu::gpu ]())
		return 0;
	return option [ OptionKeys::gpu::device ]();
}

int GPU::Init()
{
	// Check if already initialized and enabled
	if(initialized_)
		return 1;

	if(ndevice_ < 0)
		return 0;

	if((unsigned int)ndevice_ > sizeof(devices_)/sizeof(*devices_)) {
		TR.Error << "GPU device index out of bounds." << endl;
		ndevice_ = 0;
		return 0;
	}

	// Platform initialization
	if(!context_) {

		cl_uint numPlatforms;
		cl_platform_id firstPlatformId;

		// Init CL context
		errNum_ = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
		if (errNum_ != CL_SUCCESS || numPlatforms <= 0) {
			TR.Error << "Failed to find any OpenCL platforms: " << errstr(errNum_) << endl;
			return 0;
		}

		cl_context_properties contextProperties[] =	{
			CL_CONTEXT_PLATFORM,
			(cl_context_properties)firstPlatformId,
			0
		};

		context_ = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU, NULL, NULL, &errNum_);
		if (!context_ || errNum_ != CL_SUCCESS) {
			TR.Warning << "Could not create GPU context: " << errstr(errNum_) << ". Trying CPU..." << endl;
			context_ = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU, NULL, NULL, &errNum_);
			if (!context_ || errNum_ != CL_SUCCESS) {
				TR.Error << "Failed to create an OpenCL GPU or CPU context: " << errstr(errNum_) << endl;
				return 0;
			}
		}
	}

	// Device initialization: 1) locate GPU
	if(!this_device().initialized) {

		// Find the requested device
		cl_device_id *devices;
		size_t deviceBufferSize = -1;

		errNum_ = clGetContextInfo(context_, CL_CONTEXT_DEVICES, 0, NULL, &deviceBufferSize);
		if(errNum_ != CL_SUCCESS) {
			TR.Error << "Failed to enumerate CL devices: " << errstr(errNum_) << endl;
			return 0;
		}

		if(deviceBufferSize <= 0) {
			TR.Error << "No CL devices found." << endl;
			return 0;
		}

		int max_devices = deviceBufferSize / sizeof(cl_device_id);
		devices = new cl_device_id[max_devices];
		errNum_ = clGetContextInfo(context_, CL_CONTEXT_DEVICES, deviceBufferSize, devices, NULL);
		if (errNum_ != CL_SUCCESS) {
			TR.Error << "Failed to get device IDs: " << errstr(errNum_) << endl;
			return 0;
		}

		if(ndevice_ >= max_devices) {
			TR.Error << "Device #" << (ndevice_+1) << " unavailable. Using first device." << endl;
			ndevice_ = 0;
		}

		this_device().ndevice = ndevice_;
		this_device().device = devices[ndevice_];
		delete [] devices;
	}

	// Device initialization: 1) init GPU
	GPU_DEV &dev = this_device();

	if(!dev.initialized) {

		// Create device command queue
		dev.commandQueue = clCreateCommandQueue(context_, dev.device, profiling_ ? CL_QUEUE_PROFILING_ENABLE : 0, NULL);
		if(!dev.commandQueue) {
			TR.Error << "Failed to create commandQueue for device: " << errstr(errNum_) << endl;
			return 0;
		}

		// Get device info
		errNum_ = clGetDeviceInfo(dev.device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(dev.threads), &dev.threads, NULL);
		if(errNum_ != CL_SUCCESS) {
			TR.Error << "Failed to obtain max device threads: " << errstr(errNum_) << endl;
			return 0;
		}

		// These are non-critical, so failure is OK
		clGetDeviceInfo(dev.device, CL_DEVICE_NAME, sizeof(dev.name), dev.name, NULL);
		clGetDeviceInfo(dev.device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(dev.clockRate), &dev.clockRate, NULL);
		clGetDeviceInfo(dev.device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(dev.multiProcessorCount), &dev.multiProcessorCount, NULL);

		TR << "GPU #" << (ndevice_+1) << ": " << (dev.name) <<
			" [" << (dev.clockRate) << " MHz] with " <<
			(dev.multiProcessorCount) << " processors, " <<
			(dev.threads) << " threads." << endl;

		using namespace basic::options;
		unsigned int max_threads = option [ OptionKeys::gpu::threads ]();

		if(max_threads < dev.threads) {
			TR << "Using " << max_threads << " threads" << endl;
			dev.threads = max_threads;
		}

	}

	// All done
	this_device().initialized = 1;
	this_device().usage++;
	initialized_ = 1;
	return 1;
}

int GPU::Release()
{
	GPU_DEV &dev = this_device();

	if(initialized_) {
		dev.usage--;
		initialized_ = 0;
	}

	if(dev.usage > 0) {
		TR.Error << "GPU instances in use; cannot release GPU #" << ndevice_ << "." << endl;
		return 0;
	}

	/// BELOW objects are GPU-dependent (command queue)

	// Release command queue
	if(dev.commandQueue) {
		clReleaseCommandQueue(dev.commandQueue);
		dev.commandQueue = NULL;
	}

	dev.initialized = 0;
	return 1;
}

int GPU::Free()
{
	for(unsigned int i = 0; i < sizeof(devices_)/sizeof(*devices_); ++i) {
		GPU g(i);
		if(g.Release())	// Release all other instances if unused
			return 0;
	}

	/// BELOW objects are GPU-independent (kernel, program, context, memory)

	// Release all shared memory objects
	for(map<string,cl_mem>::iterator i = sharedMemoryObjects_.begin();
		i != sharedMemoryObjects_.end();
		i++)
		clReleaseMemObject(i->second);
	sharedMemoryObjects_.clear();

	// Release all kernels memory objects
	for(map<string,cl_kernel>::iterator i = kernels_.begin();
		i != kernels_.end();
		i++)
		clReleaseKernel(i->second);
	kernels_.clear();

	// Release all programs
	for(map<string,cl_program>::iterator i = programs_.begin();
		i != programs_.end();
		i++)
		clReleaseProgram(i->second);
	programs_.clear();

	// Release CL context
	if(context_) {
		clReleaseContext(context_);
		context_ = NULL;
	}

	return 1;
}

int GPU::RegisterProgram(const char *filename)
{
	std::string file(filename);
	return RegisterProgram(file);
}

int GPU::RegisterProgram(string& filename)
{
	std::vector<std::string> files;
	files.push_back(filename);
	return RegisterProgram(files);
}

int GPU::RegisterProgram(std::vector<std::string> & files)
{
	cl_program program;

	if(!Init()) return 0;
	if(!files.size()) return 0;

	std::string & filename(files[0]);
	if(programs_[filename])
		return 1;

	std::vector<std::string> sources;

	for(std::vector<std::string>::const_iterator file = files.begin();
		file < files.end();
		++file) {

		// Try location provided location first
		ifstream kernelFile(file->c_str(), ios::in);

		// Next try to find the file in the database
		if (!kernelFile.is_open()) {
			std::string database_fn = basic::database::full_name(file->c_str());
			kernelFile.open(database_fn.c_str(), ios::in);
		}

		if (!kernelFile.is_open()) {
			TR.Error << "Failed to open kernel file " << *file << endl;
			return 0;
		}

		ostringstream oss;
		oss << kernelFile.rdbuf();
		string srcStdStr = oss.str();
		sources.push_back(srcStdStr);
	}

	const char **psources = (const char**) &(*sources.begin());
	program = clCreateProgramWithSource(context_, sources.size(), psources, NULL, &errNum_);
	if(!program) {
		TR.Error << "Failed to create CL program from " << filename << ": " << errstr(errNum_) << endl;
		return 0;
	}

	errNum_ = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if (errNum_ != CL_SUCCESS) {
		TR.Error << "Failed to compile CL program from " << filename << ": " << errstr(errNum_) << endl;
		// Determine the reason for the error
		char buildLog[16384];
		std::cout << "buildlog: " << sizeof(buildLog) << std::endl;
		platform::Size actual_size(0);
		clGetProgramBuildInfo(program, this_device().device, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, &actual_size);
		TR.Error << "Error log size: " << actual_size << std::endl;
		TR.Error << "Error in kernel: " << endl;
		std::string build_log_string( buildLog, actual_size );
		TR.Error << build_log_string << std::endl;
		TR.Error << buildLog << std::endl;;
		clReleaseProgram(program);
		return 0;
	}

	programs_[filename] = program;
	TR.Trace << "OpenCL program " << filename << " loaded" << std::endl;
	return 1;
}

cl_kernel GPU::BuildKernel(const char *kernel_name)
{
	// Do we already have this kernel built?
	cl_kernel kernel = kernels_[kernel_name];
	if(kernel)
		return kernel;

	// Build it from available programs
	for(map<string,cl_program>::const_iterator i = programs_.begin();
		i != programs_.end();
		i++) {
		kernel = clCreateKernel(i->second, kernel_name, NULL);
		if(kernel)
			break;
	}
	if(!kernel)
		return NULL;

	kernels_[kernel_name] = kernel;
	return kernel;
}

cl_mem GPU::AllocateMemory(unsigned int size, void *data, int flags, const char *name)
{
	cl_mem r;

	if(!Init()) return NULL;

	if(!flags)
		flags = data ? (CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR) : CL_MEM_READ_WRITE;

	r = clCreateBuffer(context_, flags, size, data, &errNum_);
	if(!r || errNum_ != CL_SUCCESS) {
		TR.Error << "Failed to allocate GPU memory, " << size << " bytes";
		if(name)
			TR.Error << ", shared with name " << name;
		TR.Error << ": " << errstr(errNum_) << endl;
	}

	if(r) {
		TR.Trace << "Allocated " << size << " bytes = " << hex << r << dec;
		if(name)
			TR.Trace << ", with shared name " << name;
		if(data)
			TR.Trace << ", initialized with data " << hex << data << dec;
		TR.Trace << std::endl;

		if(name)
			sharedMemoryObjects_[name] = r;
		else
			privateMemoryObjects_.push_back(r);
	}

	return r;
}

cl_mem GPU::GetSharedMemory(const char *name)
{
	return sharedMemoryObjects_[name];
}

cl_mem GPU::AllocateMemoryReuse(cl_mem &old_mem, unsigned int &old_size, unsigned int new_size, int flags)
{
	if(old_mem && old_size < new_size) {
		Free(old_mem);
		old_mem = NULL;
	}

	cl_mem mem = old_mem;

	if(!mem)
		mem = AllocateMemory(new_size, NULL, flags);

	if(mem) {
		old_mem = mem;
		old_size = new_size;
	}
	return mem;
}

void GPU::Free(cl_mem h)
{
	{
	int found = 0;
	vector<cl_mem>::iterator i;
	for(i = privateMemoryObjects_.begin();
		i != privateMemoryObjects_.end();
		i++) {
		if(*i == h) {
			found = 1;
			break;
		}
	}
	if(found)
		privateMemoryObjects_.erase(i);
	}

	{
	int found = 0;
	map<string,cl_mem>::iterator i;
	for(i = sharedMemoryObjects_.begin();
		i != sharedMemoryObjects_.end();
		i++) {
		if(i->second == h) {
			found = 1;
			break;
		}
	}
	if(found)
		sharedMemoryObjects_.erase(i);
	}

	clReleaseMemObject(h);
	TR.Trace << "Released memory " << hex << h << dec << std::endl;
}

// Execute kernel on GPU
int GPU::_ExecuteKernel(const char *kernel_name, int total_threads, int max_conc_threads_high, int max_conc_threads_low, GPU_KERNEL_ARG *args, int async)
{
	cl_kernel kernel = BuildKernel(kernel_name);
	if(!kernel) {
		TR.Error << "Unknown kernel " << kernel_name << ": " << errstr(errNum_) << endl;
		return 0;
	}

	TR.Trace << "Execute kernel " << kernel_name << ": " << total_threads << " total, " << max_conc_threads_high << " concurrent threads" << std::endl;

	if(args) {
		for(int i =0; args[i].type; ++i) {
			errNum_	= clSetKernelArg(kernel, i, args[i].k_size, &args[i].k_p);
			if (errNum_ != CL_SUCCESS) {
				TR.Error << "Error setting kernel argument " << i << ": " << errstr(errNum_) << endl;
				return 0;
			}
		}
	}

	// Set supported thread limits
	int use_threads = this_device().threads;
	if(use_threads < max_conc_threads_high)
		max_conc_threads_high = use_threads;
	if(use_threads < max_conc_threads_low)
		max_conc_threads_low = use_threads;

	// Set requested upper limit
	if(max_conc_threads_high < use_threads)
		use_threads = max_conc_threads_high;

	cl_event kernelEvent;
	cl_ulong start, end;
	kernelRuntime_ = 0.0;

	do {
		size_t globalWorkSize = { UPPER_MULTIPLE(total_threads, use_threads) };
		size_t localWorkSize = { use_threads };

		errNum_ = clEnqueueNDRangeKernel(this_device().commandQueue, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, &kernelEvent );

		if(errNum_ == CL_OUT_OF_RESOURCES) {
			// Try fewer threads
			use_threads /= 2;
			if(use_threads < max_conc_threads_low)
				// Fail
				break;
		}

	} while(errNum_ != CL_SUCCESS);

	if(errNum_ != CL_SUCCESS) {
		TR.Error << "Failed to launch kernel " << kernel_name << ": " << errstr(errNum_) << endl;
		return 0;
	}

	if(async)
		return 1;

	clFinish(this_device().commandQueue);

	if( profiling_) {
		clFinish(this_device().commandQueue);
		if(	clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL) == CL_SUCCESS &&
			clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL) == CL_SUCCESS
		) {
			cl_ulong t = end - start;
			kernelRuntime_ = t * 1.0e-6;
		} else
			kernelRuntime_ = 0.0;
	}

	return 1;
}

// Execute kernel with custom kernel arguments
int GPU::ExecuteKernel(const char *kernel_name, int total_threads, int max_conc_threads_high, int max_conc_threads_low, ...)
{
	if(!Init()) return 0;

	va_list vl;
	int need_finish = 0;
	GPU_KERNEL_ARG args[32], *parg;

	// Allocate memory and copy data to GPU
	va_start(vl, max_conc_threads_low);
	parg = args;
	do {
		GPU_KERNEL_ARG &arg = *parg;
		arg.type = va_arg(vl, int);
		if(!arg.type)
			break;

		arg.mem = NULL;
		switch(arg.type & GPU_ARG_TYPE) {
			case GPU_INT:
				// Direct numeric integer argument
				arg.i = va_arg(vl, int);
				arg.k_size = sizeof(arg.i);
				arg.k_p = (void*)arg.i;
				break;
			case GPU_FLOAT:
				// warning: float is promoted to double when passed through ... (so you should pass double not float to va_arg)
			case GPU_DOUBLE:
				// Direct numeric float argument
				arg.f = va_arg(vl, double);
				arg.k_size = sizeof(arg.f);
				arg.k_p = (void*)arg.i;
				break;
			case GPU_DEVMEM:
				// Device memory pointer
				arg.k_size = sizeof(cl_mem);
				arg.k_p = va_arg(vl, cl_mem);
				break;
			case GPU_MEM:
			default:
				// Host memory pointer
				arg.size = va_arg(vl, int);
				arg.ptr = va_arg(vl, void*);
				arg.mem = AllocateMemory(arg.size, (arg.type & GPU_IN) ? arg.ptr : NULL);
				arg.k_size = sizeof(arg.mem);
				arg.k_p = arg.mem;
				need_finish = 1;
				break;
		}

		parg++;

	} while( parg < (args + sizeof(args)/sizeof(*args)) );

	va_end(vl);

	// Run kernel
	if(_ExecuteKernel(kernel_name, total_threads, max_conc_threads_high, max_conc_threads_low, args)) {
		for(GPU_KERNEL_ARG *arg = args; arg && arg->type; ++arg) {
			if(arg->type & GPU_OUT && arg->ptr && arg->mem && arg->size) {
				if(!ReadData(arg->ptr, arg->mem, arg->size, false)) {
					need_finish = 1;
					return 0;
				}
			}
		}
	}

	if(need_finish)
		clFinish(this_device().commandQueue);

	// Free memory
	for(GPU_KERNEL_ARG *arg = args; arg && arg->type; ++arg)
		if(arg->mem)
			Free(arg->mem);

	return 1;
}

int GPU::ReadData(void *dst, cl_mem src, unsigned int size, int blocking)
{
	errNum_ = clEnqueueReadBuffer(this_device().commandQueue, src, blocking, 0, size, dst, 0, NULL, NULL);
	if(errNum_ != CL_SUCCESS) {
		TR.Error << "Failed to read data from GPU (" << size << " bytes, " << hex << src << " -> " << dst << dec << "): " << errstr(errNum_) << endl;
		return 0;
	} else
		TR.Trace << "Read data from GPU: " << size << " bytes, " << hex << src << " -> " << dst << dec << endl;
	return 1;
}

int GPU::WriteData(cl_mem dst, void *src, unsigned int size, int blocking)
{
	errNum_ = clEnqueueWriteBuffer(this_device().commandQueue, dst, blocking, 0, size, src, 0, NULL, NULL);
	if(errNum_ != CL_SUCCESS) {
		TR.Error << "Failed to write data to GPU (" << size << " bytes, " << hex << src << " -> " << dst << dec << "): " << errstr(errNum_) << endl;
		return 0;
	} else
		TR.Trace << "Wrote data to GPU: " << size << " bytes, " << hex << src << " -> " << dst << dec << endl;
	return 1;
}

int GPU::setKernelArg(cl_kernel kernel, cl_uint i, size_t size, const void *p)
{
	errNum_	= clSetKernelArg(kernel, i, size, p);

	if(errNum_ != CL_SUCCESS) {
		TR.Error << "Failed to set kernel argument " << i << ": " << errstr(errNum_) << endl;
		return 0;
	}
	return 1;
}

} // gpu
} // basic

#endif // INCLUDED_basic_gpu_GPU_cc
#endif // USEOPENCL
