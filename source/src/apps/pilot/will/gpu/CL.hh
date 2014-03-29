#ifdef MAC
#include <OpenCL/cl_platform.h>
#include <OpenCL/opencl.h>
#include <OpenCL/cl.h>
#else
#include <CL/cl_platform.h>
#include <CL/opencl.h>
#include <CL/cl.h>
#endif

#include <map>
#include <string>

#include <time.h>
#include <sys/time.h>
#include <stdio.h>

#ifdef __MACH__
#include <mach/mach_time.h>
#endif


std::string zeromemfunc("__kernel void CL_HH__AUTO_GEN__zero_mem( __global float* output ) { output[get_global_id(0)] = 0.0f; } ");

double time_highres() {
#ifdef __MACH__
  mach_timebase_info_data_t info;
  mach_timebase_info(&info);
  return mach_absolute_time() / 1000000000.0;
  //uint64_t duration = mach_absolute_time();
  //duration *= info.numer;
  //duration /= info.denom;
#else
  timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  return tp.tv_sec + tp.tv_nsec/1000000000.0;
#endif
  return 0;
}


inline float const native_recip(float const x) { return 1.0f/x; }

struct CL {
  int err_;                            // error code returned from api calls
  Size wgsize_;                       // wgsize domain size for our calculation
  cl_platform_id platform_;
  cl_device_id device_id_;             // compute device id
  cl_context context_;                 // compute context
  cl_command_queue queue_;          // compute command queue
  cl_program program_;                 // compute program
  cl_mem zmem_size_;
  bool exitonerr_;
  std::ostream & out_;
  std::map<std::string,cl_kernel> kernels_;
  CL(std::ostream & out) : exitonerr_(true),out_(out) {
    cl_uint numPlatforms;
    err_ = clGetPlatformIDs(0, NULL, &numPlatforms);
    cl_platform_id* platforms = new cl_platform_id[numPlatforms];
    err_ = clGetPlatformIDs(numPlatforms, platforms, NULL);
    platform_ = platforms[0];

    out << "clGetDeviceIDs" << endl;
    err_ = clGetDeviceIDs(platform_, CL_DEVICE_TYPE_GPU, 2, &device_id_, NULL);
    if(err_ != CL_SUCCESS) {
      out_ << "No GPU 2, falling back to GPU 1!!!!!!!!!!" << std::endl;
      err_ = clGetDeviceIDs(platform_, CL_DEVICE_TYPE_GPU, 1, &device_id_, NULL);
    }
    if(err_ != CL_SUCCESS) {
      out_ << "No GPU, falling back to CPU!!!!!!!!!!" << std::endl;
      err_ = clGetDeviceIDs(platform_, CL_DEVICE_TYPE_CPU, 1, &device_id_, NULL);
      if(err_ != CL_SUCCESS) handle_error("Error: Failed to create a device group! ERR="+errstr(err_));
    }
    out << "clCreateContext(0, 1, &device_id_, NULL, NULL, &err_);" << endl;
    context_ = clCreateContext(0, 1, &device_id_, NULL, NULL, &err_);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to create a compute context: "+errstr(err_));

    queue_ = clCreateCommandQueue(context_, device_id_, NULL, &err_);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to create a command queue: "+errstr(err_));

		out << "reading kernels from: '" << basic::options::option[basic::options::OptionKeys::gpu::kernel]() << "'" << endl;
    utility::io::izstream kin(basic::options::option[basic::options::OptionKeys::gpu::kernel]());
    std::string stmp((std::istreambuf_iterator<char>(kin)), std::istreambuf_iterator<char>());
    stmp += zeromemfunc;
    const char * kernelsource = stmp.c_str();
    //out_ << kernelsource << std::endl;
    program_ = clCreateProgramWithSource(context_, 1, (const char **) &kernelsource, NULL, &err_);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to create compute program! "+errstr(err_));
    out << "clBuildProgram" << endl;
    double t = time_highres();
    err_ = clBuildProgram(program_, 0, NULL, "-I src/apps/pilot/will -cl-single-precision-constant -cl-mad-enable -cl-no-signed-zeros -cl-fast-relaxed-math -w", NULL, NULL);
    //err_ = clBuildProgram(program_, 0, NULL, "-Isrc/apps/pilot/will -cl-single-precision-constant -w", NULL, NULL);
    out << "build took " << time_highres()-t << endl;
    if(err_ != CL_SUCCESS) {
      Size len_status, len_options, len_log;
      char buffer_options[81920];
      char buffer_log[81920];
      cout << "Error: Failed to build program executable!" << endl;

      cl_build_status status;
      clGetProgramBuildInfo(program_, device_id_, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status , &len_status);
      cout << "CL_PROGRAM_BUILD_STATUS: '" << status << "'" << endl;

      clGetProgramBuildInfo(program_, device_id_, CL_PROGRAM_BUILD_OPTIONS, sizeof(buffer_options), buffer_options, &len_options);
      cout << "CL_PROGRAM_BUILD_OPTIONS: '";
      for(uint i=0;i<len_options;++i) cout << buffer_options[i];
      cout << "'" << endl;

      clGetProgramBuildInfo(program_, device_id_, CL_PROGRAM_BUILD_LOG, sizeof(buffer_log    ), buffer_log    , &len_log);
      cout << "---- start build log ----" << endl;
      for(uint i=0;i<len_log;++i) cout << buffer_log[i];
      cout << endl << "----  end build log  ----" << endl;

      handle_error( errstr(err_) );
    }
    make_kernel("CL_HH__AUTO_GEN__zero_mem");
    zmem_size_ = makeROmem(sizeof(uint));
    err_ = clGetKernelWorkGroupInfo(kernels_["CL_HH__AUTO_GEN__zero_mem"], device_id_, CL_KERNEL_WORK_GROUP_SIZE, sizeof(wgsize_), &wgsize_, NULL);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to retrieve kernel work group info for kernel CL_HH__AUTO_GEN__zero_mem! ERR="+errstr(err_));
  }
  ~CL(){
    for(std::map<std::string,cl_kernel>::iterator i = kernels_.begin(); i != kernels_.end(); ++i) clReleaseKernel(i->second);
    clReleaseProgram(program_);
    clReleaseCommandQueue(queue_);
    clReleaseContext(context_);
  }
  void make_kernel(std::string kname) {
    kernels_[kname] = clCreateKernel(program_,kname.c_str(),&err_);
    if(!(kernels_[kname]) || err_ != CL_SUCCESS) handle_error("Error: Failed to create compute kernel "+kname+"! "+errstr(err_));
  }
  cl_kernel & operator[](std::string kname) {
    return kernels_[kname];
  }
  void handle_error(std::string msg) {
    if(exitonerr_) utility_exit_with_message(msg);
    else           out_ << "ERROR:" << msg << std::endl;
  }
  void exit_on_error() { exitonerr_ = true; }
  void dont_exit_on_error() { exitonerr_ = false; }
  void zero_gpu_mem(cl_mem buf, uint size) {
    err_ = 0;
    err_  = clSetKernelArg(kernels_["CL_HH__AUTO_GEN__zero_mem"], 0, sizeof(cl_mem), &buf);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to set kernel arguments to CL_HH__AUTO_GEN__zero_mem! "+errstr(err_));
    q1d("CL_HH__AUTO_GEN__zero_mem",size,min(64u,size));
    clFinish(queue_);
  }
  cl_mem makeROmem(size_t size) {
    cl_mem tmp = clCreateBuffer(context_, CL_MEM_READ_ONLY , size, NULL, NULL);
    if(!tmp) handle_error("Error: Failed to allocate device memory (RO)!");
    return(tmp);
  }
  cl_mem makeWOmem(size_t size) {
    cl_mem tmp = clCreateBuffer(context_, CL_MEM_WRITE_ONLY , size, NULL, NULL);
    if(!tmp) handle_error("Error: Failed to allocate device memory (WO)!");
    return(tmp);
  }
  cl_mem makeRWmem(size_t size) {
    cl_mem tmp = clCreateBuffer(context_, CL_MEM_READ_WRITE , size, NULL, NULL);
    if(!tmp) handle_error("Error: Failed to allocate device memory (RW)!");
    return(tmp);
  }
  void cpu2gpu(void const * const dat, cl_mem buf, std::size_t size) {
    err_ = clEnqueueWriteBuffer(queue_, buf, CL_TRUE, 0, size, dat, 0, NULL, NULL);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to write to source array! "+errstr(err_));
  }
  void gpu2cpu(cl_mem buf, void* dat, std::size_t size) {
    err_ = clEnqueueReadBuffer( queue_, buf, CL_TRUE, 0, size, dat, 0, NULL, NULL );
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to read output array! "+errstr(err_));
  }
  bool havekernel( std::string kname ) {
    return kernels_.find(kname) != kernels_.end();
  }
  void q1d( std::string kname, std::size_t gdim0, std::size_t ldim0 = 0 ) {
    if(!havekernel(kname)) handle_error("no kernel "+kname);
    if(ldim0==0) ldim0 = wgsize_;
    err_ = clEnqueueNDRangeKernel(queue_, kernels_[kname], 1, NULL, &gdim0, &ldim0, 0, NULL, NULL);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to q1d kernel! ERR="+errstr(err_)+"; max_wg_size: "+str(wgsize_)+" local_dim: ("+str(ldim0)+") global_dim: ("+str(gdim0)+")");
  }
  void q2d( std::string kname, std::size_t gdim0, sstd::ize_t gdim1, std::size_t ldim0 = 0, std::size_t ldim1 = 1 ) {
    if(!havekernel(kname)) handle_error("no kernel "+kname);
    if(ldim0==0) {
      ldim0 = min(gdim0,wgsize_);
      ldim1 = min(gdim1,wgsize_/ldim0);
    }
    std::size_t const gdim[2] = { gdim0,gdim1 };
    std::size_t const ldim[2] = { ldim0,ldim1 };
    err_ = clEnqueueNDRangeKernel(queue_, kernels_[kname], 2, NULL, gdim, ldim, 0, NULL, NULL);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to q2d kernel! ERR="+errstr(err_)+"; max_wg_size: "+str(wgsize_)+" local_dim: ("+str(ldim0)+","+str(ldim1)+") global_dim: ("+str(gdim0)+","+str(gdim1)+")");
  }
  void q3d( std::string kname, std::size_t gdim0, std::size_t gdim1, std::size_t gdim2, std::size_t ldim0 = 0, std::size_t ldim1 = 1, std::size_t ldim2 = 1 ) {
    if(!havekernel(kname)) handle_error("no kernel "+kname);
    if(ldim0==0) {
      ldim0 = min(gdim0,wgsize_);
      ldim1 = min(gdim1,wgsize_/ldim0);
      ldim2 = min(gdim1,wgsize_/ldim0/ldim1);
    }
    std::size_t gdim[3] = { gdim0,gdim1,gdim2 };
    std::size_t ldim[3] = { ldim0,ldim1,ldim2 };
    err_ = clEnqueueNDRangeKernel(queue_, kernels_[kname], 3, NULL, gdim, ldim, 0, NULL, NULL);
    if(err_ != CL_SUCCESS) handle_error("Error: Failed to q3d kernel! ERR="+errstr(err_)+"; max_wg_size: "+str(wgsize_)+" local_dim: ("+str(ldim0)+","+str(ldim1)+","+str(ldim2)+") global_dim: ("+str(gdim0)+","+str(gdim1)+","+str(gdim2)+")");
  }
  void finish() { clFinish(queue_); }
  std::string errstr(int input)
  {
    int errorCode = (int)input;
    switch(errorCode)
      {
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
        //case CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR:
        //        return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
        //case CL_PLATFORM_NOT_FOUND_KHR:
        //return "CL_PLATFORM_NOT_FOUND_KHR";
        //case CL_INVALID_PROPERTY_EXT:
        //    return "CL_INVALID_PROPERTY_EXT";
        //case CL_DEVICE_PARTITION_FAILED_EXT:
        //return "CL_DEVICE_PARTITION_FAILED_EXT";
        //case CL_INVALID_PARTITION_COUNT_EXT:
        //return "CL_INVALID_PARTITION_COUNT_EXT";
      default:
        return "unknown error code";
      }

    return "unknown error code";
  }
  ///////////////// arg setting convenience funcs ////////////////////
  void setargs(std::string kname, cl_mem a1, cl_mem a2, cl_mem a3, cl_mem a4, cl_mem a5) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 5 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 5 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 5 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 5 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],4,sizeof(cl_mem),&a5); if(err_!=CL_SUCCESS) handle_error(s+"5 of 5 type cl_mem for "+kname+"! ERR="+str(err_));
  }


#include <apps/pilot/will/gpu/CL_args.hh>

  // void setargs(std::string kname, cl_mem m1) { std::string const s("Error: Failed to set kernel argument "); err_ = 0;
  //   err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&m1);if(err_!=CL_SUCCESS)handle_error(s+"1 of 1 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  // }
  // void setargs(std::string kname, cl_mem m1, cl_mem m2) { std::string const s("Error: Failed to set kernel argument "); err_ = 0;
  //   err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&m1);if(err_!=CL_SUCCESS)handle_error(s+"1 of 2 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  //   err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&m2);if(err_!=CL_SUCCESS)handle_error(s+"2 of 2 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  // }
  // void setargs(std::string kname, cl_mem m1, cl_mem m2, cl_mem m3) { std::string const s("Error: Failed to set kernel argument "); err_ = 0;
  //   err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&m1);if(err_!=CL_SUCCESS)handle_error(s+"1 of 3 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  //   err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&m2);if(err_!=CL_SUCCESS)handle_error(s+"2 of 3 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  //   err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&m3);if(err_!=CL_SUCCESS)handle_error(s+"3 of 3 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  // }
  // void setargs(std::string kname, cl_mem m1, cl_mem m2, cl_mem m3, cl_mem m4) { std::string const s("Error: Failed to set kernel argument "); err_ = 0;
  //   err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&m1);if(err_!=CL_SUCCESS)handle_error(s+"1 of 4 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  //   err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&m2);if(err_!=CL_SUCCESS)handle_error(s+"2 of 4 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  //   err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&m3);if(err_!=CL_SUCCESS)handle_error(s+"3 of 4 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  //   err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&m4);if(err_!=CL_SUCCESS)handle_error(s+"4 of 4 (cl_mem) for "+kname+"! ERR="+errstr(err_));
  // }



};


uint gid_,lid_,gsz_,lsz_;
uint get_global_id  (int) { return gid_; }
uint get_global_size(int) { return gsz_; }
uint get_local_id   (int) { return lid_; }
uint get_local_size (int) { return lsz_; }

