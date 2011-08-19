  void setargs(std::string kname, cl_mem a1, cl_mem a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1, cl_mem a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,    int a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,   uint a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1,  float a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname, cl_mem a1) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(cl_mem),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 1 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1, cl_mem a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,    int a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,   uint a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1,  float a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,    int a1) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(   int),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 1 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1, cl_mem a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,    int a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,   uint a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1,  float a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,   uint a1) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof(  uint),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 1 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1, cl_mem a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(cl_mem),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,    int a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(   int),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,   uint a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof(  uint),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2, cl_mem a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2, cl_mem a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2, cl_mem a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2, cl_mem a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2, cl_mem a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(cl_mem),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,    int a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,    int a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,    int a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,    int a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type    int for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,    int a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(   int),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,   uint a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,   uint a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,   uint a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,   uint a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type   uint for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,   uint a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof(  uint),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,  float a3, cl_mem a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(cl_mem),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type cl_mem for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,  float a3,    int a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(   int),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type    int for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,  float a3,   uint a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof(  uint),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type   uint for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,  float a3,  float a4) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 4 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],3,sizeof( float),&a4); if(err_!=CL_SUCCESS) handle_error(s+"4 of 4 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2,  float a3) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 3 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],2,sizeof( float),&a3); if(err_!=CL_SUCCESS) handle_error(s+"3 of 3 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1,  float a2) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 2 type  float for "+kname+"! ERR="+str(err_));
    err_=clSetKernelArg(kernels_[kname],1,sizeof( float),&a2); if(err_!=CL_SUCCESS) handle_error(s+"2 of 2 type  float for "+kname+"! ERR="+str(err_));
  }
  void setargs(std::string kname,  float a1) {
    std::string const s("Error: Failed to set kernel argument "); err_ = 0;
    err_=clSetKernelArg(kernels_[kname],0,sizeof( float),&a1); if(err_!=CL_SUCCESS) handle_error(s+"1 of 1 type  float for "+kname+"! ERR="+str(err_));
  }
