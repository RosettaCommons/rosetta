void test_k_square() {
  std::size_t count(1024);
  float results[count],data[count];
  for(int i = 0; i < count; i++) data[i] = rand() / (float)RAND_MAX;

  CL cl((std::ofstream &)TR);
  cl.make_kernel("k_square");
  cl.dont_exit_on_error();

  cl_mem input  = cl.makeROmem(sizeof(float)*count); // TODO: sanity check on GPU
  cl_mem output = cl.makeWOmem(sizeof(float)*count); // CPU doesn't respect RO vs. WO!
  cl.cpu2gpu(data,input,sizeof(float)*count);
  cl.setargs("k_square",input,output);

  bool fail11=0,fail21=0,fail31=0,fail12=0,fail22=0,fail32=0;
  int correct = 0;
  cl.q1d("k_square",count,       64    ); clFinish(cl.queue_); cl.gpu2cpu(output,results,sizeof(float)*count); correct=0; for(int i = 0; i < count; i++) if(results[i] != data[i] * data[i]) fail11=true; cl.zero_gpu_mem(output,count);
  cl.q2d("k_square",count/16,16, 64,1  ); clFinish(cl.queue_); cl.gpu2cpu(output,results,sizeof(float)*count); correct=0; for(int i = 0; i < count; i++) if(results[i] != data[i] * data[i]) fail21=true; cl.zero_gpu_mem(output,count);
  cl.q3d("k_square",count/16,4,4,64,1,1); clFinish(cl.queue_); cl.gpu2cpu(output,results,sizeof(float)*count); correct=0; for(int i = 0; i < count; i++) if(results[i] != data[i] * data[i]) fail31=true; cl.zero_gpu_mem(output,count);
  cl.q1d("k_square",count,        4    ); clFinish(cl.queue_); cl.gpu2cpu(output,results,sizeof(float)*count); correct=0; for(int i = 0; i < count; i++) if(results[i] != data[i] * data[i]) fail12=true; cl.zero_gpu_mem(output,count);
  cl.q2d("k_square",count/16,16,  4,4  ); clFinish(cl.queue_); cl.gpu2cpu(output,results,sizeof(float)*count); correct=0; for(int i = 0; i < count; i++) if(results[i] != data[i] * data[i]) fail22=true; cl.zero_gpu_mem(output,count);
  cl.q3d("k_square",count/16,4,4, 4,4,4); clFinish(cl.queue_); cl.gpu2cpu(output,results,sizeof(float)*count); correct=0; for(int i = 0; i < count; i++) if(results[i] != data[i] * data[i]) fail32=true; cl.zero_gpu_mem(output,count);

  TR << "test11 " << (fail11?"fail":"pass") << std::endl;
  TR << "test21 " << (fail21?"fail":"pass") << std::endl;
  TR << "test31 " << (fail31?"fail":"pass") << std::endl;
  TR << "test12 " << (fail12?"fail":"pass") << std::endl;
  TR << "test22 " << (fail22?"fail":"pass") << std::endl;
  TR << "test32 " << (fail32?"fail":"pass") << std::endl;
  TR << std::endl << "fails last 2 on my laptop CPU, workunit dim seems to be only linear... why? will it work on gpu?" << std::endl;

  clReleaseMemObject(input); // cleanup
  clReleaseMemObject(output);
}
