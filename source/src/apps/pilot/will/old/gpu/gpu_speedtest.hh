void test_gpu_speed() {
  CL cl(std::cout);
  cl.make_kernel("speed_test");
  cl_mem output = cl.makeWOmem(sizeof(struct MAT)*1);
  cl.setargs("speed_test",output);
  uint N = 10;
  MAT buf[N];
  for(int i = 0; i < N; ++i) {
    std::cerr << i << " ";
    cl.q1d("speed_test",12800*10,256);
    clFinish(cl.queue_);
    cl.gpu2cpu(output,buf+i,sizeof(struct MAT)*1);
    cl.zero_gpu_mem(output,sizeof(struct MAT)*1);
  }
  std::cerr << std::endl;
}
