float test_gpu_speed() {
  CL cl(std::cout);
  cl.make_kernel("speed_test");
  cl_mem output = cl.makeWOmem(sizeof(struct MAT)*1);
  cl.setargs("speed_test",output);
  for(int i = 0; i < 10; ++i) {
    std::cerr << i << " ";
    cl.q1d("test_funcs",128*20*50,128);
    clFinish(cl.queue_);
  }
  std::cerr << std::endl;
  return r;
}
