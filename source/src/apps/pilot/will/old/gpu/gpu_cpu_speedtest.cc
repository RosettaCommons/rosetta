#include <iostream>


float test_gpu_inlining() {
  // CL cl(std::cout);
  // cl.make_kernel("speed_test");
  // cl_mem output = cl.makeWOmem(sizeof(MAT)*1);
  // cl.setargs("speed_test",output);
  // for(int i = 0; i < 10; ++i) {
  //   std::cerr << i << " ";
  //   cl.q1d("test_funcs",128*20*50,128); 
  //   clFinish(cl.queue_);
  // }
  // std::cerr << std::endl;
  float r;
  for(int j = 0; j <= 12800; ++j) {
    float xx = 1024.04f;
    float xy = 504.1f;
    float xz = 6436.04f;
    float yx = 16024.04f;
    float yy = 5034.1f;
    float yz = 646.04f;
    float zx = 10724.04f;
    float zy = 5.4f;
    float zz = 64736.04f;
  int i = 0;
  while( i < 40000u ) {
    {
      float const t1=xx+xy;
      float const t2=xz+yx;
      float const t3=yy+yz; 
      float const t4=zy+zz;
      float const t5=zy+xx;
      float const t6=t1+t2;
      float const t7=t3+t4;
      float const t8=t1+t3;
      float const t9=t2+t4;
      float const t0=t2+t3;
      xx = t1+t2;
      xy = t3+t4;
      xz = t5+t6;
      yx = t7+t8;
      yy = t9+t0;
    }
    {
      float const t1=xx+xy;
      float const t2=xz+yx;
      float const t3=yy+yz; 
      float const t4=zy+zz;
      float const t5=zy+xx;
      float const t6=t1+t2;
      float const t7=t3+t4;
      float const t8=t1+t3;
      float const t9=t2+t4;
      float const t0=t2+t3;
      xx = t1+t2;
      xy = t3+t4;
      xz = t5+t6;
      yx = t7+t8;
      yy = t9+t0;
    }
    {
      float const t1=xx+xy;
      float const t2=xz+yx;
      float const t3=yy+yz; 
      float const t4=zy+zz;
      float const t5=zy+xx;
      float const t6=t1+t2;
      float const t7=t3+t4;
      float const t8=t1+t3;
      float const t9=t2+t4;
      float const t0=t2+t3;
      xx = t1+t2;
      xy = t3+t4;
      xz = t5+t6;
      yx = t7+t8;
      yy = t9+t0;
    }
    {
      float const t1=xx+xy;
      float const t2=xz+yx;
      float const t3=yy+yz; 
      float const t4=zy+zz;
      float const t5=zy+xx;
      float const t6=t1+t2;
      float const t7=t3+t4;
      float const t8=t1+t3;
      float const t9=t2+t4;
      float const t0=t2+t3;
      xx = t1+t2;
      xy = t3+t4;
      xz = t5+t6;
      yx = t7+t8;
      yy = t9+t0;
    }
    {
      float const t1=xx+xy;
      float const t2=xz+yx;
      float const t3=yy+yz; 
      float const t4=zy+zz;
      float const t5=zy+xx;
      float const t6=t1+t2;
      float const t7=t3+t4;
      float const t8=t1+t3;
      float const t9=t2+t4;
      float const t0=t2+t3;
      xx = t1+t2;
      xy = t3+t4;
      xz = t5+t6;
      yx = t7+t8;
      yy = t9+t0;
    }
    //m = multmm(m,n);
    i++;
  }
   r = xx+yy+zz;
}
  return r;
}


int main (int argc, char *argv[]) {
  std::cout << "time this... HD5870 does 100X the work in 2.9 sec."  << std::endl;
  std::cout << test_gpu_inlining() << std::endl;
  std::cout << "DONE." << std::endl;

  //test_k_square();
}
