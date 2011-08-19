#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

struct OctreeIdx {
  uint b1,b2,b3,b4,b5,b6,b7,b8,b9;
  uint e1,e2,e3,e4,e5,e6,e7,e8,e9;
};


// packing pattern:
// 4 4 4 10 10
// 4 4 4 10 10
// 4 4 4 10 10
// 10 10 10

// cl_uint4 pack14(OctreeIdx const & idx) {
//   cl_uint const i1 = (idx.e1-idx.b1) + 16*(idx.e2-idx.b2) + 256*(idx.e3-idx.b3) + 4096*idx.b1 + 4194304*idx.b2;
//   cl_uint const i2 = (idx.e4-idx.b4) + 16*(idx.e5-idx.b5) + 256*(idx.e6-idx.b6) + 4096*idx.b4 + 4194304*idx.b5;
//   cl_uint const i3 = (idx.e7-idx.b7) + 16*(idx.e8-idx.b8) + 256*(idx.e9-idx.b9) + 4096*idx.b7 + 4194304*idx.b8;
//   cl_uint const i4 = idx.b3 + 1024*idx.b6 + 1048576*idx.b9;
//   cl_uint4 r;
//   r.x = i1;
//   r.y = i2;
//   r.z = i3;
//   r.w = i4;
//   return r;
// }

// OctreeIdx pack14(cl_uint4 const & idx) {
//   OctreeIdx r;
//   r.b1 = (idx.x/4096u) % 1024;
//   r.b2 = (idx.x/4194304u);
//   r.b3 = (idx.w%1024u);
//   r.b4 = (idx.y/4096u) % 1024;
//   r.b5 = (idx.y/4194304u);
//   r.b6 = (idx.w/1024u) % 1024;
//   r.b7 = (idx.z/4096u) % 1024;
//   r.b8 = (idx.z/4194304u);
//   r.b9 = (idx.w/1048576u);
//   r.e1 = r.b1 + (idx.x    ) % 16u;
//   r.e2 = r.b2 + (idx.x/ 16) % 16u;
//   r.e3 = r.b3 + (idx.x/256) % 16u;
//   r.e4 = r.b4 + (idx.y    ) % 16u;
//   r.e5 = r.b5 + (idx.y/ 16) % 16u;
//   r.e6 = r.b6 + (idx.y/256) % 16u;
//   r.e7 = r.b7 + (idx.z    ) % 16u;
//   r.e8 = r.b8 + (idx.z/ 16) % 16u;
//   r.e9 = r.b9 + (idx.z/256) % 16u;
//   return r;
// }

cl_uint4 pack14(OctreeIdx const & idx) {
  cl_uint const i1 = (idx.e1-idx.b1) + ((idx.e2-idx.b2)<<4) + ((idx.e3-idx.b3)<<8) + (idx.b1<<12) + (idx.b2<<22);
  cl_uint const i2 = (idx.e4-idx.b4) + ((idx.e5-idx.b5)<<4) + ((idx.e6-idx.b6)<<8) + (idx.b4<<12) + (idx.b5<<22);
  cl_uint const i3 = (idx.e7-idx.b7) + ((idx.e8-idx.b8)<<4) + ((idx.e9-idx.b9)<<8) + (idx.b7<<12) + (idx.b8<<22);
  cl_uint const i4 = idx.b3 + (idx.b6<<10) + (idx.b9<<20);
  cl_uint4 r;
  r.x = i1;
  r.y = i2;
  r.z = i3;
  r.w = i4;
  return r;
}

OctreeIdx pack14(cl_uint4 const & idx) {
  OctreeIdx r;
  r.b1 = (idx.x>>12) & 1023u;
  r.b2 = (idx.x>>22);
  r.b3 = (idx.w    ) & 1023u;
  r.b4 = (idx.y>>12) & 1023u;
  r.b5 = (idx.y>>22);
  r.b6 = (idx.w>>10) & 1023u;
  r.b7 = (idx.z>>12) & 1023u;
  r.b8 = (idx.z>>22);
  r.b9 = (idx.w>>20);
  r.e1 = r.b1 + ((idx.x   ) & 15u);
  r.e2 = r.b2 + ((idx.x>>4) & 15u);
  r.e3 = r.b3 + ((idx.x>>8) & 15u);
  r.e4 = r.b4 + ((idx.y   ) & 15u);
  r.e5 = r.b5 + ((idx.y>>4) & 15u);
  r.e6 = r.b6 + ((idx.y>>8) & 15u);
  r.e7 = r.b7 + ((idx.z   ) & 15u);
  r.e8 = r.b8 + ((idx.z>>4) & 15u);
  r.e9 = r.b9 + ((idx.z>>8) & 15u);
  return r;
}

void test_pack14(int N) {
    OctreeIdx i;
    i.b1 = rand()%1024; i.e1 = i.b1 + min(15,rand()%16);
    i.b2 = rand()%1024; i.e2 = i.b2 + min(15,rand()%16);
    i.b3 = rand()%1024; i.e3 = i.b3 + min(15,rand()%16);
    i.b4 = rand()%1024; i.e4 = i.b4 + min(15,rand()%16);
    i.b5 = rand()%1024; i.e5 = i.b5 + min(15,rand()%16);
    i.b6 = rand()%1024; i.e6 = i.b6 + min(15,rand()%16);
    i.b7 = rand()%1024; i.e7 = i.b7 + min(15,rand()%16);
    i.b8 = rand()%1024; i.e8 = i.b8 + min(15,rand()%16);
    i.b9 = rand()%1024; i.e9 = i.b9 + min(15,rand()%16);

    for(int ii = 0; ii < N; ++ii) {
      pack14(pack14(i));

    // if( i.b1 != j.b1 || i.e1 != j.e1 ||
    //     i.b2 != j.b2 || i.e2 != j.e2 ||
    //     i.b3 != j.b3 || i.e3 != j.e3 ||
    //     i.b4 != j.b4 || i.e4 != j.e4 ||
    //     i.b5 != j.b5 || i.e5 != j.e5 ||
    //     i.b6 != j.b6 || i.e6 != j.e6 ||
    //     i.b7 != j.b7 || i.e7 != j.e7 ||
    //     i.b8 != j.b8 || i.e8 != j.e8 ||
    //     i.b9 != j.b9 || i.e9 != j.e9
    //     ) {
    //   std::cout << i.b1 << " " << j.b1 << "   " << i.e1 << " " << j.e1 << std::endl;
    //   std::cout << i.b2 << " " << j.b2 << "   " << i.e2 << " " << j.e2 << std::endl;
    //   std::cout << i.b3 << " " << j.b3 << "   " << i.e3 << " " << j.e3 << std::endl;
    //   std::cout << i.b4 << " " << j.b4 << "   " << i.e4 << " " << j.e4 << std::endl;
    //   std::cout << i.b5 << " " << j.b5 << "   " << i.e5 << " " << j.e5 << std::endl;
    //   std::cout << i.b6 << " " << j.b6 << "   " << i.e6 << " " << j.e6 << std::endl;
    //   std::cout << i.b7 << " " << j.b7 << "   " << i.e7 << " " << j.e7 << std::endl;
    //   std::cout << i.b8 << " " << j.b8 << "   " << i.e8 << " " << j.e8 << std::endl;
    //   std::cout << i.b9 << " " << j.b9 << "   " << i.e9 << " " << j.e9 << std::endl;
    // }
  }
}
