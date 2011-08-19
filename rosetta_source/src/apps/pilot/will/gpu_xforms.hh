#include <apps/pilot/will/gpu_mat_vec.hh>
#include <apps/pilot/will/gpu_bb_struct.hh>

inline struct MAT bb2bb(BB const bb, uint const rsd1, uint const rsd2) {
  float16 const r1 = bb.xyz[rsd1];
  float16 const r2 = bb.xyz[rsd2];
  struct VEC const x1 = subvv(get_BB__N(r1),get_BB_CA(r1));
  struct VEC const y1 = subvv(get_BB__C(r1),get_BB_CA(r1));
  struct VEC const x2 = subvv(get_BB__N(r2),get_BB_CA(r2));
  struct VEC const y2 = subvv(get_BB__C(r2),get_BB_CA(r2));
  return tform(x1,y1,x2,y2);
}

inline struct MAT his_b2m(float const chi1, float const chi2, float const chi3) {
  // return xform from bb in canonical to metal FOR
}

inline struct MAT his_m2b(float const chi1, float const chi2, float const chi3) {
  // return xform from metal in canonical to bb FOR
}
