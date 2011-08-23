#include <apps/pilot/will/gpu_mat_vec.hh>
#include <apps/pilot/will/gpu_bb_struct.hh>


#define CYS_SG  vec( 0.00000000f, 0.0f, 1.80880300f)
#define CYS_HG  vec(-1.32211517f, 0.0f, 1.94748754f)
#define ASP_CG  vec( 0.00000000f, 0.0f, 1.51724700f)
#define ASP_OD1 vec( 1.06283315f, 0.0f, 2.12095157f)
#define ASP_OD2 vec(-1.06272609f, 0.0f, 2.09109568f)
#define HIS_CG  vec( 0.00000000f, 0.0f, 1.51724700f)
#define HIS_ND1 vec( 1.16150173f, 0.0f, 2.26090072f)
#define HIS_CD2 vec(-1.01837843f, 0.0f, 2.40902742f)
#define HIS_CE1 vec( 0.85638485f, 0.0f, 3.54714067f)
#define HIS_NE2 vec(-0.45888741f, 0.0f, 3.66308773f)
#define HIS_MD  vec( 3.12034018f, 0.0f, 1.50396484f)
#define HIS_ME  vec(-1.53298549f, 0.0f, 5.46761352f)
#define HIS_MDN vec( 0.93278021f, 0.0f,-0.360445659)
#define HIS_MEN vec(-0.511475f  , 0.0f, 0.859298f  )

#define COS_CA_CB_CG 0.389893133175 
#define SIN_CA_CB_CG 0.920860111365

inline struct XFORM const
bb2bb(struct VEC const n1, struct VEC const ca1, struct VEC const c1,
      struct VEC const n2, struct VEC const ca2, struct VEC const c2)
{
  struct VEC const x1 = subvv(n1,ca1);
  struct VEC const y1 = subvv(c1,ca1);
  struct VEC const x2 = subvv(n2,ca2);
  struct VEC const y2 = subvv(c2,ca2);
  return vvcxform(x1,y1,x2,y2,ca1,ca2);
}

inline struct XFORM const
his_bb2m(float const chi1, float const chi2, float const chi3, bool DorE){
  // return xform from bb in canonical to metal FOR
  float const cos1 = native_cos(chi1);
  float const cos2 = native_cos(chi2);
  float const cos3 = native_cos(chi3);
  float const sin1 = native_sin(chi1);
  float const sin2 = native_sin(chi2);
  float const sin3 = native_sin(chi3);
  VEC X = (DorE) ? HIS_MDN : HIS_MEN;
  VEC Y = (DorE) ? vec(-0.3577517*(-sin3),cos3,-0.9338167*(-sin3)) : vec(0.8593772*(-sin3),cos3, 0.511342*(-sin3));
  VEC M = (DorE) ? HIS_MD  : HIS_ME;
  // rot chi2
  { float const x=X.x, y=X.y;
    X.x = cos2*x - sin2*y;
    X.y = sin2*x + cos2*y;  }
  { float const x=Y.x, y=Y.y;
    Y.x = cos2*x - sin2*y;
    Y.y = sin2*x + cos2*y;  }
  { float const x=M.x, y=M.y;
    M.x = cos2*x - sin2*y;
    M.y = sin2*x + cos2*y;  }
  // align chi1 Z
  { float const x=X.x, z=X.z;
    X.x =  COS_CA_CB_CG*x + SIN_CA_CB_CG*z;
    X.z = -SIN_CA_CB_CG*x + COS_CA_CB_CG*z;  }
  { float const x=Y.x, z=Y.z;
    Y.x =  COS_CA_CB_CG*x + SIN_CA_CB_CG*z;
    Y.z = -SIN_CA_CB_CG*x + COS_CA_CB_CG*z;  }
  { float const x=M.x, z=M.z;
    M.x =  COS_CA_CB_CG*x + SIN_CA_CB_CG*z;
    M.z = -SIN_CA_CB_CG*x + COS_CA_CB_CG*z;  }
  // rot chi1
  { float const x=X.x, y=X.y;
    X.x = cos1*x - sin1*y;
    X.y = sin1*x + cos1*y;  }
  { float const x=Y.x, y=Y.y;
    Y.x = cos1*x - sin1*y;
    Y.y = sin1*x + cos1*y;  }
  { float const x=M.x, y=M.y;
    M.x = cos1*x - sin1*y;
    M.y = sin1*x + cos1*y;  }
  // reset chi1 from Z
  { float const x=X.x, z=X.z;
    X.x =  COS_CA_CB_CG*x - SIN_CA_CB_CG*z;
    X.z =  SIN_CA_CB_CG*x + COS_CA_CB_CG*z;  }
  { float const x=Y.x, z=Y.z;
    Y.x =  COS_CA_CB_CG*x - SIN_CA_CB_CG*z;
    Y.z =  SIN_CA_CB_CG*x + COS_CA_CB_CG*z;  }
  { float const x=M.x, z=M.z;
    M.x =  COS_CA_CB_CG*x - SIN_CA_CB_CG*z;
    M.z =  SIN_CA_CB_CG*x + COS_CA_CB_CG*z;  }
  struct VEC const Z = crossvv(X,Y);
  std::cout << "hisd_m2bb: " << X << std::endl;
  std::cout << "hisd_m2bb: " << Y << std::endl;
  // std::cout << "hisd_m2bb: " << Z << std::endl;
  std::cout << "hisd_m2bb: " << M << std::endl;

  struct MAT const R = multmm( cols(X,Y,Z), rows(vec(0.6876640316337250,0.0,0.7260290487282525),vec(0.1124134291172747,0.7885885444652520,-0.6045587882847050),vec(-0.5725381907761041,0.4973487487177556,0.5422839777871462)) );
  return xform(R,multmv(R,M));
}

inline struct XFORM const
hisd_bb2m(float const chi1, float const chi2, float const chi3){
  return his_bb2m(chi1,chi2,chi3,true);
}
inline struct XFORM const
hise_bb2m(float const chi1, float const chi2, float const chi3){
  return his_bb2m(chi1,chi2,chi3,false);
}
inline struct XFORM const
hisd_m2bb(float const chi1, float const chi2, float const chi3){
  return xrev( his_bb2m(chi1,chi2,chi3,true) );
}
inline struct XFORM const
hise_m2bb(float const chi1, float const chi2, float const chi3){
  return xrev( his_bb2m(chi1,chi2,chi3,false) );
}

