#include <apps/pilot/will/gpu_mat_vec.hh>
#include <apps/pilot/will/gpu_bb_struct.hh>

#define C1_SG           vec( 0.70523987, 1.66565453, 0.0)
#define C1_SG_M         vec(-0.45817111, 2.30884722, 0.0) 
#define D1_CG           vec( 0.59156419, 1.39717224, 0.0)
#define D2_OD1          vec( 0.60370457, 1.06283315, 0.0)
#define D2_OD1_M1       vec( 0.04702106, 2.01157099, 0.0)
#define D2_OD1_M2       vec( 1.70367738, 1.07056630, 0.0)
#define D2_OD2          vec( 0.57384868,-1.06272609, 0.0)
#define D21_rot001_COS  0.389893133f
#define D21_rot001_SIN -0.920860111f
#define D21_mov100     -1.517247000f
#define E1_CG           vec( 0.59156419, 1.39717224, 0.0)
#define E2_CD           vec( 0.56433087, 1.40081729, 0.0)
#define E3_OE1          vec( 0.60445958, 1.05973045, 0.0)
#define E3_OE1_M1       vec( 0.04947919, 2.00946555, 0.0)
#define E3_OE1_M2       vec( 1.70444450, 1.06549001, 0.0)
#define E3_OE2          vec( 0.57411109,-1.06347084, 0.0)
#define E21_rot001_COS  0.389893133f
#define E21_rot001_SIN -0.920860111f
#define E21_mov100     -1.517247000f
#define E32_rot001_COS  0.373675107f
#define E32_rot001_SIN -0.927559656f
#define E32_mov100     -1.510218000f
#define H1_CG           vec( 0.59156419, 1.39717224, 0.0)
#define H2_ND1          vec( 0.74366239, 1.16149143, 0.0)
#define H2_ND1_MD       vec(-0.01328216 ,3.12034011, 0.0)
#define H2_CD2          vec( 0.89177249,-1.01838273, 0.0)
#define H2_CE1          vec( 2.02989334, 0.85637562, 0.0)
#define H2_NE2          vec( 2.14584052,-0.45890166, 0.0)
#define H2_NE2_ME       vec( 3.95036650,-1.53298545, 0.0)
#define H21_rot001_COS  0.389893133f
#define H21_rot001_SIN  0.920860111f
#define H21_mov100      1.517247000f


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
  struct VEC Y  = (DorE) ? vec(-0.9338167*(-sin3),-0.3577517*(-sin3),cos3) : vec(0.511342*(-sin3),0.8593772*(-sin3),cos3);
  struct VEC ML = (DorE) ? H2_ND1_MD : H2_NE2_ME;
  struct VEC ND = H2_ND1;
  struct VEC NE = H2_NE2;
  struct VEC CE = H2_CE1;
  struct VEC CD = H2_CD2;
  struct VEC CG = H1_CG;

  rotx(& Y,sin2,cos2);  // rot chi2
  rotx(&ML,sin2,cos2);
  rotx(&ND,sin2,cos2);
  rotx(&NE,sin2,cos2);
  rotx(&CE,sin2,cos2);
  rotx(&CD,sin2,cos2);
  ML.x += H21_mov100;  // move up
  ND.x += H21_mov100;
  NE.x += H21_mov100;
  CE.x += H21_mov100;
  CD.x += H21_mov100;
  rotz(& Y,H21_rot001_SIN,H21_rot001_COS);  // tip so chi1 along x
  rotz(&ML,H21_rot001_SIN,H21_rot001_COS);
  rotz(&ND,H21_rot001_SIN,H21_rot001_COS);
  rotz(&NE,H21_rot001_SIN,H21_rot001_COS);
  rotz(&CE,H21_rot001_SIN,H21_rot001_COS);
  rotz(&CD,H21_rot001_SIN,H21_rot001_COS);
  rotx(& Y,sin1,cos1); // rot chi1
  rotx(&ML,sin1,cos1);
  rotx(&ND,sin1,cos1);
  rotx(&NE,sin1,cos1);
  rotx(&CE,sin1,cos1);
  rotx(&CD,sin1,cos1);
  rotx(&CG,sin1,cos1); // CG here only!!!

  // clash check?? here or other func?

  struct VEC const X = normalizedv(subvv(ML,(DorE)?ND:NE));
  struct VEC const Z = crossvv(X,Y);
  // std::cout << "hisd_bb2m ML: " << ML << std::endl;
  // std::cout << "hisd_bb2m ND: " << ND << std::endl;
  // std::cout << "hisd_bb2m NE: " << NE << std::endl;
  // std::cout << "hisd_bb2m CE: " << CE << std::endl;
  // std::cout << "hisd_bb2m CD: " << CD << std::endl;
  // std::cout << "hisd_bb2m CG: " << CG << std::endl;

  struct MAT const R = cols(X,Y,Z);
  return xform(R,ML);
}
inline struct XFORM const hisd_bb2m(float const chi1, float const chi2, float const chi3){ return       his_bb2m(chi1,chi2,chi3,true )  ; }
inline struct XFORM const hise_bb2m(float const chi1, float const chi2, float const chi3){ return       his_bb2m(chi1,chi2,chi3,false)  ; }
inline struct XFORM const hisd_m2bb(float const chi1, float const chi2, float const chi3){ return xrev( his_bb2m(chi1,chi2,chi3,true ) ); }
inline struct XFORM const hise_m2bb(float const chi1, float const chi2, float const chi3){ return xrev( his_bb2m(chi1,chi2,chi3,false) ); }

inline struct XFORM const
cys_bb2m(float const chi1, float const chi2){
  // return xform from bb in canonical to metal FOR
  float const cos1 = native_cos(chi1);
  float const cos2 = native_cos(chi2);
  float const cos3 = native_cos(chi3);
  float const sin1 = native_sin(chi1);
  float const sin2 = native_sin(chi2);
  float const sin3 = native_sin(chi3);
  TR << chi1 << " " << chi2 << endl;
  struct VEC Y  = vec(0.920937120054*(sin3),0.389711201926*sin3,cos3);
  struct VEC ML = C1_SG_M;
  struct VEC SG = C1_SG;

  // rotx(& Y,sin1,cos1);  // rot chi2
  // rotx(&ML,sin1,cos1);
  // rotx(&SG,sin1,cos1);

  // clash check?? here or other func?

  struct VEC const X = normalizedv(subvv(ML,SG));
  struct VEC const Z = crossvv(X,Y);

  struct MAT const R = cols(X,Y,Z);
  return xform(R,ML);
}
inline struct XFORM const cys_m2bb(float const chi1, float const chi2){ return xrev( cys_bb2m(chi1,chi2) ); }


