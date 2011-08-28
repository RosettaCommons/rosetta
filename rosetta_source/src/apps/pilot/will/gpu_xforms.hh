#include <apps/pilot/will/gpu_mat_vec.hh>
#include <apps/pilot/will/gpu_bb_struct.hh>

// C[CYS_M]
#define C1_SG           vec( 0.70523987, 1.66565453, 0.00000000)  // 17
#define C2_M            vec( 0.10432359, 0.99454341, 0.00000000)  // 53
#define C3_MY           vec( 0.00000000, 1.00000000, 0.00000000)  // 53
#define C21_rot001_COS  0.389893133f
#define C21_rot001_SIN -0.920860111f
#define C21_movx -1.80880300
#define C32_rot001_COS  0.104323588f
#define C32_rot001_SIN -0.994543407f
#define C32_movx -1.00000000
// D[ASP_M1]
#define D1_CG           vec( 0.59156419, 1.39717224, 0.00000000)  // 2
#define D2_OD1          vec( 0.61116150, 1.05856277, 0.00000000)  // 15
#define D2_OD2          vec( 0.60388100,-1.04595257, 0.00000000)  // 15
#define D21_rot001_COS  0.389893133f
#define D21_rot001_SIN -0.920860111f
#define D21_movx -1.51724700
#define D2_M1           vec(-0.43883850, 2.87721612, 0.00000000)  // 53
#define D3_M1Y          vec( 0.00000000, 1.00000000, 0.00000000)  // 53
#define D32M1_rot001_COS -0.500000000f
#define D32M1_rot001_SIN -0.866025404f
#define D32M1_movx -2.71116150
#define D32M1_movy  1.05856277
#define D2_M2           vec( 2.71116150, 1.05856277, 0.00000000)  // 53
#define D3_M2Y          vec( 0.00000000, 1.00000000, 0.00000000)  // 53
#define D32M2_rot001_COS  1.000000000f
#define D32M2_rot001_SIN -0.000000000f
#define D32M2_movx -2.71116150
#define D32M2_movy  1.05856277
// E[GLU_M1]
#define E1_CG           vec( 0.59156419, 1.39717224, 0.00000000)  // 4
#define E2_CD           vec( 0.56433087, 1.40081729, 0.00000000)  // 2
#define E3_OE1          vec( 0.61000000, 1.05655099, 0.00000000)  // 15
#define E3_OE2          vec( 0.60427100,-1.04662807, 0.00000000)  // 15
#define E21_rot001_COS  0.389893133f
#define E21_rot001_SIN -0.920860111f
#define E21_movx -1.51724700
#define E32_rot001_COS  0.373675107f
#define E32_rot001_SIN -0.927559656f
#define E32_movx -1.51021800
#define E3_M1           vec(-0.44000000, 2.87520434, 0.00000000)  // 53
#define E4_M1Y          vec( 0.00000000, 1.00000000, 0.00000000)  // 53
#define E43M1_rot001_COS -0.500000000f
#define E43M1_rot001_SIN -0.866025404f
#define E43M1_movx -2.71000000
#define E43M1_movy  1.05655099
#define E3_M2           vec( 2.71000000, 1.05655099, 0.00000000)  // 53
#define E4_M2Y          vec( 0.00000000, 1.00000000, 0.00000000)  // 53
#define E43M2_rot001_COS  1.000000000f
#define E43M2_rot001_SIN -0.000000000f
#define E43M2_movx -2.71000000
#define E43M2_movy  1.05655099
// H[HIS_MD]
#define H1_CG           vec( 0.59156419, 1.39717224, 0.00000000)  // 6
#define H2_ND1          vec( 0.74365372, 1.16150173, 0.00000000)  // 7
#define H2_CD2          vec( 0.89178042,-1.01837843, 0.00000000)  // 6
#define H2_CE1          vec( 2.02989367, 0.85638485, 0.00000000)  // 6
#define H2_NE2          vec( 2.14584073,-0.45888741, 0.00000000)  // 8
#define H21_rot001_COS  0.389893133f
#define H21_rot001_SIN  0.920860111f
#define H21_movx    1.51724700
#define H2_MD           vec(-0.01328216, 3.12034018, 0.00000000)  // 53
#define H3_MDY          vec( 0.00000000, 1.00000000, 0.00000000)  // 53
#define H32MD_rot001_COS -0.360445659f
#define H32MD_rot001_SIN  0.932780214f
#define H32MD_movx  2.91537908
#define H32MD_movy -1.11232374
#define H2_ME           vec( 3.95036652,-1.53298549, 0.00000000)  // 53
#define H3_MEY          vec( 0.00000000, 1.00000000, 0.00000000)  // 53
#define H32ME_rot001_COS  0.859297995f
#define H32ME_rot001_SIN -0.511475273f
#define H32ME_movx  4.17862620
#define H32ME_movy  0.70322344


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
  struct VEC ML = (DorE) ? H2_MD  : H2_ME;
  struct VEC MY = (DorE) ? H3_MDY : H3_MEY;
  struct VEC ND = H2_ND1;
  struct VEC NE = H2_NE2;
  struct VEC CE = H2_CE1;
  struct VEC CD = H2_CD2;
  struct VEC CG = H1_CG;
  
  rotx(&MY,sin3,cos3);  // rot chi3
  MY.x += (DorE) ? H32MD_movx : H32ME_movx;  // move up
  MY.y += (DorE) ? H32MD_movy : H32ME_movy; 
  rotz(&MY,(DorE)?H32MD_rot001_SIN:H32ME_rot001_SIN,(DorE)?H32MD_rot001_COS:H32ME_rot001_COS); // tip so chi1 along x
  rotx(&MY,sin2,cos2);  // rot chi2
  rotx(&ML,sin2,cos2);
  rotx(&ND,sin2,cos2);
  rotx(&NE,sin2,cos2);
  rotx(&CE,sin2,cos2);
  rotx(&CD,sin2,cos2);
  MY.x += H21_movx;  // move up
  ML.x += H21_movx;  
  ND.x += H21_movx;
  NE.x += H21_movx;
  CE.x += H21_movx;
  CD.x += H21_movx;
  rotz(&MY,H21_rot001_SIN,H21_rot001_COS);  // tip so chi1 along x
  rotz(&ML,H21_rot001_SIN,H21_rot001_COS);
  rotz(&ND,H21_rot001_SIN,H21_rot001_COS);
  rotz(&NE,H21_rot001_SIN,H21_rot001_COS);
  rotz(&CE,H21_rot001_SIN,H21_rot001_COS);
  rotz(&CD,H21_rot001_SIN,H21_rot001_COS);
  rotx(&MY,sin1,cos1); // rot chi1
  rotx(&ML,sin1,cos1);
  rotx(&ND,sin1,cos1);
  rotx(&NE,sin1,cos1);
  rotx(&CE,sin1,cos1);
  rotx(&CD,sin1,cos1);
  rotx(&CG,sin1,cos1); // CG here only!!!

  // clash check?? here or other func?

  struct VEC const X = normalizedv(subvv(ML,(DorE)?ND:NE));
  struct VEC const Y = normalizedv(subvv(MY,ML));
  struct VEC const Z = crossvv(X,Y);
  std::cout << "hisd_bb2m MY: " << MY << std::endl;
  std::cout << "hisd_bb2m ML: " << ML << std::endl;
  std::cout << "hisd_bb2m ND: " << ND << std::endl;
  std::cout << "hisd_bb2m NE: " << NE << std::endl;
  std::cout << "hisd_bb2m CE: " << CE << std::endl;
  std::cout << "hisd_bb2m CD: " << CD << std::endl;
  std::cout << "hisd_bb2m CG: " << CG << std::endl;

  struct MAT const R = cols(X,Y,Z);
  return xform(R,ML);
}
inline struct XFORM const hisd_bb2m(float const chi1, float const chi2, float const chi3){ return       his_bb2m(chi1,chi2,chi3,true )  ; }
inline struct XFORM const hise_bb2m(float const chi1, float const chi2, float const chi3){ return       his_bb2m(chi1,chi2,chi3,false)  ; }
inline struct XFORM const hisd_m2bb(float const chi1, float const chi2, float const chi3){ return xrev( his_bb2m(chi1,chi2,chi3,true ) ); }
inline struct XFORM const hise_m2bb(float const chi1, float const chi2, float const chi3){ return xrev( his_bb2m(chi1,chi2,chi3,false) ); }

