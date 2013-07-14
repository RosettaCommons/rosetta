#include </Users/sheffler/mini/src/apps/pilot/will/gpu/gpu_mat_vec.cl>
//#include </Users/sheffler/mini/src/apps/pilot/will/gpu/gpu_refold.cl>

#define _N 0u
#define CA 3u
#define _C 6u
#define Xi 0u
#define Yi 1u
#define Zi 2u
#define COS_C_N_CA 0.5254717f
#define SIN_C_N_CA 0.8508111f
#define COS_N_CA_C 0.3616246f
#define SIN_N_CA_C 0.9323238f
#define COS_CA_C_N 0.4415059f
#define SIN_CA_C_N 0.8972584f
#define DIS_N_CA 1.458001f
#define DIS_CA_C 1.523258f
#define DIS_C_N  1.328685f


void refold(
  __global const float *tor, 
  uint N,
  __local float *coord6
){
  __local float *N__xyz = coord6 + 0u*3u*N;
  __local float *CA_xyz = coord6 + 1u*3u*N;
  __local float *C__xyz = coord6 + 2u*3u*N;
  __local float *O__xyz = coord6 + 3u*3u*N;
  __local float *CB_xyz = coord6 + 4u*3u*N;
  __local float *H__xyz = coord6 + 5u*3u*N;

  // first phase uses O,CB,H as temp storage
  // could be done more efficiently
  for(uint ichunk = 0; ichunk < N; ichunk+=get_local_size(0)) {
    uint const i = get_local_id(0) + ichunk;
    uint const i3 = 3u*i;
    if( i < N-1u ) {
      __local float *N__tmp = ((i==0u)?N__xyz:O__xyz);
      __local float *CA_tmp = ((i==0u)?CA_xyz:CB_xyz);
      __local float *C__tmp = ((i==0u)?C__xyz:H__xyz);
      uint const in = 3u*(i+1u);
      N__tmp[i3+Xi] =  0.0f;
      N__tmp[i3+Yi] =  SIN_N_CA_C * DIS_N_CA;
      N__tmp[i3+Zi] = -DIS_CA_C - COS_N_CA_C * DIS_N_CA;
      CA_tmp[i3+Xi] =  0.0f;
      CA_tmp[i3+Yi] =  0.0f;
      CA_tmp[i3+Zi] = -DIS_CA_C;
      C__tmp[i3+Xi] =  0.0f;
      C__tmp[i3+Yi] =  0.0f;
      C__tmp[i3+Zi] =  0.0f;
      N__xyz[in+Xi] =  0.0f;
      N__xyz[in+Yi] =  0.0f;
      N__xyz[in+Zi] =  0.0f;
      CA_xyz[in+Xi] =  0.0f;
      CA_xyz[in+Yi] =  0.0f;
      CA_xyz[in+Zi] =  0.0f;
      C__xyz[in+Xi] =  0.0f;
      C__xyz[in+Yi] =  0.0f;
      C__xyz[in+Zi] =  0.0f;

      // this first block could be merged into init!!!
      { // rot CA-C-N-CA, rot CA-N-C and drop
        float const COS_CA_C = cos(tor[3u*i+1u]);
        float const SIN_CA_C = sin(tor[3u*i+1u]);
        float const Nx = N__tmp[i3+Xi];
        float const Ny = N__tmp[i3+Yi];
        N__tmp[i3+Xi] = mad(  COS_CA_C , Nx, + SIN_CA_C * Ny );
        N__tmp[i3+Yi] = mad( -SIN_CA_C , Nx, + COS_CA_C * Ny );
        // CA x/y are 0 at start, no need to rot!
      }
      { // rot CA_C_N bond angle
        float const Ny  = N__tmp[i3+Yi];
        float const Nz  = N__tmp[i3+Zi];
        float const CAy = CA_tmp[i3+Yi];
        float const CAz = CA_tmp[i3+Zi];
        N__tmp[i3+Yi] = mad( COS_CA_C_N ,  Ny, -SIN_CA_C_N *  Nz );
        N__tmp[i3+Zi] = mad( SIN_CA_C_N ,  Ny,  COS_CA_C_N *  Nz );
        CA_tmp[i3+Yi] = mad( COS_CA_C_N , CAy, -SIN_CA_C_N * CAz );
        CA_tmp[i3+Zi] = mad( SIN_CA_C_N , CAy,  COS_CA_C_N * CAz );
      }
      N__tmp[i3+Zi] -= DIS_C_N;
      CA_tmp[i3+Zi] -= DIS_C_N;
      C__tmp[i3+Zi] -= DIS_C_N;
      { // rot omega2
        float const COS_C_N = cos(tor[3u*i+2u]);
        float const SIN_C_N = sin(tor[3u*i+2u]);
        float const  Nx = N__tmp[i3+Xi];
        float const  Ny = N__tmp[i3+Yi];
        float const CAx = CA_tmp[i3+Xi];
        float const CAy = CA_tmp[i3+Yi];
        // float const  Cx = C__tmp[i3+Xi];
        // float const  Cy = C__tmp[i3+Yi];
        N__tmp[i3+Xi] = mad(  COS_C_N ,  Nx, + SIN_C_N *  Ny );
        N__tmp[i3+Yi] = mad( -SIN_C_N ,  Nx, + COS_C_N *  Ny );
        CA_tmp[i3+Xi] = mad(  COS_C_N , CAx, + SIN_C_N * CAy );
        CA_tmp[i3+Yi] = mad( -SIN_C_N , CAx, + COS_C_N * CAy );
      }
      { // rot C_N_CA angle
        float const  Ny = N__tmp[i3+Yi];
        float const  Nz = N__tmp[i3+Zi];
        float const CAy = CA_tmp[i3+Yi];
        float const CAz = CA_tmp[i3+Zi];
        float const  Cy = C__tmp[i3+Yi];
        float const  Cz = C__tmp[i3+Zi];
        N__tmp[i3+Yi] = mad( COS_C_N_CA ,  Ny, -SIN_C_N_CA *  Nz );
        N__tmp[i3+Zi] = mad( SIN_C_N_CA ,  Ny,  COS_C_N_CA *  Nz );
        CA_tmp[i3+Yi] = mad( COS_C_N_CA , CAy, -SIN_C_N_CA * CAz );
        CA_tmp[i3+Zi] = mad( SIN_C_N_CA , CAy,  COS_C_N_CA * CAz );
        C__tmp[i3+Yi] = mad( COS_C_N_CA ,  Cy, -SIN_C_N_CA *  Cz );
        C__tmp[i3+Zi] = mad( SIN_C_N_CA ,  Cy,  COS_C_N_CA *  Cz );
      }
      N__tmp[i3+Zi] -= DIS_N_CA;
      CA_tmp[i3+Zi] -= DIS_N_CA;
      C__tmp[i3+Zi] -= DIS_N_CA;
      N__xyz[in+Zi] -= DIS_N_CA;
      { // rot phi2
        float const COS_N_CA = cos(tor[3u*i+3u]);
        float const SIN_N_CA = sin(tor[3u*i+3u]);
        float const  Nx = N__tmp[i3+Xi];
        float const  Ny = N__tmp[i3+Yi];
        float const CAx = CA_tmp[i3+Xi];
        float const CAy = CA_tmp[i3+Yi];
        float const  Cx = C__tmp[i3+Xi];
        float const  Cy = C__tmp[i3+Yi];
        float const N2x = N__xyz[in+Xi];
        float const N2y = N__xyz[in+Yi];
        N__tmp[i3+Xi] = mad(  COS_N_CA ,  Nx, SIN_N_CA *  Ny );
        N__tmp[i3+Yi] = mad( -SIN_N_CA ,  Nx, COS_N_CA *  Ny );
        CA_tmp[i3+Xi] = mad(  COS_N_CA , CAx, SIN_N_CA * CAy );
        CA_tmp[i3+Yi] = mad( -SIN_N_CA , CAx, COS_N_CA * CAy );
        C__tmp[i3+Xi] = mad(  COS_N_CA ,  Cx, SIN_N_CA *  Cy );
        C__tmp[i3+Yi] = mad( -SIN_N_CA ,  Cx, COS_N_CA *  Cy );
      }
      { // rot C_CA_N angle
        float const  Ny = N__tmp[i3+Yi];
        float const  Nz = N__tmp[i3+Zi];
        float const CAy = CA_tmp[i3+Yi];
        float const CAz = CA_tmp[i3+Zi];
        float const  Cy = C__tmp[i3+Yi];
        float const  Cz = C__tmp[i3+Zi];
        float const N2y = N__xyz[in+Yi];
        float const N2z = N__xyz[in+Zi];
        N__tmp[i3+Yi] = mad( COS_N_CA_C ,  Ny, -SIN_N_CA_C *  Nz );
        N__tmp[i3+Zi] = mad( SIN_N_CA_C ,  Ny,  COS_N_CA_C *  Nz );
        CA_tmp[i3+Yi] = mad( COS_N_CA_C , CAy, -SIN_N_CA_C * CAz );
        CA_tmp[i3+Zi] = mad( SIN_N_CA_C , CAy,  COS_N_CA_C * CAz );
        C__tmp[i3+Yi] = mad( COS_N_CA_C ,  Cy, -SIN_N_CA_C *  Cz );
        C__tmp[i3+Zi] = mad( SIN_N_CA_C ,  Cy,  COS_N_CA_C *  Cz );
        N__xyz[in+Yi] = mad( COS_N_CA_C , N2y, -SIN_N_CA_C * N2z );
        N__xyz[in+Zi] = mad( SIN_N_CA_C , N2y,  COS_N_CA_C * N2z );
      }
      //TR<<F(5,2,N__tmp[i3+Xi])<<" "<<F(5,2,N__tmp[i3+Yi])<<" "<<F(5,2,CA_tmp[i3+Xi])<<" "<<F(5,2,CA_tmp[i3+Yi])<<" "<<F(5,2,C__tmp[i3+Xi])<<" "<<F(5,2,C__tmp[i3+Yi])<<" "<<F(5,2,C__tmp[i3+Xi])<<" "<<F(5,2,C__tmp[i3+Yi])<<endl;
      N__tmp[i3+Zi] -= DIS_CA_C;
      CA_tmp[i3+Zi] -= DIS_CA_C;
      C__tmp[i3+Zi] -= DIS_CA_C;
      N__xyz[in+Zi] -= DIS_CA_C;
      CA_xyz[in+Zi] -= DIS_CA_C;
    }

  }

  for(uint c = 2u; c < N*2u-3u; c=2u*c) {
    for(uint ichunk = 0; ichunk < N; ichunk+=get_local_size(0)) {
      //uint const i3 = 3u*(get_local_id(0)+ichunk);
      uint j = get_local_id(0)+ichunk;
      uint i = (max(0,((int)j))/c)*c+c/2;
      uint const i3 = 3u*(i);
      struct MAT R;
      barrier(CLK_LOCAL_MEM_FENCE);
      if( j < N && !(j > i || i > N-2u)) {
        if(c<4u) { // skip "to" if 1st iter
          struct VEC const az = normalizedv(         vec(H__xyz[i3+Xi]-CB_xyz[i3+Xi],H__xyz[i3+Yi]-CB_xyz[i3+Yi],H__xyz[i3+Zi]-CB_xyz[i3+Zi]) );
          struct VEC const ay = normalizedv(pproj(az,vec(O__xyz[i3+Xi]-CB_xyz[i3+Xi],O__xyz[i3+Yi]-CB_xyz[i3+Yi],O__xyz[i3+Zi]-CB_xyz[i3+Zi])));
          R = cols(crossvv(ay,az),ay,az);
        } else {
          struct VEC       az = normalizedv(         vec(H__xyz[i3+Xi]-CB_xyz[i3+Xi],H__xyz[i3+Yi]-CB_xyz[i3+Yi],H__xyz[i3+Zi]-CB_xyz[i3+Zi]) );
          struct VEC       ay = normalizedv(pproj(az,vec(O__xyz[i3+Xi]-CB_xyz[i3+Xi],O__xyz[i3+Yi]-CB_xyz[i3+Yi],O__xyz[i3+Zi]-CB_xyz[i3+Zi])));
          struct MAT const to = cols(crossvv(ay,az),ay,az);
          az = normalizedv(                   vec(C__xyz[i3+Xi]-CA_xyz[i3+Xi],C__xyz[i3+Yi]-CA_xyz[i3+Yi],C__xyz[i3+Zi]-CA_xyz[i3+Zi]) );
          ay = normalizedv(pproj(az,          vec(N__xyz[i3+Xi]-CA_xyz[i3+Xi],N__xyz[i3+Yi]-CA_xyz[i3+Yi],N__xyz[i3+Zi]-CA_xyz[i3+Zi])));
          R = multmm(to,rows(crossvv(ay,az),ay,az));
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      if( j < N && !(j > i || i > N-2u)) {
        struct VEC T = vec(N__xyz[i3+Xi],N__xyz[i3+Yi],N__xyz[i3+Zi]);
        T = multmv(R,T);  T.x = O__xyz[i3+0u]-T.x;  T.y = O__xyz[i3+1u]-T.y;  T.z = O__xyz[i3+2u]-T.z;
        uint const start = max((int)i-(int)c/2,0);
        uint const j3 = 3u*j;
        __local float *N__tmp = ((j==i-c/2u&&j!=0u) ? O__xyz : N__xyz);
        __local float *CA_tmp = ((j==i-c/2u&&j!=0u) ? CB_xyz : CA_xyz);
        __local float *C__tmp = ((j==i-c/2u&&j!=0u) ? H__xyz : C__xyz);                
        struct VEC v1 = multmv(R,vec(N__tmp[j3+Xi],N__tmp[j3+Yi],N__tmp[j3+Zi]));
        struct VEC v2 = multmv(R,vec(CA_tmp[j3+Xi],CA_tmp[j3+Yi],CA_tmp[j3+Zi]));
        struct VEC v3 = multmv(R,vec(C__tmp[j3+Xi],C__tmp[j3+Yi],C__tmp[j3+Zi]));
        N__tmp[j3+Xi]=v1.x+T.x; N__tmp[j3+Yi]=v1.y+T.y; N__tmp[j3+Zi]=v1.z+T.z;
        CA_tmp[j3+Xi]=v2.x+T.x; CA_tmp[j3+Yi]=v2.y+T.y; CA_tmp[j3+Zi]=v2.z+T.z;
        C__tmp[j3+Xi]=v3.x+T.x; C__tmp[j3+Yi]=v3.y+T.y; C__tmp[j3+Zi]=v3.z+T.z;
      }
      barrier(CLK_LOCAL_MEM_FENCE);      
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);


}


__kernel void refold_test(
                     __global const float *tors,
                     __global const uint  *nres,
                     __global       float *out,
                     __local        float *coord6
                     ){
//  __local float  coord6[7u*128u*3u];

  uint const N = nres[0];
  refold(tors,N,coord6);

  __local float *N__xyz = coord6 + 0u*3u*N;
  __local float *CA_xyz = coord6 + 1u*3u*N;
  __local float *C__xyz = coord6 + 2u*3u*N;
  if(get_global_id(0)==get_local_id(0)) {
    for(uint ichunk = 0; ichunk < N; ichunk+=get_local_size(0)) {
      uint const i = get_local_id(0) + ichunk;
      if( i < N) {
        uint const i9  =  9u*i;
        uint const i3  =  3u*i;    
        out[i9+0] = N__xyz[i3+0];
        out[i9+1] = N__xyz[i3+1];
        out[i9+2] = N__xyz[i3+2];
        out[i9+3] = CA_xyz[i3+0];
        out[i9+4] = CA_xyz[i3+1];
        out[i9+5] = CA_xyz[i3+2];
        out[i9+6] = C__xyz[i3+0];
        out[i9+7] = C__xyz[i3+1];
        out[i9+8] = C__xyz[i3+2];
      }
    }
  }
	
	// if(get_local_id(0)==0) {
	// 	for(int i = 0; i < N*3; ++i) {
	// 		out[i] = tors[i];
	// 	}
	// }


}
