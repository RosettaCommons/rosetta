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


enum SCORETYPE {
  VDW,
  RG
};

void refold(
  __global const float *tor, 
  uint N,
  __local float *coord7,
  __constant uint *aas,
  __constant struct VEC * CB_LCOR,  
  __constant struct VEC *CEN_LCOR,
  __local float *scratch,
  __local float *scratchN,
  __global float *scores
){
  __local float *N__xyz = coord7 + 0u*3u*N;
  __local float *CA_xyz = coord7 + 1u*3u*N;
  __local float *C__xyz = coord7 + 2u*3u*N;
  __local float *O__xyz = coord7 + 3u*3u*N;
  __local float *CB_xyz = coord7 + 4u*3u*N;
  __local float *H__xyz = coord7 + 5u*3u*N;
  __local float *CENxyz = coord7 + 6u*3u*N;

  // first phase uses O,CB,H as temp storage
  // could be done more efficiently
  for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
    uint const i = get_global_id(0) + ichunk;
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
    for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
      //uint const i3 = 3u*(get_global_id(0)+ichunk);
      uint j = get_global_id(0)+ichunk;
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

  for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
    uint const i = get_global_id(0u) + ichunk;
    if(i >= N) break;
    uint const  i3 = 3u*(i                              );
    uint const ip3 = 3u*(i==0 ? 0 : get_global_id(0u)-1u);
    uint const in3 = 3u*(i==N ? N : get_global_id(0u)+1u);
    O__xyz[i3+0u] = C__xyz[i3+0u] + 1.6410f*(C__xyz[i3+0u]-(CA_xyz[i3+0u]+1.147f*N__xyz[in3+0u])/2.147f);
    O__xyz[i3+1u] = C__xyz[i3+1u] + 1.6410f*(C__xyz[i3+1u]-(CA_xyz[i3+1u]+1.147f*N__xyz[in3+1u])/2.147f);
    O__xyz[i3+2u] = C__xyz[i3+2u] + 1.6410f*(C__xyz[i3+2u]-(CA_xyz[i3+2u]+1.147f*N__xyz[in3+2u])/2.147f);
    H__xyz[i3+0u] = N__xyz[i3+0u] + 1.4915f*(N__xyz[i3+0u]-(CA_xyz[i3+0u]+1.100f*C__xyz[ip3+0u])/2.100f);
    H__xyz[i3+1u] = N__xyz[i3+1u] + 1.4915f*(N__xyz[i3+1u]-(CA_xyz[i3+1u]+1.100f*C__xyz[ip3+1u])/2.100f);
    H__xyz[i3+2u] = N__xyz[i3+2u] + 1.4915f*(N__xyz[i3+2u]-(CA_xyz[i3+2u]+1.100f*C__xyz[ip3+2u])/2.100f);
    struct XFORM bbstub = stub(vec(CA_xyz[i3+0u],CA_xyz[i3+1u],CA_xyz[i3+2u]),
                               vec(N__xyz[i3+0u],N__xyz[i3+1u],N__xyz[i3+2u]),
                               vec(C__xyz[i3+0u],C__xyz[i3+1u],C__xyz[i3+2u]));
    const struct VEC cb = multxv(bbstub, CB_LCOR[aas[get_global_id(0)]]);
    const struct VEC cn = multxv(bbstub,CEN_LCOR[aas[get_global_id(0)]]);
    CB_xyz[i3+0u] = cb.x;
    CB_xyz[i3+1u] = cb.y;
    CB_xyz[i3+2u] = cb.z;
    CENxyz[i3+0u] = cn.x;
    CENxyz[i3+1u] = cn.y;
    CENxyz[i3+2u] = cn.z;
  }


  // // COMPUTE RG // ASSUMES LOCAL DIM >= NRES!!!!!!!!

  // // blows away CEN coords
  // for(uint c=get_local_size(0)/2;c>0;c/=2) {
  //   barrier(CLK_LOCAL_MEM_FENCE);
  //   if( c > get_local_id(0) && c < N ) {
  //     uint i  = 3u* get_local_id(0)   ;
  //     uint ic = 3u*(get_local_id(0)+c);
  //     CENxyz[i+0u] += CENxyz[ic+0u];
  //     CENxyz[i+1u] += CENxyz[ic+1u];
  //     CENxyz[i+2u] += CENxyz[ic+2u];            
  //   }
  // }
  // barrier(CLK_LOCAL_MEM_FENCE);
  // if(get_local_id(0) < 3) {
  //   scratch[get_local_id(0)] = CENxyz[get_local_id(0)] / (float)N;
  // }

  // // // restore CEN coords ... top method is clever, but slower
  // // for(uint c = 1; c<get_local_size(0);c*=2) {
  // //   barrier(CLK_LOCAL_MEM_FENCE);
  // //   if( c > get_local_id(0) && c < N ) {
  // //     uint i  = 3u* get_local_id(0)   ;
  // //     uint ic = 3u*(get_local_id(0)+c);
  // //     CENxyz[i+0u] -= CENxyz[ic+0u];
  // //     CENxyz[i+1u] -= CENxyz[ic+1u];
  // //     CENxyz[i+2u] -= CENxyz[ic+2u];            
  // //   }
  // // }
  // // barrier(CLK_LOCAL_MEM_FENCE);  
  // for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
  //   uint const i = get_global_id(0u) + ichunk;
  //   uint const  i3 = 3u*i;
  //   struct XFORM bbstub = stub(vec(CA_xyz[i3+0u],CA_xyz[i3+1u],CA_xyz[i3+2u]),
  //                              vec(N__xyz[i3+0u],N__xyz[i3+1u],N__xyz[i3+2u]),
  //                              vec(C__xyz[i3+0u],C__xyz[i3+1u],C__xyz[i3+2u]));
  //   const struct VEC cn = multxv(bbstub,CEN_LCOR[aas[get_global_id(0)]]);
  //   CENxyz[i3+0u] = cn.x;
  //   CENxyz[i3+1u] = cn.y;
  //   CENxyz[i3+2u] = cn.z;
  // }

  // for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
  //   uint const i = get_global_id(0u) + ichunk;
  //   uint const  i3 = 3u*i;
  //   float const d0 = CENxyz[i3+0u]-scratch[0];
  //   float const d1 = CENxyz[i3+1u]-scratch[1];
  //   float const d2 = CENxyz[i3+2u]-scratch[2];
  //   scratchN[i] = d0*d0+d1*d1+d2*d2;
  // }

  // efficiently loop over ir,jr
  float localRG = 0.0;
  float localVDW = 0.0;

  uint const Nwork = (N*(N+1u))/2u-1;
  for(uint ichunk = N; ichunk <= Nwork; ichunk+=get_local_size(0) ) {
    if( ichunk+get_local_id(0) <= Nwork) {            
      uint const i = Nwork - (ichunk+get_local_id(0));
      uint const n = floor((native_sqrt(1.0f+8.0f*(float)i)-1.0f)/2.0f);
      uint const jr = i - (n*(n+1u))/2u;
      uint const ir = n-jr;
      {
        float const d0 = CENxyz[3u*ir+0u]-CENxyz[3u*jr+0];
        float const d1 = CENxyz[3u*ir+1u]-CENxyz[3u*jr+1];
        float const d2 = CENxyz[3u*ir+2u]-CENxyz[3u*jr+2];
        localRG += d0*d0+d1*d1+d2*d2;
      }

      {
        float const d0 = N__xyz[3u*ir+0u]-N__xyz[3u*jr+0];
        float const d1 = N__xyz[3u*ir+1u]-N__xyz[3u*jr+1];
        float const d2 = N__xyz[3u*ir+2u]-N__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
        if(d0*d0+d1*d1+d2*d2 > 100.0) continue;
      }
      {
        float const d0 = N__xyz[3u*ir+0u]-CA_xyz[3u*jr+0];
        float const d1 = N__xyz[3u*ir+1u]-CA_xyz[3u*jr+1];
        float const d2 = N__xyz[3u*ir+2u]-CA_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = N__xyz[3u*ir+0u]-C__xyz[3u*jr+0];
        float const d1 = N__xyz[3u*ir+1u]-C__xyz[3u*jr+1];
        float const d2 = N__xyz[3u*ir+2u]-C__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = N__xyz[3u*ir+0u]-O__xyz[3u*jr+0];
        float const d1 = N__xyz[3u*ir+1u]-O__xyz[3u*jr+1];
        float const d2 = N__xyz[3u*ir+2u]-O__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = N__xyz[3u*ir+0u]-CB_xyz[3u*jr+0];
        float const d1 = N__xyz[3u*ir+1u]-CB_xyz[3u*jr+1];
        float const d2 = N__xyz[3u*ir+2u]-CB_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CA_xyz[3u*ir+0u]-N__xyz[3u*jr+0];
        float const d1 = CA_xyz[3u*ir+1u]-N__xyz[3u*jr+1];
        float const d2 = CA_xyz[3u*ir+2u]-N__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CA_xyz[3u*ir+0u]-CA_xyz[3u*jr+0];
        float const d1 = CA_xyz[3u*ir+1u]-CA_xyz[3u*jr+1];
        float const d2 = CA_xyz[3u*ir+2u]-CA_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CA_xyz[3u*ir+0u]-C__xyz[3u*jr+0];
        float const d1 = CA_xyz[3u*ir+1u]-C__xyz[3u*jr+1];
        float const d2 = CA_xyz[3u*ir+2u]-C__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CA_xyz[3u*ir+0u]-O__xyz[3u*jr+0];
        float const d1 = CA_xyz[3u*ir+1u]-O__xyz[3u*jr+1];
        float const d2 = CA_xyz[3u*ir+2u]-O__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CA_xyz[3u*ir+0u]-CB_xyz[3u*jr+0];
        float const d1 = CA_xyz[3u*ir+1u]-CB_xyz[3u*jr+1];
        float const d2 = CA_xyz[3u*ir+2u]-CB_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = C__xyz[3u*ir+0u]-N__xyz[3u*jr+0];
        float const d1 = C__xyz[3u*ir+1u]-N__xyz[3u*jr+1];
        float const d2 = C__xyz[3u*ir+2u]-N__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = C__xyz[3u*ir+0u]-CA_xyz[3u*jr+0];
        float const d1 = C__xyz[3u*ir+1u]-CA_xyz[3u*jr+1];
        float const d2 = C__xyz[3u*ir+2u]-CA_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = C__xyz[3u*ir+0u]-C__xyz[3u*jr+0];
        float const d1 = C__xyz[3u*ir+1u]-C__xyz[3u*jr+1];
        float const d2 = C__xyz[3u*ir+2u]-C__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = C__xyz[3u*ir+0u]-O__xyz[3u*jr+0];
        float const d1 = C__xyz[3u*ir+1u]-O__xyz[3u*jr+1];
        float const d2 = C__xyz[3u*ir+2u]-O__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = C__xyz[3u*ir+0u]-CB_xyz[3u*jr+0];
        float const d1 = C__xyz[3u*ir+1u]-CB_xyz[3u*jr+1];
        float const d2 = C__xyz[3u*ir+2u]-CB_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = O__xyz[3u*ir+0u]-N__xyz[3u*jr+0];
        float const d1 = O__xyz[3u*ir+1u]-N__xyz[3u*jr+1];
        float const d2 = O__xyz[3u*ir+2u]-N__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = O__xyz[3u*ir+0u]-CA_xyz[3u*jr+0];
        float const d1 = O__xyz[3u*ir+1u]-CA_xyz[3u*jr+1];
        float const d2 = O__xyz[3u*ir+2u]-CA_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = O__xyz[3u*ir+0u]-C__xyz[3u*jr+0];
        float const d1 = O__xyz[3u*ir+1u]-C__xyz[3u*jr+1];
        float const d2 = O__xyz[3u*ir+2u]-C__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = O__xyz[3u*ir+0u]-O__xyz[3u*jr+0];
        float const d1 = O__xyz[3u*ir+1u]-O__xyz[3u*jr+1];
        float const d2 = O__xyz[3u*ir+2u]-O__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = O__xyz[3u*ir+0u]-CB_xyz[3u*jr+0];
        float const d1 = O__xyz[3u*ir+1u]-CB_xyz[3u*jr+1];
        float const d2 = O__xyz[3u*ir+2u]-CB_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CB_xyz[3u*ir+0u]-N__xyz[3u*jr+0];
        float const d1 = CB_xyz[3u*ir+1u]-N__xyz[3u*jr+1];
        float const d2 = CB_xyz[3u*ir+2u]-N__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CB_xyz[3u*ir+0u]-CA_xyz[3u*jr+0];
        float const d1 = CB_xyz[3u*ir+1u]-CA_xyz[3u*jr+1];
        float const d2 = CB_xyz[3u*ir+2u]-CA_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CB_xyz[3u*ir+0u]-C__xyz[3u*jr+0];
        float const d1 = CB_xyz[3u*ir+1u]-C__xyz[3u*jr+1];
        float const d2 = CB_xyz[3u*ir+2u]-C__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CB_xyz[3u*ir+0u]-O__xyz[3u*jr+0];
        float const d1 = CB_xyz[3u*ir+1u]-O__xyz[3u*jr+1];
        float const d2 = CB_xyz[3u*ir+2u]-O__xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
      {
        float const d0 = CB_xyz[3u*ir+0u]-CB_xyz[3u*jr+0];
        float const d1 = CB_xyz[3u*ir+1u]-CB_xyz[3u*jr+1];
        float const d2 = CB_xyz[3u*ir+2u]-CB_xyz[3u*jr+2];
        localVDW += d0*d0+d1*d1+d2*d2;
      }
    }
  }

  // sum RG
  scratchN[get_local_id(0)] = localRG;
  for(uint c=get_local_size(0)/2;c>0;c/=2) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if(c>get_local_id(0) ) scratchN[get_local_id(0)] += scratchN[get_local_id(0)+c];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if(get_local_id(0)==0) scores[RG] = native_sqrt( scratchN[0] / (float)(N*(N-1)) );

  // sum VDW
  scratchN[get_local_id(0)] = localVDW;
  for(uint c=get_local_size(0)/2;c>0;c/=2) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if(c>get_local_id(0) ) scratchN[get_local_id(0)] += scratchN[get_local_id(0)+c];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if(get_local_id(0)==0) scores[VDW] = native_sqrt( scratchN[0] / (float)(N*(N-1)) );

}


__kernel void abinitio(
                     __global       float *tors,
                     __global const uint  *nres,
                     __global       float *out,
                     __global       float *scores,
                     __local float *coord7,
                     __local float *scratchN,
                     __constant uint *aas,
                     __constant struct VEC * CB_LCOR,
                     __constant struct VEC *CEN_LCOR        
                     ){
//  __local float  coord7[7u*128u*3u];

  uint N = nres[0];
  __local float scratch[3u];
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
//  refold(torsions,N,bb,b2);

  __local float *N__xyz = coord7 + 0u*3u*N;
  __local float *CA_xyz = coord7 + 1u*3u*N;
  __local float *C__xyz = coord7 + 2u*3u*N;
  __local float *O__xyz = coord7 + 3u*3u*N;
  __local float *CB_xyz = coord7 + 4u*3u*N;
  __local float *H__xyz = coord7 + 5u*3u*N;
  __local float *CENxyz = coord7 + 6u*3u*N;
  if(get_global_id(0)==get_local_id(0)) {
    for(uint ichunk = 0; ichunk < N; ichunk+=get_local_size(0)) {
      uint const i = get_local_id(0) + ichunk;
      if( i < N) {
        uint const i21 = 21u*i;
        uint const i3  =  3u*i;    
        out[21u*N*get_local_id(1)+i21+_N+Xi] = N__xyz[i3+Xi];
        out[21u*N*get_local_id(1)+i21+_N+Yi] = N__xyz[i3+Yi];
        out[21u*N*get_local_id(1)+i21+_N+Zi] = N__xyz[i3+Zi];
        out[21u*N*get_local_id(1)+i21+CA+Xi] = CA_xyz[i3+Xi];
        out[21u*N*get_local_id(1)+i21+CA+Yi] = CA_xyz[i3+Yi];
        out[21u*N*get_local_id(1)+i21+CA+Zi] = CA_xyz[i3+Zi];
        out[21u*N*get_local_id(1)+i21+_C+Xi] = C__xyz[i3+Xi];
        out[21u*N*get_local_id(1)+i21+_C+Yi] = C__xyz[i3+Yi];
        out[21u*N*get_local_id(1)+i21+_C+Zi] = C__xyz[i3+Zi];
        out[21u*N*get_local_id(1)+i21+ 9+Xi] = O__xyz[i3+Xi];
        out[21u*N*get_local_id(1)+i21+ 9+Yi] = O__xyz[i3+Yi];
        out[21u*N*get_local_id(1)+i21+ 9+Zi] = O__xyz[i3+Zi];
        out[21u*N*get_local_id(1)+i21+12+Xi] = CB_xyz[i3+Xi];
        out[21u*N*get_local_id(1)+i21+12+Yi] = CB_xyz[i3+Yi];
        out[21u*N*get_local_id(1)+i21+12+Zi] = CB_xyz[i3+Zi];
        out[21u*N*get_local_id(1)+i21+15+Xi] = CENxyz[i3+Xi];
        out[21u*N*get_local_id(1)+i21+15+Yi] = CENxyz[i3+Yi];
        out[21u*N*get_local_id(1)+i21+15+Zi] = CENxyz[i3+Zi];
        out[21u*N*get_local_id(1)+i21+18+Xi] = H__xyz[i3+Xi];
        out[21u*N*get_local_id(1)+i21+18+Yi] = H__xyz[i3+Yi];
        out[21u*N*get_local_id(1)+i21+18+Zi] = H__xyz[i3+Zi];
      }
    }
  }
}
