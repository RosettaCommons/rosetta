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
  ALL,
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

  for(uint ichunk = 0; ichunk < N; ichunk+=get_local_size(0)) {
    uint const i = get_local_id(0u) + ichunk;
    if(i >= N) break;
    uint const  i3 = 3u*(i                              );
    uint const ip3 = 3u*(i==0 ? 0 : get_local_id(0u)-1u);
    uint const in3 = 3u*(i==N ? N : get_local_id(0u)+1u);
    O__xyz[i3+0u] = C__xyz[i3+0u] + 1.6410f*(C__xyz[i3+0u]-(CA_xyz[i3+0u]+1.147f*N__xyz[in3+0u])/2.147f);
    O__xyz[i3+1u] = C__xyz[i3+1u] + 1.6410f*(C__xyz[i3+1u]-(CA_xyz[i3+1u]+1.147f*N__xyz[in3+1u])/2.147f);
    O__xyz[i3+2u] = C__xyz[i3+2u] + 1.6410f*(C__xyz[i3+2u]-(CA_xyz[i3+2u]+1.147f*N__xyz[in3+2u])/2.147f);
    H__xyz[i3+0u] = N__xyz[i3+0u] + 1.4915f*(N__xyz[i3+0u]-(CA_xyz[i3+0u]+1.100f*C__xyz[ip3+0u])/2.100f);
    H__xyz[i3+1u] = N__xyz[i3+1u] + 1.4915f*(N__xyz[i3+1u]-(CA_xyz[i3+1u]+1.100f*C__xyz[ip3+1u])/2.100f);
    H__xyz[i3+2u] = N__xyz[i3+2u] + 1.4915f*(N__xyz[i3+2u]-(CA_xyz[i3+2u]+1.100f*C__xyz[ip3+2u])/2.100f);
    struct XFORM bbstub = stub(vec(CA_xyz[i3+0u],CA_xyz[i3+1u],CA_xyz[i3+2u]),
                               vec(N__xyz[i3+0u],N__xyz[i3+1u],N__xyz[i3+2u]),
                               vec(C__xyz[i3+0u],C__xyz[i3+1u],C__xyz[i3+2u]));
    const struct VEC cb = multxv(bbstub, CB_LCOR[aas[get_local_id(0)]]);
    const struct VEC cn = multxv(bbstub,CEN_LCOR[aas[get_local_id(0)]]);
    CB_xyz[i3+0u] = cb.x;
    CB_xyz[i3+1u] = cb.y;
    CB_xyz[i3+2u] = cb.z;
    CENxyz[i3+0u] = cn.x;
    CENxyz[i3+1u] = cn.y;
    CENxyz[i3+2u] = cn.z;
  }


  // COMPUTE RG // ASSUMES LOCAL DIM >= NRES!!!!!!!!
  // blows away CEN coords
  for(uint c=get_local_size(0)/2;c>0;c/=2) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if( c > get_local_id(0) && get_local_id(0)+c < N ) {
      uint i  = 3u* get_local_id(0)   ;
      uint ic = 3u*(get_local_id(0)+c);
      CENxyz[i+0u] += CENxyz[ic+0u];
      CENxyz[i+1u] += CENxyz[ic+1u];
      CENxyz[i+2u] += CENxyz[ic+2u];            
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if(get_local_id(0) < 3) {
    scratch[get_local_id(0)] = CENxyz[get_local_id(0)] / (float)N;
  }
  // // restore CEN coords ... top method is clever, but slower
  // for(uint c = 1; c<get_local_size(0);c*=2) {
  //   barrier(CLK_LOCAL_MEM_FENCE);
  //   if( c > get_local_id(0) && c < N ) {
  //     uint i  = 3u* get_local_id(0)   ;
  //     uint ic = 3u*(get_local_id(0)+c);
  //     CENxyz[i+0u] -= CENxyz[ic+0u];
  //     CENxyz[i+1u] -= CENxyz[ic+1u];
  //     CENxyz[i+2u] -= CENxyz[ic+2u];            
  //   }
  // }
  // barrier(CLK_LOCAL_MEM_FENCE);  
  for(uint ichunk = 0; ichunk < N; ichunk+=get_local_size(0)) {
    uint const i = get_local_id(0u) + ichunk;
    if(i >= N) break;
    uint const  i3 = 3u*i;
    struct XFORM bbstub = stub(vec(CA_xyz[i3+0u],CA_xyz[i3+1u],CA_xyz[i3+2u]),
                               vec(N__xyz[i3+0u],N__xyz[i3+1u],N__xyz[i3+2u]),
                               vec(C__xyz[i3+0u],C__xyz[i3+1u],C__xyz[i3+2u]));
    const struct VEC cn = multxv(bbstub,CEN_LCOR[aas[get_local_id(0)]]);
    CENxyz[i3+0u] = cn.x;
    CENxyz[i3+1u] = cn.y;
    CENxyz[i3+2u] = cn.z;
  }
  barrier(CLK_LOCAL_MEM_FENCE);  
  scratchN[get_local_id(0)] = 0.0f;
  for(uint ichunk = 0; ichunk < N; ichunk+=get_local_size(0)) {
    uint const i = get_local_id(0u) + ichunk;
    if(i >= N) break;
    float const d0 = CENxyz[i*3u+0u]-scratch[0u];
    float const d1 = CENxyz[i*3u+1u]-scratch[1u];
    float const d2 = CENxyz[i*3u+2u]-scratch[2u];
    scratchN[get_local_id(0)] += mad(d0,d0,mad(d1,d1,d2*d2));
  }
  // sum RG
  for(uint c=get_local_size(0)/2;c>0;c/=2) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if(c>get_local_id(0) ) scratchN[get_local_id(0)] += scratchN[get_local_id(0)+c];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if(get_local_id(0)==0) {
    scores[RG] = native_sqrt( scratchN[0] / (float)(N-1) );
    scores[ALL] += scores[RG];
  }

  uint nbc20 = 0;
  uint nbc16 = 0;
  uint nbc12 = 0;
  uint nbc10 = 0;
  uint nbc8  = 0;
  uint nbc6  = 0;            
  if(get_local_id(0) < N) {
    for(uint jr = 0; jr < N; ++jr) {
      float const d0 = CENxyz[3u*get_local_id(0)+0u]-CENxyz[3u*jr+0];
      float const d1 = CENxyz[3u*get_local_id(0)+1u]-CENxyz[3u*jr+1];
      float const d2 = CENxyz[3u*get_local_id(0)+2u]-CENxyz[3u*jr+2];
      float const dsq = mad(d0,d0,mad(d1,d1,d2*d2));
      nbc20 += (dsq < 400.0f);
      nbc16 += (dsq < 256.0f);      
      nbc12 += (dsq < 144.0f);
      nbc10 += (dsq < 100.0f);
      nbc8  += (dsq <  64.0f);
      nbc6  += (dsq <  36.0f);                
    }
  }
  //scratchN[get_local_id(0)] = nbc20+nbc16+nbc12+nbc10+nbc8+nbc6;

  // loop over ir,jr triangle, involves SQRT... not worth it? ~25% faster on 76 res pose
  // float localRG = 0.0;
  // float localVDW = 0.0;
  scratchN[get_local_id(0)] = 0.0f;
  uint const Nwork = (N*(N+1u))/2u-1;
  for(uint ichunk = N; ichunk <= Nwork; ichunk+=get_local_size(0) ) {
    if( ichunk+get_local_id(0) > Nwork) break;
    {
      uint const iwork = Nwork - (ichunk+get_local_id(0));
      uint const n = floor((native_sqrt(1.0f+8.0f*(float)iwork)-1.0f)/2.0f);
      uint const tmp = iwork - (n*(n+1u))/2u;
      uint const ir = min(tmp,n-tmp);
      uint const jr = max(tmp,n-tmp);
      // {
      //   float const d0 = CENxyz[3u*ir+0u]-CENxyz[3u*jr+0];
      //   float const d1 = CENxyz[3u*ir+1u]-CENxyz[3u*jr+1];
      //   float const d2 = CENxyz[3u*ir+2u]-CENxyz[3u*jr+2];
      //   if( 200.0 < mad(d0,d0,mad(d1,d1,d2*d2)) ) continue;
      // }
        
      {
        //for(uint ia = 0; ia < 5; ia++) {
            float const d0aa = coord7[(N*0+ir)*3u+0u]-coord7[(N*0u+jr)*3u+0];
            float const d1aa = coord7[(N*0+ir)*3u+1u]-coord7[(N*0u+jr)*3u+1];
            float const d2aa = coord7[(N*0+ir)*3u+2u]-coord7[(N*0u+jr)*3u+2];
            float const d0ab = coord7[(N*0+ir)*3u+0u]-coord7[(N*1u+jr)*3u+0];
            float const d1ab = coord7[(N*0+ir)*3u+1u]-coord7[(N*1u+jr)*3u+1];
            float const d2ab = coord7[(N*0+ir)*3u+2u]-coord7[(N*1u+jr)*3u+2];
            float const d0ac = coord7[(N*0+ir)*3u+0u]-coord7[(N*2u+jr)*3u+0];
            float const d1ac = coord7[(N*0+ir)*3u+1u]-coord7[(N*2u+jr)*3u+1];
            float const d2ac = coord7[(N*0+ir)*3u+2u]-coord7[(N*2u+jr)*3u+2];
            float const d0ad = coord7[(N*0+ir)*3u+0u]-coord7[(N*3u+jr)*3u+0];
            float const d1ad = coord7[(N*0+ir)*3u+1u]-coord7[(N*3u+jr)*3u+1];
            float const d2ad = coord7[(N*0+ir)*3u+2u]-coord7[(N*3u+jr)*3u+2];
            float const d0ae = coord7[(N*0+ir)*3u+0u]-coord7[(N*4u+jr)*3u+0];
            float const d1ae = coord7[(N*0+ir)*3u+1u]-coord7[(N*4u+jr)*3u+1];
            float const d2ae = coord7[(N*0+ir)*3u+2u]-coord7[(N*4u+jr)*3u+2];

            float const d0ba = coord7[(N*1+ir)*3u+0u]-coord7[(N*0u+jr)*3u+0];
            float const d1ba = coord7[(N*1+ir)*3u+1u]-coord7[(N*0u+jr)*3u+1];
            float const d2ba = coord7[(N*1+ir)*3u+2u]-coord7[(N*0u+jr)*3u+2];
            float const d0bb = coord7[(N*1+ir)*3u+0u]-coord7[(N*1u+jr)*3u+0];
            float const d1bb = coord7[(N*1+ir)*3u+1u]-coord7[(N*1u+jr)*3u+1];
            float const d2bb = coord7[(N*1+ir)*3u+2u]-coord7[(N*1u+jr)*3u+2];
            float const d0bc = coord7[(N*1+ir)*3u+0u]-coord7[(N*2u+jr)*3u+0];
            float const d1bc = coord7[(N*1+ir)*3u+1u]-coord7[(N*2u+jr)*3u+1];
            float const d2bc = coord7[(N*1+ir)*3u+2u]-coord7[(N*2u+jr)*3u+2];
            float const d0bd = coord7[(N*1+ir)*3u+0u]-coord7[(N*3u+jr)*3u+0];
            float const d1bd = coord7[(N*1+ir)*3u+1u]-coord7[(N*3u+jr)*3u+1];
            float const d2bd = coord7[(N*1+ir)*3u+2u]-coord7[(N*3u+jr)*3u+2];
            float const d0be = coord7[(N*1+ir)*3u+0u]-coord7[(N*4u+jr)*3u+0];
            float const d1be = coord7[(N*1+ir)*3u+1u]-coord7[(N*4u+jr)*3u+1];
            float const d2be = coord7[(N*1+ir)*3u+2u]-coord7[(N*4u+jr)*3u+2];

            float const d0ca = coord7[(N*2+ir)*3u+0u]-coord7[(N*0u+jr)*3u+0];
            float const d1ca = coord7[(N*2+ir)*3u+1u]-coord7[(N*0u+jr)*3u+1];
            float const d2ca = coord7[(N*2+ir)*3u+2u]-coord7[(N*0u+jr)*3u+2];
            float const d0cb = coord7[(N*2+ir)*3u+0u]-coord7[(N*1u+jr)*3u+0];
            float const d1cb = coord7[(N*2+ir)*3u+1u]-coord7[(N*1u+jr)*3u+1];
            float const d2cb = coord7[(N*2+ir)*3u+2u]-coord7[(N*1u+jr)*3u+2];
            float const d0cc = coord7[(N*2+ir)*3u+0u]-coord7[(N*2u+jr)*3u+0];
            float const d1cc = coord7[(N*2+ir)*3u+1u]-coord7[(N*2u+jr)*3u+1];
            float const d2cc = coord7[(N*2+ir)*3u+2u]-coord7[(N*2u+jr)*3u+2];
            float const d0cd = coord7[(N*2+ir)*3u+0u]-coord7[(N*3u+jr)*3u+0];
            float const d1cd = coord7[(N*2+ir)*3u+1u]-coord7[(N*3u+jr)*3u+1];
            float const d2cd = coord7[(N*2+ir)*3u+2u]-coord7[(N*3u+jr)*3u+2];
            float const d0ce = coord7[(N*2+ir)*3u+0u]-coord7[(N*4u+jr)*3u+0];
            float const d1ce = coord7[(N*2+ir)*3u+1u]-coord7[(N*4u+jr)*3u+1];
            float const d2ce = coord7[(N*2+ir)*3u+2u]-coord7[(N*4u+jr)*3u+2];

            float const d0da = coord7[(N*3+ir)*3u+0u]-coord7[(N*0u+jr)*3u+0];
            float const d1da = coord7[(N*3+ir)*3u+1u]-coord7[(N*0u+jr)*3u+1];
            float const d2da = coord7[(N*3+ir)*3u+2u]-coord7[(N*0u+jr)*3u+2];
            float const d0db = coord7[(N*3+ir)*3u+0u]-coord7[(N*1u+jr)*3u+0];
            float const d1db = coord7[(N*3+ir)*3u+1u]-coord7[(N*1u+jr)*3u+1];
            float const d2db = coord7[(N*3+ir)*3u+2u]-coord7[(N*1u+jr)*3u+2];
            float const d0dc = coord7[(N*3+ir)*3u+0u]-coord7[(N*2u+jr)*3u+0];
            float const d1dc = coord7[(N*3+ir)*3u+1u]-coord7[(N*2u+jr)*3u+1];
            float const d2dc = coord7[(N*3+ir)*3u+2u]-coord7[(N*2u+jr)*3u+2];
            float const d0dd = coord7[(N*3+ir)*3u+0u]-coord7[(N*3u+jr)*3u+0];
            float const d1dd = coord7[(N*3+ir)*3u+1u]-coord7[(N*3u+jr)*3u+1];
            float const d2dd = coord7[(N*3+ir)*3u+2u]-coord7[(N*3u+jr)*3u+2];
            float const d0de = coord7[(N*3+ir)*3u+0u]-coord7[(N*4u+jr)*3u+0];
            float const d1de = coord7[(N*3+ir)*3u+1u]-coord7[(N*4u+jr)*3u+1];
            float const d2de = coord7[(N*3+ir)*3u+2u]-coord7[(N*4u+jr)*3u+2];

            float const d0ea = coord7[(N*4+ir)*3u+0u]-coord7[(N*0u+jr)*3u+0];
            float const d1ea = coord7[(N*4+ir)*3u+1u]-coord7[(N*0u+jr)*3u+1];
            float const d2ea = coord7[(N*4+ir)*3u+2u]-coord7[(N*0u+jr)*3u+2];
            float const d0eb = coord7[(N*4+ir)*3u+0u]-coord7[(N*1u+jr)*3u+0];
            float const d1eb = coord7[(N*4+ir)*3u+1u]-coord7[(N*1u+jr)*3u+1];
            float const d2eb = coord7[(N*4+ir)*3u+2u]-coord7[(N*1u+jr)*3u+2];
            float const d0ec = coord7[(N*4+ir)*3u+0u]-coord7[(N*2u+jr)*3u+0];
            float const d1ec = coord7[(N*4+ir)*3u+1u]-coord7[(N*2u+jr)*3u+1];
            float const d2ec = coord7[(N*4+ir)*3u+2u]-coord7[(N*2u+jr)*3u+2];
            float const d0ed = coord7[(N*4+ir)*3u+0u]-coord7[(N*3u+jr)*3u+0];
            float const d1ed = coord7[(N*4+ir)*3u+1u]-coord7[(N*3u+jr)*3u+1];
            float const d2ed = coord7[(N*4+ir)*3u+2u]-coord7[(N*3u+jr)*3u+2];
            float const d0ee = coord7[(N*4+ir)*3u+0u]-coord7[(N*4u+jr)*3u+0];
            float const d1ee = coord7[(N*4+ir)*3u+1u]-coord7[(N*4u+jr)*3u+1];
            float const d2ee = coord7[(N*4+ir)*3u+2u]-coord7[(N*4u+jr)*3u+2];

            float const daa = ((jr-ir)<2) ? 0.0f : max( 0.0f, (5.518f-mad(d0aa,d0aa,mad(d1aa,d1aa,d2aa*d2aa))) * 0.4257054f );
            float const dab = ((jr-ir)<2) ? 0.0f : max( 0.0f, (5.803f-mad(d0ab,d0ab,mad(d1ab,d1ab,d2ab*d2ab))) * 0.4151201f );
            float const dac = ((jr-ir)<1) ? 0.0f : max( 0.0f, (7.508f-mad(d0ac,d0ac,mad(d1ac,d1ac,d2ac*d2ac))) * 0.3649538f );
            float const dad = ((jr-ir)<1) ? 0.0f : max( 0.0f, (4.584f-mad(d0ad,d0ad,mad(d1ad,d1ad,d2ad*d2ad))) * 0.4670654f );
            float const dae = ((jr-ir)<1) ? 0.0f : max( 0.0f, (9.175f-mad(d0ae,d0ae,mad(d1ae,d1ae,d2ae*d2ae))) * 0.3301391f );
            float const dba = ((jr-ir)<2) ? 0.0f : max( 0.0f, (5.803f-mad(d0ba,d0ba,mad(d1ba,d1ba,d2ba*d2ba))) * 0.4151201f );
            float const dbb = ((jr-ir)<2) ? 0.0f : max( 0.0f, (7.670f-mad(d0bb,d0bb,mad(d1bb,d1bb,d2bb*d2bb))) * 0.3610791f );
            float const dbc = ((jr-ir)<2) ? 0.0f : max( 0.0f, (5.871f-mad(d0bc,d0bc,mad(d1bc,d1bc,d2bc*d2bc))) * 0.4127090f );
            float const dbd = ((jr-ir)<1) ? 0.0f : max( 0.0f, (6.744f-mad(d0bd,d0bd,mad(d1bd,d1bd,d2bd*d2bd))) * 0.3850714f );
            float const dbe = ((jr-ir)<2) ? 0.0f : max( 0.0f, (9.797f-mad(d0be,d0be,mad(d1be,d1be,d2be*d2be))) * 0.3194872f );
            float const dca = ((jr-ir)<3) ? 0.0f : max( 0.0f, (7.508f-mad(d0ca,d0ca,mad(d1ca,d1ca,d2ca*d2ca))) * 0.3649538f );
            float const dcb = ((jr-ir)<2) ? 0.0f : max( 0.0f, (5.871f-mad(d0cb,d0cb,mad(d1cb,d1cb,d2cb*d2cb))) * 0.4127090f );
            float const dcc = ((jr-ir)<2) ? 0.0f : max( 0.0f, (6.880f-mad(d0cc,d0cc,mad(d1cc,d1cc,d2cc*d2cc))) * 0.3812464f );
            float const dcd = ((jr-ir)<2) ? 0.0f : max( 0.0f, (6.600f-mad(d0cd,d0cd,mad(d1cd,d1cd,d2cd*d2cd))) * 0.3892495f );
            float const dce = ((jr-ir)<2) ? 0.0f : max( 0.0f, (9.778f-mad(d0ce,d0ce,mad(d1ce,d1ce,d2ce*d2ce))) * 0.3197974f );
            float const dda = ((jr-ir)<2) ? 0.0f : max( 0.0f, (4.584f-mad(d0da,d0da,mad(d1da,d1da,d2da*d2da))) * 0.4670654f );
            float const ddb = ((jr-ir)<2) ? 0.0f : max( 0.0f, (6.744f-mad(d0db,d0db,mad(d1db,d1db,d2db*d2db))) * 0.3850714f );
            float const ddc = ((jr-ir)<2) ? 0.0f : max( 0.0f, (6.600f-mad(d0dc,d0dc,mad(d1dc,d1dc,d2dc*d2dc))) * 0.3892495f );
            float const ddd = ((jr-ir)<1) ? 0.0f : max( 0.0f, (4.889f-mad(d0dd,d0dd,mad(d1dd,d1dd,d2dd*d2dd))) * 0.4522619f );
            float const dde = ((jr-ir)<2) ? 0.0f : max( 0.0f, (6.589f-mad(d0de,d0de,mad(d1de,d1de,d2de*d2de))) * 0.3895743f );
            float const dea = ((jr-ir)<2) ? 0.0f : max( 0.0f, (9.175f-mad(d0ea,d0ea,mad(d1ea,d1ea,d2ea*d2ea))) * 0.3301391f );
            float const deb = ((jr-ir)<2) ? 0.0f : max( 0.0f, (9.797f-mad(d0eb,d0eb,mad(d1eb,d1eb,d2eb*d2eb))) * 0.3194872f );
            float const dec = ((jr-ir)<1) ? 0.0f : max( 0.0f, (9.778f-mad(d0ec,d0ec,mad(d1ec,d1ec,d2ec*d2ec))) * 0.3197974f );
            float const ded = ((jr-ir)<1) ? 0.0f : max( 0.0f, (6.589f-mad(d0ed,d0ed,mad(d1ed,d1ed,d2ed*d2ed))) * 0.3895743f );
            float const dee = ((jr-ir)<1) ? 0.0f : max( 0.0f, (9.653f-mad(d0ee,d0ee,mad(d1ee,d1ee,d2ee*d2ee))) * 0.3218614f );

            scratchN[get_local_id(0)] += mad(daa,daa,mad(dab,dab,mad(dac,dac,mad(dad,dad,dae*dae))))
                                      +  mad(dba,dba,mad(dbb,dbb,mad(dbc,dbc,mad(dbd,dbd,dbe*dbe))))
                                      +  mad(dca,dca,mad(dcb,dcb,mad(dcc,dcc,mad(dcd,dcd,dce*dce))))
                                      +  mad(dda,dda,mad(ddb,ddb,mad(ddc,ddc,mad(ddd,ddd,dde*dde))))
                                      +  mad(dea,dea,mad(deb,deb,mad(dec,dec,mad(ded,ded,dee*dee))));

          //}
        //}
      }
    }
  }


  // sum VDW
  // scratchN[get_local_id(0)] = localVDW;
  for(uint c=get_local_size(0)/2;c>0;c/=2) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if(c>get_local_id(0) ) scratchN[get_local_id(0)] += scratchN[get_local_id(0)+c];
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if(get_local_id(0)==0) {
    scores[VDW] = native_sqrt( scratchN[0] / (float)(N*(N-1)) );
    scores[ALL] += scores[VDW];
  }

  // // sum VDW
  // // scratchN[get_local_id(0)] = localVDW;
  // for(uint c=get_local_size(0)/2;c>0;c/=2) {
  //   barrier(CLK_LOCAL_MEM_FENCE);
  //   if(c>get_local_id(0) ) scratchN[get_local_id(0)] += scratchN[get_local_id(0)+c];
  // }
  // barrier(CLK_LOCAL_MEM_FENCE);
  // if(get_local_id(0)==0) {
  //   scores[VDW] = native_sqrt( scratchN[0] / (float)(N*(N-1)) );
  //   scores[ALL] += scores[VDW];
  // }

  // // sum VDW
  // // scratchN[get_local_id(0)] = localVDW;
  // for(uint c=get_local_size(0)/2;c>0;c/=2) {
  //   barrier(CLK_LOCAL_MEM_FENCE);
  //   if(c>get_local_id(0) ) scratchN[get_local_id(0)] += scratchN[get_local_id(0)+c];
  // }
  // barrier(CLK_LOCAL_MEM_FENCE);
  // if(get_local_id(0)==0) {
  //   scores[VDW] = native_sqrt( scratchN[0] / (float)(N*(N-1)) );
  //   scores[ALL] += scores[VDW];
  // }

}

void copy_torsions(uint N, __global float *a, __global float const *b) {
  for(uint ichunk = 0u; ichunk < 3u*N; ichunk+=get_local_size(0)) {
    uint i = ichunk + get_local_id(0);
    if(i < 3u*N) a[i] = b[i];
  }
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
  printf("initial refold %i",N);
  refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);

  __local float      score;
  __local float best_score;
  __local float  low_score;    
  if(get_local_id(0)==0)      score = scores[ALL];
  if(get_local_id(0)==0) best_score = scores[ALL];
  if(get_local_id(0)==0)  low_score = scores[ALL];
  __global float *best_tors = tors + 3u*N;
  __global float * low_tors = tors + 6u*N;
  copy_torsions(N,best_tors,tors);
  copy_torsions(N, low_tors,tors);

  for(uint iter = 1; iter < 10; ++iter) {
    //tors[43] += 10.0f;
    refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
    if(get_local_id(0)==0) score = scores[ALL];
    if(score < best_score) copy_torsions(N,best_tors,tors);
    if(score <  low_score) copy_torsions(N, low_tors,tors);

    // if(score > best_score) {
    //   copy_torsions(N,tors,best_tors);
    //   refold(tors,N,coord7,aas,CB_LCOR,CEN_LCOR,scratch,scratchN,scores);
    //   if(get_local_id(0)==0) score = scores[ALL];
    // }
  }

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
