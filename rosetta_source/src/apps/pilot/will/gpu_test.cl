#include </Users/sheffler/mini/src/apps/pilot/will/gpu_mat_vec.cl>
//#include </Users/sheffler/mini/src/apps/pilot/will/gpu_refold.cl>


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
  __local const float *tor, 
  uint N,
  __local float *bb,
  __local float *b2
  // __local float *N__xyz,
  // __local float *CA_xyz,
  // __local float *C__xyz,
  // __local float *O__xyz,
  // __local float *CB_xyz,
  // __local float *H__xyz,
  // __local float *CENxyz
){

  for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
    uint const i = get_global_id(0) + ichunk;
    uint const i9 = 9u*i;
    if( i < N-1u ) {
      __local float *XX = ((i==0u)?b2:bb);
      uint const i2 = 9u*(i+1u);
      XX[i9+_N+Xi] =  0.0f;
      XX[i9+_N+Yi] =  SIN_N_CA_C * DIS_N_CA;
      XX[i9+_N+Zi] = -DIS_CA_C - COS_N_CA_C * DIS_N_CA;
      XX[i9+CA+Xi] =  0.0f;
      XX[i9+CA+Yi] =  0.0f;
      XX[i9+CA+Zi] = -DIS_CA_C;
      XX[i9+_C+Xi] =  0.0f;
      XX[i9+_C+Yi] =  0.0f;
      XX[i9+_C+Zi] =  0.0f;
      b2[i2+_N+Xi] =  0.0f;
      b2[i2+_N+Yi] =  0.0f;
      b2[i2+_N+Zi] =  0.0f;
      b2[i2+CA+Xi] =  0.0f;
      b2[i2+CA+Yi] =  0.0f;
      b2[i2+CA+Zi] =  0.0f;
      b2[i2+_C+Xi] =  0.0f;
      b2[i2+_C+Yi] =  0.0f;
      b2[i2+_C+Zi] =  0.0f;

      float ofst = ((float)get_global_id(1))/100000.0;

      // this first block could be merged into init!!!
      { // rot CA-C-N-CA, rot CA-N-C and drop
        float const COS_CA_C = cos(tor[3u*i+1u]+ofst);
        float const SIN_CA_C = sin(tor[3u*i+1u]+ofst);
        float const Nx = XX[i9+_N+Xi];
        float const Ny = XX[i9+_N+Yi];
        XX[i9+_N+Xi] = mad(  COS_CA_C , Nx, + SIN_CA_C * Ny );
        XX[i9+_N+Yi] = mad( -SIN_CA_C , Nx, + COS_CA_C * Ny );
        // CA x/y are 0 at start, no need to rot!
      }
      { // rot CA_C_N bond angle
        float const Ny  = XX[i9+_N+Yi];
        float const Nz  = XX[i9+_N+Zi];
        float const CAy = XX[i9+CA+Yi];
        float const CAz = XX[i9+CA+Zi];
        XX[i9+_N+Yi] = mad( COS_CA_C_N ,  Ny, -SIN_CA_C_N *  Nz );
        XX[i9+_N+Zi] = mad( SIN_CA_C_N ,  Ny,  COS_CA_C_N *  Nz );
        XX[i9+CA+Yi] = mad( COS_CA_C_N , CAy, -SIN_CA_C_N * CAz );
        XX[i9+CA+Zi] = mad( SIN_CA_C_N , CAy,  COS_CA_C_N * CAz );
      }
      XX[i9+_N+Zi] -= DIS_C_N;
      XX[i9+CA+Zi] -= DIS_C_N;
      XX[i9+_C+Zi] -= DIS_C_N;
      { // rot omega2
        float const COS_C_N = cos(tor[3u*i+2u]+ofst);
        float const SIN_C_N = sin(tor[3u*i+2u]+ofst);
        float const  Nx = XX[i9+_N+Xi];
        float const  Ny = XX[i9+_N+Yi];
        float const CAx = XX[i9+CA+Xi];
        float const CAy = XX[i9+CA+Yi];
        // float const  Cx = XX[i9+_C+Xi];
        // float const  Cy = XX[i9+_C+Yi];
        XX[i9+_N+Xi] = mad(  COS_C_N ,  Nx, + SIN_C_N *  Ny );
        XX[i9+_N+Yi] = mad( -SIN_C_N ,  Nx, + COS_C_N *  Ny );
        XX[i9+CA+Xi] = mad(  COS_C_N , CAx, + SIN_C_N * CAy );
        XX[i9+CA+Yi] = mad( -SIN_C_N , CAx, + COS_C_N * CAy );
      }
      { // rot C_N_CA angle
        float const  Ny = XX[i9+_N+Yi];
        float const  Nz = XX[i9+_N+Zi];
        float const CAy = XX[i9+CA+Yi];
        float const CAz = XX[i9+CA+Zi];
        float const  Cy = XX[i9+_C+Yi];
        float const  Cz = XX[i9+_C+Zi];
        XX[i9+_N+Yi] = mad( COS_C_N_CA ,  Ny, -SIN_C_N_CA *  Nz );
        XX[i9+_N+Zi] = mad( SIN_C_N_CA ,  Ny,  COS_C_N_CA *  Nz );
        XX[i9+CA+Yi] = mad( COS_C_N_CA , CAy, -SIN_C_N_CA * CAz );
        XX[i9+CA+Zi] = mad( SIN_C_N_CA , CAy,  COS_C_N_CA * CAz );
        XX[i9+_C+Yi] = mad( COS_C_N_CA ,  Cy, -SIN_C_N_CA *  Cz );
        XX[i9+_C+Zi] = mad( SIN_C_N_CA ,  Cy,  COS_C_N_CA *  Cz );
      }
      XX[i9+_N+Zi] -= DIS_N_CA;
      XX[i9+CA+Zi] -= DIS_N_CA;
      XX[i9+_C+Zi] -= DIS_N_CA;
      b2[i2+_N+Zi] -= DIS_N_CA;
      { // rot phi2
        float const COS_N_CA = cos(tor[3u*i+3u]+ofst);
        float const SIN_N_CA = sin(tor[3u*i+3u]+ofst);
        float const  Nx = XX[i9+_N+Xi];
        float const  Ny = XX[i9+_N+Yi];
        float const CAx = XX[i9+CA+Xi];
        float const CAy = XX[i9+CA+Yi];
        float const  Cx = XX[i9+_C+Xi];
        float const  Cy = XX[i9+_C+Yi];
        float const N2x = b2[i2+_N+Xi];
        float const N2y = b2[i2+_N+Yi];
        XX[i9+_N+Xi] = mad(  COS_N_CA ,  Nx, SIN_N_CA *  Ny );
        XX[i9+_N+Yi] = mad( -SIN_N_CA ,  Nx, COS_N_CA *  Ny );
        XX[i9+CA+Xi] = mad(  COS_N_CA , CAx, SIN_N_CA * CAy );
        XX[i9+CA+Yi] = mad( -SIN_N_CA , CAx, COS_N_CA * CAy );
        XX[i9+_C+Xi] = mad(  COS_N_CA ,  Cx, SIN_N_CA *  Cy );
        XX[i9+_C+Yi] = mad( -SIN_N_CA ,  Cx, COS_N_CA *  Cy );
      }
      { // rot C_CA_N angle
        float const  Ny = XX[i9+_N+Yi];
        float const  Nz = XX[i9+_N+Zi];
        float const CAy = XX[i9+CA+Yi];
        float const CAz = XX[i9+CA+Zi];
        float const  Cy = XX[i9+_C+Yi];
        float const  Cz = XX[i9+_C+Zi];
        float const N2y = b2[i2+_N+Yi];
        float const N2z = b2[i2+_N+Zi];
        XX[i9+_N+Yi] = mad( COS_N_CA_C ,  Ny, -SIN_N_CA_C *  Nz );
        XX[i9+_N+Zi] = mad( SIN_N_CA_C ,  Ny,  COS_N_CA_C *  Nz );
        XX[i9+CA+Yi] = mad( COS_N_CA_C , CAy, -SIN_N_CA_C * CAz );
        XX[i9+CA+Zi] = mad( SIN_N_CA_C , CAy,  COS_N_CA_C * CAz );
        XX[i9+_C+Yi] = mad( COS_N_CA_C ,  Cy, -SIN_N_CA_C *  Cz );
        XX[i9+_C+Zi] = mad( SIN_N_CA_C ,  Cy,  COS_N_CA_C *  Cz );
        b2[i2+_N+Yi] = mad( COS_N_CA_C , N2y, -SIN_N_CA_C * N2z );
        b2[i2+_N+Zi] = mad( SIN_N_CA_C , N2y,  COS_N_CA_C * N2z );
      }
      //TR<<F(5,2,XX[i9+_N+Xi])<<" "<<F(5,2,XX[i9+_N+Yi])<<" "<<F(5,2,XX[i9+CA+Xi])<<" "<<F(5,2,XX[i9+CA+Yi])<<" "<<F(5,2,XX[i9+_C+Xi])<<" "<<F(5,2,XX[i9+_C+Yi])<<" "<<F(5,2,XX[i9+_C+Xi])<<" "<<F(5,2,XX[i9+_C+Yi])<<endl;
      XX[i9+_N+Zi] -= DIS_CA_C;
      XX[i9+CA+Zi] -= DIS_CA_C;
      XX[i9+_C+Zi] -= DIS_CA_C;
      b2[i2+_N+Zi] -= DIS_CA_C;
      b2[i2+CA+Zi] -= DIS_CA_C;
    }
  }

  for(uint c = 2u; c < N*2u-3u; c=2u*c) {
    for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
      uint j = get_global_id(0)+ichunk;
      uint i = (max(0,((int)j))/c)*c+c/2;
      uint const i9 = 9u*(i);
      struct MAT R;
      barrier(CLK_LOCAL_MEM_FENCE);
      if( j < N && !(j > i || i > N-2u)) {
        if(c<4u) { // skip "to" if 1st iter
          struct VEC const az = normalizedv(         vec(bb[i9+6u]-bb[i9+3u],bb[i9+7u]-bb[i9+4u],bb[i9+8u]-bb[i9+5u]) );
          struct VEC const ay = normalizedv(pproj(az,vec(bb[i9+0u]-bb[i9+3u],bb[i9+1u]-bb[i9+4u],bb[i9+2u]-bb[i9+5u])));
          R = cols(crossvv(ay,az),ay,az);
        } else {
          struct VEC       az = normalizedv(         vec(bb[i9+6u]-bb[i9+3u],bb[i9+7u]-bb[i9+4u],bb[i9+8u]-bb[i9+5u]) );
          struct VEC       ay = normalizedv(pproj(az,vec(bb[i9+0u]-bb[i9+3u],bb[i9+1u]-bb[i9+4u],bb[i9+2u]-bb[i9+5u])));
          struct MAT const to = cols(crossvv(ay,az),ay,az);
          az = normalizedv(                   vec(b2[i9+6u]-b2[i9+3u],b2[i9+7u]-b2[i9+4u],b2[i9+8u]-b2[i9+5u]) );
          ay = normalizedv(pproj(az,          vec(b2[i9+0u]-b2[i9+3u],b2[i9+1u]-b2[i9+4u],b2[i9+2u]-b2[i9+5u])));
          R = multmm(to,rows(crossvv(ay,az),ay,az));
        }
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      if( j < N && !(j > i || i > N-2u)) {
        struct VEC T = vec(b2[i9+0u],b2[i9+1u],b2[i9+2u]);
        T = multmv(R,T);  T.x = bb[i9+0u]-T.x;  T.y = bb[i9+1u]-T.y;  T.z = bb[i9+2u]-T.z;
        uint const start = max((int)i-(int)c/2,0);
        uint const j9 = 9u*j;
        __local float *XX = ((j==i-c/2u&&j!=0u) ? bb : b2);
        struct VEC v1 = multmv(R,vec(XX[j9+0u],XX[j9+1u],XX[j9+2u]));
        struct VEC v2 = multmv(R,vec(XX[j9+3u],XX[j9+4u],XX[j9+5u]));
        struct VEC v3 = multmv(R,vec(XX[j9+6u],XX[j9+7u],XX[j9+8u]));
        XX[j9+0u]=v1.x+T.x; XX[j9+1u]=v1.y+T.y; XX[j9+2u]=v1.z+T.z;
        XX[j9+3u]=v2.x+T.x; XX[j9+4u]=v2.y+T.y; XX[j9+5u]=v2.z+T.z;
        XX[j9+6u]=v3.x+T.x; XX[j9+7u]=v3.y+T.y; XX[j9+8u]=v3.z+T.z;
      }
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
}


__kernel void abinitio(
                     __global const float *tor_in,
                     __global const uint  *nres,
                     __global       float *out
                     ){
  uint N = nres[0];
  __local float torsions[341*3];
  __local float       bb[341*9];
  __local float       b2[341*9];
  __local float      cen[341*3];
  // __local float N__xyz[128*3];
  // __local float CA_xyz[128*3];
  // __local float C__xyz[128*3];
  // __local float O__xyz[128*3];
  // __local float CB_xyz[128*3];
  // __local float H__xyz[128*3];
  // __local float CENxyz[128*3];

  for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
    int i = 3*(get_global_id(0)+ichunk);
    if(i < 3*N) {
      torsions[i+0] = tor_in[i+0];
      torsions[i+1] = tor_in[i+1];
      torsions[i+2] = tor_in[i+2];            
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);

//  refold(tor,N,out,N__xyz,CA_xyz,C__xyz,O__xyz,CB_xyz,H__xyz,CENxyz);
  refold(torsions,N,bb,b2);

  for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
    uint const i = get_global_id(0) + ichunk;
    uint const i9 = 9u*i;
    out[9u*N*get_global_id(1)+i9+0u] = b2[i9+0u];
    out[9u*N*get_global_id(1)+i9+1u] = b2[i9+1u];
    out[9u*N*get_global_id(1)+i9+2u] = b2[i9+2u];
    out[9u*N*get_global_id(1)+i9+3u] = b2[i9+3u];
    out[9u*N*get_global_id(1)+i9+4u] = b2[i9+4u];
    out[9u*N*get_global_id(1)+i9+5u] = b2[i9+5u];
    out[9u*N*get_global_id(1)+i9+6u] = b2[i9+6u];
    out[9u*N*get_global_id(1)+i9+7u] = b2[i9+7u];
    out[9u*N*get_global_id(1)+i9+8u] = b2[i9+8u];
  }
}
