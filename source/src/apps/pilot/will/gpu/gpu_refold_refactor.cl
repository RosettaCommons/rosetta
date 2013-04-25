void refold(
  __local const float *tor, 
  uint N,
  __local float *N__xyz,
  __local float *CA_xyz,
  __local float *C__xyz,
  __local float *O__xyz,
  __local float *CB_xyz,
  __local float *H__xyz
){

  // first fold up each individual residue... this is the only time trig functions are needed
  // geometry is hard-coded for canonical bond lengths / angles... this would neeed to be made
  // more flexible
  for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
    // first phase uses O,CB,H as temp storage
    // could be done more efficiently
    uint const i = get_global_id(0) + ichunk;
    uint const i3 = 3u*i;
    if( i < N-1u ) {
      __local float *N__tmp = ((i==0u)?N__xyz:O__xyz);
      __local float *CA_tmp = ((i==0u)?CA_xyz:CB_xyz);
      __local float *C__tmp = ((i==0u)?C__xyz:H__xyz);
      uint const in = i3+3;
      N__tmp[i3+Xi] =  0.0f;
      N__tmp[i3+Yi] =  SIN_N_CA_C * DIS_N_CA;
      N__tmp[i3+Zi] = -DIS_CA_C - COS_N_CA_C * DIS_N_CA;
      CA_tmp[i3+Xi] =  0.0f;
      CA_tmp[i3+Yi] =  0.0f;
      CA_tmp[i3+Zi] = -DIS_CA_C;
      C__tmp[i3+Xi] =  0.0f;
      C__tmp[i3+Yi] =  0.0f;
      C__tmp[i3+Zi] =  0.0f;
      N__xyz[i3+Xi] =  0.0f;
      N__xyz[i3+Yi] =  0.0f;
      N__xyz[i3+Zi] =  0.0f;
      CA_xyz[i3+Xi] =  0.0f;
      CA_xyz[i3+Yi] =  0.0f;
      CA_xyz[i3+Zi] =  0.0f;
      C__xyz[i3+Xi] =  0.0f;
      C__xyz[i3+Yi] =  0.0f;
      C__xyz[i3+Zi] =  0.0f;

      float ofst = ((float)get_global_id(1))/100000.0;

      // this first block could be merged into init!!!
      { // rot CA-C-N-CA, rot CA-N-C and drop
        float const COS_CA_C = cos(tor[3u*i+1u]+ofst);
        float const SIN_CA_C = sin(tor[3u*i+1u]+ofst);
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
        float const COS_C_N = cos(tor[3u*i+2u]+ofst);
        float const SIN_C_N = sin(tor[3u*i+2u]+ofst);
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
      N__xyz[i3+Zi] -= DIS_N_CA;
      { // rot phi2
        float const COS_N_CA = cos(tor[3u*i+3u]+ofst);
        float const SIN_N_CA = sin(tor[3u*i+3u]+ofst);
        float const  Nx = N__tmp[i3+Xi];
        float const  Ny = N__tmp[i3+Yi];
        float const CAx = CA_tmp[i3+Xi];
        float const CAy = CA_tmp[i3+Yi];
        float const  Cx = C__tmp[i3+Xi];
        float const  Cy = C__tmp[i3+Yi];
        float const N2x = N__xyz[i3+Xi];
        float const N2y = N__xyz[i3+Yi];
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
        float const N2y = N__xyz[i3+Yi];
        float const N2z = N__xyz[i3+Zi];
        N__tmp[i3+Yi] = mad( COS_N_CA_C ,  Ny, -SIN_N_CA_C *  Nz );
        N__tmp[i3+Zi] = mad( SIN_N_CA_C ,  Ny,  COS_N_CA_C *  Nz );
        CA_tmp[i3+Yi] = mad( COS_N_CA_C , CAy, -SIN_N_CA_C * CAz );
        CA_tmp[i3+Zi] = mad( SIN_N_CA_C , CAy,  COS_N_CA_C * CAz );
        C__tmp[i3+Yi] = mad( COS_N_CA_C ,  Cy, -SIN_N_CA_C *  Cz );
        C__tmp[i3+Zi] = mad( SIN_N_CA_C ,  Cy,  COS_N_CA_C *  Cz );
        N__xyz[i3+Yi] = mad( COS_N_CA_C , N2y, -SIN_N_CA_C * N2z );
        N__xyz[i3+Zi] = mad( SIN_N_CA_C , N2y,  COS_N_CA_C * N2z );
      }
      //TR<<F(5,2,N__tmp[i3+Xi])<<" "<<F(5,2,N__tmp[i3+Yi])<<" "<<F(5,2,CA_tmp[i3+Xi])<<" "<<F(5,2,CA_tmp[i3+Yi])<<" "<<F(5,2,C__tmp[i3+Xi])<<" "<<F(5,2,C__tmp[i3+Yi])<<" "<<F(5,2,C__tmp[i3+Xi])<<" "<<F(5,2,C__tmp[i3+Yi])<<endl;
      N__tmp[i3+Zi] -= DIS_CA_C;
      CA_tmp[i3+Zi] -= DIS_CA_C;
      C__tmp[i3+Zi] -= DIS_CA_C;
      N__xyz[i3+Zi] -= DIS_CA_C;
      CA_xyz[i3+Zi] -= DIS_CA_C;
    }
  }

  // now join segments of length 2 / 4 / 8 / 16 / 32 / 64 / 128 / .....
  for(uint c = 2u; c < N*2u-3u; c=2u*c) {
    for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
      uint const i3 = 3u*(get_global_id(0)+ichunk);
      uint j = get_global_id(0)+ichunk;
      barrier(CLK_LOCAL_MEM_FENCE);
      if( j < N ) {
        uint i = (max(0,((int)j))/c)*c+c/2;
        if(j > i || i > N-2u) continue;;
        uint const i9 = 9u*(i);
        struct MAT R;
        if(c<4u) { // skip "to" if 1st iter
          struct VEC const az = normalizedv(         vec(H__xyz[i3+0u]-CB_xyz[i3+0u],H__xyz[i3+1u]-CB_xyz[i3+1u],H__xyz[i3+2u]-CB_xyz[i3+2u]) );
          struct VEC const ay = normalizedv(pproj(az,vec(O__xyz[i3+0u]-CB_xyz[i3+0u],O__xyz[i3+1u]-CB_xyz[i3+1u],O__xyz[i3+2u]-CB_xyz[i3+2u])));
          R = cols(crossvv(ay,az),ay,az);
        } else {
          struct VEC       az = normalizedv(         vec(H__xyz[i3+0u]-CB_xyz[i3+0u],H__xyz[i3+1u]-CB_xyz[i3+1u],H__xyz[i3+2u]-CB_xyz[i3+2u]) );
          struct VEC       ay = normalizedv(pproj(az,vec(O__xyz[i3+0u]-CB_xyz[i3+0u],O__xyz[i3+1u]-CB_xyz[i3+1u],O__xyz[i3+2u]-CB_xyz[i3+2u])));
          struct MAT const to = cols(crossvv(ay,az),ay,az);
          az =                  normalizedv(         vec(C__xyz[i3+0u]-CA_xyz[i3+0u],C__xyz[i3+1u]-CA_xyz[i3+1u],C__xyz[i3+2u]-CA_xyz[i3+2u]) );
          ay =                  normalizedv(pproj(az,vec(N__xyz[i3+0u]-CA_xyz[i3+0u],N__xyz[i3+1u]-CA_xyz[i3+1u],N__xyz[i3+2u]-CA_xyz[i3+2u])));
          R = multmm(to,rows(crossvv(ay,az),ay,az));
        }
        struct VEC T = vec(N__xyz[i3+0u],N__xyz[i3+1u],N__xyz[i3+2u]);
        T = multmv(R,T);  T.x = O__xyz[i3+0u]-T.x;  T.y = O__xyz[i3+1u]-T.y;  T.z = O__xyz[i3+2u]-T.z;
        uint const start = max((int)i-(int)c/2,0);
        uint const j3 = 3u*j;
        __local float *N__tmp = ((j==i-c/2u&&j!=0u) ? O__xyz : N__xyz);
        __local float *CA_tmp = ((j==i-c/2u&&j!=0u) ? CB_xyz : CA_xyz);
        __local float *C__tmp = ((j==i-c/2u&&j!=0u) ? H__xyz : C__xyz);                
        struct VEC v1 = multmv(R,vec(N__tmp[i3+0u],N__tmp[i3+1u],N__tmp[i3+2u]));
        struct VEC v2 = multmv(R,vec(CA_tmp[i3+0u],CA_tmp[i3+1u],CA_tmp[i3+2u]));
        struct VEC v3 = multmv(R,vec(C__tmp[i3+0u],C__tmp[i3+1u],C__tmp[i3+2u]));
        N__tmp[i3+0u]=v1.x+T.x; N__tmp[i3+1u]=v1.y+T.y; N__tmp[i3+2u]=v1.z+T.z;
        CA_tmp[i3+0u]=v2.x+T.x; CA_tmp[i3+1u]=v2.y+T.y; CA_tmp[i3+2u]=v2.z+T.z;
        C__tmp[i3+0u]=v3.x+T.x; C__tmp[i3+1u]=v3.y+T.y; C__tmp[i3+2u]=v3.z+T.z;
      }
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
}
