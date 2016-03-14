typedef float Gfloat;

uint  const _N         = 0u;
uint  const CA         = 3u;
uint  const _C         = 6u;
uint  const Xi         = 0u;
uint  const Yi         = 1u;
uint  const Zi         = 2u;
Gfloat const COS_C_N_CA = 0.5254717f;
Gfloat const SIN_C_N_CA = 0.8508111f;
Gfloat const COS_N_CA_C = 0.3616246f;
Gfloat const SIN_N_CA_C = 0.9323238f;
Gfloat const COS_CA_C_N = 0.4415059f;
Gfloat const SIN_CA_C_N = 0.8972584f;
Gfloat const DIS_N_CA   = 1.4580010f;
Gfloat const DIS_CA_C   = 1.5232580f;
Gfloat const DIS_C_N    = 1.3286850f;

// void refold_gpu_test(uint N, float const *tor, float *bb, float *b2) {
//   Gfloat *SIN = new Gfloat[3*N];
//   Gfloat *COS = new Gfloat[3*N];

//   double t = time_highres();
//   for(uint i = 0; i < 3*N; ++i) {
//     SIN[i] = sin(-tor[i]);
//     COS[i] = cos( tor[i]);
//   }
//   for(uint i = 0; i < N-1u; ++i) {
//     //TR<<"init build "<<i<<std::endl;
//     float *XX = ((i==0u)?b2:bb);
//     uint const i9 = 9u*(i   );
//     uint const i2 = 9u*(i+1u);
//     XX[i9+_N+Xi] =  0.0f;
//     XX[i9+_N+Yi] =  SIN_N_CA_C * DIS_N_CA;
//     XX[i9+_N+Zi] = -DIS_CA_C - COS_N_CA_C * DIS_N_CA;
//     XX[i9+CA+Xi] =  0.0f;
//     XX[i9+CA+Yi] =  0.0f;
//     XX[i9+CA+Zi] = -DIS_CA_C;
//     XX[i9+_C+Xi] =  0.0f;
//     XX[i9+_C+Yi] =  0.0f;
//     XX[i9+_C+Zi] =  0.0f;
//     b2[i2+_N+Xi] =  0.0f;
//     b2[i2+_N+Yi] =  0.0f;
//     b2[i2+_N+Zi] =  0.0f;
//     b2[i2+CA+Xi] =  0.0f;
//     b2[i2+CA+Yi] =  0.0f;
//     b2[i2+CA+Zi] =  0.0f;
//     b2[i2+_C+Xi] =  0.0f;
//     b2[i2+_C+Yi] =  0.0f;
//     b2[i2+_C+Zi] =  0.0f;

//     // this first block could be merged into init!!!
//     { // rot CA-C-N-CA, rot CA-N-C and drop
//       Gfloat const COS_CA_C = COS[3u*i+1u];
//       Gfloat const SIN_CA_C = SIN[3u*i+1u];
//       Gfloat const Nx = XX[i9+_N+Xi];
//       Gfloat const Ny = XX[i9+_N+Yi];
//       XX[i9+_N+Xi] = mad( COS_CA_C , Nx, - SIN_CA_C * Ny );
//       XX[i9+_N+Yi] = mad( SIN_CA_C , Nx, + COS_CA_C * Ny );
//       // CA x/y are 0 at start, no need to rot!
//     }
//     { // rot CA_C_N bond angle
//       Gfloat const Ny  = XX[i9+_N+Yi];
//       Gfloat const Nz  = XX[i9+_N+Zi];
//       Gfloat const CAy = XX[i9+CA+Yi];
//       Gfloat const CAz = XX[i9+CA+Zi];
//       XX[i9+_N+Yi] = mad( COS_CA_C_N ,  Ny, - SIN_CA_C_N *  Nz );
//       XX[i9+_N+Zi] = mad( SIN_CA_C_N ,  Ny, + COS_CA_C_N *  Nz );
//       XX[i9+CA+Yi] = mad( COS_CA_C_N , CAy, - SIN_CA_C_N * CAz );
//       XX[i9+CA+Zi] = mad( SIN_CA_C_N , CAy, + COS_CA_C_N * CAz );
//     }
//     XX[i9+_N+Zi] -= DIS_C_N;
//     XX[i9+CA+Zi] -= DIS_C_N;
//     XX[i9+_C+Zi] -= DIS_C_N;
//     { // rot omega2
//       Gfloat const COS_C_N = COS[3u*i+2u];
//       Gfloat const SIN_C_N = SIN[3u*i+2u];
//       Gfloat const  Nx = XX[i9+_N+Xi];
//       Gfloat const  Ny = XX[i9+_N+Yi];
//       Gfloat const CAx = XX[i9+CA+Xi];
//       Gfloat const CAy = XX[i9+CA+Yi];
//       // Gfloat const  Cx = XX[i9+_C+Xi];
//       // Gfloat const  Cy = XX[i9+_C+Yi];
//       XX[i9+_N+Xi] = mad( COS_C_N ,  Nx, - SIN_C_N *  Ny );
//       XX[i9+_N+Yi] = mad( SIN_C_N ,  Nx, + COS_C_N *  Ny );
//       XX[i9+CA+Xi] = mad( COS_C_N , CAx, - SIN_C_N * CAy );
//       XX[i9+CA+Yi] = mad( SIN_C_N , CAx, + COS_C_N * CAy );
//     }
//     { // rot C_N_CA angle
//       Gfloat const  Ny = XX[i9+_N+Yi];
//       Gfloat const  Nz = XX[i9+_N+Zi];
//       Gfloat const CAy = XX[i9+CA+Yi];
//       Gfloat const CAz = XX[i9+CA+Zi];
//       Gfloat const  Cy = XX[i9+_C+Yi];
//       Gfloat const  Cz = XX[i9+_C+Zi];
//       XX[i9+_N+Yi] = mad( COS_C_N_CA ,  Ny, - SIN_C_N_CA *  Nz );
//       XX[i9+_N+Zi] = mad( SIN_C_N_CA ,  Ny, + COS_C_N_CA *  Nz );
//       XX[i9+CA+Yi] = mad( COS_C_N_CA , CAy, - SIN_C_N_CA * CAz );
//       XX[i9+CA+Zi] = mad( SIN_C_N_CA , CAy, + COS_C_N_CA * CAz );
//       XX[i9+_C+Yi] = mad( COS_C_N_CA ,  Cy, - SIN_C_N_CA *  Cz );
//       XX[i9+_C+Zi] = mad( SIN_C_N_CA ,  Cy, + COS_C_N_CA *  Cz );
//     }
//     XX[i9+_N+Zi] -= DIS_N_CA;
//     XX[i9+CA+Zi] -= DIS_N_CA;
//     XX[i9+_C+Zi] -= DIS_N_CA;
//     b2[i2+_N+Zi] -= DIS_N_CA;
//     { // rot phi2
//       Gfloat const COS_N_CA = COS[3u*i+3u];
//       Gfloat const SIN_N_CA = SIN[3u*i+3u];
//       Gfloat const  Nx = XX[i9+_N+Xi];
//       Gfloat const  Ny = XX[i9+_N+Yi];
//       Gfloat const CAx = XX[i9+CA+Xi];
//       Gfloat const CAy = XX[i9+CA+Yi];
//       Gfloat const  Cx = XX[i9+_C+Xi];
//       Gfloat const  Cy = XX[i9+_C+Yi];
//       Gfloat const N2x = b2[i2+_N+Xi];
//       Gfloat const N2y = b2[i2+_N+Yi];
//       XX[i9+_N+Xi] = mad( COS_N_CA ,  Nx, - SIN_N_CA *  Ny );
//       XX[i9+_N+Yi] = mad( SIN_N_CA ,  Nx, + COS_N_CA *  Ny );
//       XX[i9+CA+Xi] = mad( COS_N_CA , CAx, - SIN_N_CA * CAy );
//       XX[i9+CA+Yi] = mad( SIN_N_CA , CAx, + COS_N_CA * CAy );
//       XX[i9+_C+Xi] = mad( COS_N_CA ,  Cx, - SIN_N_CA *  Cy );
//       XX[i9+_C+Yi] = mad( SIN_N_CA ,  Cx, + COS_N_CA *  Cy );
//     }
//     { // rot C_CA_N angle
//       Gfloat const  Ny = XX[i9+_N+Yi];
//       Gfloat const  Nz = XX[i9+_N+Zi];
//       Gfloat const CAy = XX[i9+CA+Yi];
//       Gfloat const CAz = XX[i9+CA+Zi];
//       Gfloat const  Cy = XX[i9+_C+Yi];
//       Gfloat const  Cz = XX[i9+_C+Zi];
//       Gfloat const N2y = b2[i2+_N+Yi];
//       Gfloat const N2z = b2[i2+_N+Zi];
//       XX[i9+_N+Yi] = mad( COS_N_CA_C ,  Ny, - SIN_N_CA_C *  Nz );
//       XX[i9+_N+Zi] = mad( SIN_N_CA_C ,  Ny, + COS_N_CA_C *  Nz );
//       XX[i9+CA+Yi] = mad( COS_N_CA_C , CAy, - SIN_N_CA_C * CAz );
//       XX[i9+CA+Zi] = mad( SIN_N_CA_C , CAy, + COS_N_CA_C * CAz );
//       XX[i9+_C+Yi] = mad( COS_N_CA_C ,  Cy, - SIN_N_CA_C *  Cz );
//       XX[i9+_C+Zi] = mad( SIN_N_CA_C ,  Cy, + COS_N_CA_C *  Cz );
//       b2[i2+_N+Yi] = mad( COS_N_CA_C , N2y, - SIN_N_CA_C * N2z );
//       b2[i2+_N+Zi] = mad( SIN_N_CA_C , N2y, + COS_N_CA_C * N2z );
//     }
//     //TR<<F(5,2,XX[i9+_N+Xi])<<" "<<F(5,2,XX[i9+_N+Yi])<<" "<<F(5,2,XX[i9+CA+Xi])<<" "<<F(5,2,XX[i9+CA+Yi])<<" "<<F(5,2,XX[i9+_C+Xi])<<" "<<F(5,2,XX[i9+_C+Yi])<<" "<<F(5,2,XX[i9+_C+Xi])<<" "<<F(5,2,XX[i9+_C+Yi])<<endl;
//     XX[i9+_N+Zi] -= DIS_CA_C;
//     XX[i9+CA+Zi] -= DIS_CA_C;
//     XX[i9+_C+Zi] -= DIS_CA_C;
//     b2[i2+_N+Zi] -= DIS_CA_C;
//     b2[i2+CA+Zi] -= DIS_CA_C;

//   }
  
//   // for(int i = 0; i < p.n_residue(); ++i) {
//   //   Pose tmp;
//   //   core::pose::make_pose_from_sequence(tmp,"G",*crs,false);
//   //   xyzStripeHashPoseWithMeta(tmp);
//   //   tmp.set_xyz(AtomID(1,1),Vec(bb[9*i+0],bb[9*i+1],bb[9*i+2]));
//   //   tmp.set_xyz(AtomID(2,1),Vec(bb[9*i+3],bb[9*i+4],bb[9*i+5]));
//   //   tmp.set_xyz(AtomID(3,1),Vec(bb[9*i+6],bb[9*i+7],bb[9*i+8]));
//   //   tmp.dump_pdb("bb"+str(i)+".pdb");
//   //   tmp.set_xyz(AtomID(1,1),Vec(b2[9*i+0],b2[9*i+1],b2[9*i+2]));
//   //   tmp.set_xyz(AtomID(2,1),Vec(b2[9*i+3],b2[9*i+4],b2[9*i+5]));
//   //   tmp.set_xyz(AtomID(3,1),Vec(b2[9*i+6],b2[9*i+7],b2[9*i+8]));
//   //   tmp.dump_pdb("b2"+str(i)+".pdb");
//   // }

//   for(uint c = 2u; c < N*2u-3u; c=2u*c) {
//     //    TR << c << "----------------------" << endl;
//     //for(uint i = c/2u; i < N-1u; i+=c) {
//     for(uint j = 0; j < N; ++j) {
//       uint i = (max(0,((int)j))/c)*c+c/2;
//       if(j > i || i > N-2u) continue;
//       uint const i9 = 9u*(i);
//       MAT R;
//       if(c<4u) { // skip "to" if 1st iter
//         VEC const az = normalizedv(         vec(bb[i9+6u]-bb[i9+3u],bb[i9+7u]-bb[i9+4u],bb[i9+8u]-bb[i9+5u]) );
//         VEC const ay = normalizedv(pproj(az,vec(bb[i9+0u]-bb[i9+3u],bb[i9+1u]-bb[i9+4u],bb[i9+2u]-bb[i9+5u])));
//         R = cols(crossvv(ay,az),ay,az);
//       } else {
//         VEC       az = normalizedv(         vec(bb[i9+6u]-bb[i9+3u],bb[i9+7u]-bb[i9+4u],bb[i9+8u]-bb[i9+5u]) );
//         VEC       ay = normalizedv(pproj(az,vec(bb[i9+0u]-bb[i9+3u],bb[i9+1u]-bb[i9+4u],bb[i9+2u]-bb[i9+5u])));
//         MAT const to = cols(crossvv(ay,az),ay,az);
//         az = normalizedv(                   vec(b2[i9+6u]-b2[i9+3u],b2[i9+7u]-b2[i9+4u],b2[i9+8u]-b2[i9+5u]) );
//         ay = normalizedv(pproj(az,          vec(b2[i9+0u]-b2[i9+3u],b2[i9+1u]-b2[i9+4u],b2[i9+2u]-b2[i9+5u])));
//         R = multmm(to,rows(crossvv(ay,az),ay,az));
//       }
//       VEC T(b2[i9+0u],b2[i9+1u],b2[i9+2u]);
//       T = multmv(R,T);  T.x = bb[i9+0u]-T.x;  T.y = bb[i9+1u]-T.y;  T.z = bb[i9+2u]-T.z;
//       uint const start = max((int)i-(int)c/2,0);
//       //TR << "join " << c << " at " << i << " xform " << i <<"-"<< stop << endl;
//       //TR << c << " " << stop << endl;
//       //for(uint j = start; j <= i; ++j) {
//       //TR << c << " " << j << " " << i << endl;
//       uint const j9 = 9u*j;
//       float *XX = ((j==i-c/2u&&j!=0u) ? bb : b2);
//       //TR << "join " << c << " at " << i << " xform " << start <<"-"<< i <<" @ "<< j << " " << ((j==i-c/2u&&j!=0u)?"bb":"b2") << " " << endl;;
//       VEC v1 = multmv(R,vec(XX[j9+0u],XX[j9+1u],XX[j9+2u]));
//       VEC v2 = multmv(R,vec(XX[j9+3u],XX[j9+4u],XX[j9+5u]));
//       VEC v3 = multmv(R,vec(XX[j9+6u],XX[j9+7u],XX[j9+8u]));
//       XX[j9+0u]=v1.x+T.x; XX[j9+1u]=v1.y+T.y; XX[j9+2u]=v1.z+T.z;
//       XX[j9+3u]=v2.x+T.x; XX[j9+4u]=v2.y+T.y; XX[j9+5u]=v2.z+T.z;
//       XX[j9+6u]=v3.x+T.x; XX[j9+7u]=v3.y+T.y; XX[j9+8u]=v3.z+T.z;
//       //TR << v1 << " " << v2 << " " << v3;
//       //TR << endl;
//       //}
//     }

//   }
// }


//__kernel
void refold_first(
                  /*__global*/ const float *tor,
                  /*__global*/ const uint  *nres,
                  float *N__xyz, float *CA_xyz, float *C__xyz, float *O__xyz, float *CB_xyz, float *H__xyz, float *CENxyz
                  ){
  //  __local Gfloat bb[64*9];
  //  __local Gfloat b2[64*9];
  uint N = nres[0];
  //  for(uint i = 0; i < N-1u; ++i) {
  uint const i = get_global_id(0);
  uint const i3 = 3u*i;

  if( i < N-1u ) {
    //TR<<"init build "<<i<<std::endl;
    /*__local*/ float *N__tmp = ((i==0u)?N__xyz:O__xyz);
    /*__local*/ float *CA_tmp = ((i==0u)?CA_xyz:CB_xyz);
    /*__local*/ float *C__tmp = ((i==0u)?C__xyz:H__xyz);        
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
      Gfloat const COS_CA_C = cos(tor[3u*i+1u]);
      Gfloat const SIN_CA_C = sin(tor[3u*i+1u]);
      Gfloat const Nx = N__tmp[i3+Xi];
      Gfloat const Ny = N__tmp[i3+Yi];
      N__tmp[i3+Xi] = mad(  COS_CA_C , Nx, + SIN_CA_C * Ny );
      N__tmp[i3+Yi] = mad( -SIN_CA_C , Nx, + COS_CA_C * Ny );
      // CA x/y are 0 at start, no need to rot!
    }
    { // rot CA_C_N bond angle
      Gfloat const Ny  = N__tmp[i3+Yi];
      Gfloat const Nz  = N__tmp[i3+Zi];
      Gfloat const CAy = CA_tmp[i3+Yi];
      Gfloat const CAz = CA_tmp[i3+Zi];
      N__tmp[i3+Yi] = mad( COS_CA_C_N ,  Ny, -SIN_CA_C_N *  Nz );
      N__tmp[i3+Zi] = mad( SIN_CA_C_N ,  Ny,  COS_CA_C_N *  Nz );
      CA_tmp[i3+Yi] = mad( COS_CA_C_N , CAy, -SIN_CA_C_N * CAz );
      CA_tmp[i3+Zi] = mad( SIN_CA_C_N , CAy,  COS_CA_C_N * CAz );
    }
    N__tmp[i3+Zi] -= DIS_C_N;
    CA_tmp[i3+Zi] -= DIS_C_N;
    C__tmp[i3+Zi] -= DIS_C_N;
    { // rot omega2
      Gfloat const COS_C_N = cos(tor[3u*i+2u]);
      Gfloat const SIN_C_N = sin(tor[3u*i+2u]);
      Gfloat const  Nx = N__tmp[i3+Xi];
      Gfloat const  Ny = N__tmp[i3+Yi];
      Gfloat const CAx = CA_tmp[i3+Xi];
      Gfloat const CAy = CA_tmp[i3+Yi];
      // Gfloat const  Cx = C__tmp[i3+Xi];
      // Gfloat const  Cy = C__tmp[i3+Yi];
      N__tmp[i3+Xi] = mad(  COS_C_N ,  Nx, + SIN_C_N *  Ny );
      N__tmp[i3+Yi] = mad( -SIN_C_N ,  Nx, + COS_C_N *  Ny );
      CA_tmp[i3+Xi] = mad(  COS_C_N , CAx, + SIN_C_N * CAy );
      CA_tmp[i3+Yi] = mad( -SIN_C_N , CAx, + COS_C_N * CAy );
    }
    { // rot C_N_CA angle
      Gfloat const  Ny = N__tmp[i3+Yi];
      Gfloat const  Nz = N__tmp[i3+Zi];
      Gfloat const CAy = CA_tmp[i3+Yi];
      Gfloat const CAz = CA_tmp[i3+Zi];
      Gfloat const  Cy = C__tmp[i3+Yi];
      Gfloat const  Cz = C__tmp[i3+Zi];
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
      Gfloat const COS_N_CA = cos(tor[3u*i+3u]);
      Gfloat const SIN_N_CA = sin(tor[3u*i+3u]);
      Gfloat const  Nx = N__tmp[i3+Xi];
      Gfloat const  Ny = N__tmp[i3+Yi];
      Gfloat const CAx = CA_tmp[i3+Xi];
      Gfloat const CAy = CA_tmp[i3+Yi];
      Gfloat const  Cx = C__tmp[i3+Xi];
      Gfloat const  Cy = C__tmp[i3+Yi];
      Gfloat const N2x = N__xyz[in+Xi];
      Gfloat const N2y = N__xyz[in+Yi];
      N__tmp[i3+Xi] = mad(  COS_N_CA ,  Nx, SIN_N_CA *  Ny );
      N__tmp[i3+Yi] = mad( -SIN_N_CA ,  Nx, COS_N_CA *  Ny );
      CA_tmp[i3+Xi] = mad(  COS_N_CA , CAx, SIN_N_CA * CAy );
      CA_tmp[i3+Yi] = mad( -SIN_N_CA , CAx, COS_N_CA * CAy );
      C__tmp[i3+Xi] = mad(  COS_N_CA ,  Cx, SIN_N_CA *  Cy );
      C__tmp[i3+Yi] = mad( -SIN_N_CA ,  Cx, COS_N_CA *  Cy );
    }
    { // rot C_CA_N angle
      Gfloat const  Ny = N__tmp[i3+Yi];
      Gfloat const  Nz = N__tmp[i3+Zi];
      Gfloat const CAy = CA_tmp[i3+Yi];
      Gfloat const CAz = CA_tmp[i3+Zi];
      Gfloat const  Cy = C__tmp[i3+Yi];
      Gfloat const  Cz = C__tmp[i3+Zi];
      Gfloat const N2y = N__xyz[in+Yi];
      Gfloat const N2z = N__xyz[in+Zi];
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


void refold_second(
                   /*__global*/ const float *tor,
                   /*__global*/ const uint  *nres,
                   uint c,
                   float *N__xyz, float *CA_xyz, float *C__xyz, float *O__xyz, float *CB_xyz, float *H__xyz, float *CENxyz
                   ){
  uint N = nres[0];
  //for(uint c = 2u; c < N*2u-3u; c=2u*c) {
  uint j = get_global_id(0);
  if( j < N ) {
    uint i = (max(0,((int)j))/c)*c+c/2;
    if(j > i || i > N-2u) return; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    uint const i3 = 3u*(i);
    struct MAT R;
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
    struct VEC T = vec(N__xyz[i3+Xi],N__xyz[i3+Yi],N__xyz[i3+Zi]);
    T = multmv(R,T); 
    T.x = O__xyz[i3+Xi]-T.x; T.y = O__xyz[i3+Yi]-T.y; T.z = O__xyz[i3+Zi]-T.z;
    uint const start = max((int)i-(int)c/2,0);
    uint const j3 = 3u*j;
    /*__local*/ float *N__tmp = ((j==i-c/2u&&j!=0u) ? O__xyz : N__xyz );
    /*__local*/ float *CA_tmp = ((j==i-c/2u&&j!=0u) ? CB_xyz : CA_xyz );
    /*__local*/ float *C__tmp = ((j==i-c/2u&&j!=0u) ? H__xyz : C__xyz );        
    //TR << "join " << c << " at " << i << " xform " << start <<"-"<< i <<" @ "<< j << " " << ((j==i-c/2u&&j!=0u)?"bb":"b2") << " " << endl;;
    struct VEC v1 = multmv(R,vec(N__tmp[j3+Xi],N__tmp[j3+Yi],N__tmp[j3+Zi]));
    struct VEC v2 = multmv(R,vec(CA_tmp[j3+Xi],CA_tmp[j3+Yi],CA_tmp[j3+Zi]));
    struct VEC v3 = multmv(R,vec(C__tmp[j3+Xi],C__tmp[j3+Yi],C__tmp[j3+Zi]));
    N__tmp[j3+Xi]=v1.x+T.x; N__tmp[j3+Yi]=v1.y+T.y; N__tmp[j3+Zi]=v1.z+T.z;
    CA_tmp[j3+Xi]=v2.x+T.x; CA_tmp[j3+Yi]=v2.y+T.y; CA_tmp[j3+Zi]=v2.z+T.z;
    C__tmp[j3+Xi]=v3.x+T.x; C__tmp[j3+Yi]=v3.y+T.y; C__tmp[j3+Zi]=v3.z+T.z;
  }
}

void refold_third(
                  /*__global*/ const float *tor,
                  /*__global*/ const uint  *nres,
                  float *N__xyz, float *CA_xyz, float *C__xyz, float *O__xyz, float *CB_xyz, float *H__xyz, float *CENxyz,
                  const struct VEC *CB_LCOR, const struct VEC *CEN_LCOR,
                  const uint * aas
){
  uint const  i3 = 3u*(get_global_id(0u)                                              );
  uint const ip3 = 3u*(get_global_id(0u) ==      0u  ?      0u  : get_global_id(0u)-1u);
  uint const in3 = 3u*(get_global_id(0u) == nres[0u] ? nres[0u] : get_global_id(0u)+1u);
  O__xyz[i3+0u] = C__xyz[i3+0u] + 1.641f*(C__xyz[i3+0u]-(CA_xyz[i3+0u]+1.147f*N__xyz[in3+0u])/2.147f);
  O__xyz[i3+1u] = C__xyz[i3+1u] + 1.641f*(C__xyz[i3+1u]-(CA_xyz[i3+1u]+1.147f*N__xyz[in3+1u])/2.147f);
  O__xyz[i3+2u] = C__xyz[i3+2u] + 1.641f*(C__xyz[i3+2u]-(CA_xyz[i3+2u]+1.147f*N__xyz[in3+2u])/2.147f);
  H__xyz[i3+0u] = N__xyz[i3+0u] + 1.4915f*(N__xyz[i3+0u]-(CA_xyz[i3+0u]+1.1f*C__xyz[ip3+0u])/2.1f);
  H__xyz[i3+1u] = N__xyz[i3+1u] + 1.4915f*(N__xyz[i3+1u]-(CA_xyz[i3+1u]+1.1f*C__xyz[ip3+1u])/2.1f);
  H__xyz[i3+2u] = N__xyz[i3+2u] + 1.4915f*(N__xyz[i3+2u]-(CA_xyz[i3+2u]+1.1f*C__xyz[ip3+2u])/2.1f);
  struct XFORM bbstub = stub(vec(CA_xyz[i3+0u],CA_xyz[i3+1u],CA_xyz[i3+2u]),
                             vec(N__xyz[i3+0u],N__xyz[i3+1u],N__xyz[i3+2u]),
                             vec(C__xyz[i3+0u],C__xyz[i3+1u],C__xyz[i3+2u]));
  const struct VEC cb = multxv(bbstub,CB_LCOR [aas[get_global_id(0)]]);
  const struct VEC cn = multxv(bbstub,CEN_LCOR[aas[get_global_id(0)]]);
  CB_xyz[i3+0u] = cb.x;
  CB_xyz[i3+1u] = cb.y;
  CB_xyz[i3+2u] = cb.z;
  CENxyz[i3+0u] = cn.x;
  CENxyz[i3+1u] = cn.y;
  CENxyz[i3+2u] = cn.z;

//   uint N = nres[0];
//   float scratch[3];
//   float scratchN[N];
//   float scores[50];
//   uint RG = 1;


//   // COMPUTE RG // ASSUMES LOCAL DIM >= NRES!!!!!!!!

//   // blows away CEN coords
//   for(uint c=get_local_size(0)/2;c>0;c/=2) {
// //    barrier(CLK_LOCAL_MEM_FENCE);
//     if( c > get_local_id(0) && c < N ) {
//       uint i  = 3u* get_local_id(0)   ;
//       uint ic = 3u*(get_local_id(0)+c);
//       CENxyz[i+0u] += CENxyz[ic+0u];
//       CENxyz[i+1u] += CENxyz[ic+1u];
//       CENxyz[i+2u] += CENxyz[ic+2u];            
//     }
//   }
// //  barrier(CLK_LOCAL_MEM_FENCE);
//   if(get_local_id(0) < 3) {
//     scratch[get_local_id(0)] = CENxyz[get_local_id(0)] / (float)N;
//   }

//   // // restore CEN coords ... top method is clever, but slower
//   // for(uint c = 1; c<get_local_size(0);c*=2) {
//   //   barrier(CLK_LOCAL_MEM_FENCE);
//   //   if( c > get_local_id(0) && c < N ) {
//   //     uint i  = 3u* get_local_id(0)   ;
//   //     uint ic = 3u*(get_local_id(0)+c);
//   //     CENxyz[i+0u] -= CENxyz[ic+0u];
//   //     CENxyz[i+1u] -= CENxyz[ic+1u];
//   //     CENxyz[i+2u] -= CENxyz[ic+2u];            
//   //   }
//   // }
//   // barrier(CLK_LOCAL_MEM_FENCE);  
//   for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
//     uint const i = get_global_id(0u) + ichunk;
//     uint const  i3 = 3u*i;
//     struct XFORM bbstub = stub(vec(CA_xyz[i3+0u],CA_xyz[i3+1u],CA_xyz[i3+2u]),
//                                vec(N__xyz[i3+0u],N__xyz[i3+1u],N__xyz[i3+2u]),
//                                vec(C__xyz[i3+0u],C__xyz[i3+1u],C__xyz[i3+2u]));
//     const struct VEC cn = multxv(bbstub,CEN_LCOR[aas[get_global_id(0)]]);
//     CENxyz[i3+0u] = cn.x;
//     CENxyz[i3+1u] = cn.y;
//     CENxyz[i3+2u] = cn.z;
//   }

//   for(uint ichunk = 0; ichunk < N; ichunk+=get_global_size(0)) {
//     uint const i = get_global_id(0u) + ichunk;
//     uint const  i3 = 3u*i;
//     float const d0 = CENxyz[i3+0u]-scratch[0];
//     float const d1 = CENxyz[i3+1u]-scratch[1];
//     float const d2 = CENxyz[i3+2u]-scratch[2];
//     scratchN[i] = d0*d0+d1*d1+d2*d2;
//   }
//   for(uint c=get_local_size(0)/2;c>0;c/=2) {
// //    barrier(CLK_LOCAL_MEM_FENCE);
//     if(c>get_local_id(0) && c < N) scratchN[get_local_id(0)] += scratchN[get_local_id(0)+c];
//   }
// //  barrier(CLK_LOCAL_MEM_FENCE);
//   if(get_local_id(0)==0) scores[RG] = native_sqrt( scratchN[0] / (float)(N-1) );

}


