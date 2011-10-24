#include <apps/pilot/will/gpu_mat_vec.hh>
#include <apps/pilot/will/gpu_refold_test_cpu.hh>
#include <numeric/xyz.io.hh>

void set_pose_to_ideal(core::pose::Pose & p) {
  float const COS_C_N_CA = 0.5254717f;
  float const SIN_C_N_CA = 0.8508111f;
  float const COS_N_CA_C = 0.3616246f;
  float const SIN_N_CA_C = 0.9323238f;
  float const COS_CA_C_N = 0.4415059f;
  float const SIN_CA_C_N = 0.8972584f;
  float const DIS_N_CA   = 1.4580010f;
  float const DIS_CA_C   = 1.5232580f;
  float const DIS_C_N    = 1.3286850f;
  Real ang_n = /*3.1415926535897932384626433832795-*/(acos(COS_CA_C_N)+asin(SIN_CA_C_N))/2.0;
  Real angca = /*3.1415926535897932384626433832795-*/(acos(COS_C_N_CA)+asin(SIN_C_N_CA))/2.0;
  Real ang_c = /*3.1415926535897932384626433832795-*/(acos(COS_N_CA_C)+asin(SIN_N_CA_C))/2.0;  
  for(Size i = 1; i <= p.n_residue(); ++i) {
    ;        p.set_dof(DOF_ID(AtomID(1,i),core::id::THETA),ang_n);
    if(i!=1) p.set_dof(DOF_ID(AtomID(2,i),core::id::THETA),angca);
    ;        p.set_dof(DOF_ID(AtomID(3,i),core::id::THETA),ang_c);
    ;        p.set_dof(DOF_ID(AtomID(1,i),core::id::D),DIS_C_N);
    if(i!=1) p.set_dof(DOF_ID(AtomID(2,i),core::id::D),DIS_N_CA);
    ;        p.set_dof(DOF_ID(AtomID(3,i),core::id::D),DIS_CA_C);
  }
}

void gpu_refold_test(uint const NITER) {

  core::chemical::ResidueTypeSetCAP crs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
  core::chemical::ResidueTypeSetCAP frs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
  using basic::options::option;
  using namespace basic::options::OptionKeys;
  using core::pose::Pose;

  int err;
  struct VEC *CB_LCOR = new VEC[20];
  CB_LCOR[ 0u] = vec(0.5333336334f,-0.7751260775f,1.195999688f); // A
  CB_LCOR[ 1u] = vec(0.5385020267f,-0.571040073f,1.312000289f); // C
  CB_LCOR[ 2u] = vec(0.5256167631f,-0.7793338556f,1.207998884f); // D
  CB_LCOR[ 3u] = vec(0.525755529f,-0.7785836494f,1.208000359f); // E
  CB_LCOR[ 4u] = vec(0.5262055372f,-0.7772538616f,1.207999598f); // F
  CB_LCOR[ 5u] = vec(0.0f,0.0f,0.0f); // G
  CB_LCOR[ 6u] = vec(0.5234921973f,-0.7836302636f,1.207997917f); // H
  CB_LCOR[ 7u] = vec(0.5041208954f,-0.8015844639f,1.213999277f); // I
  CB_LCOR[ 8u] = vec(0.5258157482f,-0.7785053493f,1.206997076f); // K
  CB_LCOR[ 9u] = vec(0.5336851968f,-0.7754772515f,1.210999657f); // L
  CB_LCOR[10u] = vec(0.5284624554f,-0.7709555162f,1.208000858f); // M
  CB_LCOR[11u] = vec(0.5246473761f,-0.7987818032f,1.178999268f); // N
  CB_LCOR[12u] = vec(0.5787360676f,-0.6275511491f,1.279085868f); // P
  CB_LCOR[13u] = vec(0.5166017522f,-0.7752828156f,1.214998881f); // Q
  CB_LCOR[14u] = vec(0.5569299709f,-0.8316861918f,1.145999356f); // R
  CB_LCOR[15u] = vec(0.5281901647f,-0.7647611169f,1.198002116f); // S
  CB_LCOR[16u] = vec(0.5029925471f,-0.801331108f,1.215000707f); // T
  CB_LCOR[17u] = vec(0.5035780983f,-0.801583241f,1.2150076f); // V
  CB_LCOR[18u] = vec(0.5244585036f,-0.7753769855f,1.209996785f); // W
  CB_LCOR[19u] = vec(0.5248130885f,-0.7777290247f,1.209000815f); // Y
  struct VEC *CEN_LCOR = new VEC[20];
  CEN_LCOR[ 0u] = vec(0.5324949443f,-0.7748415953f,1.195342832f); // A
  CEN_LCOR[ 1u] = vec(1.265400934f,-1.276711201f,1.47384148f); // C
  CEN_LCOR[ 2u] = vec(1.256405064f,-1.444627167f,1.454677027f); // D
  CEN_LCOR[ 3u] = vec(1.624075779f,-1.893351417f,1.882064572f); // E
  CEN_LCOR[ 4u] = vec(1.675429398f,-1.824215788f,1.539187475f); // F
  CEN_LCOR[ 5u] = vec(-0.000378412968f,-0.0001466455924f,-0.0003408792981f); // G
  CEN_LCOR[ 6u] = vec(1.581532088f,-1.68863299f,1.511134449f); // H
  CEN_LCOR[ 7u] = vec(1.314107304f,-1.420581615f,1.600644937f); // I
  CEN_LCOR[ 8u] = vec(1.847861549f,-2.235910231f,1.984348918f); // K
  CEN_LCOR[ 9u] = vec(1.616501126f,-1.860619718f,1.363157218f); // L
  CEN_LCOR[10u] = vec(1.618967482f,-2.066666501f,1.666680953f); // M
  CEN_LCOR[11u] = vec(1.322834211f,-1.446019727f,1.391972756f); // N
  CEN_LCOR[12u] = vec(1.280126925f,0.9946088254f,1.528740065f); // P
  CEN_LCOR[13u] = vec(1.690375247f,-1.948444963f,1.756352754f); // Q
  CEN_LCOR[14u] = vec(1.986485326f,-2.417875528f,2.408043756f); // R
  CEN_LCOR[15u] = vec(0.6995892515f,-0.8583423264f,1.705672815f); // S
  CEN_LCOR[16u] = vec(0.9334680122f,-1.029926495f,1.537327961f); // T
  CEN_LCOR[17u] = vec(0.9210670044f,-1.31522194f,1.406024671f); // V
  CEN_LCOR[18u] = vec(1.497323704f,-2.166713044f,1.648385057f); // W
  CEN_LCOR[19u] = vec(1.752540328f,-1.956289875f,1.624625239f); // Y
  // {
  //   Vec dummy(0,0,0);
  //   Pose tmp;
  //   cout << std::setprecision(10);
  //   for(Size i = 1; i <= 20; ++i) {
  //     make_pose_from_sequence(tmp,"A"+str(core::chemical::oneletter_code_from_aa((core::chemical::AA)i))+"A",*crs,false);
  //     core::kinematics::Stub s(tmp.xyz(AtomID(2,2)),tmp.xyz(AtomID(3,2)),tmp.xyz(AtomID(1,2)));
  //     Vec v = tmp.residue(2).has("CB") ? s.global2local(tmp.residue(2).xyz("CB")) : dummy;
  //     cout << " CB_LCOR["<<lzs(i-1,2)<<"u] = vec(" << v.x() << "f," << v.y() << "f," << v.z() << "f); // " << core::chemical::oneletter_code_from_aa((core::chemical::AA)i) << endl;
  //   }
  //   for(Size i = 1; i <= 20; ++i) {
  //     make_pose_from_sequence(tmp,"A"+str(core::chemical::oneletter_code_from_aa((core::chemical::AA)i))+"A",*crs,false);
  //     core::kinematics::Stub s(tmp.xyz(AtomID(2,2)),tmp.xyz(AtomID(3,2)),tmp.xyz(AtomID(1,2)));
  //     Vec v = tmp.residue(2).has("CEN") ? s.global2local(tmp.residue(2).xyz("CEN")) : dummy;
  //     cout << " CEN_LCOR["<<lzs(i-1,2)<<"u] = vec(" << v.x() << "f," << v.y() << "f," << v.z() << "f); // " << core::chemical::oneletter_code_from_aa((core::chemical::AA)i) << endl;
  //   }
  // }
  
  CL cl(TR);

  Pose p;
  core::import_pose::pose_from_pdb(p,option[in::file::s]()[1],crs);
  for(Size i = 1; i <= p.n_residue(); ++i) {
    if(p.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(p,i);
    if(p.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(p,i);
  }
  uint N = p.n_residue();
  vector1<xyzVector<Real> > nat_crd(3*N);
  for(Size i=1;i<=N;++i){nat_crd[3*i-2]=p.xyz(AtomID(1,i));nat_crd[3*i-1]=p.xyz(AtomID(2,i));nat_crd[3*i-0]=p.xyz(AtomID(3,i));}
  p.dump_pdb("refold_natv.pdb");
  Pose nat(p);
  uint * aas = new uint[N];
  for(Size i = 1; i <= N; ++i) aas[i-1u] = nat.residue(i).aa();

  float *tor = new float[3*N];
  float *degrees_tor = new float[3*N];
  for(uint i = 0; i < N; ++i) {
    tor[3u*i+0u] = dihedral_radians(safe_xyz(p,3,i  ),safe_xyz(p,1,i+1),safe_xyz(p,2,i+1),safe_xyz(p,3,i+1));
    tor[3u*i+1u] = dihedral_radians(safe_xyz(p,1,i+1),safe_xyz(p,2,i+1),safe_xyz(p,3,i+1),safe_xyz(p,1,i+2));
    tor[3u*i+2u] = dihedral_radians(safe_xyz(p,2,i+1),safe_xyz(p,3,i+1),safe_xyz(p,1,i+2),safe_xyz(p,2,i+2));
    degrees_tor[3u*i+0u] = dihedral_degrees(safe_xyz(p,3,i  ),safe_xyz(p,1,i+1),safe_xyz(p,2,i+1),safe_xyz(p,3,i+1));
    degrees_tor[3u*i+1u] = dihedral_degrees(safe_xyz(p,1,i+1),safe_xyz(p,2,i+1),safe_xyz(p,3,i+1),safe_xyz(p,1,i+2));
    degrees_tor[3u*i+2u] = dihedral_degrees(safe_xyz(p,2,i+1),safe_xyz(p,3,i+1),safe_xyz(p,1,i+2),safe_xyz(p,2,i+2));
  }

  float *N__xyz = new float[3*N];
  float *CA_xyz = new float[3*N];
  float *C__xyz = new float[3*N];
  float *O__xyz = new float[3*N];
  float *CB_xyz = new float[3*N];
  float *H__xyz = new float[3*N];
  float *CENxyz = new float[3*N];

  double tpar = 0.0;
  for(uint iter = 0; iter < 1*NITER; ++iter) {
    double t = time_highres();
    ;                               for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_first (tor,&N,  N__xyz,CA_xyz,C__xyz,O__xyz,CB_xyz,H__xyz,CENxyz); }
    for(uint c=2u;c<N*2u-3u;c=2u*c) for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_second(tor,&N,c,N__xyz,CA_xyz,C__xyz,O__xyz,CB_xyz,H__xyz,CENxyz); }
    ;                               for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_third (tor,&N,  N__xyz,CA_xyz,C__xyz,O__xyz,CB_xyz,H__xyz,CENxyz,CB_LCOR,CEN_LCOR,aas); }
    tpar += time_highres() - t;
  }

  Pose tmp;
  string seq = "";
  for(int i = 0; i < N; ++i) seq += "A";
  core::pose::make_pose_from_sequence(tmp,seq,*crs,false);
  for(Size i = 1; i <= tmp.n_residue(); ++i) {
    if(tmp.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(tmp,i);
    if(tmp.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(tmp,i);
  }
  for(uint i = 0u; i < N; ++i) {
    tmp.set_xyz(AtomID(1u,i+1u),Vec(N__xyz[3u*i+0u],N__xyz[3u*i+1u],N__xyz[3u*i+2u]));
    tmp.set_xyz(AtomID(2u,i+1u),Vec(CA_xyz[3u*i+0u],CA_xyz[3u*i+1u],CA_xyz[3u*i+2u]));
    tmp.set_xyz(AtomID(3u,i+1u),Vec(C__xyz[3u*i+0u],C__xyz[3u*i+1u],C__xyz[3u*i+2u]));
    tmp.set_xyz(AtomID(4u,i+1u),Vec(O__xyz[3u*i+0u],O__xyz[3u*i+1u],O__xyz[3u*i+2u]));
    tmp.set_xyz(AtomID(5u,i+1u),Vec(CB_xyz[3u*i+0u],CB_xyz[3u*i+1u],CB_xyz[3u*i+2u]));
    tmp.set_xyz(AtomID(6u,i+1u),Vec(CENxyz[3u*i+0u],CENxyz[3u*i+1u],CENxyz[3u*i+2u]));
    tmp.set_xyz(AtomID(7u,i+1u),Vec(H__xyz[3u*i+0u],H__xyz[3u*i+1u],H__xyz[3u*i+2u]));
  }
  cout << "GEOM H " << tmp.xyz(AtomID(1,2)).distance(tmp.xyz(AtomID(7,2))) << " " 
       << angle_degrees(tmp.xyz(AtomID(3,1)),tmp.xyz(AtomID(1,2)),tmp.xyz(AtomID(7,2))) << " "
       << angle_degrees(tmp.xyz(AtomID(2,2)),tmp.xyz(AtomID(1,2)),tmp.xyz(AtomID(7,2))) << endl;       
  cout << "GEOM O " << tmp.xyz(AtomID(4,2)).distance(tmp.xyz(AtomID(3,2))) << " " 
       << angle_degrees(tmp.xyz(AtomID(2,2)),tmp.xyz(AtomID(3,2)),tmp.xyz(AtomID(4,2))) << " "
       << angle_degrees(tmp.xyz(AtomID(1,3)),tmp.xyz(AtomID(3,2)),tmp.xyz(AtomID(4,2))) << endl;       
  vector1<xyzVector<Real> > cpu_crd(3u*N);
  for(Size i=1;i<=N;++i){cpu_crd[3*i-2]=tmp.xyz(AtomID(1,i));cpu_crd[3*i-1]=tmp.xyz(AtomID(2,i));cpu_crd[3*i-0]=tmp.xyz(AtomID(3,i));}
  core::scoring::calpha_superimpose_pose(tmp,nat);
  tmp.dump_pdb("refold_cpu.pdb");
  
  utility_exit_with_message("airst");
  
  vector1<xyzVector<Real> > ros_crd(3*N), idl_crd(3*N);
  double tros = 0.0;
  {
    Pose ideal_bla(nat);
    set_pose_to_ideal(ideal_bla);
    for(Size i = 1; i <= N; ++i) {
      idl_crd[3*i-2] = ideal_bla.xyz(AtomID(1,i));
      idl_crd[3*i-1] = ideal_bla.xyz(AtomID(2,i));
      idl_crd[3*i-0] = ideal_bla.xyz(AtomID(3,i));
    }
    Vec t1,t2,t3,t4;
    for(uint iter = 0; iter < 1*NITER; ++iter) {
      for(Size i = 1; i <= N; ++i) {
        ideal_bla.set_phi  (i,0.0);
        ideal_bla.set_psi  (i,0.0);
        ideal_bla.set_omega(i,0.0);
      }
      t1 += ideal_bla.xyz(AtomID(1,1));
      t2 += ideal_bla.xyz(AtomID(3,N));
      double t = time_highres();
      for(Size i = 1; i <= N; ++i) {
        ideal_bla.set_phi  (i,degrees_tor[3u*i-3u]);
        ideal_bla.set_psi  (i,degrees_tor[3u*i-2u]);
        ideal_bla.set_omega(i,degrees_tor[3u*i-1u]);
      }
      t3 += ideal_bla.xyz(AtomID(1,1));
      t4 += ideal_bla.xyz(AtomID(3,N));
      tros += time_highres() - t;
    }
    for(Size i = 1; i <= N; ++i) {
      ros_crd[3u*i-2u] = ideal_bla.xyz(AtomID(1u,i));
      ros_crd[3u*i-1u] = ideal_bla.xyz(AtomID(2u,i));
      ros_crd[3u*i-0u] = ideal_bla.xyz(AtomID(3u,i));
    }
    core::scoring::calpha_superimpose_pose(ideal_bla,nat);
    ideal_bla.dump_pdb("refold_rosetta.pdb");
  }


  string Kname = "abinitio";
  cl.make_kernel(Kname);
  vector1<cl_mem> outs;
  cl_mem clmtor = cl.makeROmem(sizeof(cl_float)*3u*N);
  cl_mem clmN   = cl.makeROmem(sizeof(cl_uint)      );
  cl.cpu2gpu( &N , clmN      , sizeof(cl_float)     );
  cl_mem clmout = cl.makeWOmem(sizeof(cl_float)*100u*21u*N);
  cl.setargs(Kname,clmtor,clmN,clmout);
  float *out = new float[21u*N*100u];

  double tcl = time_highres();
  cl.cpu2gpu( tor , clmtor   , sizeof(cl_float)*3u*N);
  for(Size iter = 0; iter < NITER; ++iter){
    cl.q2d(Kname,256u,100u,256u,1u);
  }
  cl.finish();
  tcl = time_highres() - tcl;
  for(Size i = 0u; i < 21u*N; ++i) out[i] = 0.0f;
  cl.gpu2cpu( clmout, out, sizeof(cl_float)*100u*N*21u);

  vector1<xyzVector<Real> > gpu_crd(3*N); {
    Pose tmp;
    string seq = "";
    for(int i = 0; i < N; ++i) seq += "A";
    core::pose::make_pose_from_sequence(tmp,seq,*crs,false);
    for(Size i = 1; i <= tmp.n_residue(); ++i) {
      if(tmp.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(tmp,i);
      if(tmp.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(tmp,i);
    }
    for(int i = 0; i < N; ++i) {
      tmp.set_xyz(AtomID(1,i+1),Vec(out[21*i+ 0],out[21*i+ 1],out[21*i+ 2]));
      tmp.set_xyz(AtomID(2,i+1),Vec(out[21*i+ 3],out[21*i+ 4],out[21*i+ 5]));
      tmp.set_xyz(AtomID(3,i+1),Vec(out[21*i+ 6],out[21*i+ 7],out[21*i+ 8]));
      tmp.set_xyz(AtomID(4,i+1),Vec(out[21*i+ 9],out[21*i+10],out[21*i+11]));
      tmp.set_xyz(AtomID(5,i+1),Vec(out[21*i+12],out[21*i+13],out[21*i+14]));
      tmp.set_xyz(AtomID(6,i+1),Vec(out[21*i+15],out[21*i+16],out[21*i+17]));
      tmp.set_xyz(AtomID(7,i+1),Vec(out[21*i+18],out[21*i+19],out[21*i+20]));
    }
    core::scoring::calpha_superimpose_pose(tmp,nat);
    tmp.dump_pdb("refold_gpu.pdb");
    for(Size i = 1; i <= N; ++i) {
      gpu_crd[3u*i-2u] = tmp.xyz(AtomID(1u,i));
      gpu_crd[3u*i-1u] = tmp.xyz(AtomID(2u,i));
      gpu_crd[3u*i-0u] = tmp.xyz(AtomID(3u,i));
    }
  }
  TR << "CA RMSD: ros/gpu " << F(11,6,numeric::model_quality::calc_rms(ros_crd,gpu_crd)) 
     <<       "   cpu/gpu " << F(11,6,numeric::model_quality::calc_rms(cpu_crd,gpu_crd)) 
     <<       "   nat/gpu " << F(11,6,numeric::model_quality::calc_rms(nat_crd,gpu_crd)) << endl;
  TR << "RUNTIME: ros/gpu " << F(10,5,100.0*tros/tcl) << "x   cpu/gpu " << F(10,5,100.0*tpar/tcl) << "x" << endl;
  TR << endl;
  //TR << "tros: " << tros << "                " << t3.x() << t4.x() << std::endl;
  TR << "IDL/CPU RMS: " << numeric::model_quality::calc_rms(idl_crd,cpu_crd) << endl;
  TR << "ROS/CPU RMS: " << numeric::model_quality::calc_rms(cpu_crd,ros_crd) << endl;
  TR << "ROS/NAT RMS: " << numeric::model_quality::calc_rms(nat_crd,ros_crd) << endl;
  TR << "IDL/NAT RMS: " << numeric::model_quality::calc_rms(nat_crd,idl_crd) << endl;  
  TR << "CPU/NAT RMS: " << numeric::model_quality::calc_rms(nat_crd,cpu_crd) << endl;
  // TR << "ROSFLD/GPU RMS: " << numeric::model_quality::calc_rms(ros_crd,gpu_crd) << endl;  


}
