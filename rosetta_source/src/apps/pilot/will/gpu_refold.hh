#include <apps/pilot/will/gpu_mat_vec.hh>
#include <apps/pilot/will/gpu_refold_test_cpu.hh>

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

  CL cl(TR);

  Pose p;
  core::import_pose::pose_from_pdb(p,option[in::file::s]()[1]);
  for(Size i = 1; i <= p.n_residue(); ++i) {
    if(p.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(p,i);
    if(p.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(p,i);
  }
  uint N = p.n_residue();
  vector1<xyzVector<Real> > nat_crd(3*N);
  for(Size i=1;i<=N;++i){nat_crd[3*i-2]=p.xyz(AtomID(1,i));nat_crd[3*i-1]=p.xyz(AtomID(2,i));nat_crd[3*i-0]=p.xyz(AtomID(3,i));}
  p.dump_pdb("refold_natv.pdb");
  Pose nat(p);

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

  float *bb     = new float[9*N];
  float *b2     = new float[9*N];
  float *out    = new float[1000*9*N];
  for(uint i = 0u; i < 9*N; ++i) {bb[i]=0.0f; b2[i]=0.0f; out[i]=0.0f; }

  double tpar = 0.0;
  for(uint iter = 0; iter < 10*NITER; ++iter) {
    double t = time_highres();
    ;                               for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_first (tor,&N,out,  bb,b2); }
    for(uint c=2u;c<N*2u-3u;c=2u*c) for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_second(tor,&N,out,c,bb,b2); }
    ;                               for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_third (tor,&N,out,  bb,b2); }
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
  for(int i = 0; i < N; ++i) {
    tmp.set_xyz(AtomID(1,i+1),Vec(b2[9*i+0],b2[9*i+1],b2[9*i+2]));
    tmp.set_xyz(AtomID(2,i+1),Vec(b2[9*i+3],b2[9*i+4],b2[9*i+5]));
    tmp.set_xyz(AtomID(3,i+1),Vec(b2[9*i+6],b2[9*i+7],b2[9*i+8]));
  }
  vector1<xyzVector<Real> > cpu_crd(3*N);
  for(Size i=1;i<=N;++i){cpu_crd[3*i-2]=tmp.xyz(AtomID(1,i));cpu_crd[3*i-1]=tmp.xyz(AtomID(2,i));cpu_crd[3*i-0]=tmp.xyz(AtomID(3,i));}
  core::scoring::calpha_superimpose_pose(tmp,nat);
  tmp.dump_pdb("refold_cpu.pdb");


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
    for(uint iter = 0; iter < 10*NITER; ++iter) {
      for(Size i = 1; i <= N; ++i) {
        ideal_bla.set_phi  (i,0.0);
        ideal_bla.set_psi  (i,0.0);
        ideal_bla.set_omega(i,0.0);
      }
      t1 += ideal_bla.xyz(AtomID(1,1));
      t2 += ideal_bla.xyz(AtomID(3,N));
      double t = time_highres();
      for(Size i = 1; i <= N; ++i) {
        ideal_bla.set_phi  (i,degrees_tor[3*i-3]);
        ideal_bla.set_psi  (i,degrees_tor[3*i-2]);
        ideal_bla.set_omega(i,degrees_tor[3*i-1]);
      }
      t3 += ideal_bla.xyz(AtomID(1,1));
      t4 += ideal_bla.xyz(AtomID(3,N));
      tros += time_highres() - t;
    }
    for(Size i = 1; i <= N; ++i) {
      ros_crd[3*i-2] = ideal_bla.xyz(AtomID(1,i));
      ros_crd[3*i-1] = ideal_bla.xyz(AtomID(2,i));
      ros_crd[3*i-0] = ideal_bla.xyz(AtomID(3,i));
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
  cl_mem clmout = cl.makeWOmem(sizeof(cl_float)*1000*9u*N);
  cl.setargs(Kname,clmtor,clmN,clmout);
  float *result = new float[9u*N];

  double tcl = time_highres();
  cl.cpu2gpu( tor , clmtor   , sizeof(cl_float)*3u*N);
  for(Size iter = 0; iter < NITER; ++iter){
    cl.q2d(Kname,256,1000,256,1);
  }
  cl.finish();
  tcl = time_highres() - tcl;
  for(Size i = 0; i < 9*N; ++i) out[i] = 0.0f;
  cl.gpu2cpu( clmout, out, sizeof(cl_float)*1000*N*9u);

  vector1<xyzVector<Real> > gpu_crd(3*N); {
    Pose tmp;
    string seq = "";
    for(int i = 0; i < N; ++i) seq += "G";
    core::pose::make_pose_from_sequence(tmp,seq,*crs,false);
    for(Size i = 1; i <= tmp.n_residue(); ++i) {
      if(tmp.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(tmp,i);
      if(tmp.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(tmp,i);
    }
    for(int i = 0; i < N; ++i) {
      tmp.set_xyz(AtomID(1,i+1),Vec(out[9*i+0],out[9*i+1],out[9*i+2]));
      tmp.set_xyz(AtomID(2,i+1),Vec(out[9*i+3],out[9*i+4],out[9*i+5]));
      tmp.set_xyz(AtomID(3,i+1),Vec(out[9*i+6],out[9*i+7],out[9*i+8]));
    }
    core::scoring::calpha_superimpose_pose(tmp,nat);
    tmp.dump_pdb("refold_gpu.pdb");
    for(Size i = 1; i <= N; ++i) {
      gpu_crd[3*i-2] = tmp.xyz(AtomID(1,i));
      gpu_crd[3*i-1] = tmp.xyz(AtomID(2,i));
      gpu_crd[3*i-0] = tmp.xyz(AtomID(3,i));
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
