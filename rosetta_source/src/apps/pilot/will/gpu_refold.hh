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
    ;                               for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_third (tor,&N,  N__xyz,CA_xyz,C__xyz,O__xyz,CB_xyz,H__xyz,CENxyz); }
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
