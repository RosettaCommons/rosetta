#include <apps/pilot/will/gpu/gpu_mat_vec.hh>
#include <apps/pilot/will/gpu/gpu_refold_test_cpu.hh>
#include <numeric/xyz.io.hh>

#include <apps/pilot/will/gpu/set_pose_to_ideal.ihh>

void gpu_refold_test(uint const NITER) {
  uint LOCAL_DIM = option[basic::options::OptionKeys::gpu::numthreads_per_workunit]();

  TR << "gpu_refold_test" << endl;

  core::chemical::ResidueTypeSetCAP crs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
  core::chemical::ResidueTypeSetCAP frs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
  using basic::options::option;
  using namespace basic::options::OptionKeys;
  using core::pose::Pose;

  int err;
  struct VEC *CB_LCOR = new VEC[20];
  CB_LCOR[ 0u] = vec(0.5298019858f,-0.7775442472f,-1.195999688f); // A
  CB_LCOR[ 1u] = vec(0.3376587247f,-0.7085603598f,-1.312000289f); // C
  CB_LCOR[ 2u] = vec(0.5365158091f,-0.7718711205f,-1.207998884f); // D
  CB_LCOR[ 3u] = vec(0.5357659908f,-0.7717293427f,-1.208000359f); // E
  CB_LCOR[ 4u] = vec(0.5343635718f,-0.7716679376f,-1.207999598f); // F
  CB_LCOR[ 5u] = vec(0.0f,0.0f,0.0f); // G
  CB_LCOR[ 6u] = vec(0.5412897477f,-0.7714440224f,-1.207997917f); // H
  CB_LCOR[ 7u] = vec(0.565034011f,-0.7598763691f,-1.213999277f); // I
  CB_LCOR[ 8u] = vec(0.5356713209f,-0.7717570965f,-1.206997076f); // K
  CB_LCOR[ 9u] = vec(0.5300015401f,-0.7779995015f,-1.210999657f); // L
  CB_LCOR[10u] = vec(0.5276755193f,-0.7714943428f,-1.208000858f); // M
  CB_LCOR[11u] = vec(0.5549986657f,-0.7779998198f,-1.178999268f); // N
  CB_LCOR[12u] = vec(0.3757952229f,-0.7665075545f,-1.279085868f); // P
  CB_LCOR[13u] = vec(0.5359990809f,-0.762001181f,-1.214998881f); // Q
  CB_LCOR[14u] = vec(0.5740011989f,-0.8199972791f,-1.145999356f); // R
  CB_LCOR[15u] = vec(0.5219984905f,-0.7690006449f,-1.198002116f); // S
  CB_LCOR[16u] = vec(0.5652062499f,-0.7587324576f,-1.215000707f); // T
  CB_LCOR[17u] = vec(0.5652287619f,-0.7593701601f,-1.2150076f); // V
  CB_LCOR[18u] = vec(0.5332451765f,-0.7693606264f,-1.209996785f); // W
  CB_LCOR[19u] = vec(0.5353099059f,-0.7705417045f,-1.209000815f); // Y
  struct VEC *CEN_LCOR = new VEC[20];
  CEN_LCOR[ 0u] = vec(0.5298400468f,-0.7766594416f,-1.195342832f); // A
  CEN_LCOR[ 1u] = vec(0.7327082569f,-1.641453509f,-1.47384148f); // C
  CEN_LCOR[ 2u] = vec(0.8925138823f,-1.693788743f,-1.454677027f); // D
  CEN_LCOR[ 3u] = vec(1.177911f,-2.198846835f,-1.882064572f); // E
  CEN_LCOR[ 4u] = vec(1.094883789f,-2.221723745f,-1.539187475f); // F
  CEN_LCOR[ 5u] = vec(0.0002735646292f,0.0002997727436f,0.0003408792981f); // G
  CEN_LCOR[ 6u] = vec(1.002432522f,-2.085150872f,-1.511134449f); // H
  CEN_LCOR[ 7u] = vec(0.8492291086f,-1.738890466f,-1.600644937f); // I
  CEN_LCOR[ 8u] = vec(1.416360674f,-2.531365108f,-1.984348918f); // K
  CEN_LCOR[ 9u] = vec(1.150131621f,-2.179949284f,-1.363157218f); // L
  CEN_LCOR[10u] = vec(1.341344973f,-2.256758693f,-1.666680953f); // M
  CEN_LCOR[11u] = vec(0.8697909727f,-1.756225232f,-1.391972756f); // N
  CEN_LCOR[12u] = vec(-1.390223339f,-0.8338169625f,-1.528740065f); // P
  CEN_LCOR[13u] = vec(1.205301427f,-2.280582101f,-1.756352754f); // Q
  CEN_LCOR[14u] = vec(1.535880763f,-2.726410883f,-2.408043756f); // R
  CEN_LCOR[15u] = vec(0.5472643686f,-0.9626413564f,-1.705672815f); // S
  CEN_LCOR[16u] = vec(0.6226610786f,-1.242740639f,-1.537327961f); // T
  CEN_LCOR[17u] = vec(0.8931319791f,-1.334349446f,-1.406024671f); // V
  CEN_LCOR[18u] = vec(1.478608645f,-2.17952751f,-1.648385057f); // W
  CEN_LCOR[19u] = vec(1.190133765f,-2.341377649f,-1.624625239f); // Y

  uint const NSCPRM = 1000u;
  float CL_SCPRM[NSCPRM];
  for(Size i = 0; i < NSCPRM; ++i) CL_SCPRM[i] = 0.0f;

  // {
  //   Vec dummy(0,0,0);
  //   Pose tmp;
  //   cout << std::setprecision(10);
  //   for(Size i = 1; i <= 20; ++i) {
  //     make_pose_from_sequence(tmp,"A"+str(core::chemical::oneletter_code_from_aa((core::chemical::AA)i))+"A",*crs,false);
  //     core::kinematics::Stub s(tmp.xyz(AtomID(2,2)),tmp.xyz(AtomID(1,2)),tmp.xyz(AtomID(3,2)));
  //     Vec v = tmp.residue(2).has("CB") ? s.global2local(tmp.residue(2).xyz("CB")) : dummy;
  //     cout << " CB_LCOR["<<lzs(i-1,2)<<"u] = vec(" << v.x() << "f," << v.y() << "f," << v.z() << "f); // " << core::chemical::oneletter_code_from_aa((core::chemical::AA)i) << endl;
  //   }
  //   for(Size i = 1; i <= 20; ++i) {
  //     make_pose_from_sequence(tmp,"A"+str(core::chemical::oneletter_code_from_aa((core::chemical::AA)i))+"A",*crs,false);
  //     core::kinematics::Stub s(tmp.xyz(AtomID(2,2)),tmp.xyz(AtomID(1,2)),tmp.xyz(AtomID(3,2)));
  //     Vec v = tmp.residue(2).has("CEN") ? s.global2local(tmp.residue(2).xyz("CEN")) : dummy;
  //     cout << " CEN_LCOR["<<lzs(i-1,2)<<"u] = vec(" << v.x() << "f," << v.y() << "f," << v.z() << "f); // " << core::chemical::oneletter_code_from_aa((core::chemical::AA)i) << endl;
  //   }
  // }
  
  TR << "init CL" << endl;
  CL cl(TR);

  TR << "setup poses" << endl;
  Pose p( *core::import_pose::pose_from_file(*crs,option[in::file::s]()[1],false) , core::import_pose::PDB_file); 
  for(Size i = 1; i <= p.n_residue(); ++i) {
    if(p.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(p,i);
    if(p.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(p,i);
  }
  uint N = p.n_residue();
  vector1<xyzVector<Real> > nat_crd(3*N);
  for(Size i=1;i<=N;++i){nat_crd[3*i-2]=p.xyz(AtomID(1,i));nat_crd[3*i-1]=p.xyz(AtomID(2,i));nat_crd[3*i-0]=p.xyz(AtomID(3,i));}
  // p.dump_pdb("refold_natv.pdb");
  Pose nat(p);
  uint * aas = new uint[N];
  for(Size i = 1; i <= N; ++i) aas[i-1u] = (uint)nat.residue(i).aa()-1;

  TR << "get torsions" << endl;
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

  TR << "CPU start" << endl;
  double tpar = time_highres();
  for(uint iter = 0; iter < 10*NITER; ++iter) {
    ;                               for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_first (tor,&N,  N__xyz,CA_xyz,C__xyz,O__xyz,CB_xyz,H__xyz,CENxyz); }
    for(uint c=2u;c<N*2u-3u;c=2u*c) for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_second(tor,&N,c,N__xyz,CA_xyz,C__xyz,O__xyz,CB_xyz,H__xyz,CENxyz); }
    ;                               for(uint gi = 0; gi < N; ++gi) { gid_ = gi; refold_third (tor,&N,  N__xyz,CA_xyz,C__xyz,O__xyz,CB_xyz,H__xyz,CENxyz,CB_LCOR,CEN_LCOR,aas); }
  }
  tpar = time_highres() - tpar;
  TR << "CPU stop: " << tpar << endl;

  Pose tmp;
  core::pose::make_pose_from_sequence(tmp,nat.sequence(),*crs,false);
  for(Size i = 1; i <= tmp.n_residue(); ++i) {
    if(tmp.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(tmp,i);
    if(tmp.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(tmp,i);
  }
  for(uint i = 0u; i < N; ++i) {
    if(tmp.residue(i+1u).has( "N" )) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "N" ),i+1u),Vec(N__xyz[3u*i+0u],N__xyz[3u*i+1u],N__xyz[3u*i+2u]));
    if(tmp.residue(i+1u).has( "CA")) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "CA"),i+1u),Vec(CA_xyz[3u*i+0u],CA_xyz[3u*i+1u],CA_xyz[3u*i+2u]));
    if(tmp.residue(i+1u).has( "C" )) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "C" ),i+1u),Vec(C__xyz[3u*i+0u],C__xyz[3u*i+1u],C__xyz[3u*i+2u]));
    if(tmp.residue(i+1u).has( "O" )) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "O" ),i+1u),Vec(O__xyz[3u*i+0u],O__xyz[3u*i+1u],O__xyz[3u*i+2u]));
    if(tmp.residue(i+1u).has( "CB")) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "CB"),i+1u),Vec(CB_xyz[3u*i+0u],CB_xyz[3u*i+1u],CB_xyz[3u*i+2u]));
    if(tmp.residue(i+1u).has("CEN")) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index("CEN"),i+1u),Vec(CENxyz[3u*i+0u],CENxyz[3u*i+1u],CENxyz[3u*i+2u]));
    if(tmp.residue(i+1u).has( "H" )) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "H" ),i+1u),Vec(H__xyz[3u*i+0u],H__xyz[3u*i+1u],H__xyz[3u*i+2u]));
  }
  // cout << "GEOM H " << tmp.xyz(AtomID(1,2)).distance(tmp.xyz(AtomID(7,2))) << " " 
  //      << angle_degrees(tmp.xyz(AtomID(3,1)),tmp.xyz(AtomID(1,2)),tmp.xyz(AtomID(7,2))) << " "
  //      << angle_degrees(tmp.xyz(AtomID(2,2)),tmp.xyz(AtomID(1,2)),tmp.xyz(AtomID(7,2))) << endl;       
  // cout << "GEOM O " << tmp.xyz(AtomID(4,2)).distance(tmp.xyz(AtomID(3,2))) << " " 
  //      << angle_degrees(tmp.xyz(AtomID(2,2)),tmp.xyz(AtomID(3,2)),tmp.xyz(AtomID(4,2))) << " "
  //      << angle_degrees(tmp.xyz(AtomID(1,3)),tmp.xyz(AtomID(3,2)),tmp.xyz(AtomID(4,2))) << endl;       
  vector1<xyzVector<Real> > cpu_crd(3u*N);
  for(Size i=1;i<=N;++i){cpu_crd[3*i-2]=tmp.xyz(AtomID(1,i));cpu_crd[3*i-1]=tmp.xyz(AtomID(2,i));cpu_crd[3*i-0]=tmp.xyz(AtomID(3,i));}
  core::scoring::calpha_superimpose_pose(tmp,nat);
  // tmp.dump_pdb("refold_cpu.pdb");

  ScoreFunctionOP sf = new ScoreFunction();
  sf->set_weight(core::scoring::rg     ,1.0);
  sf->set_weight(core::scoring::vdw    ,1.0);
  sf->set_weight(core::scoring::pair   ,1.0);
  sf->set_weight(core::scoring::env    ,1.0);
  sf->set_weight(core::scoring::cbeta  ,1.0);
  sf->set_weight(core::scoring::cenpack,1.0);

  vector1<xyzVector<Real> > ros_crd(3*N), idl_crd(3*N);
  double tros = 0.0;
  {
    Pose ideal_bla(nat);
    set_pose_to_ideal(ideal_bla);
    sf->show(ideal_bla);
    for(Size i = 1; i <= N; ++i) {
      idl_crd[3*i-2] = ideal_bla.xyz(AtomID(1,i));
      idl_crd[3*i-1] = ideal_bla.xyz(AtomID(2,i));
      idl_crd[3*i-0] = ideal_bla.xyz(AtomID(3,i));
    }
    Vec t1,t2,t3,t4;
    TR << "ROS start" << endl;
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
        ideal_bla.set_phi  (i,degrees_tor[3u*i-3u]);
        ideal_bla.set_psi  (i,degrees_tor[3u*i-2u]);
        ideal_bla.set_omega(i,degrees_tor[3u*i-1u]);
      }
      sf->score(ideal_bla);
      // t3 += ideal_bla.xyz(AtomID(1,1));
      // t4 += ideal_bla.xyz(AtomID(3,N));
      tros += time_highres() - t;
    }
    TR << "ROS stop: " << tros << endl;    
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
  cl_mem clmtor   = cl.makeRWmem(sizeof(cl_float)*3u*3u*N      );
  cl_mem clmN     = cl.makeROmem(sizeof(cl_uint)            );
  cl_mem clmout   = cl.makeWOmem(sizeof(cl_float)*100u*21u*N);
  // cl_mem clmtst1 = cl.makeWOmem(sizeof(cl_uint)*128*127);
  // cl_mem clmtst2 = cl.makeWOmem(sizeof(cl_uint)*128*127);  
  cl_mem clmsc    = cl.makeWOmem(sizeof(cl_float)*50u)      ; // 50 scores????
  cl_mem clmaas   = cl.makeROmem(sizeof(cl_uint)*N          );
  cl_mem clmCB    = cl.makeROmem(sizeof(struct VEC)*20u     );
  cl_mem clmCEN   = cl.makeROmem(sizeof(struct VEC)*20u     );
  cl_mem clmscprm = cl.makeROmem(sizeof(cl_float)*NSCPRM    );
  cl.cpu2gpu( tor     , clmtor, sizeof(cl_float)*3u*N      );
  cl.cpu2gpu( &N      , clmN  , sizeof(cl_float)           );
  //          out               sizeof(cl_float)*100u*21u*N);  
  cl.cpu2gpu( aas     , clmaas, sizeof(cl_uint)*N          );
  cl.cpu2gpu(  CB_LCOR, clmCB , sizeof(struct VEC)*20u     );  
  cl.cpu2gpu( CEN_LCOR, clmCEN, sizeof(struct VEC)*20u     );    
  cl.cpu2gpu( CL_SCPRM, clmscprm, sizeof(cl_float)*NSCPRM     );  

  err = clSetKernelArg(cl.kernels_[Kname],0,sizeof(cl_mem), &clmtor   ); if(err!=CL_SUCCESS) cout << "ERR " << 0 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],1,sizeof(cl_mem), &clmN     ); if(err!=CL_SUCCESS) cout << "ERR " << 1 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],2,sizeof(cl_mem), &clmout   ); if(err!=CL_SUCCESS) cout << "ERR " << 2 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],3,sizeof(cl_mem), &clmsc    ); if(err!=CL_SUCCESS) cout << "ERR " << 3 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],4,sizeof(cl_float)*21u*N, NULL ); if(err!=CL_SUCCESS) cout << "ERR " << 4 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],5,sizeof(cl_float)*LOCAL_DIM, NULL ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],6,sizeof(cl_mem), &clmaas   ); if(err!=CL_SUCCESS) cout << "ERR " << 6 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],7,sizeof(cl_mem), &clmCB    ); if(err!=CL_SUCCESS) cout << "ERR " << 7 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],8,sizeof(cl_mem), &clmCEN   ); if(err!=CL_SUCCESS) cout << "ERR " << 8 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname],9,sizeof(cl_mem), &clmscprm ); if(err!=CL_SUCCESS) cout << "ERR " << 9 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname],8,sizeof(cl_mem), &clmtst1 ); if(err!=CL_SUCCESS) cout << "ERR " << 8 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname],9,sizeof(cl_mem), &clmtst2 ); if(err!=CL_SUCCESS) cout << "ERR " << 9 << " " << cl.errstr(err) << endl;  

  // cl_mem clmtorsions = cl.makeRWmem(sizeof(cl_float)*3u*N);
  // cl_mem clmN__xyz   = cl.makeRWmem(sizeof(cl_float)*3u*N);
  // cl_mem clmCA_xyz   = cl.makeRWmem(sizeof(cl_float)*3u*N);
  // cl_mem clmC__xyz   = cl.makeRWmem(sizeof(cl_float)*3u*N);
  // cl_mem clmO__xyz   = cl.makeRWmem(sizeof(cl_float)*3u*N);
  // cl_mem clmCB_xyz   = cl.makeRWmem(sizeof(cl_float)*3u*N);
  // cl_mem clmH__xyz   = cl.makeRWmem(sizeof(cl_float)*3u*N);
  // cl_mem clmCENxyz   = cl.makeRWmem(sizeof(cl_float)*3u*N);
  // err = clSetKernelArg(cl.kernels_[Kname], 6,sizeof(cl_mem), &clmtorsions ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname], 7,sizeof(cl_mem), &clmN__xyz   ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname], 8,sizeof(cl_mem), &clmCA_xyz   ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname], 9,sizeof(cl_mem), &clmC__xyz   ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname],10,sizeof(cl_mem), &clmO__xyz   ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname],11,sizeof(cl_mem), &clmCB_xyz   ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname],12,sizeof(cl_mem), &clmH__xyz   ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;
  // err = clSetKernelArg(cl.kernels_[Kname],13,sizeof(cl_mem), &clmCENxyz   ); if(err!=CL_SUCCESS) cout << "ERR " << 5 << " " << cl.errstr(err) << endl;

  float *out = new float[21u*N*100u];
  float *sc  = new float[50u];

  TR << "GPU start" << endl;
  double tcl = time_highres();
  for(Size iter = 0; iter < NITER; ++iter){
    cl.q2d(Kname,LOCAL_DIM,100u,LOCAL_DIM,1u);
  }
  cl.finish();
  tcl = time_highres() - tcl;
  TR << "GPU stop: " << tcl << endl;
  for(Size i = 0u; i < 21u*N; ++i) out[i] = 0.0f;
  cl.gpu2cpu( clmout, out, sizeof(cl_float)*100u*N*21u);
  cl.gpu2cpu( clmsc , sc , sizeof(cl_float)*50u);

  // uint test1[128*127];
  // uint test2[128*127];  
  // cl.gpu2cpu( clmtst1 , test1 , sizeof(cl_float)*128*127);
  // cl.gpu2cpu( clmtst2 , test2 , sizeof(cl_float)*128*127);
  // for(uint i = 0; i <= 127*128; ++i) {
  //   cout << test1[i] << " " << test2[i] << " " << i << endl;

  // }

  TR << "GPU VDW:  " << sc[1] << endl;
  TR << "GPU RG :  " << sc[2] << endl;

  vector1<xyzVector<Real> > gpu_crd(3*N); {
    Pose tmp;
    core::pose::make_pose_from_sequence(tmp,nat.sequence(),*crs,false);
    for(Size i = 1; i <= tmp.n_residue(); ++i) {
      if(tmp.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(tmp,i);
      if(tmp.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(tmp,i);
    }
    for(int i = 0; i < N; ++i) {
      if(tmp.residue(i+1u).has( "N" )) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "N" ),i+1u),Vec(out[21*i+ 0],out[21*i+ 1],out[21*i+ 2]));
      if(tmp.residue(i+1u).has( "CA")) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "CA"),i+1u),Vec(out[21*i+ 3],out[21*i+ 4],out[21*i+ 5]));
      if(tmp.residue(i+1u).has( "C" )) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "C" ),i+1u),Vec(out[21*i+ 6],out[21*i+ 7],out[21*i+ 8]));
      if(tmp.residue(i+1u).has( "O" )) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "O" ),i+1u),Vec(out[21*i+ 9],out[21*i+10],out[21*i+11]));
      if(tmp.residue(i+1u).has( "CB")) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "CB"),i+1u),Vec(out[21*i+12],out[21*i+13],out[21*i+14]));
      if(tmp.residue(i+1u).has("CEN")) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index("CEN"),i+1u),Vec(out[21*i+15],out[21*i+16],out[21*i+17]));
      if(tmp.residue(i+1u).has( "H" )) tmp.set_xyz(AtomID(tmp.residue(i+1u).atom_index( "H" ),i+1u),Vec(out[21*i+18],out[21*i+19],out[21*i+20]));
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
  TR << "RUNTIME: ros/gpu " << F(10,5,100.0*tros/tcl) << "x   gpu/cpu " << F(10,5,100.0*tpar/tcl) << "x" << endl;
  TR << endl;
  //TR << "tros: " << tros << "                " << t3.x() << t4.x() << std::endl;
  TR << "IDL/CPU RMS: " << numeric::model_quality::calc_rms(idl_crd,cpu_crd) << endl;
  TR << "ROS/CPU RMS: " << numeric::model_quality::calc_rms(cpu_crd,ros_crd) << endl;
  TR << "ROS/NAT RMS: " << numeric::model_quality::calc_rms(nat_crd,ros_crd) << endl;
  TR << "IDL/NAT RMS: " << numeric::model_quality::calc_rms(nat_crd,idl_crd) << endl;  
  TR << "CPU/NAT RMS: " << numeric::model_quality::calc_rms(nat_crd,cpu_crd) << endl;
  // TR << "ROSFLD/GPU RMS: " << numeric::model_quality::calc_rms(ros_crd,gpu_crd) << endl;  


}
