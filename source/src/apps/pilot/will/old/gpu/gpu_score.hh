#ifndef INCLUDED_apps_pilot_will_gpu_gpu_score_hh
#define INCLUDED_apps_pilot_will_gpu_gpu_score_hh


#include <apps/pilot/will/xyzStripeHashWithMeta.hh>
#include <apps/pilot/will/xyzStripeHashPoseWithMeta.hh>
//#include <apps/pilot/will/gpu/gpu_refold.hh>


void gpu_score_test(xyzStripeHashPoseWithMetaMode pom) {
  using namespace basic::options;
  using core::pose::Pose;


  Pose p;
  core::import_pose::pose_from_file(p,option[OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);
	for(Size ir = 1; ir <= p.size(); ++ir) {
		if( p.residue(ir).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(p,ir);
		if( p.residue(ir).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(p,ir);		
	}

  //uint natom = pose_natom(p);
  uint N = p.size();

  uint runmode = option[OptionKeys::out::nstruct]();
  uint NITER = runmode/10000;
  Real GRADIUS = Real((runmode/10)%1000)/100.0;
  runmode = runmode%10;
  TR<<"MODE: "<<runmode<<" "<<NITER<<" "<<GRADIUS<<std::endl;


  xyzStripeHashPoseWithMeta poc(GRADIUS,p,pom);
//  protocols::scoring::ImplicitFastClashCheck ifc(p,GRADIUS);


  uint dim[8];
  dim[0] = poc.xdim();
  dim[1] = poc.ydim();
  dim[2] = poc.zdim();
  dim[3] = dim[0]*dim[1]*dim[2];
  dim[4] = poc.natom();

  CL cl(cerr);
  cl.make_kernel("octree");
  cl_mem gatom   = cl.makeROmem(sizeof(cl_float4)*(poc.natom()+4));
  cl_mem gstripe = cl.makeROmem(sizeof(cl_ushort2)*dim[3]);
  cl_mem gsize   = cl.makeROmem(sizeof(cl_float)*1);
  cl_mem gdim    = cl.makeROmem(sizeof(cl_uint8)*1);
  cl_mem output  = cl.makeWOmem(sizeof(cl_float)*500u*NITER);
  cl.cpu2gpu( poc.grid_atoms() , gatom  , sizeof(cl_float4 )*poc.natom() );
  cl.cpu2gpu( poc.grid_stripe(), gstripe, sizeof(cl_ushort2)*dim[3] );
  cl.cpu2gpu( &poc.grid_size() , gsize  , sizeof(cl_float)*1 );
  cl.cpu2gpu( dim              , gdim   , sizeof(cl_uint8)*1 );
  cl.setargs("octree",gatom,gstripe,gsize,gdim,output);


  float count = 0.0f;
  double tcount = 0.0f;

  // double tom = time_highres();
  // count = 0u; tcount = 0u;
  // for(Size iter=0;iter<NITER/100;++iter) {
  float sctest = octree_test( (float4*)poc.grid_atoms(), (ushort2*)poc.grid_stripe(), &poc.grid_size(), (uint8*)dim,tcount);
  //   tcount += count;
  // }
  // tom = time_highres()-tom;
  TR<<"octree gpu mimic: "<< poc.natom() << " " << F(20,10,sctest) << "                                " <<endl;
  // //TR<<poc.nbcount( poc.grid_atoms()[0].x, poc.grid_atoms()[0].y, poc.grid_atoms()[0].z )<<std::endl;
  
  core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
  //sf.set_weight(core::scoring::fa_atr,1.0);
  //sf.set_weight(core::scoring::fa_rep,1.0);
  //sf.set_weight(core::scoring::fa_sol,1.0);
  
  double tno = 0.0;
  for(Size iter=0;iter<NITER*1;++iter) {
    trans_pose(p,Vec(0.00001,0,0));
    double t = time_highres();
    sf->score(p);
    tno += time_highres()-t;
  }
  TR << "score time: " << tno << std::endl;
  // tcount = 0u;
  // for(Size iter=0;iter<NITER/10;++iter) {
  //   count = 0u;
  //   for(Size i = 0; i < natom; ++i) {
  //     float4 const & a = *((float4*)(&poc.grid_atoms()[i]));
  //     for(Size j = i+1; j < natom; ++j) {
  //       float4 const & a1 = *((float4*)(&poc.grid_atoms()[j]));
  //       if( sqr(a.x-a1.x)+sqr(a.y-a1.y)+sqr(a.z-a1.z) <= GRADIUS*GRADIUS ) {
  //         ++count;
  //       }
  //     }
  //   }
  //   tcount += count;
  // }
  //tno = time_highres()-tno;
  // TR<<"no octree:        "<<natom<<" "<<count<<"                                "<<tcount<<std::endl;

  // double tcl = 0;
  // for(Size iter=0;iter<NITER*1;++iter) {
  //   double t = time_highres();
  //   cl.q1d("octree",256u*1000u,256u);
  //   cl.finish();
  //   tcl += time_highres()-t;
  // }

  double tcl = time_highres();
  cl.q1d("octree",256u*NITER*500u,256u);
  cl.finish();
  tcl = time_highres()-tcl;
  TR << "gpu score time: " << tcl << std::endl;
  float *result = new float[NITER*500u];
  for(uint i = 0; i<NITER*500;++i) result[i] = 0.0f;
  cl.gpu2cpu(output,result,sizeof(cl_uint)*500*NITER);
  tcount = 0; for(uint i = 0; i<NITER*500;++i) {
    tcount += result[i];
    //TR << i << " " << F(20,10,result[i]) << std::endl;
  }

  TR << "gpue:             "<<poc.natom()<<" "<< F(20,10,result[0])<<"                                "<<F(20,10,tcount/(double)NITER/(double)500)<<std::endl;
  TR << "gpu speedup: " << 500*tno / tcl << std::endl;

}

#endif
