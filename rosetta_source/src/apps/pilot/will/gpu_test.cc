#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/gpu.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <protocols/moves/MinMover.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include "apps/pilot/will/will_util.hh"
#include "mynamespaces.hh"

static basic::Tracer TR("gpu_test");

#include <apps/pilot/will/CL.hh>


typedef numeric::xyzVector<float> Vecf;
typedef numeric::xyzMatrix<float> Matf;
typedef cl_float16 float16;
typedef cl_float8  float8;
typedef cl_float4  float4;
typedef cl_float3  float3;
typedef cl_float2  float2;
typedef cl_ushort2 ushort2;
typedef cl_uint8 uint8;

#include <apps/pilot/will/gpu_xforms.hh>
#include <apps/pilot/will/gpu_bit_utils.hh>
//#include <apps/pilot/will/gpu_score.hh>
//#include <apps/pilot/will/gpu_pack14.h>
//#include <apps/pilot/will/dump_unique_atoms.hh>
//#include <apps/pilot/will/gpu_test_k_square.hh>
//#include <apps/pilot/will/gpu_speedtest.hh>




int main(int argc, char *argv[]) {
  core::init(argc,argv);
  using namespace basic::options;

  // gpu_refold_test(1);
  // gpu_score_test(HVY);
  // utility_exit_with_message("RST");

  Pose his,glu,asp,cys,p;
  pose_from_pdb(p,option[OptionKeys::in::file::s]()[1]);
  make_pose_from_sequence(cys,"C","fa_standard",false); remove_termini(cys); to_canonical_sc_frame(cys); 
  make_pose_from_sequence(asp,"D","fa_standard",false); remove_termini(asp); to_canonical_sc_frame(asp); 
  make_pose_from_sequence(glu,"E","fa_standard",false); remove_termini(glu); to_canonical_sc_frame(glu); 
  make_pose_from_sequence(his,"H[HIS_DE]","fa_standard",false); remove_termini(his); to_canonical_sc_frame(his);
  his.set_dof(DOF_ID(AtomID(his.residue(1).atom_index("HD1"),1),core::id::D),2.1);
  his.set_dof(DOF_ID(AtomID(his.residue(1).atom_index("HE2"),1),core::id::D),2.1);
  // cys.set_phi(  1,uniform()*360.0); asp.set_phi(  1,uniform()*360.0); glu.set_phi(  1,uniform()*360.0); his.set_phi(  1,uniform()*360.0); 
  // cys.set_psi(  1,uniform()*360.0); asp.set_psi(  1,uniform()*360.0); glu.set_psi(  1,uniform()*360.0); his.set_psi(  1,uniform()*360.0); 
  // cys.set_omega(1,uniform()*360.0); asp.set_omega(1,uniform()*360.0); glu.set_omega(1,uniform()*360.0); his.set_omega(1,uniform()*360.0); 
  // cys.set_chi(1,1,uniform()*360.0); asp.set_chi(1,1,uniform()*360.0); glu.set_chi(1,1,uniform()*360.0); his.set_chi(1,1,uniform()*360.0); 
  // cys.set_chi(2,1,uniform()*360.0); asp.set_chi(2,1,uniform()*360.0); glu.set_chi(2,1,uniform()*360.0); his.set_chi(2,1,uniform()*360.0); 
  // to_canonical_sc_frame_from_bb(cys);
  // to_canonical_sc_frame_from_bb(asp);
  // to_canonical_sc_frame_from_bb(glu);
  // to_canonical_sc_frame_from_bb(his);

  cout << std::setprecision(8) << std::fixed;

  for(Size ia = 6; ia <= cys.residue(1).natoms(); ++ia) {
    cout << "#define CYS_" << strip(cys.residue(1).atom_name(ia)) << " " << "vec("<<cys.residue(1).xyz(ia).x()<<","<<cys.residue(1).xyz(ia).y()<<","<<cys.residue(1).xyz(ia).z()<<")"<<endl;
  }
  for(Size ia = 6; ia <= asp.residue(1).nheavyatoms(); ++ia) {
    cout << "#define ASP_" << strip(asp.residue(1).atom_name(ia)) << " " << "vec("<<asp.residue(1).xyz(ia).x()<<","<<asp.residue(1).xyz(ia).y()<<","<<asp.residue(1).xyz(ia).z()<<")"<<endl;
  }
  // for(Size ia = 6; ia <= e.nheavyatoms(); ++ia) {
  //   cout << "#define GLU_" << strip(e.atom_name(ia)) << " " << "vec("<<e.xyz(ia).x()<<","<<e.xyz(ia).y()<<","<<e.xyz(ia).z()<<")"<<endl;
  // }
  for(Size ia = 6; ia <= his.residue(1).natoms(); ++ia) {
    cout << "#define HIS_" << strip(his.residue(1).atom_name(ia)) << " " << "vec("<<his.residue(1).xyz(ia).x()<<","<<his.residue(1).xyz(ia).y()<<","<<his.residue(1).xyz(ia).z()<<")"<<endl;
  }

  cys.dump_pdb("cys_canon.pdb");
  asp.dump_pdb("asp_canon.pdb");
  //glu.dump_pdb("glu_canon.pdb");
  his.dump_pdb("his_canon.pdb");
  his.set_chi(1,1,0);
  his.dump_pdb("his_0.pdb");
  his.set_chi(1,1,20);
  his.dump_pdb("his_20.pdb");
}
