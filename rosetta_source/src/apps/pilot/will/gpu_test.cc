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
//#include <protocols/scoring/ImplicitFastClashCheck.hh>
//#include <protocols/moves/MinMover.hh>
//#include <protocols/moves/PackRotamersMover.hh>
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


// void repack(Pose & arg) {
//   ScoreFunctionOP sf = core::scoring::getScoreFunction();
//   core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(arg);
//   task->restrict_to_repacking();
//   for(Size i=1; i<=arg.n_residue(); ++i) if(arg.residue(i).name3()=="HIS") task->nonconst_residue_task(i).prevent_repacking();
//   protocols::moves::PackRotamersMover repack( sf, task );
//   repack.apply(arg);
// }


int main(int argc, char *argv[]) {
  core::init(argc,argv);
  using namespace basic::options;

  // gpu_refold_test(1);
  // gpu_score_test(HVY);
  // utility_exit_with_message("RST");

  Size ERRTH = 0.01;

  Pose his,glu,asp,cys,p;
  pose_from_pdb(p,option[OptionKeys::in::file::s]()[1]);
  make_pose_from_sequence(cys,"C","fa_standard",false); remove_termini(cys); to_canonical_sc_frame(cys);
  make_pose_from_sequence(asp,"D","fa_standard",false); remove_termini(asp); to_canonical_sc_frame(asp);
  make_pose_from_sequence(glu,"E","fa_standard",false); remove_termini(glu); to_canonical_sc_frame(glu);
  make_pose_from_sequence(his,"H[HIS_DE]","fa_standard",false); remove_termini(his); to_canonical_sc_frame(his);
  his.set_dof(DOF_ID(AtomID(his.residue(1).atom_index("HD1"),1),core::id::D),2.1);
  his.set_dof(DOF_ID(AtomID(his.residue(1).atom_index("HE2"),1),core::id::D),2.1);
  his.set_dof(DOF_ID(AtomID(cys.residue(1).atom_index("HG"),1),core::id::D),2.1);

  cout << std::setprecision(8) << std::fixed;

  float c11 = -69.9394;
  float c21 = 166.2;
  float c31 = uniform()*360.0;
  float c12 = 0.0;//-70.0784;
  float c22 = 0.0;//165.854;
  float c32 = 0.0;//0.0f;//uniform()*360.0;
  p.replace_residue(19,his.residue(1),true);
  p.replace_residue(65,his.residue(1),true);
  p.set_chi(1,19,c11);
  p.set_chi(2,19,c21);
  p.set_chi(1,65,c12);
  p.set_chi(2,65,c22);

  // repack(p);
  // TR << p.chi(1,19) << " " << p.chi(2,19) << std::endl;
  // TR << p.chi(1,65) << " " << p.chi(2,65) << std::endl;
  // his.set_chi(1,1,p.chi(1,19));
  // his.set_chi(2,1,p.chi(2,19));
  // his.dump_pdb("his0.pdb");



  vector1<XFORM> toloc(p.n_residue());
  for(uint i = 1; i <= p.n_residue(); ++i) {
    toloc[i] = stubrev(VEC(p.xyz(AtomID(5,i))),VEC(p.xyz(AtomID(2,i))),VEC(p.xyz(AtomID(1,i))));
  }
  xform_pose(p,toloc[19].stub());
  //core::kinematics::Stub s(p.xyz(AtomID(5,19)),p.xyz(AtomID(2,19)),p.xyz(AtomID(1,19)));
  //xform_pose_rev(p,s);
  //  p.dump_pdb("p0.pdb");
  // {
  //   XFORM x1 = hisd_bb2m(radians(p.chi(1,19)),radians(p.chi(2,19)),radians(c31));
  //   Vecf X = (p.residue(19).xyz("HD1")-p.residue(19).xyz("ND1")).normalized();
  //   Vecf Y = projperp(X, p.residue(19).xyz("ND1")-p.residue(19).xyz("CE1") ).normalized();
  //   Y = (Matf)rotation_matrix_degrees((Vec)X,90.0+c31) * Y;
  //   Vecf M = p.residue(19).xyz("HD1");
  //   if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH ||
  //       fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH ||
  //       fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ) {
  //     TR << "xform from hisd_bb2m FAIL:!" << std::endl;
  //     cout << std::setprecision(8) << std::fixed;
  //     cout << "MAIN: " << X << " " << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl;
  //     cout << "MAIN: " << Y << " " << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl;
  //     cout << "MAIN: " << M << " " << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
  //   } else {
  //     TR << "xform from hisd_bb2m seems correct!" << std::endl;
  //   }
  // }
  // {
  //   XFORM x1 = hise_bb2m(radians(p.chi(1,19)),radians(p.chi(2,19)),radians(c31));
  //   Vecf X = (p.residue(19).xyz("HE2")-p.residue(19).xyz("NE2")).normalized();
  //   Vecf Y = projperp(X, p.residue(19).xyz("NE2")-p.residue(19).xyz("CD2") ).normalized();
  //   Y = (Matf)rotation_matrix_degrees((Vec)X,90.0+c31) * Y;
  //   Vecf M = p.residue(19).xyz("HE2");
  //   if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH ||
  //       fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH ||
  //       fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ) {
  //     TR << "xform from hise_bb2m FAIL:!" << std::endl;
  //     cout << std::setprecision(8) << std::fixed;
  //     cout << "MAIN: " << X << " " << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl;
  //     cout << "MAIN: " << Y << " " << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl;
  //     cout << "MAIN: " << M << " " << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
  //   } else {
  //     TR << "xform from hise_bb2m seems correct!" << std::endl;
  //   }
  // }


  c11=0.0;
  c12=0.0;
  c13=0.0;
  p.replace_residue(19,cys.residue(1),true);
  p.set_chi(1,19,c11);
  p.set_chi(2,19,c12);
  core::kinematics::Stub s(p.xyz(AtomID(5,19)),p.xyz(AtomID(2,19)),p.xyz(AtomID(1,19)));
  xform_pose_rev(p,s);
  p.dump_pdb("test.pdb");
  {
    XFORM x1 = cys_bb2m(radians(p.chi(1,19)),radians(p.chi(2,19)),radians(c13));
    Vecf X = (p.residue(19).xyz("HG")-p.residue(19).xyz("SG")).normalized();
    Vecf Y = projperp(X, p.residue(19).xyz("SG")-p.residue(19).xyz("CB") ).normalized();
    Y = (Matf)rotation_matrix_degrees((Vec)X,90.0+c13) * Y;
    Vecf M = p.residue(19).xyz("HG");
    if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH ||
        fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH ||
        fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ) {
      TR << "xform from cys_bb2m FAIL:!" << std::endl;
      cout << std::setprecision(8) << std::fixed;
      cout << "MAIN: " << X << " " << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl;
      cout << "MAIN: " << Y << " " << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl;
      cout << "MAIN: " << M << " " << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
    } else {
      TR << "xform from hise_bb2m seems correct!" << std::endl;
    }
  }


}
