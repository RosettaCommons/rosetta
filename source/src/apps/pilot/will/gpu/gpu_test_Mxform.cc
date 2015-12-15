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
#include <devel/init.hh>
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
//#include <protocols/simple_moves/MinMover.hh>
//#include <protocols/simple_moves/PackRotamersMover.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>

static THREAD_LOCAL basic::Tracer TR( "gpu_test" );

#include <apps/pilot/will/gpu/CL.hh>


typedef numeric::xyzVector<float> Vecf;
typedef numeric::xyzMatrix<float> Matf;
typedef cl_float16 float16;
typedef cl_float8  float8;
typedef cl_float4  float4;
typedef cl_float3  float3;
typedef cl_float2  float2;
typedef cl_ushort2 ushort2;
typedef cl_uint8 uint8;

#include <apps/pilot/will/gpu/gpu_xforms.hh>
#include <apps/pilot/will/gpu/gpu_bit_utils.hh>
//#include <apps/pilot/will/gpu/gpu_score.hh>
//#include <apps/pilot/will/gpu/gpu_pack14.h>
//#include <apps/pilot/will/dump_unique_atoms.hh>
//#include <apps/pilot/will/gpu/gpu_test_k_square.hh>
//#include <apps/pilot/will/gpu/gpu_speedtest.hh>


// void repack(Pose & arg) {
//   ScoreFunctionOP sf = core::scoring::get_score_function();
//   core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(arg);
//   task->restrict_to_repacking();
//   for(Size i=1; i<=arg.n_residue(); ++i) if(arg.residue(i).name3()=="HIS") task->nonconst_residue_task(i).prevent_repacking();
//   protocols::simple_moves::PackRotamersMover repack( sf, task );
//   repack.apply(arg);
// }


int main(int argc, char *argv[]) {
  devel::init(argc,argv);
  using namespace basic::options;

  // gpu_refold_test(1);
  // gpu_score_test(HVY);
  // utility_exit_with_message("RST");

  Real ERRTH = 0.0003;

  Pose hmd,hme,dm1,dm2,em1,em2,cys,p;
  pose_from_pdb(p,option[OptionKeys::in::file::s]()[1]);
  make_pose_from_sequence(cys,"C[CYS_M]" ,core::chemical::FA_STANDARD,false); remove_termini(cys); to_canonical_sc_frame(cys);
  make_pose_from_sequence(dm1,"D[ASP_M1]",core::chemical::FA_STANDARD,false); remove_termini(dm1); to_canonical_sc_frame(dm1);
  make_pose_from_sequence(dm2,"D[ASP_M2]",core::chemical::FA_STANDARD,false); remove_termini(dm2); to_canonical_sc_frame(dm2);
  make_pose_from_sequence(em1,"E[GLU_M1]",core::chemical::FA_STANDARD,false); remove_termini(em1); to_canonical_sc_frame(em1);
  make_pose_from_sequence(em2,"E[GLU_M2]",core::chemical::FA_STANDARD,false); remove_termini(em2); to_canonical_sc_frame(em2);
  make_pose_from_sequence(hmd,"H[HIS_MD]",core::chemical::FA_STANDARD,false); remove_termini(hmd); to_canonical_sc_frame(hmd);
  make_pose_from_sequence(hme,"H[HIS_ME]",core::chemical::FA_STANDARD,false); remove_termini(hme); to_canonical_sc_frame(hme);

  cout << std::setprecision(8) << std::fixed;

  for(Size iter = 0; iter < option[OptionKeys::out::nstruct](); ++iter) {

    float c11 = uniform()*360.0;
    float c21 = uniform()*360.0;
    float c31 = uniform()*360.0;
    float c41 = uniform()*360.0;
    float c12 = uniform()*360.0;
    float c22 = uniform()*360.0;
    float c32 = uniform()*360.0;
    float c42 = uniform()*360.0;
    p.replace_residue(19,hmd.residue(1),true);
    p.replace_residue(65,hme.residue(1),true);
    p.set_chi(1,19,c11);
    p.set_chi(2,19,c21);
    p.set_chi(3,19,c31);
    p.set_chi(1,65,c12);
    p.set_chi(2,65,c22);
    p.set_chi(3,65,c32);

    // repack(p);
    // TR << p.chi(1,19) << " " << p.chi(2,19) << std::endl;
    // TR << p.chi(1,65) << " " << p.chi(2,65) << std::endl;
    // his.set_chi(1,1,p.chi(1,19));
    // his.set_chi(2,1,p.chi(2,19));
    // his.dump_pdb("his0.pdb");


    
    
    // xform_pose_rev(p,core::kinematics::Stub(p.residue(19).xyz("CB"),p.residue(19).xyz("CA"),p.residue(19).xyz("N")));
    // p.dump_pdb("test0.pdb");
    // core::kinematics::Stub sc(p.residue(19).xyz("CA"),p.residue(19).xyz("N"),p.residue(19).xyz("C"));
    // cout << sc.M.transposed() << endl;
    // cout << sc.M.transposed()*-sc.v << endl;
    // cout << endl;
    // Vec ct = sc.v;
    // Mat cr = sc.M;
    // Real cw;
    // Vec ca;
    // ca = rotation_axis(cr,cw);
    // TR << ca << endl;
    // TR << ct << endl;
    // TR << cw << endl;
    // xform_pose_rev(p,sc);
    // p.dump_pdb("test1.pdb");

    vector1<XFORM> toloc(p.n_residue());
    for(uint i = 1; i <= p.n_residue(); ++i) {
      //toloc[i] = stubcrev(VEC(p.xyz(AtomID(5,i))),VEC(p.xyz(AtomID(5,i))),VEC(p.xyz(AtomID(2,i))),VEC(p.xyz(AtomID(1,i))));
      toloc[i] = stubrev(VEC(p.xyz(AtomID(2,i))),VEC(p.xyz(AtomID(1,i))),VEC(p.xyz(AtomID(3,i))));
    }
    Pose init = p;

    {
      p = init; xform_pose(p,toloc[19].stub());
      // p.dump_pdb("p0.pdb");
      // std::cout << "MAIN MY: " << p.residue(19).xyz("MDY") << std::endl;
      // std::cout << "MAIN ML: " << p.residue(19).xyz("MD" ) << std::endl;
      // std::cout << "MAIN ND: " << p.residue(19).xyz("ND1") << std::endl;
      // std::cout << "MAIN NE: " << p.residue(19).xyz("NE2") << std::endl;
      // std::cout << "MAIN CE: " << p.residue(19).xyz("CE1") << std::endl;
      // std::cout << "MAIN CD: " << p.residue(19).xyz("CD2") << std::endl;
      // std::cout << "MAIN CG: " << p.residue(19).xyz("CG" ) << std::endl;
      XFORM x1 = hisd_bb2m(radians(p.chi(1,19)),radians(p.chi(2,19)),radians(p.chi(3,19)));
      Vecf X = (p.residue(19).xyz("MD" )-p.residue(19).xyz("ND1")).normalized();
      Vecf Y = (p.residue(19).xyz("MDY")-p.residue(19).xyz("MD" )).normalized();
      Vecf M =  p.residue(19).xyz("MD");
      if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH || fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH || fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ){
        TR << "xform from hisd_bb2m FAIL:!" << std::endl; cout << "MAIN: " << endl << X << endl << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl;      cout << "MAIN: " << endl << Y << endl << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl;     cout << "MAIN: " << endl << M << endl << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
        utility_exit_with_message("YOU SUCK!!!!");
      } else TR << "xform from hisd_bb2m seems correct!" << std::endl;
    }
    //utility_exit_with_message("airoeth");
    {
      p = init; xform_pose(p,toloc[65].stub());
      XFORM x1 = hise_bb2m(radians(p.chi(1,65)),radians(p.chi(2,65)),radians(p.chi(3,65)));
      Vecf X = (p.residue(65).xyz("ME" )-p.residue(65).xyz("NE2")).normalized();
      Vecf Y = (p.residue(65).xyz("MEY")-p.residue(65).xyz("ME" )).normalized();
      Vecf M =  p.residue(65).xyz("ME");
      if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH || fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH || fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ){
        TR << "xform from hise_bb2m FAIL:!" << std::endl; cout << std::setprecision(8) << std::fixed; cout << "MAIN: " << endl << X << endl << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl; cout << "MAIN: " << endl << Y << endl << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl; cout << "MAIN: " << endl << M << endl << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
        utility_exit_with_message("YOU SUCK!!!!");
      } else TR << "xform from hise_bb2m seems correct!" << std::endl;
    }

    p = init;
    p.replace_residue(19,dm1.residue(1),true);
    p.replace_residue(65,dm2.residue(1),true);
    p.set_chi(1,19,c11);
    p.set_chi(2,19,c21);
    p.set_chi(3,19,c31);
    p.set_chi(1,65,c12);
    p.set_chi(2,65,c22);
    p.set_chi(3,65,c32);
    init = p;


    {
      //core::kinematics::Stub s(p.residue(19).xyz("CG"),p.residue(19).xyz("CB"),p.residue(19).xyz("CA"));
      //core::kinematics::Stub s(p.residue(19).xyz("MD"),p.residue(19).xyz("ND1"),p.residue(19).xyz("CG"));
      //xform_pose_rev(p,s);
      //p.dump_pdb("p0.pdb");
      p = init; xform_pose(p,toloc[19].stub());
      XFORM x1 = aspd_bb2m(radians(p.chi(1,19)),radians(p.chi(2,19)),radians(p.chi(3,19)));
      Vecf X = (p.residue(19).xyz("M1" )-p.residue(19).xyz("OD1")).normalized();
      Vecf Y = (p.residue(19).xyz("M1Y")-p.residue(19).xyz("M1" )).normalized();
      Vecf M =  p.residue(19).xyz("M1");
      if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH || fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH || fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ){
        TR << "xform from aspd_bb2m FAIL:!" << std::endl; cout << "MAIN: " << endl << X << endl << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl;      cout << "MAIN: " << endl << Y << endl << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl;     cout << "MAIN: " << endl << M << endl << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
        utility_exit_with_message("YOU SUCK!!!!");
      } else TR << "xform from aspd_bb2m seems correct!" << std::endl;
    }


    {
      //core::kinematics::Stub s(p.residue(65).xyz("CB"),p.residue(65).xyz("CA"),p.residue(65).xyz("N"));
      //core::kinematics::Stub s(p.residue(65).xyz("CG"),p.residue(65).xyz("CB"),p.residue(65).xyz("CA"));
      //core::kinematics::Stub s(p.residue(65).xyz("MD"),p.residue(65).xyz("ND1"),p.residue(65).xyz("CG"));
      //xform_pose_rev(p,s);
      // std::cout << "MAIN MY: " << p.residue(65).xyz("M2Y") << std::endl;
      // std::cout << "MAIN ML: " << p.residue(65).xyz("M2" ) << std::endl;
      // std::cout << "MAIN O1: " << p.residue(65).xyz("OD1") << std::endl;
      // std::cout << "MAIN O2: " << p.residue(65).xyz("OD2") << std::endl;
      // std::cout << "MAIN CG: " << p.residue(65).xyz("CG" ) << std::endl;
      //    p.dump_pdb("p0.pdb");
      p = init; xform_pose(p,toloc[65].stub());
      XFORM x1 = aspe_bb2m(radians(p.chi(1,65)),radians(p.chi(2,65)),radians(p.chi(3,65)));
      Vecf X = (p.residue(65).xyz("M2" )-p.residue(65).xyz("OD1")).normalized();
      Vecf Y = (p.residue(65).xyz("M2Y")-p.residue(65).xyz("M2" )).normalized();
      Vecf M =  p.residue(65).xyz("M2");
      if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH || fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH || fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ){
        TR << "xform from aspe_bb2m FAIL:!" << std::endl; cout << std::setprecision(8) << std::fixed; cout << "MAIN: " << endl << X << endl << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl; cout << "MAIN: " << endl << Y << endl << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl; cout << "MAIN: " << endl << M << endl << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
        utility_exit_with_message("YOU SUCK!!!!");
      } else TR << "xform from aspe_bb2m seems correct!" << std::endl;
    }


    p = init;
    p.replace_residue(19,em1.residue(1),true);
    p.replace_residue(65,em2.residue(1),true);
    p.set_chi(1,19,c11);
    p.set_chi(2,19,c21);
    p.set_chi(3,19,c31);
    p.set_chi(4,19,c41);
    p.set_chi(1,65,c12);
    p.set_chi(2,65,c22);
    p.set_chi(3,65,c32);
    p.set_chi(4,65,c42);
    init = p;

    {
      //core::kinematics::Stub s(p.residue(19).xyz("CG"),p.residue(19).xyz("CB"),p.residue(19).xyz("CA"));
      //core::kinematics::Stub s(p.residue(19).xyz("MD"),p.residue(19).xyz("ND1"),p.residue(19).xyz("CG"));
      //xform_pose_rev(p,s);
      //p.dump_pdb("p0.pdb");
      p = init; xform_pose(p,toloc[19].stub());
      XFORM x1 = glud_bb2m(radians(p.chi(1,19)),radians(p.chi(2,19)),radians(p.chi(3,19)),radians(p.chi(4,19)));
      Vecf X = (p.residue(19).xyz("M1" )-p.residue(19).xyz("OE1")).normalized();
      Vecf Y = (p.residue(19).xyz("M1Y")-p.residue(19).xyz("M1" )).normalized();
      Vecf M =  p.residue(19).xyz("M1");
      if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH || fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH || fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ){
        TR << "xform from glud_bb2m FAIL:!" << std::endl; cout << "MAIN: " << endl << X << endl << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl;      cout << "MAIN: " << endl << Y << endl << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl;     cout << "MAIN: " << endl << M << endl << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
        utility_exit_with_message("YOU SUCK!!!!");
      } else TR << "xform from glud_bb2m seems correct!" << std::endl;
    }
    {
      //core::kinematics::Stub s(p.residue(65).xyz("CB"),p.residue(65).xyz("CA"),p.residue(65).xyz("N"));
      //core::kinematics::Stub s(p.residue(65).xyz("CG"),p.residue(65).xyz("CB"),p.residue(65).xyz("CA"));
      //core::kinematics::Stub s(p.residue(65).xyz("MD"),p.residue(65).xyz("ND1"),p.residue(65).xyz("CG"));
      //xform_pose_rev(p,s);
      // std::cout << "MAIN MY: " << p.residue(65).xyz("M2Y") << std::endl;
      // std::cout << "MAIN ML: " << p.residue(65).xyz("M2") << std::endl;
      // std::cout << "MAIN O1: " << p.residue(65).xyz("OD1") << std::endl;
      // std::cout << "MAIN O2: " << p.residue(65).xyz("OD2") << std::endl;
      // std::cout << "MAIN CG: " << p.residue(65).xyz("CG") << std::endl;
      //    p.dump_pdb("p0.pdb");
      p = init; xform_pose(p,toloc[65].stub());
      XFORM x1 = glue_bb2m(radians(p.chi(1,65)),radians(p.chi(2,65)),radians(p.chi(3,65)),radians(p.chi(4,65)));
      Vecf X = (p.residue(65).xyz("M2" )-p.residue(65).xyz("OE1")).normalized();
      Vecf Y = (p.residue(65).xyz("M2Y")-p.residue(65).xyz("M2" )).normalized();
      Vecf M =  p.residue(65).xyz("M2");
      if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH || fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH || fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ){
        TR << "xform from glue_bb2m FAIL:!" << std::endl; cout << std::setprecision(8) << std::fixed; cout << "MAIN: " << endl << X << endl << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl; cout << "MAIN: " << endl << Y << endl << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl; cout << "MAIN: " << endl << M << endl << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
        utility_exit_with_message("YOU SUCK!!!!");
      } else TR << "xform from glue_bb2m seems correct!" << std::endl;
    }


    p = init;
    p.replace_residue(19,cys.residue(1),true);
    p.set_chi(1,19,c11);
    p.set_chi(2,19,c21);
    p.set_chi(3,19,c31);
    init = p;

    {
      //core::kinematics::Stub s(p.residue(19).xyz("CG"),p.residue(19).xyz("CB"),p.residue(19).xyz("CA"));
      //core::kinematics::Stub s(p.residue(19).xyz("MD"),p.residue(19).xyz("ND1"),p.residue(19).xyz("CG"));
      //xform_pose_rev(p,s);
      //p.dump_pdb("p0.pdb");
      p = init; xform_pose(p,toloc[19].stub());
      XFORM x1 = cys_bb2m(radians(p.chi(1,19)),radians(p.chi(2,19)),radians(p.chi(3,19)));
      Vecf X = (p.residue(19).xyz("M" )-p.residue(19).xyz("SG")).normalized();
      Vecf Y = (p.residue(19).xyz("MY")-p.residue(19).xyz("M" )).normalized();
      Vecf M =  p.residue(19).xyz("M");
      if( fabs(X.x()-x1.R.xx) > ERRTH || fabs(X.y()-x1.R.yx) > ERRTH || fabs(X.z()-x1.R.zx) > ERRTH || fabs(Y.x()-x1.R.xy) > ERRTH || fabs(Y.y()-x1.R.yy) > ERRTH || fabs(Y.z()-x1.R.zy) > ERRTH || fabs(M.x()-x1.t.x ) > ERRTH || fabs(M.y()-x1.t.y ) > ERRTH || fabs(M.z()-x1.t.z ) > ERRTH ){
        TR << "xform from cys_bb2m FAIL:!" << std::endl; cout << "MAIN: " << endl << X << endl << x1.R.xx << " " << x1.R.yx << " " << x1.R.zx << std::endl;      cout << "MAIN: " << endl << Y << endl << x1.R.xy << " " << x1.R.yy << " " << x1.R.zy << std::endl;     cout << "MAIN: " << endl << M << endl << x1.t.x  << " " << x1.t.y  << " " << x1.t.z  << std::endl;
        utility_exit_with_message("YOU SUCK!!!!");
      } else TR << "xform from cys_bb2m seems correct!" << std::endl;
    }

  }


  float c11 = uniform()*360.0; float c21 = uniform()*360.0; float c31 = uniform()*360.0; float c41 = uniform()*360.0;
  float c12 = uniform()*360.0; float c22 = uniform()*360.0; float c32 = uniform()*360.0; float c42 = uniform()*360.0;
  p.replace_residue(19,hmd.residue(1),true);
  p.replace_residue(65,hme.residue(1),true);
  p.set_chi(1,19,c11);
  p.set_chi(2,19,c21);
  p.set_chi(3,19,c31);
  p.set_chi(1,65,c12);
  p.set_chi(2,65,c22);
  p.set_chi(3,65,c32);
  Pose q(p);

  struct XFORM tol19 = stubrev(p,19);
  struct XFORM tog19 = stub   (p,19);
  struct XFORM tol65 = stubrev(p,65);
  struct XFORM tog65 = stub   (p,65);  
  struct XFORM b2m19 = hisd_bb2m(c11,c21,c31);
  struct XFORM m2b19 = hisd_m2bb(c11,c21,c31);
  struct XFORM b2m65 = hise_bb2m(c12,c22,c32);
  struct XFORM m2b65 = hise_m2bb(c12,c22,c32);

  tol19.apply(p);
  tol65.apply(q);

  TR << b2m19.t << endl;
  TR << m2b19.t << endl;
  TR << p.residue(19).xyz("MD") << endl;

  //  m2b19.apply(p);
  //  m2b65.apply(q);


  p.dump_pdb("testp.pdb");
  q.dump_pdb("testq.pdb");


}
