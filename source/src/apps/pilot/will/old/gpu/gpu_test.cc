#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
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
//#include <protocols/minimization_packing/MinMover.hh>
//#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


static basic::Tracer TR( "gpu_test" );


OPT_1GRP_KEY( Integer, gpu, numthreads_per_workunit )
OPT_1GRP_KEY( File   , gpu, kernel )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( gpu::numthreads_per_workunit ,"local work unit dim", 128 );
	NEW_OPT( gpu::kernel ,"kernel file", "" );


}

#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>
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

//#include <apps/pilot/will/gpu/gpu_xforms.hh>
//#include <apps/pilot/will/gpu/gpu_bit_utils.hh>
//#include <apps/pilot/will/gpu/gpu_score.hh>
//#include <apps/pilot/will/gpu/gpu_pack14.h>
//#include <apps/pilot/will/dump_unique_atoms.hh>
//#include <apps/pilot/will/gpu/gpu_test_k_square.hh>
//#include <apps/pilot/will/gpu/gpu_speedtest.hh>
#include <apps/pilot/will/gpu/gpu_refold.hh>

#include <core/scoring/AtomVDW.hh>

int main(int argc, char *argv[]) {
	register_options();
	devel::init(argc,argv);
	using namespace basic::options;

	//   {
	//   uint const N = 8;
	//   uint const Nwork = (N*(N+1u))/2u-1;
	//   for(uint ichunk = N; ichunk <= Nwork; ichunk+=8) {
	//     for(uint ith = 0; ith < 8; ith++) {
	//        if( ichunk+ith <= Nwork) {
	//         uint const i = Nwork - (ichunk+ith);
	//         uint const n = floor((sqrt(1.0f+8.0f*(float)i)-1.0f)/2.0f);
	//         uint const jr = i - (n*(n+1u))/2u;
	//         uint const ir = n-jr;
	//         cout << ir << " " << jr << "   " << ir+jr << "   " << ichunk << " " << ith << endl;
	//     }
	// }
	//   }
	//   utility_exit_with_message("aorstn");
	// }

	// Nbb     1 1.28109
	// CAbb    2 1.32174
	// CB      3 1.29833
	// CObb    4 1.28592
	// OCbb    5 0.843484
	// HNbb    6 1
	//  vector1<Size> at;
	//  at.push_back(1); // N
	//  at.push_back(2); // CA
	//  at.push_back(4); // C
	//  at.push_back(5); // O
	//  at.push_back(3); // CB
	// core::scoring::AtomVDW const & atom_vdw( core::scoring::ScoringManager::get_instance()->get_AtomVDW( "centroid" ) );
	// for(Size i = 1; i <= 5; ++i) {
	//  for(Size j = 1; j <= 5; ++j) {
	//   float v = atom_vdw(at[i])[at[j]];
	//   cout << v << "// " << at[i] << " " << at[j] << endl;
	//  }
	// }

	// utility_exit_with_message("arst");

	gpu_refold_test(max(1,option[OptionKeys::out::nstruct]()));


}
