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


static THREAD_LOCAL basic::Tracer TR( "gpu_score_test" );


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

#include <apps/pilot/will/gpu/gpu_score.hh>

int main(int argc, char *argv[]) {
	register_options();
  devel::init(argc,argv);
  using namespace basic::options;


	gpu_score_test(BB);


}
