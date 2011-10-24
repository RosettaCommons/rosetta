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

//#include <apps/pilot/will/gpu_xforms.hh>
//#include <apps/pilot/will/gpu_bit_utils.hh>
//#include <apps/pilot/will/gpu_score.hh>
//#include <apps/pilot/will/gpu_pack14.h>
//#include <apps/pilot/will/dump_unique_atoms.hh>
//#include <apps/pilot/will/gpu_test_k_square.hh>
//#include <apps/pilot/will/gpu_speedtest.hh>
#include <apps/pilot/will/gpu_refold.hh>


int main(int argc, char *argv[]) {
  core::init(argc,argv);
  using namespace basic::options;

  gpu_refold_test(max(1,option[OptionKeys::out::nstruct]()));


}
