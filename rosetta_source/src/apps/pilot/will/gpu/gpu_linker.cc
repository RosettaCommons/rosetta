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
#include <core/conformation/Conformation.hh>
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
//#include <protocols/simple_moves/MinMover.hh>
//#include <protocols/simple_moves/PackRotamersMover.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

static basic::Tracer TR("gpu_linker");


OPT_1GRP_KEY( Integer, gpu, numthreads_per_workunit )
OPT_1GRP_KEY( File   , gpu, kernel )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( gpu::numthreads_per_workunit ,"local work unit dim", 128 );
	NEW_OPT( gpu::kernel ,"kernel file", "NO_SOURCE_FILE" );
}

#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/gpu/CL.hh>
#include <apps/pilot/will/gpu/set_pose_to_ideal.hh>

typedef numeric::xyzVector<float> Vecf;
typedef numeric::xyzMatrix<float> Matf;
typedef cl_float16 float16;
typedef cl_float8  float8;
typedef cl_float4  float4;
typedef cl_float3  float3;
typedef cl_float2  float2;
typedef cl_ushort2 ushort2;
typedef cl_uint8 uint8;


void dump_points_pdb(float const *p, int n, string fn) {
	std::ofstream o(fn.c_str());
	for(Size i = 0; i < n; ++i) {
		o<<"ATOM  "<<I(5,3*i+1)<<' '<<"  N "<<' '<<"GLY"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,p[9*i+0])<<F(8,3,p[9*i+1])<<F(8,3,p[9*i+2])<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		o<<"ATOM  "<<I(5,3*i+2)<<' '<<" CA "<<' '<<"GLY"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,p[9*i+3])<<F(8,3,p[9*i+4])<<F(8,3,p[9*i+5])<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		o<<"ATOM  "<<I(5,3*i+3)<<' '<<"  C "<<' '<<"GLY"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,p[9*i+6])<<F(8,3,p[9*i+7])<<F(8,3,p[9*i+8])<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';				
	}
	o.close();
}

void dump_points_pdb(Pose const & p, string fn) {
	std::ofstream o(fn.c_str());
	for(Size i = 0; i < p.n_residue(); ++i) {
		o<<"ATOM  "<<I(5,3*i+1)<<' '<<"  N "<<' '<<"GLY"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,p.xyz(AtomID(1,i+1)).x())<<F(8,3,p.xyz(AtomID(1,i+1)).y())<<F(8,3,p.xyz(AtomID(1,i+1)).z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		o<<"ATOM  "<<I(5,3*i+2)<<' '<<" CA "<<' '<<"GLY"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,p.xyz(AtomID(2,i+1)).x())<<F(8,3,p.xyz(AtomID(2,i+1)).y())<<F(8,3,p.xyz(AtomID(2,i+1)).z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		o<<"ATOM  "<<I(5,3*i+3)<<' '<<"  C "<<' '<<"GLY"<<' '<<"A"<<I(4,i+1)<<"    "<<F(8,3,p.xyz(AtomID(3,i+1)).x())<<F(8,3,p.xyz(AtomID(3,i+1)).y())<<F(8,3,p.xyz(AtomID(3,i+1)).z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';				
	}
	o.close();
}

utility::vector1<numeric::xyzVector<Real> > array2vecs(float const *xyz, int N) {
	utility::vector1<numeric::xyzVector<Real> > v(N);
	for(int i = 0; i < N; ++i) {
		v[i+1].x() = xyz[3*i+0];
		v[i+1].y() = xyz[3*i+1];
		v[i+1].z() = xyz[3*i+2];				
	}
	return v;
}

 // radians!
void fere_torsions(int N, float *tor) {
	core::pose::Pose tmp;
	core::import_pose::pose_from_pdb(tmp,"input/fr52re.pdb");
	for(int i = 0; i < N; ++i) {
		tor[3*i+0] = numeric::conversions::radians(tmp.phi  (i+1));
		tor[3*i+1] = numeric::conversions::radians(tmp.psi  (i+1));
		tor[3*i+2] = numeric::conversions::radians(tmp.omega(i+1));				
	}
}

 // radians!
void random_torsions(int N, float *tor) {
	for(int i = 0; i < 3*N; ++i) {
		tor[i] = numeric::conversions::radians(numeric::random::uniform()*360.0);
	}
}

void refold_ros(core::pose::Pose & p, float const *tor, float *crd) {
	for(int i = 0; i < p.n_residue(); ++i) {
		p.set_phi  (i+1, numeric::conversions::degrees(tor[3*i+0]) );
		p.set_psi  (i+1, numeric::conversions::degrees(tor[3*i+1]) );
		p.set_omega(i+1, numeric::conversions::degrees(tor[3*i+2]) );
	}
	for(int i = 0; i < p.n_residue(); ++i) {
		crd[9*i+0] = p.xyz(AtomID(1,i+1)).x();
		crd[9*i+1] = p.xyz(AtomID(1,i+1)).y();
		crd[9*i+2] = p.xyz(AtomID(1,i+1)).z();
		crd[9*i+3] = p.xyz(AtomID(2,i+1)).x();
		crd[9*i+4] = p.xyz(AtomID(2,i+1)).y();
		crd[9*i+5] = p.xyz(AtomID(2,i+1)).z();
		crd[9*i+6] = p.xyz(AtomID(3,i+1)).x();
		crd[9*i+7] = p.xyz(AtomID(3,i+1)).y();
		crd[9*i+8] = p.xyz(AtomID(3,i+1)).z();
	}
}

void refold_gpu(float const *tor, int N, float *crd, CL & cl, cl_mem & clmtor, cl_mem & clmout) {
  cl.cpu2gpu( tor , clmtor, sizeof(cl_float)*3u*N      );
	cl.q2d("refold_test",64,1,64,1);
  cl.finish();
  cl.gpu2cpu( clmout, crd, sizeof(cl_float)*N*9u); 	
}


void test_refold() {
	using namespace basic::options;
	core::chemical::ResidueTypeSetCAP crs = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

  int N = 64;
	float tor[3*N];
	float roscrd[9*N];
	float gpucrd[9*N];	

  //TR << "init CL" << endl;
	CL cl(TR);
	string Kname = "refold_test";
  cl.make_kernel(Kname);
  vector1<cl_mem> outs;
  cl_mem clmtor   = cl.makeRWmem(sizeof(cl_float)*3u*3u*N   );
  cl_mem clmN     = cl.makeROmem(sizeof(cl_uint)            );
  cl_mem clmout   = cl.makeWOmem(sizeof(cl_float)*100u*21u*N);
  cl.cpu2gpu( &N  , clmN  , sizeof(cl_float)           ); // no doesn't change per-run
	int err;
  err = clSetKernelArg(cl.kernels_[Kname],0,sizeof(cl_mem), &clmtor   ); if(err!=CL_SUCCESS) cout << "ERR " << 0 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],1,sizeof(cl_mem), &clmN     ); if(err!=CL_SUCCESS) cout << "ERR " << 1 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],2,sizeof(cl_mem), &clmout   ); if(err!=CL_SUCCESS) cout << "ERR " << 2 << " " << cl.errstr(err) << endl;
  err = clSetKernelArg(cl.kernels_[Kname],3,sizeof(cl_float)*6*N, NULL); if(err!=CL_SUCCESS) cout << "ERR " << 3 << " " << cl.errstr(err) << endl;

  //TR << "setup pose" << endl;
	Pose p; { std::string seq = "";
		for(int i = 0; i < N; ++i) seq += "G";
		core::pose::make_pose_from_sequence(p,seq,*crs,false);
		set_pose_to_ideal(p);
	}
	
	fere_torsions(N,tor); // radians!
	refold_ros(p,tor,roscrd);
	refold_gpu(tor,N,gpucrd, cl,clmtor,clmout);
	dump_points_pdb(roscrd,N,"fere_ros.pdb");
	dump_points_pdb(gpucrd,N,"fere_gpu.pdb");	
	//cout << "fere rms: " << numeric::model_quality::calc_rms(array2vecs(roscrd,3*N),array2vecs(gpucrd,3*N)) << endl;;

	Real minrms=9e9,maxrms=-9e9;
	for(int i = 0; i < 10000; ++i) {
		random_torsions(N,tor); // radians!
		refold_ros(p,tor,roscrd);
		refold_gpu(tor,N,gpucrd, cl,clmtor,clmout);
		Real rms = numeric::model_quality::calc_rms(array2vecs(roscrd,3*N),array2vecs(gpucrd,3*N));
		minrms = min(rms,minrms);
		maxrms = max(rms,maxrms);
	}
	cout << "10,000 random 64-res refolds with " << minrms << " < rms < " << maxrms << endl;
  
}

int main(int argc, char *argv[]) {
	register_options();
  core::init(argc,argv);

	test_refold();
}
