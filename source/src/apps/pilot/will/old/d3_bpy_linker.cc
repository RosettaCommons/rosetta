// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#define CONTACT_D2 20.25
#define CONTACT_TH 30
#define NSS 1812 // 672 8192 17282
#define MIN_HELEX_RES 20
#define MAX_CYS_RES 3
#define MAX_NRES 200

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <numeric/xyzVector.hh>
#include <protocols/sic_dock/SICFast.hh>
#include <protocols/sic_dock/designability_score.hh>

#include <apps/pilot/will/will_util.ihh>

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

using core::id::AtomID;
using basic::options::option;
using core::pose::Pose;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::rotation_matrix_degrees;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_file;
using core::kinematics::Stub;

static basic::Tracer TR( "D3BPYDIG" );
static core::io::silent::SilentFileData sfd;

// OPT_1GRP_KEY( FileVector, cxdock, bench2 )
// OPT_1GRP_KEY( FileVector, cxdock, bench3 )
// OPT_1GRP_KEY( FileVector, cxdock, bench4 )
// OPT_1GRP_KEY( FileVector, cxdock, bench5 )
// OPT_1GRP_KEY( FileVector, cxdock, bench6 )
// OPT_1GRP_KEY( IntegerVector, cxdock, nfold )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// NEW_OPT( cxdock::bench2 , "benchmark dimers (must be aligned-z)" , "" );
	// NEW_OPT( cxdock::bench3 , "benchmark trimers (must be aligned-z)" , "" );
	// NEW_OPT( cxdock::bench4 , "benchmark tetra (must be aligned-z)" , "" );
	// NEW_OPT( cxdock::bench5 , "benchmark penta (must be aligned-z)" , "" );
	// NEW_OPT( cxdock::bench6 , "benchmark hexa (must be aligned-z)" , "" );
	// NEW_OPT( cxdock::nfold  , "C syms to test", 1 );
}

#ifdef USE_OPENMP
#include <omp.h>
#endif
int num_threads() {
#ifdef USE_OPENMP
	return omp_get_max_threads();
#else
	return 1;
#endif
}
int thread_num() {
#ifdef USE_OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}

// info about a "hit"
struct Hit {
	int iss,irt,sym;
	Real score;
	core::kinematics::Stub stub;
	Hit(int _iss, int _irt, Real _score, int _sym) : iss(_iss),irt(_irt),sym(_sym),score(_score) {}
};
bool cmpHit(Hit i,Hit j) { return i.score > j.score; }


class SlideTask;
typedef utility::pointer::owning_ptr< SlideTask > SlideTaskOP;
typedef utility::pointer::owning_ptr< SlideTask const > SlideTaskCOP;


class SlideTask : public utility::pointer::ReferenceCount {
public:
	SlideTask(
		core::pose::Pose const & pose,
		vector1<Vec> const & ssamp
	): pose_(pose),
		ssamp_(ssamp)
	{}
	core::pose::Pose const & pose () const { return pose_; }
	vector1<Vec>     const & ssamp() const { return ssamp_; }
private:
	core::pose::Pose const & pose_;
	vector1<Vec> const & ssamp_;
	std::vector<protocols::sic_dock::SICFastCOP> sics_;

};

void
make_nmer(
	SlideTaskOP stask
){

}


// dock pose against rotated version of itself
// rotate all different ways by looping over rotation axes on unit sphere point
// samples in ssamp, then loopint over rotation magnitudes from 1-180°
// picks cases with most CB contacts in helices
void
dock(
	Pose const & init_pose,
	std::string const & fn,
	vector1<Vec> const & ssamp
){
	using namespace basic::options;
	using namespace OptionKeys;

	// set up one SIC per thread
	std::vector<protocols::sic_dock::SICFast*> sics_;
	sics_.resize(num_threads());
	for ( int i = 0; i < num_threads(); ++i ) sics_[i] = new protocols::sic_dock::SICFast;
	for ( int i = 0; i < num_threads(); ++i ) sics_[i]->init(init_pose);

	//store initial coords & angle samples 1-180 (other distribution of samps would be better...)
	vector1<double> asamp; for ( Real i = 1.0; i < 180.0; i += 5.0 ) asamp.push_back(i);

	Size nfold = 3;
	Mat R2 = rotation_matrix_degrees(Vec(1,0,0),180.0);
	Mat R3 = rotation_matrix_degrees(Vec(1,0,0),120.0);
	Mat RswapXZ = rotation_matrix_degrees(Vec(1,0,1),180.0);
	vector1<Hit> hits;

	// loop over axes on sphere points & rotaton amounts 1-180
#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
#endif
	for ( int iss = 1; iss <= (int)ssamp.size(); ++iss ) {
		if ( iss%1000==0 ) TR << iss << " of " << NSS << std::endl;
		for ( int irt = 1; irt <= (int)asamp.size(); ++irt ) {
			Mat R = rotation_matrix_degrees(ssamp[iss],(Real)asamp[irt]);
			Stub x1(    R , Vec(0,0,0) );
			Stub x2( R3*R , Vec(0,0,0) );
			Real score;
			Real t = sics_[thread_num()]->slide_into_contact(x1,x2,Vec(0,0,1),score);
			if ( score >= CONTACT_TH ) {
				Hit h(iss,irt,score,3);
				h.stub = x1;
				h.stub.v += Vec(0,fabs(t)/2.0/tan(numeric::constants::d::pi/(Real)h.sym),t/2.0);
				h.stub.M = RswapXZ * h.stub.M;
				h.stub.v = RswapXZ * h.stub.v;
				h.score = score;
#ifdef USE_OPENMP
				#pragma omp critical
#endif
				hits.push_back(h);
			}
		}
	}

	Size nstruct = basic::options::option[basic::options::OptionKeys::out::nstruct]();
	std::string outdir = basic::options::option[basic::options::OptionKeys::out::file::o]()+"/";

	// dump top results
	std::sort(hits.begin(),hits.end(),cmpHit);
	for ( int i = 1; i <= (int)min(nstruct,hits.size()); ++i ) {
		Hit & h(hits[i]);
		std::string tag = utility::file_basename(fn);
		if ( tag.substr(tag.size()-3)==".gz" ) tag = tag.substr(0,tag.size()-3);
		if ( tag.substr(tag.size()-4)==".pdb" ) tag = tag.substr(0,tag.size()-4);
		tag += "_C"+ObjexxFCL::string_of(h.sym)+"_"+ObjexxFCL::string_of(i)+".pdb";
		cout << "RESULT " << h.sym << " " << h.iss << " " << NSS  << " " << h.irt << " " << h.score << " " << tag << endl;
		Pose tmp(init_pose);
		xform_pose(tmp,h.stub);
		tmp.dump_pdb( option[out::file::o]()+"/"+tag);
	}

}


void
read_sphere(
	vector1<Vec> & ssamp
){
	izstream is;
	if ( !basic::database::open(is,"sampling/spheres/sphere_"+str(NSS)+".dat") ) {
		utility_exit_with_message("can't open sphere data");
	}
	for ( int i = 1; i <= NSS; ++i ) {
		double x,y,z;
		is >> x >> y >> z;
		ssamp[i] = Vec(x,y,z);
	}
	is.close();
}


int main(int argc, char *argv[]) {

	try {

		register_options();
		devel::init(argc,argv);
		using namespace basic::options::OptionKeys;

		vector1<Vec> ssamp(NSS); read_sphere(ssamp);

		string fn = option[in::file::s]()[1];
		TR << "searching d3_bpy_linker " << fn << std::endl;

		Pose pnat;
		core::import_pose::pose_from_file(pnat,fn, core::import_pose::PDB_file);
		core::scoring::dssp::Dssp dssp(pnat);
		dssp.insert_ss_into_pose(pnat);

		Vec cen = center_of_geom(pnat,1,pnat.size());
		trans_pose(pnat,-cen);

		dock(pnat,fn,ssamp);


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

