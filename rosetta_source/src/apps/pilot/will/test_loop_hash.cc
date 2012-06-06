#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/HomogeneousTransform.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <protocols/sic_dock/designability_score.hh>

//loophash
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashMap.hh>

#include <apps/pilot/will/will_util.ihh>


typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

using core::id::AtomID;
using basic::options::option;
using namespace basic::options::OptionKeys;
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
using ObjexxFCL::fmt::F;
using ObjexxFCL::fmt::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;
using core::kinematics::Stub;

#define MAX_CYS_RES 0
#define MAX_NRES 1000

OPT_1GRP_KEY( Integer , dstat, min_seq_sep )
OPT_1GRP_KEY( File    , dstat, make_hist   )
OPT_1GRP_KEY( File   , dstat, test_score  )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( dstat::min_seq_sep    ,"",  20 );
	NEW_OPT( dstat::make_hist      ,"",  "" );
	NEW_OPT( dstat::test_score     ,"",  "" );
}

int main(int argc, char *argv[]) {
	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using namespace core::kinematics;

	Size lstart = 2, lstop = 5;


	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, option[in::file::s]()[1] );


	Real6 loop_transform;
	get_rt_over_leap_fast( pose, lstart, lstop, loop_transform );
	// std::cout << loop_transform[1] << " ";
	// std::cout << loop_transform[2] << " ";
	// std::cout << loop_transform[3] << " ";
	// std::cout << loop_transform[4] << " ";
	// std::cout << loop_transform[5] << " ";
	// std::cout << loop_transform[6] << std::endl;

	// Pose tmp;
	// tmp.append_residue_by_jump(pose.residue(lstart),1);
	// tmp.append_residue_by_jump(pose.residue(lstop),1);
	// RT rt = tmp.jump(1).rt();
	// numeric::HomogeneousTransform< Real > ht( rt.get_rotation() , rt.get_translation() );
	// numeric::xyzVector < Real > euler_angles =  ht.euler_angles_rad();

	// std::cout << rt.get_translation().x() << " ";
	// std::cout << rt.get_translation().y() << " ";
	// std::cout << rt.get_translation().z() << " ";
	// std::cout << euler_angles.x()*180/3.1415 << " ";
	// std::cout << euler_angles.y()*180/3.1415 << " ";
	// std::cout << euler_angles.z()*180/3.1415 << std::endl;


	// utility_exit_with_message("FOO");


	//======= Define the loop size database you need to collect ==========
	utility::vector1<core::Size> loopsizes;
	loopsizes.push_back(3);
	loopsizes.push_back(4);
	loopsizes.push_back(5);
	loopsizes.push_back(6);
	loopsizes.push_back(7);
	loopsizes.push_back(8);

	//======= load the library with the sizes you need
	LoopHashLibraryOP loop_hash_library = new LoopHashLibrary ( loopsizes, 1, 0 );

	loop_hash_library->load_mergeddb();

	//===========Extrac the backbone torsions from the library==========
	//===========BackboneSegment is a container of torsions ======

	const BackboneDB & bbdb_ = loop_hash_library->backbone_database();

	for(utility::vector1<core::Size>::const_iterator i = loopsizes.begin(); i != loopsizes.end(); ++i){
		Size loopsize = *i;

		BackboneSegment backbone_;

		LoopHashMap &hashmap = loop_hash_library->gethash(loopsize);

		//hash indices are stored in this vector for later access
		std::vector < core::Size > leap_index_list;

		std::cout << "radius = ";
		for (Size radius = 0; radius < 2 ; radius++ ){
			hashmap.radial_lookup( radius, loop_transform, leap_index_list );
			std::cout << radius << "... ";
			if (leap_index_list.size() < 1000){ //making sure at least harvest 1000 segment to build.
				continue;
			} else {
				break;
			}
		}

		std::cout << "collected " << leap_index_list.size() << " fragments at loopsize: " << loopsize << std::endl;

		//here you load all the fragments to a vector
		std::vector < BackboneSegment > bs_vec_;
		for( std::vector < core::Size >::const_iterator itx = leap_index_list.begin(); itx != leap_index_list.end(); ++itx ){
			core::Size bb_index = *itx;
			LeapIndex cp = hashmap.get_peptide( bb_index );
			bbdb_.get_backbone_segment( cp.index, cp.offset , loopsize , backbone_ );
			bs_vec_.push_back( backbone_ );
		}

		//=============Apply the fragments =======

		int count = 0;
		for ( std::vector< BackboneSegment >::iterator i = bs_vec_.begin(), ie = bs_vec_.end(); i != ie; ++i) {
			std::vector<core::Real> phi = (*i).phi();
			std::vector<core::Real> psi = (*i).psi();
			std::vector<core::Real> omega = (*i).omega();
			Size seg_length = (*i).length();

			for ( Size i = 0; i < seg_length; i++){
				Size ires = lstart+i;  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
				if (ires > pose.total_residue() ) break;
				pose.set_phi  ( ires, phi[i]);
				pose.set_psi  ( ires, psi[i]);
				pose.set_omega( ires, omega[i]);
			}
			pose.dump_pdb("test_lh_"+str(loopsize)+"_"+str(++count)+".pdb");
		}
	}
}

