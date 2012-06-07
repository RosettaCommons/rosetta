#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/HomogeneousTransform.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//loophash
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashMap.hh>

#include <apps/pilot/will/will_util.ihh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/tree/Atom_.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <core/kinematics/tree/JumpAtom.hh>

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

void
testfunc( core::pose::Pose const & pose,
	int a1, int r1, int a2, int r2, int a3, int r3,
	int a4, int r4, int a5, int r5, int a6, int r6, int a7, int r7,
	numeric::geometry::hashing::Real6 &rt_6 
){
		using namespace core::kinematics;
		core::id::StubID id1( AtomID(a1,r1), AtomID(a2,r2), AtomID(a3,r3) );
		core::id::StubID id2( AtomID(a7,r7), AtomID(a4,r4), AtomID(a5,r5), AtomID(a6,r6) );	
		RT rt2 = pose.conformation().get_stub_transform(id1,id2);
		numeric::HomogeneousTransform< Real > ht2( rt2.get_rotation(), rt2.get_translation() );
		numeric::xyzVector < Real > euler_angles2 =  ht2.euler_angles_rad();
		if( fabs(rt_6[1]-rt2.get_translation().x()) < 0.01 &&
		    fabs(rt_6[2]-rt2.get_translation().y()) < 0.01 &&
		    fabs(rt_6[3]-rt2.get_translation().z()) < 0.01 && (
		    	((fabs(rt_6[4]-euler_angles2.x()*180.0/numeric::constants::d::pi) < 0.01 || fabs(rt_6[4]-euler_angles2.x()*180.0/numeric::constants::d::pi-180.0) < 0.01) &&
		    	 (fabs(rt_6[5]-euler_angles2.y()*180.0/numeric::constants::d::pi) < 0.01 || fabs(rt_6[5]-euler_angles2.y()*180.0/numeric::constants::d::pi-180.0) < 0.01) &&
		    	 (fabs(rt_6[6]-euler_angles2.z()*180.0/numeric::constants::d::pi) < 0.01 || fabs(rt_6[6]-euler_angles2.z()*180.0/numeric::constants::d::pi-180.0) < 0.01) ))
		    )
		{
			std::cout << "HIT ";

			std::cout << rt2.get_translation().x() << " ";
			std::cout << rt2.get_translation().y() << " ";
			std::cout << rt2.get_translation().z() << " ";
			std::cout << euler_angles2.x()*180.0/numeric::constants::d::pi << " ";
			std::cout << euler_angles2.y()*180.0/numeric::constants::d::pi << " ";
			std::cout << euler_angles2.z()*180.0/numeric::constants::d::pi << "     ";

			std::cout << a1 << "," << r1 << "  "
			          << a2 << "," << r2 << "  "
			          << a3 << "," << r3 << "  "
			          << a4 << "," << r4 << "  "
			          << a5 << "," << r5 << "  "
			          << a6 << "," << r6 << "  "
			          << a7 << "," << r7 << std::endl;
		}
}

void
test_replicate_xform(){
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using namespace core::kinematics;

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, option[in::file::s]()[1] );

	for(Size lstart =        2; lstart < pose.n_residue()-1; ++lstart ){
	for(Size lstop  = lstart+1; lstop  < pose.n_residue()  ; ++lstop  ){

		std::cout << "testing... " << lstart << " " << lstop << std::endl;

		Real6 lxf_std, lxf_mine;
		get_rt_over_leap_fast( pose, lstart, lstop, lxf_std );
		get_rt_over_leap_without_foldtree_bs( pose, lstart, lstop, lxf_mine );

		for(int i = 1; i <= 6; ++i){
			if( fabs(lxf_mine[i]-lxf_std[i]) > 0.0001 ){
				std::cout << "ERROR! " << lstart << " " << lstop << std::endl;
				std::cout << "std:  " <<  lxf_std[1] << " ";
				std::cout << "std:  " <<  lxf_std[2] << " ";
				std::cout << "std:  " <<  lxf_std[3] << " ";
				std::cout << "std:  " <<  lxf_std[4] << " ";
				std::cout << "std:  " <<  lxf_std[5] << " ";
				std::cout << "std:  " <<  lxf_std[6] << std::endl;;
				std::cout << "mine: " << lxf_mine[1] << " ";
				std::cout << "mine: " << lxf_mine[2] << " ";
				std::cout << "mine: " << lxf_mine[3] << " ";
				std::cout << "mine: " << lxf_mine[4] << " ";
				std::cout << "mine: " << lxf_mine[5] << " ";
				std::cout << "mine: " << lxf_mine[6] << std::endl;;
			}
		}


	}}

}

int main(int argc, char *argv[]) {
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using namespace core::kinematics;

	register_options();
	devel::init(argc,argv);

	test_replicate_xform();

	std::string fname = option[in::file::s]()[1];
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, fname );



	Real6 loop_transform;
	get_rt_over_leap_fast( pose, 4, 10, loop_transform );

	assert(false);

	cout << "DONE" << std::endl;

/*

	//======= Define the loop size database you need to collect ==========
	utility::vector1<core::Size> loopsizes;
	// loopsizes.push_back(3);
	// loopsizes.push_back(4);
	loopsizes.push_back(5);
	// loopsizes.push_back(6);
	// loopsizes.push_back(7);
	// loopsizes.push_back(8);
	// loopsizes.push_back(9);
	// loopsizes.push_back(10);

	//======= load the library with the sizes you need
	LoopHashLibraryOP loop_hash_library = new LoopHashLibrary ( loopsizes, 1, 0 );
	loop_hash_library->load_mergeddb();


	//===========Extrac the backbone torsions from the library==========
	//===========BackboneSegment is a container of torsions ======
	const BackboneDB & bbdb_ = loop_hash_library->backbone_database();


	std::cout << "start counts" << std::endl;
	int ncalls = 0;
	for(int i = 1; i <= 10; ++i){
		for(Size lstart =        2; lstart < pose.n_residue()-2; ++lstart ){
		for(Size lstop  = lstart+2; lstop  < pose.n_residue()  ; ++lstop  ){
			Real6 loop_transform;
			get_rt_over_leap_fast( pose, lstart, lstop, loop_transform );
			ncalls++;
			Size count = 0;
			for(utility::vector1<core::Size>::const_iterator i = loopsizes.begin(); i != loopsizes.end(); ++i){
				count += loop_hash_library->gethash(*i).radial_count(1,loop_transform);
			}
			// if(count) std::cout << lstart << " " << lstop << " " << count << std::endl;
		}}
	}
	std::cout << "done counts " << ncalls << std::endl;

	utility_exit_with_message("exit");

	for(utility::vector1<core::Size>::const_iterator i = loopsizes.begin(); i != loopsizes.end(); ++i){
		Size loopsize = *i;

		for(Size lstart = 2; lstart < pose.n_residue()-loopsize; ++lstart ){
			Size lstop = lstart+loopsize;
				Real6 loop_transform;
				get_rt_over_leap_fast( pose, lstart, lstop, loop_transform );

			BackboneSegment backbone_;

			LoopHashMap &hashmap = loop_hash_library->gethash(loopsize);

			//hash indices are stored in this vector for later access
			std::vector < core::Size > leap_index_list;

			hashmap.radial_lookup( 0, loop_transform, leap_index_list );
			if( leap_index_list.size() == 0 ) {
				std::cout << "FAIL " << fname << " " << lstart << " " << lstop << std::endl;
			}

			//here you load all the fragments to a vector
			std::vector < BackboneSegment > bs_vec_;
			for( std::vector < core::Size >::const_iterator itx = leap_index_list.begin(); itx != leap_index_list.end(); ++itx ){
				core::Size bb_index = *itx;
				LeapIndex cp = hashmap.get_peptide( bb_index );
				bbdb_.get_backbone_segment( cp.index, cp.offset , loopsize , backbone_ );
				bs_vec_.push_back( backbone_ );
			}
			//=============Apply the fragments =======
			Pose tmp(pose);
			// FoldTree f;
			// ObjexxFCL::FArray2D<int> jump(2,1); jump(1,1)=1; jump(2,1)=lstop;
			// ObjexxFCL::FArray1D<int> cuts(1,1); cuts(1)=lstop-1;
			// f.tree_from_jumps_and_cuts( tmp.n_residue(), 1, jump, cuts, 1 );
			// f.reorder(1);
			// tmp.fold_tree(f);

			for ( std::vector< BackboneSegment >::iterator i = bs_vec_.begin(), ie = bs_vec_.end(); i != ie; ++i) {
				std::vector<core::Real> phi = (*i).phi();
				std::vector<core::Real> psi = (*i).psi();
				std::vector<core::Real> omega = (*i).omega();
				Size seg_length = (*i).length();
				for ( Size i = 0; i < seg_length; i++){
					Size ires = lstart+i;  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
					if (ires > pose.total_residue() ) break;
					tmp.set_phi  ( ires, phi[i]);
					tmp.set_psi  ( ires, psi[i]);
					tmp.set_omega( ires, omega[i]);
				}
				tmp.dump_pdb("test_lh_"+str(lstart)+"_"+str(lstop)+"_"+str(0)+".pdb");
			}
		}

	}
*/
}

