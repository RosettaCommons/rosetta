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

#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>


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
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_file;
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
test_replicate_xform(){
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using namespace core::kinematics;

	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, option[in::file::s]()[1] , core::import_pose::PDB_file);

	for(Size lstart =        2; lstart < pose.n_residue()-1; ++lstart ){
	for(Size lstop  = lstart+1; lstop  < pose.n_residue()  ; ++lstop  ){

		// std::cout << "testing... " << lstart << " " << lstop << std::endl;

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
	std::cout << "all pairs OK" << std::endl;

}

void
test_lh_counts(core::pose::Pose const & pose){
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using namespace core::kinematics;

	utility::vector1<core::Size> loopsizes;
	loopsizes.push_back(3);
	loopsizes.push_back(4);
	loopsizes.push_back(5);
	loopsizes.push_back(6);
	loopsizes.push_back(7);
	loopsizes.push_back(8);
	loopsizes.push_back(9);
	loopsizes.push_back(10);

	//======= load the library with the sizes you need
	LoopHashLibraryOP loop_hash_library = new LoopHashLibrary ( loopsizes, 1, 0 );
	loop_hash_library->load_mergeddb();

	std::cout << "start counts" << std::endl;
	int ncalls=0, hits0=0, hits1=0, hits2=0, hits3=0, hits4=0;
	for(int i = 1; i <= 100; ++i){
		hits0=0, hits1=0, hits2=0, hits3=0, hits4=0;
		for(Size lstart =        2; lstart < pose.n_residue()-1; ++lstart ){
		for(Size lstop  = lstart+1; lstop  < pose.n_residue()  ; ++lstop  ){
			Real6 loop_transform;
			get_rt_over_leap_without_foldtree_bs( pose, lstart, lstop, loop_transform );
			ncalls++;
			Size count0=0, count1=0, count2=0, count3=0, count4=0;
			for(utility::vector1<core::Size>::const_iterator i = loopsizes.begin(); i != loopsizes.end(); ++i){
				count0 += loop_hash_library->gethash(*i).radial_count(0,loop_transform);
				count1 += loop_hash_library->gethash(*i).radial_count(1,loop_transform);
				count2 += loop_hash_library->gethash(*i).radial_count(2,loop_transform);
				// count3 += loop_hash_library->gethash(*i).radial_count(3,loop_transform);
				// count4 += loop_hash_library->gethash(*i).radial_count(4,loop_transform);
			}
			hits0 += count0;
			hits1 += count1;
			hits2 += count2;
			hits3 += count3;
			hits4 += count4;
			// if(count) std::cout << lstart << " " << lstop << " " << count << std::endl;
		}}
	}
	std::cout << "done counts " << ncalls << " " << hits0 << " " << hits1 << " " << hits2 << " " << hits3 << " " << hits4 << std::endl;

	utility_exit_with_message("exit");

}

int main(int argc, char *argv[]) {

	try {

	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using namespace core::kinematics;

	register_options();
	devel::init(argc,argv);

	// test_replicate_xform();
	// test_lh_counts(pose);

	std::string fname = option[in::file::s]()[1];
	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, fname , core::import_pose::PDB_file);
	core::pose::remove_upper_terminus_type_from_pose_residue(pose,2);
	core::pose::remove_lower_terminus_type_from_pose_residue(pose,3);


	//======= Define the loop size database you need to collect ==========
	utility::vector1<core::Size> loopsizes;
	// loopsizes.push_back(3);
	// loopsizes.push_back(4);
	// loopsizes.push_back(5);
	loopsizes.push_back(6);
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

	for(utility::vector1<core::Size>::const_iterator i = loopsizes.begin(); i != loopsizes.end(); ++i){
		Size loopsize = *i;

		for(Size lstart = 2; lstart < pose.n_residue()-loopsize; ++lstart ){
			Size lstop = lstart+loopsize;

			Real6 loop_transform;
			get_rt_over_leap_without_foldtree_bs( pose, lstart, lstop, loop_transform );

			BackboneSegment backbone_;

			LoopHashMap &hashmap = loop_hash_library->gethash(loopsize);

			//hash indices are stored in this vector for later access
			std::vector < core::Size > leap_index_list;

			Size radius;
			for(radius=0; radius < 5; ++radius){
				hashmap.radial_lookup( radius, loop_transform, leap_index_list );
				if(leap_index_list.size() > 0) break;
			}
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
			Pose tmp;
			{
				core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
				tmp.append_residue_by_jump(pose.residue(lstart-1),1);
				tmp.append_residue_by_bond(pose.residue(lstart    ));
				for(Size i = lstart+1; i <= lstop; ++i){
					core::conformation::ResidueOP new_rsd( NULL );
					new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map("ALA") );
					cout << "apending residue " << new_rsd->name() << std::endl;
					tmp.append_residue_by_bond( *new_rsd, true );
					tmp.set_phi  ( tmp.n_residue(), 180.0 );
					tmp.set_psi  ( tmp.n_residue(), 180.0 );
					tmp.set_omega( tmp.n_residue()-1, 180.0 );
				}
				tmp.dump_pdb("test.pdb");
			}

			for ( std::vector< BackboneSegment >::iterator i = bs_vec_.begin(), ie = bs_vec_.end(); i != ie; ++i) {
				std::vector<core::Real> phi   = i->phi();
				std::vector<core::Real> psi   = i->psi();
				std::vector<core::Real> omega = i->omega();
				Size seg_length = (*i).length();
				for ( Size i = 0; i < seg_length; i++){
					Size ires = lstart+i;  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
					if (ires > pose.total_residue() ) break;
					tmp.set_phi  ( ires, phi[i]  );
					tmp.set_psi  ( ires, psi[i]  );
					tmp.set_omega( ires, omega[i]);
				}
				tmp.dump_pdb("test_lh_"+str(lstart)+"_"+str(lstop)+"_"+str(radius)+".pdb");
			}
		}

	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

