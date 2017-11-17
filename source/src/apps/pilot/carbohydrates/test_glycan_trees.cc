// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/moves/Mover.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/GlycanTree.hh>
#include <core/conformation/carbohydrates/GlycanNode.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/ResidueFactory.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/string_util.hh>


#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/carbohydrates/util.hh>

// Protocol Headers
#include <basic/Tracer.hh>


static basic::Tracer TR("glycan_info");


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::jd2;
		using namespace core::conformation::carbohydrates;
		using namespace core::pose::carbohydrates;
		using namespace core::import_pose;
		using namespace core::io::pdb;
		using namespace core::conformation;
		using namespace core::chemical;
		using namespace core::import_pose;

		devel::init( argc, argv );
		//register_options();




		core::pose::Pose pose_;
		core::pose::PoseOP man9_op_;


		pose_from_file(pose_, "/Users/jadolfbr/Documents/Rosetta/main/source/test/core/chemical/carbohydrates/gp120_2glycans_man5.pdb", PDB_file);

		pose_.delete_residue_slow(593);

		GlycanTreeSetCOP tree_set = pose_.glycan_tree_set();
		utility::vector1< GlycanTreeCOP > const trees =  tree_set->get_all_trees() ;

		std::string const man9_s( "a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
			"[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-b-D-GlcpNAc-(1->4)-b-D-GlcpNAc" );
		man9_op_ = core::pose::pose_from_saccharide_sequence( man9_s, "fa_standard", true, false ); //No need to idealize.

		std::cout << std::endl << "Deleting from 8" << std::endl;
		core::pose::PoseOP man9_copy = man9_op_->clone();
		//TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 8); //Deletes 9
		//TS_ASSERT_EQUALS(man9_copy->size(), 10);
		man9_copy->dump_pdb("man9_trim_at_8.pdb");

		std::cout << std::endl << "Deleting from 10" << std::endl;
		man9_copy = man9_op_->clone();
		//TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 10); //Deletes 11
		//TS_ASSERT_EQUALS(man9_copy->size(), 10);
		man9_copy->dump_pdb("man9_trim_at_10.pdb");

		std::cout << std::endl << "Deleting from 4" << std::endl;
		man9_copy = man9_op_->clone();
		//TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 4); //Deletes 6,5
		//TS_ASSERT_EQUALS(man9_copy->size(), 9);
		man9_copy->dump_pdb("man9_trim_at_4.pdb");

		std::cout << std::endl << "Deleting from 7" << std::endl;
		man9_copy = man9_op_->clone();
		//TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 7); //Deletes 9,8,11,10
		//TS_ASSERT_EQUALS(man9_copy->size(), 7);
		man9_copy->dump_pdb("man9_trim_at_7.pdb");

		std::cout << std::endl << "Deleting from 3" << std::endl;
		man9_copy = man9_op_->clone();
		//TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 3); //Deletes 9,8,11,10,6,5,4
		//TS_ASSERT_EQUALS(man9_copy->size(), 3);
		man9_copy->dump_pdb("man9_trim_at_3.pdb");

		std::cout << std::endl << "Deleting from 2" << std::endl;
		man9_copy = man9_op_->clone();
		//TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 2); ///Deletes 9,8,11,10,6,5,4,3
		//TS_ASSERT_EQUALS(man9_copy->size(), 2);
		man9_copy->dump_pdb("man9_trim_at_2.pdb");

		std::cout << std::endl << "Deleting from 1" << std::endl;
		man9_copy = man9_op_->clone();
		//TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 1); ///Deletes 9,8,11,10,6,5,4,3, 2
		//TS_ASSERT_EQUALS(man9_copy->size(), 1);
		man9_copy->dump_pdb("man9_trim_at_1.pdb");

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
