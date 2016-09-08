#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/AtomType.hh>

#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/import_pose/import_pose.hh>  //Read pdb in functionality
#include <core/pose/Pose.hh>        //Creation of pose object
#include <core/pose/util.hh>

#include <core/scoring/dna/setup.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <core/graph/Graph.hh>

#include <devel/dna/util_motif_loop.hh>
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>

static THREAD_LOCAL basic::Tracer tr( "devel.dna.util_motif_loop" );

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////Utility functions//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

bool hard_sphere_check( core::pose::Pose & pose, core::conformation::Residue & resi)
{
	using namespace numeric;
	for ( core::Size res(1); res <= pose.size(); ++res ) {
		if ( ! pose.residue(res).is_DNA() ) continue;

		core::conformation::Residue this_resi(pose.residue(res));

		for ( core::Size atm2(1); atm2 <= this_resi.nbr_atom(); ++atm2 ) {
			for ( core::Size atm1(1); atm1 <= resi.nbr_atom(); ++atm1 ) {
				core::Real collision_dist( resi.atom_type(atm1).lj_radius() + this_resi.atom_type(atm2).lj_radius());
				core::Real dist( resi.xyz(atm1).distance(this_resi.xyz(atm2)) );

				if ( dist < collision_dist ) return false;



			}
		}
	}
	return true;
}

core::conformation::Residue new_residue( std::string threeLetterAA )
{
	///Give it a 3 letter AA code and get back a new residue object of that type
	core::conformation::ResidueOP new_resi = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map(threeLetterAA) );
	return *new_resi;
}

core::pack::rotamer_set::RotamerSetOP get_aa_rotset( std::string & AA)
{
	using namespace core::pack::rotamer_set;
	using namespace core::pack::task;
	using namespace core::scoring;

	core::pose::Pose rotamer_building_site;
	core::conformation::Residue resi( new_residue(AA) );
	rotamer_building_site.append_residue_by_jump(resi, (core::Size) 0);

	core::pose::add_upper_terminus_type_to_pose_residue(rotamer_building_site,1);  //to allow bb independent rotamer building
	core::pose::add_lower_terminus_type_to_pose_residue(rotamer_building_site,1);

	ScoreFunction scorefxn;
	scorefxn(rotamer_building_site);

	core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task(rotamer_building_site);

	dummy_task->nonconst_residue_task(1).restrict_to_repacking();
	dummy_task->nonconst_residue_task(1).or_include_current(false);
	dummy_task->nonconst_residue_task(1).or_fix_his_tautomer(true);

	core::graph::GraphOP dummy_png = core::pack::create_packer_graph(rotamer_building_site, scorefxn, dummy_task);


	RotamerSetFactory rsf;
	RotamerSetOP rotset ( rsf.create_rotamer_set( rotamer_building_site.residue(1) ) );
	rotset->set_resid( 1 );
	rotset->build_rotamers(rotamer_building_site, scorefxn, *dummy_task, dummy_png);

	return  rotset;
}


std::map< std::string, core::pack::rotamer_set::RotamerSetOP >  build_all_aa_rotamers()
{
	//Builds the rotamers of each amino acid

	using namespace core::pack::rotamer_set;
	using namespace core::pack::task;
	using namespace core::scoring;

	std::cout << "BUILDING ROTAMERS FOR EACH AMINO ACID\n";

	std::map< std::string, RotamerSetOP > aa_rotamers;

	for ( core::Size i(1); i <= core::chemical::num_canonical_aas; ++i ) {
		std::string name(core::chemical::name_from_aa(core::chemical::AA(i)));
		if ( name == "GLY" ) continue;


		core::pose::Pose rotamer_building_site;
		core::conformation::Residue resi( new_residue(name) );
		rotamer_building_site.append_residue_by_jump(resi, (core::Size) 0);

		core::pose::add_upper_terminus_type_to_pose_residue(rotamer_building_site,1);  //to allow bb independent rotamer building
		core::pose::add_lower_terminus_type_to_pose_residue(rotamer_building_site,1);

		ScoreFunction scorefxn;
		scorefxn(rotamer_building_site);

		core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task(rotamer_building_site);

		dummy_task->nonconst_residue_task(1).restrict_to_repacking();
		dummy_task->nonconst_residue_task(1).or_include_current(false);
		dummy_task->nonconst_residue_task(1).or_fix_his_tautomer(true);

		core::graph::GraphOP dummy_png = core::pack::create_packer_graph(rotamer_building_site, scorefxn, dummy_task);


		RotamerSetFactory rsf;
		RotamerSetOP rotset ( rsf.create_rotamer_set( rotamer_building_site.residue(1) ) );
		rotset->set_resid( 1 );
		rotset->build_rotamers(rotamer_building_site, scorefxn, *dummy_task, dummy_png);

		aa_rotamers[name]= rotset;

	} //End AA Loop

	/* core::pose::Pose outpose;
	RotamerSetOP tmp = aa_rotamers["LYS"];
	for ( core::Size i(1); i != tmp->num_rotamers(); ++i)
	{
	outpose.append_residue_by_bond( *(tmp->nonconst_rotamer(i)) );
	}
	outpose.dump_pdb("built.pdb");
	*/  return aa_rotamers;
}

