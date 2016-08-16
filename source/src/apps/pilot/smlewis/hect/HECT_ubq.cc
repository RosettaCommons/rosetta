// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/smlewis/HECT/HECT_ubq.cc
/// @brief This app is intended for modeling the binding of a E2/E3/UBQ complex.  This is the "first app" in the series; it takes the starting structure and makes a few mutations and modifications.  The starting structures (PDBs 3JVZ and 3JW0) have the catalytic cysteines mutated out for stability; this protocol re-mutates to cysteine and then relaxes the ubiquitin tail/chemical linkage to account for the new bond parameters.  It uses constraints to maintain the position of the UBQ relative to E3.
/// @author Steven Lewis

// Unit Headers
#include <apps/pilot/smlewis/hect/HECT.hh>

//The headers are mostly transcluded from HECT.hh
#include <core/conformation/ResidueFactory.hh>

#include <basic/prof.hh>

//Auto Headers
#include <core/pose/util.hh>


//static basic::Tracer TR("apps.pilot.smlewis.HECT");


/// @brief HECT mover
class HECTUBQMover : public HECTMover {

	typedef HECTMover parent;

public:
	HECTUBQMover() : parent()
	{}

	/// @brief init_on_new_input system allows for initializing these details the first time apply() is called.  the job distributor will reinitialize the whole mover when the input changes (a freshly constructed mover, which will re-run this on first apply().
	virtual
	void
	init_on_new_input(core::pose::Pose & pose) {
		init_for_input_yet_ = true;

		//determine where the flexible regions are
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		bool const debug = option[basic::options::OptionKeys::debug].value();

		parse_options(pose); //held in parent class; sets up points of interest data

		//set up starting pose
		//first: mutate catalytic residues back to cysteine
		if(debug) pose.dump_pdb("test0.pdb");
		core::Real const chi1(pose.residue(e2_catalytic_res_).chi(1));

		using namespace core::chemical;
		using namespace core::conformation;
		ResidueTypeSetCAP fa_standard(ChemicalManager::get_instance()->residue_type_set(FA_STANDARD));
		core::conformation::Residue cyx_rsd( *(ResidueFactory::create_residue(fa_standard->name_map("CYX"))) );
		core::conformation::Residue cyz_rsd( *(ResidueFactory::create_residue(fa_standard->name_map("CYZ"))) );
		pose.replace_residue(e2_catalytic_res_, cyx_rsd, true); //we need to fix the fold tree to account for this
		pose.replace_residue(e3_catalytic_res_, cyz_rsd, true); //not a trivial change, but not fold-tree related
		if(debug) pose.dump_pdb("test1.pdb");

		//next, reset chi angle on new cysteine to match original serine as best as possible (this is mostly cosmetic, since this is variable later, but it makes the "test" structures much easier to examine since it prevents huge eclipsing)
		pose.set_chi(1, e2_catalytic_res_, chi1);
		if(debug) pose.dump_pdb("test2.pdb");

		//Prepare C-terminal glycine for modification by removing terminus variant (and extra oxygen!)
		core::pose::remove_upper_terminus_type_from_pose_residue(pose, ubqchain_end_);
		if(debug) pose.dump_pdb("test3.pdb");

		//set up the chemical bond - this handles the conn_id details
		pose.conformation().declare_chemical_bond(e2_catalytic_res_, "SG", ubqchain_end_, "C");
		if(debug) pose.dump_pdb("test4.pdb");

		//set up fold tree
		set_up_foldtree(pose);
		if(debug) pose.dump_pdb("test5.pdb");

		//set ideal geometry for new bond
		core::id::AtomID const e2_cys_S(pose.residue_type(e2_catalytic_res_).atom_index("SG"), e2_catalytic_res_);
		core::id::AtomID const UBQ_G_C(pose.residue_type(ubqchain_end_).atom_index("C"), ubqchain_end_);
		pose.conformation().insert_ideal_geometry_at_residue_connection( e2_catalytic_res_, 3 );
		if(debug) pose.dump_pdb("test6.pdb");

		//ok, there is an improper dihedral remaining - the O atom.  thanks to Andrew for helping with the math here!
		//the xyz of the carbonyl carbon and its bonded atoms (except O)
		core::Vector const & SGxyz(pose.residue(e2_catalytic_res_).atom("SG").xyz()),
			Cxyz(pose.residue(ubqchain_end_).atom("C").xyz()),
			CAxyz(pose.residue(ubqchain_end_).atom("CA").xyz());

		//the vector representing the two bonds from the C (except C=O)
		core::Vector const C_SG( (Cxyz-SGxyz).normalize() ),
			C_CA( (Cxyz-CAxyz).normalize());

		//copy the current distance to normalize the new vector
		core::id::AtomID const UBQ_G_CO(pose.residue_type(ubqchain_end_).atom_index("O"), ubqchain_end_);
		core::Real const COdist(pose.conformation().bond_length(UBQ_G_C, UBQ_G_CO));

		//new bond is in direction of sum of other two bond vectors, normalized, times old distance
		core::Vector const newCOpos( ((C_SG + C_CA).normalize() * COdist ) + Cxyz );
		pose.set_xyz(UBQ_G_CO, newCOpos);
		if(debug) pose.dump_pdb("test7.pdb");

		//repair the thioester dihedral to trans (it was not for ester; esters have little double bond character)
		core::id::AtomID const e2_cys_CB(pose.residue_type(e2_catalytic_res_).atom_index("CB"), e2_catalytic_res_);
		core::id::AtomID const UBQ_G_CA(pose.residue_type(ubqchain_end_).atom_index("CA"), ubqchain_end_);
		pose.conformation().set_torsion_angle( e2_cys_CB, e2_cys_S, UBQ_G_C, UBQ_G_CA, numeric::conversions::radians(180.0) );
		if(debug) pose.dump_pdb("test8.pdb");

		//setup of TaskFactory
		set_up_taskfactory();

		//setup MoveMap
		movemap_ = new core::kinematics::MoveMap;

		//initial state: all false
		//set UBQ tail flexible - last residue is not part of this; scorefunction mishandles (not a terminus, wrong rama)
		for( core::Size i(utail_start_); i<=ubqchain_end_-1; ++i) movemap_->set_bb(i, true);
		//do not set UBQ chemical linkage flexible - scorefunction mishandles
		//do not set hinge flexible for this processing
		//set jump NOT flexible
		movemap_->set_jump(false);

		//set up constraints
		ubq_constraints(pose);

		//make fragments
		std::string const seq(pose.sequence());
		fragset3mer_ubq_tail_ = make_frags(utail_start_, ubqchain_end_-1, seq);
		TR << "recovering memory?" << std::endl;
		// WARNING WARNING WARNING! THREAD UNSAFE!  WHY WOULD YOU THINK THIS IS A GOOD IDEA?
		core::fragment::picking_old::FragmentLibraryManager::get_instance()->clear_Vall();

		//build internal AtomID vector for TorsionDOFs later
		build_AtomID_vec();

		//save pose
		fixed_starting_pose_ = new core::pose::Pose(pose);
	}

	virtual ~HECTUBQMover(){};

	//no apply function needed - defined by parent

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new HECTUBQMover;
	}

};

typedef utility::pointer::owning_ptr< HECTUBQMover > HECTUBQMoverOP;

int main( int argc, char* argv[] )
{
	basic::prof_reset();
	register_options(); //in HECT header
	devel::init(argc, argv);
	protocols::jd2::JobDistributor::get_instance()->go(new HECTUBQMover);
	TR << "************************d**o**n**e**************************************" << std::endl;
	return 0;
}
