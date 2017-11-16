// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/smlewis/HECT/HECT_all.cc
/// @brief This app is intended for modeling the binding of a E2/E3/UBQ complex.  This is the "second app" in the series - see also HECT_ubq.  This application will take a starting structure prepared by HECT_all and remodel the complex.  The remodeling includes: backbone remodeling of the E3 linker region, the UBQ tail, the chemical linkage, and sidechain repacking at the variable interfaces.  This code constrains the UBQ/E3 interface based on mutational data in (Kamadurai HB, Souphron J, Scott DC, Duda DM, Miller DJ, Stringer D, Piper RC, Schulman BA.  Insights into ubiquitin transfer cascades from a structure of a UbcH5B approximately ubiquitin-HECT(NEDD4L) complex.  Mol Cell. 2009 Dec 25;36(6):1095-102) and the UBQtail/E3 catalytic C to attempt to generate catalytic-like conformations.
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
class HECTAllMover : public HECTMover {

	typedef HECTMover parent;

public:
	HECTAllMover() : parent() {}

	/// @brief init_on_new_input system allows for initializing these details the first time apply() is called.  the job distributor will reinitialize the whole mover when the input changes (a freshly constructed mover, which will re-run this on first apply().
	virtual
	void
	init_on_new_input(core::pose::Pose & pose) {
		init_for_input_yet_ = true;

		//determine where the flexible regions are
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		bool const debug = option[basic::options::OptionKeys::debug].value();

		parse_options(pose); //held in parent class; sets up points-of-interest data

		//set up starting pose
		//first: mutate catalytic residues back to cysteine
		if ( debug ) pose.dump_pdb("test0.pdb");

		using namespace core::chemical;
		using namespace core::conformation;
		ResidueTypeSetCAP fa_standard(ChemicalManager::get_instance()->residue_type_set(FA_STANDARD));
		core::conformation::Residue cyz_rsd( *(ResidueFactory::create_residue(fa_standard->name_map("CYZ"))) );

		//this should not be necessary, but it seems to be!
		pose.replace_residue(e3_catalytic_res_, cyz_rsd, true); //not a trivial change, but not fold-tree related
		if ( debug ) pose.dump_pdb("test1.pdb");

		//Mini will put this terminus on automatically, so we need to kill it
		core::pose::remove_upper_terminus_type_from_pose_residue(pose, ubqchain_end_);
		if ( debug ) pose.dump_pdb("test2.pdb");

		//make the bond?
		pose.conformation().declare_chemical_bond(e2_catalytic_res_, "SG", ubqchain_end_, "C");
		if ( debug ) pose.dump_pdb("test3.pdb");

		//set up fold tree
		set_up_foldtree(pose);
		if ( debug ) pose.dump_pdb("test4.pdb");

		//setup of TaskFactory
		set_up_taskfactory();

		//setup MoveMap
		movemap_ = new core::kinematics::MoveMap;

		//initial state: all false
		//set UBQ tail flexible - last residue is not part of this; scorefunction mishandles (not a terminus, wrong rama)
		for ( core::Size i(utail_start_);    i<=ubqchain_end_-1;    ++i ) movemap_->set_bb(i, true);
		//do not set UBQ chemical linkage flexible - scorefunction mishandles
		//set hinge flexible
		for ( core::Size i(e3_hinge_start_); i<=e3_hinge_stop_;  ++i ) movemap_->set_bb(i, true);
		//set jump NOT flexible
		movemap_->set_jump(false);

		//set up constraints
		ubq_constraints(pose);

		//the catalytic constraint we want is a weak constraint between atom C of the carbonyl and SG of the cysteine

		core::id::AtomID const UBQ_G_C(pose.residue_type(ubqchain_end_).atom_index("C"), ubqchain_end_);
		core::id::AtomID const E3_cat_SG(pose.residue_type(e3_catalytic_res_).atom_index("SG"), e3_catalytic_res_);
		using namespace core::scoring::constraints;
		using core::id::AtomID;
		pose.add_constraint( new AtomPairConstraint( E3_cat_SG, UBQ_G_C, new BoundFunc(0, 2, option[catalytic_cst_sd].value(), /*0.5,*/ "catalyze")));

		//make fragments
		std::string const seq(pose.sequence());
		fragset3mer_e3_hinge_ = make_frags(e3_hinge_start_, e3_hinge_stop_, seq);
		fragset3mer_ubq_tail_ = make_frags(utail_start_, ubqchain_end_-1, seq);

		TR << "recovering memory?" << std::endl;
		// WARNING WARNING WARNING! THREAD UNSAFE!  WHY WOULD YOU THINK THIS IS A GOOD IDEA?
		core::fragment::picking_old::FragmentLibraryManager::get_instance()->clear_Vall();

		//build internal AtomID vector for TorsionDOFs later
		build_AtomID_vec();

		//save pose
		fixed_starting_pose_ = new core::pose::Pose(pose);
	}

	virtual ~HECTAllMover(){};

	//no apply function needed - defined by parent

	virtual
	void
	add_frag_mover(protocols::moves::RandomMoverOP random_mover) {

		using protocols::abinitio::ClassicFragmentMover;
		using protocols::abinitio::ClassicFragmentMoverOP;
		ClassicFragmentMoverOP frag_mover_hinge = new ClassicFragmentMover(fragset3mer_e3_hinge_, movemap_);
		frag_mover_hinge->enable_end_bias_check(false);
		random_mover->add_mover(frag_mover_hinge, 1.0);
		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new HECTAllMover;
	}

};

typedef utility::pointer::owning_ptr< HECTAllMover > HECTAllMoverOP;

int main( int argc, char* argv[] )
{
	basic::prof_reset();
	register_options(); //in HECT header
	devel::init(argc, argv);
	protocols::jd2::JobDistributor::get_instance()->go(new HECTAllMover);
	TR << "************************d**o**n**e**************************************" << std::endl;
	return 0;
}
