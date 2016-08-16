// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/sarel/hotspot_stub_constraint_test.cc
/// @brief For running an integration test
/// @author Lei Shi

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <devel/init.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>

// Unit Headers
//#include <protocols/hotspot_hashing/HotspotHashingConstraints.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>

//options
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

//patchdock
#include <protocols/protein_interface_design/read_patchdock.hh>

#include <iostream>
#include <iomanip>
#include <ObjexxFCL/string.functions.hh>

#include <basic/options/option_macros.hh>

//job distribution
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

//viewer
#include <protocols/viewer/viewers.hh>

#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/forge/methods/pose_mod.hh>

//mc
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/RepeatMover.hh>

#include <protocols/docking/DockingLowRes.hh>

// C++ headers

static THREAD_LOCAL basic::Tracer TR( "apps.pilot_apps.shilei.score_hotspot_cst" );

typedef core::Size Size;
typedef core::Real Real;
typedef core::pose::Pose Pose;

using namespace core;
using namespace pose;
using namespace basic::options;
using namespace core::scoring;
using namespace protocols::moves;
using namespace std;

OPT_1GRP_KEY(StringVector,score_hotspot_cst,hotspot_names)
OPT_1GRP_KEY(Real,score_hotspot_cst,hotspot_score_weight)
OPT_1GRP_KEY(RealVector,score_hotspot_cst,hotspot_distcb_weight)

class run_score_hotspot : public protocols::moves::Mover {
	public:
	run_score_hotspot() { }

		void apply( pose::Pose & pose) {

		Pose archive_pose=pose;

	protocols::viewer::add_conformation_viewer( pose.conformation(), "pose" );

	core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function("interchain_cen") );
	core::scoring::ScoreFunctionOP scorefxn_emp( ScoreFunctionFactory::create_score_function("empty") );
	core::scoring::ScoreFunctionOP scorefxn_cen( ScoreFunctionFactory::create_score_function("interchain_cen") );
	//core::scoring::ScoreFunctionOP scorefxn12( get_score_function() );

	//(*scorefxn)(pose);
	//TR << "Patchdock Rosetta score is: " << pose.energies().total_energy() << std::endl;

	//scorefxn->reset();
	//(*scorefxn)(pose);
	//TR << "Pre-stub score is: " << pose.energies().total_energy() << std::endl;

	//read hotspot_name
	Size num_hotspot_name=0;
	num_hotspot_name=basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_names ].size();

	//print out input information
	/*
	for (Size i=1; i <= num_hotspot_name; i++) {
	TR << "hotspot_names["<<i<<"]: " << basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_names]()[i] << std::endl;
	TR << "hotspot_distcb_weight["<<i<<"]: " << basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_distcb_weight]()[i] << std::endl;
	}
	*/

	// Assign a fixed residue (for the constraints)
	//core::Size fixed_res(1);  // unused ~Labonte
	core::Size const chain_to_redesign = 2;
	//if ( chain_to_redesign == 1 ) fixed_res = pose.total_residue();  // unused ~Labonte
	//core::id::AtomID fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );
	core::Real const worst_allowed_stub_bonus(-1.);
	bool const apply_self_energies(false);
	core::Real const bump_cutoff(10.);
	bool const apply_ambiguous_constraints(true);

	//core::Real hotspot_distcb_weight=basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_distcb_weight];
	core::Real hotspot_score_weight=basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_score_weight];

	//read file and assign constraint to pose
	for ( Size i=1; i <= num_hotspot_name; i++ ) {
		//std::string hotspot_namei="hotspot_name"+ObjexxFCL::string_of(i);
		std::string const hotspot_name= basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_names]()[i];
		//   TR << "Reading and generate cst from: " << hotspot_name << std::endl;
		protocols::hotspot_hashing::HotspotStubSetOP hotspot_stub_setOP( new protocols::hotspot_hashing::HotspotStubSet );
		hotspot_stub_setOP->read_data( hotspot_name );
		hotspot_stub_setOP->add_hotspot_constraints_to_wholepose( pose, chain_to_redesign, hotspot_stub_setOP,
			basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_distcb_weight][i],
			worst_allowed_stub_bonus, apply_self_energies, bump_cutoff, apply_ambiguous_constraints );
		//hotspot_stub_setOP->clear();
	}

	//convert to centroid (Can only do it after setting hotspot since it uses the packer task)
	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
	scorefxn->set_weight( core::scoring::backbone_stub_linear_constraint, hotspot_score_weight );
	scorefxn_emp->set_weight( core::scoring::backbone_stub_linear_constraint, hotspot_score_weight );

	(*scorefxn_cen)(pose);
	core::Real CenScore= pose.energies().total_energies().dot( scorefxn_cen->weights() );

	//switch back to FA_STANDARD
	protocols::forge::methods::restore_residues( archive_pose, pose );
	core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD);

	//scorefxn12->set_weight( core::scoring::backbone_stub_linear_constraint, hotspot_score_weight );
	//(*scorefxn12)(pose);
	(*scorefxn_emp)(pose);
	//core::Real CstScore= pose.energies().total_energies().dot( scorefxn_emp->weights() );  // unused ~Labonte

	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	job->add_string_real_pair("Centroid_score", CenScore);
	//protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag()

}//end of apply

virtual std::string get_name() const {
	return "run_score_hotspot";
}

};//end of run_score_hotspot


///////////////////////////////////////////////////////////////////////////////
void* my_main( void* ) {
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( MoverOP( new run_score_hotspot() ) );

	if ( basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_names ].user() && basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_distcb_weight].user()  &&
			basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_names ].size()==basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_distcb_weight].size() ) {
		//          TR << "Will read: " << basic::options::option[ basic::options::OptionKeys::score_hotspot_cst::hotspot_names ].size() << " hotspot files" << std::endl;
	} else {
		throw( utility::excn::EXCN_BadInput("expected hostspot_filename and hotspot_distcb_weight this app and their size should be equal") );
	}

	try{
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
int main( int argc, char * argv [] )
{
	try {

		NEW_OPT(score_hotspot_cst::hotspot_names, "hotspot_names","hotspot_name.files");
		NEW_OPT(score_hotspot_cst::hotspot_score_weight, "weight for hotspot_score_weight",10.0);
		NEW_OPT(score_hotspot_cst::hotspot_distcb_weight, "weight for Cb distance",utility::vector1<core::Real>());

		// setup random numbers and options
		devel::init(argc, argv);

		//see David Kim app dekim/score_nonlocal_frags.cc

		// run the test
		// score_hotspot();
		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

