// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/sarel/hotspot_stub_constraint_test.cc
/// @brief For running an integration test
/// @author Lei Shi

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

// Unit Headers
//#include <protocols/hotspot_hashing/HotspotHashingConstraints.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
// AUTO-REMOVED #include <protocols/docking/DockingProtocol.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>

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
#include <protocols/docking/metrics.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>

// C++ headers

static basic::Tracer TR( "apps.pilot_apps.shilei.lowresdock_patchdock_hotspot_cst" );

typedef core::Size Size;
typedef core::Real Real;
typedef core::pose::Pose Pose;

using namespace core;
using namespace pose;
using namespace basic::options;
using namespace core::scoring;
using namespace protocols::moves;
using namespace std;

OPT_1GRP_KEY(StringVector,lowresdock_patchdock_hotspot_cst,hotspot_names)
OPT_1GRP_KEY(Real,lowresdock_patchdock_hotspot_cst,hotspot_score_weight)
OPT_1GRP_KEY(RealVector,lowresdock_patchdock_hotspot_cst,hotspot_distcb_weight)
OPT_1GRP_KEY(Real,lowresdock_patchdock_hotspot_cst,centroidscore_filter)
OPT_1GRP_KEY(Real,lowresdock_patchdock_hotspot_cst,hotspotcst_filter)
//OPT_1GRP_KEY(Real,lowresdock_patchdock_hotspot_cst,dockinglowres_rot_mag)
//OPT_1GRP_KEY(Real,lowresdock_patchdock_hotspot_cst,dockinglowres_tran_mag)

class run_score_patchdock_hotspot : public protocols::moves::Mover {
public:
        run_score_patchdock_hotspot() { }

        void apply( pose::Pose & pose) {

				//can you convert chain B to be almost alanine?
  			//pose.dump_pdb( "start.pdb" );
        Pose archive_pose=pose;
				protocols::protein_interface_design::movers::BuildAlaPose toAla( 1,2,20);
				toAla.apply( pose );
  			//pose.dump_pdb( "check_alanine.pdb" );


	//read in reference pose
	/*
	pose::Pose pose,original_pose;
  std::string pdbname;
  if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].user() ) {
  			pdbname=basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1];
  			core::import_pose::pose_from_pdb( pose, pdbname.c_str() );
  			original_pose=pose;
  } else {
	      throw( utility::excn::EXCN_BadInput("expected -s for this app") );
  }

 //output patchdock
  //std::string outpdbname0="input_"+ObjexxFCL::string_of(basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1]);
  //pose.dump_pdb( outpdbname0.c_str() );

	//read in patchdock transformation (read a random one)
	//read a random one first and then figure out how to read them all
  protocols::protein_interface_design::PatchdockReader pd_reader;
  pd_reader.read_poses( pose, original_pose, pdbname, pdbname);
  //pd_reader.read_patchdock_entry();
	std::cout << "pd size: " << pd_reader.number_of_patchdock_entries() << std::endl;
	*/

	//should try to read all
	//output patchdock
  //std::string outpdbname1="output_pd_"+ObjexxFCL::string_of(basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1]);
  //pose.dump_pdb( outpdbname1.c_str() );
  //pose.dump_pdb( "output_pd_3m17.pdb" );

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
	num_hotspot_name=basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_names ].size();

	//print out input information
  /*
	for (Size i=1; i <= num_hotspot_name; i++) {
		TR << "hotspot_names["<<i<<"]: " << basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_names]()[i] << std::endl;
		TR << "hotspot_distcb_weight["<<i<<"]: " << basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_distcb_weight]()[i] << std::endl;
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

	//core::Real hotspot_distcb_weight=basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_distcb_weight];
	core::Real hotspot_score_weight=basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_score_weight];

	//read file and assign constraint to pose
	for (Size i=1; i <= num_hotspot_name; i++) {
  		//std::string hotspot_namei="hotspot_name"+ObjexxFCL::string_of(i);
  		std::string const hotspot_name= basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_names]()[i];
//			TR << "Reading and generate cst from: " << hotspot_name << std::endl;
	    protocols::hotspot_hashing::HotspotStubSetOP hotspot_stub_setOP = new protocols::hotspot_hashing::HotspotStubSet;
			hotspot_stub_setOP->read_data( hotspot_name );
		  hotspot_stub_setOP->add_hotspot_constraints_to_wholepose( pose, chain_to_redesign, hotspot_stub_setOP, basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_distcb_weight][i], worst_allowed_stub_bonus, apply_self_energies, bump_cutoff, apply_ambiguous_constraints );
		  //hotspot_stub_setOP->clear();
	}

  //convert to centroid (Can only do it after setting hotspot since it uses the packer task)
  core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
	scorefxn->set_weight( core::scoring::backbone_stub_linear_constraint, hotspot_score_weight );
	scorefxn_emp->set_weight( core::scoring::backbone_stub_linear_constraint, hotspot_score_weight );

	// rb_pert mover
	core::Size const rb_move_jump = 1; // use the first jump as the one between partners

  //just call a docking_highres_mover_
  protocols::docking::DockingLowResOP docking_lowres_mover = new protocols::docking::DockingLowRes( scorefxn, rb_move_jump );
  //docking_lowres_mover->set_trans_magnitude(basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::dockinglowres_tran_mag]);
  //docking_lowres_mover->set_rot_magnitude(basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::dockinglowres_rot_mag]);

	docking_lowres_mover->apply(pose);

	//this censcore of interaction is the one with Alanine at the interface
	(*scorefxn_cen)(pose);
  //core::Real CenScore= pose.energies().total_energies().dot( scorefxn_cen->weights() );
  utility::vector1< core::Size > movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
  core::Real CenScore=protocols::docking::calc_interaction_energy(pose,scorefxn_cen,movable_jumps_);

	//switch back to FA_STANDARD
  protocols::forge::methods::restore_residues( archive_pose, pose );
  core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD);

	//scorefxn12->set_weight( core::scoring::backbone_stub_linear_constraint, hotspot_score_weight );
	//(*scorefxn12)(pose);
	(*scorefxn_emp)(pose);
  core::Real CstScore= pose.energies().total_energies().dot( scorefxn_emp->weights() );

   protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
   job->add_string_real_pair("Centroid_score", CenScore);
  //protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag()

   // job distributor iterates ...
    if ( CenScore > basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::centroidscore_filter]
      || CstScore> basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspotcst_filter] ) {
    	set_last_move_status( protocols::moves::FAIL_RETRY );
    	//set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
      TR<< protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag() << " Did not pass filter CentoridScore ("<< basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::centroidscore_filter] <<"): " << CenScore << " HotspotCstScore ("<< basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspotcst_filter] <<"):" << CstScore <<std::endl;
      return;
    } else {
      TR<< protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag() << " Succeed Centorid Interaction: "<< CenScore << " HotspotCstScore: " << CstScore <<std::endl;
    	set_last_move_status( protocols::moves::MS_SUCCESS);
   }

}//end of apply

  virtual std::string get_name() const {
          return "run_score_patchdock_hotspot";
  }

};//end of run_score_patchdock_hotspot


///////////////////////////////////////////////////////////////////////////////
void* my_main( void* ) {
        using namespace protocols::moves;

        SequenceMoverOP seq( new SequenceMover() );
        seq->add_mover( new run_score_patchdock_hotspot() );

        if ( basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_names ].user() && basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_distcb_weight].user()  &&
             basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_names ].size()==basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_distcb_weight].size()) {
// 	 	  		  TR << "Will read: " << basic::options::option[ basic::options::OptionKeys::lowresdock_patchdock_hotspot_cst::hotspot_names ].size() << " hotspot files" << std::endl;
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

  NEW_OPT(lowresdock_patchdock_hotspot_cst::hotspot_names, "hotspot_names","hotspot_name.files");
  NEW_OPT(lowresdock_patchdock_hotspot_cst::hotspot_score_weight, "weight for hotspot_score_weight",10);
  NEW_OPT(lowresdock_patchdock_hotspot_cst::hotspot_distcb_weight, "weight for Cb distance",utility::vector1<core::Real>());
  NEW_OPT(lowresdock_patchdock_hotspot_cst::centroidscore_filter, "centroidscore_filter",0.0);
  NEW_OPT(lowresdock_patchdock_hotspot_cst::hotspotcst_filter, "hotspotcst_filter",20);
  //NEW_OPT(lowresdock_patchdock_hotspot_cst::dockinglowres_tran_mag, "translation magnitude",3.0);
  //NEW_OPT(lowresdock_patchdock_hotspot_cst::dockinglowres_rot_mag, "rotation maginitude",8.0);

	// setup random numbers and options
	devel::init(argc, argv);

	//see David Kim app dekim/score_nonlocal_frags.cc

	// run the test
	// score_patchdock_hotspot();
  protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

