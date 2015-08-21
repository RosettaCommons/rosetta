// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  /src/apps/public/scenarios/beta_strand_homodimer_design/homodimer_maker.cc
/// @brief Takes a found exposed beta strand and tries to make the starting file for a symmetric homodimer
/// @author Ben Stranges

// Unit headers
#include <devel/init.hh>

// Project Headers
//#include <core/chemical/AA.hh>
#include <core/io/pdb/pose_io.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyzVector.hh>
//#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
//#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>

#include <core/conformation/Conformation.hh>


//#include <core/scoring/rms_util.tmpl.hh>
//#include <core/scoring/ScoreType.hh>
//#include <core/scoring/TenANeighborGraph.hh>

#include <core/id/AtomID_Map.hh>
//#include <core/id/AtomID_Map.Pose.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>

#include <core/scoring/EnergyMap.hh>
//#include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>
//#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/EtableEnergy.hh>
//#include <core/scoring/etable/CoarseEtableEnergy.hh>

//protocols
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RollMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/StructureRestrictor.hh>

//options
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Job distributor
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>


// C++ headers
#include <sstream>
#include <iostream>
#include <string>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

static thread_local basic::Tracer TR( "apps.public.beta_strand_homodimer_design.homodimer_maker" );


using namespace core;
using namespace utility;
using namespace protocols;
using namespace protocols::moves;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//specific options

IntegerOptionKey const window_size( "window_size" );
IntegerOptionKey const sheet_start( "sheet_start" );
IntegerOptionKey const sheet_stop( "sheet_stop" );
RealOptionKey const E_cutoff( "E_cutoff" );
StringOptionKey const struct_file( "struct_file" ); //filename of structure restrictor

//helper function because friend functions from numeric/xyzVector.hh are bad.
numeric::xyzVector< core::Real >
xyz_center_vector( numeric::xyzVector<core::Real> const & a, numeric::xyzVector<core::Real> const & b )
{
	return numeric::xyzVector<core::Real>(
		0.5 * ( a.x() + b.x() ),
		0.5 * ( a.y() + b.y() ),
		0.5 * ( a.z() + b.z() )
	);
}

// mover deffinition
class HDmakerMover : public Mover {
public:

	HDmakerMover();

	virtual void apply( core::pose::Pose& pose );

	virtual std::string get_name() const{
		return "HDmakerMover";
	}

	virtual numeric::xyzVector <core::Real> find_midpoint( numeric::xyzVector<core::Real> & xpoints,
		numeric::xyzVector<core::Real> & ypoints );

	virtual Real bb_score(pose::Pose & pose, core::Size aligned_chain_num, core::scoring::ScoreFunctionOP & scorefxn);


	virtual MoverOP clone() const {
		return MoverOP( new HDmakerMover( *this ) );
	}

	virtual MoverOP fresh_instance() const {
		return clone();
	}

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	int sheet_start_, sheet_stop_;
	std::string pdb_chain_;
	Real maxE_;

};

HDmakerMover::HDmakerMover() {
	// variable definitions
	scorefxn_ = core::scoring::get_score_function();
	sheet_start_ = option[ sheet_start ];
	sheet_stop_ = option[ sheet_stop ];
	pdb_chain_ = option[run::chain].def( "A");
	maxE_=option[ E_cutoff ].def(30.0);
	//window_size_= option[ window_size ].def(5)
}


//helper function to find midpoint
numeric::xyzVector <core::Real>
HDmakerMover::find_midpoint( numeric::xyzVector<core::Real> & xpoints, numeric::xyzVector<core::Real> & ypoints ){
	numeric::xyzVector<core::Real> midpoint;
	midpoint = (xpoints+ypoints);
	midpoint /= 2;
	return midpoint;
}


///////////////////////////////////////
// bb score
///////////////////////////////////////
core::Real
HDmakerMover::bb_score(pose::Pose & pose, core::Size aligned_chain_num, core::scoring::ScoreFunctionOP & scorefxn){

	// score the bb-bb energy between chains
	// This part written by P.Doug Renfrew
	// make vectors of bb atoms in each chain individually
	// the master pose will always be chain 1.
	// need to make a vector of all atoms in the chain you are comparing too

	utility::vector1<core::conformation::Atom> chain1_bb_atoms;
	utility::vector1<core::conformation::Atom> chain2_bb_atoms;
	utility::vector1<core::conformation::Atom> all_bb_atoms;

	for ( Size j = 1; j <= pose.total_residue(); ++j ) {
		core::conformation::Residue const & res( pose.residue(j) );
		core::chemical::AtomIndices bb_ai( res.mainchain_atoms() );
		//assert( bb_ai.size() == 4 );
		core::Size chain_num( res.chain() );
		for ( Size jj = 1; jj <= bb_ai.size(); ++jj ) {
			if ( chain_num == 1 ) {
				chain1_bb_atoms.push_back( res.atom(jj) );
			} else if ( chain_num == aligned_chain_num ) {
				chain2_bb_atoms.push_back( res.atom(jj) );
			}
			//optional get all the atoms not in allinged chain,
			//only need to do if more than two chains in pose
			if ( pose.conformation().num_chains() >= 3 && chain_num != 1 ) {
				all_bb_atoms.push_back( res.atom(jj) );
			}
			//end optional
		}
	}

	//NOW SCORE!
	// get instance of etable energy method
	core::scoring::methods::EnergyMethodOptions const & emo(scorefxn->energy_method_options());
	core::scoring::etable::Etable const & et(*(core::scoring::ScoringManager::get_instance()->etable(emo).lock()));

	if ( basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ] || emo.analytic_etable_evaluation() ) {
		utility_exit_with_message("homodimer_maker is incompatible with analytic_etable_evaluation. Please either fix the code or add '-analytic_etable_evalution 0' to the command line.");
	}
	core::scoring::etable::TableLookupEtableEnergy ete( et, emo );

	// iterate over both sets of atom and add into one emapvector
	//core::scoring::TwoBodyEMapVector tbemv;
	core::scoring::EMapVector tbemv;
	core::Real atr_wt( (*scorefxn).get_weight(core::scoring::fa_atr) );
	core::Real rep_wt( (*scorefxn).get_weight(core::scoring::fa_rep) );
	for ( Size ii = 1; ii <= chain1_bb_atoms.size(); ++ii ) {
		for ( Size jj = 1; jj <= chain2_bb_atoms.size(); ++jj ) {
			//calc distr squared
			Real d2( chain1_bb_atoms[ii].xyz().distance_squared( chain2_bb_atoms[jj].xyz() ) );
			ete.atom_pair_energy( chain1_bb_atoms[ii], chain2_bb_atoms[jj], 1, tbemv, d2 );
		}
	}
	core::Real bb_energy (rep_wt * tbemv[core::scoring::fa_rep] + atr_wt * tbemv[core::scoring::fa_atr] );

	// begin optional  ie skip if not needed
	core::Real all_energy;
	if ( pose.conformation().num_chains() >= 3 ) {
		//core::scoring::TwoBodyEMapVector tbemv_all;
		core::scoring::EMapVector tbemv_all;
		core::scoring::etable::TableLookupEtableEnergy ete_all( et, emo );
		for ( Size ii = 1; ii <= chain1_bb_atoms.size(); ++ii ) {
			for ( Size jj = 1; jj <= all_bb_atoms.size(); ++jj ) {
				//calc distr squared
				Real d2_all( chain1_bb_atoms[ii].xyz().distance_squared(  all_bb_atoms[jj].xyz() ) );
				ete.atom_pair_energy( chain1_bb_atoms[ii], all_bb_atoms[jj], 1, tbemv_all, d2_all );
			}
		}
		all_energy = (rep_wt * tbemv_all[core::scoring::fa_rep] + atr_wt * tbemv_all[core::scoring::fa_atr] );
	} else { //end optional for many chains
		all_energy = bb_energy ;
	}
	TR<< "Number of chains: " <<pose.conformation().num_chains()
		<<"    Backbone-backbone score: " << all_energy << std::endl;
	return all_energy;
	//return bb_energy;
}//end bb_score

/////////////////////////////////////////////////
// Actual mover apply
/////////////////////////////////////////////////
void HDmakerMover::apply (pose::Pose & pose ) {

	//job info
	protocols::jd2::JobOP const job_me( protocols::jd2::JobDistributor::get_instance()->current_job() );
	utility::file::FileName pdb_file_name (pose.pdb_info()->name());

	//restrict structures if needed
	if ( basic::options::option[ struct_file ].active() ) {
		std::string struct_filename = option[ struct_file ];
		TR << "Deleting according to structure file..."<< std::endl;
		moves::StructureRestrictorOP restrictor( new moves::StructureRestrictor( struct_filename ) );
		restrictor->apply(pose);
	}

	if ( pose.conformation().num_chains() > 1 ) {
		TR << "pose is not a monomer, skipping..."<< std::endl;
		set_last_move_status(protocols::moves::FAIL_BAD_INPUT);
		return;
	}
	//score this pose & fill hbond set
	(*scorefxn_)( pose );
	core::scoring::hbonds::HBondSet hbond_set;
	core::scoring::hbonds::fill_hbond_set(
		pose,
		false /*calc_deriv*/,
		hbond_set,
		false /*bb only*/ );
	//call to try to resize bb_don/accept arrays
	//need this for everything to work right
	hbond_set.setup_for_residue_pair_energies(pose);

	//translate pdb numbering to pose numbering
	char chain( pdb_chain_[0] );
	Size sheet_start( pose.pdb_info()->pdb2pose(chain, sheet_start_) );
	Size sheet_end( pose.pdb_info()->pdb2pose(chain, sheet_stop_) );

	//get some info about the pose
	Size n_residues (pose.total_residue());

	//dirty pose dupication
	pose::Pose copy_pose (pose);
	pose.append_residue_by_jump(copy_pose.residue( 1 ), pose.total_residue(), "" , "",  true /*start new chain*/);
	for ( Size n = 2; n <= copy_pose.total_residue(); ++n ) {
		pose.append_residue_by_bond( copy_pose.residue ( n ) );
	}

	//decision if to use the WHOLE defined sheet or a sliding window
	//assume that a sheet is locally linear and itterate through it based on some window size
	Size sheet_length = sheet_end - sheet_start + 1;
	TR << "Sheet Length: " << sheet_length;
	int window;
	if ( !option[ window_size ].active() ) {
		window = sheet_length;
	} else {
		window = option[ window_size ] ;
	}
	TR << "  Window size: "<< window << std::endl;
	Size last_start( sheet_end - window + 1 );
	Size window_num (1);
	pose::Pose const saved_pose (pose);
	for ( Size jj = sheet_start; jj <= last_start; ++jj ) {

		//skip this itteration if it is in a bb/bb hbond already
		if ( hbond_set.acc_bbg_in_bb_bb_hbond(jj) && hbond_set.don_bbg_in_bb_bb_hbond(jj) ) {
			continue;
		}

		Size this_end( jj + window -1 );
		//find center of the sheet
		Size center_residue(0); //should fail if nothing changes about this residue
		if ( window % 2 == 0 ) {
			center_residue = jj + (window / 2);
		} else { center_residue = jj + ((window - 1) / 2); }
		TR << "Sheet from: "<< jj << " to "<<  this_end << "    Window: " << window_num
			<< "  Center residue: " << center_residue <<std::endl;

		//make axis for parallel sheet making
		std::string const atom_to_use( "CA" );
		numeric::xyzVector< core::Real >  start_xyz ( pose.residue(jj).atom(atom_to_use).xyz() );
		numeric::xyzVector< core::Real >  end_xyz ( pose.residue(this_end).atom(atom_to_use).xyz() );
		numeric::xyzVector< core::Real >  center_xyz ( pose.residue(center_residue).atom(atom_to_use).xyz() );
		numeric::xyzVector<core::Real> const zero_vector (0,0,0);
		numeric::xyzVector<core::Real> midpoint( xyz_center_vector( start_xyz, end_xyz ) );
		numeric::xyzVector< core::Real > parl_vector( end_xyz - start_xyz );

#ifndef NDEBUG
		TR << "Parl axis is from pose residue: " <<jj << " to " << this_end <<"\n"
			<< "The vector between these two residues is: " << parl_vector << "\n"
			<< "The midpoint is : "<< midpoint << std::endl;
#endif

		//make axis for antiparallel sheet making
		//need another vector to define the plane that that sheet is on
		//choose the C=O  on the center residue
		numeric::xyzVector<core::Real> center_C_xyz ( pose.residue(center_residue).atom( "C" ).xyz() );
		numeric::xyzVector<core::Real> center_O_xyz ( pose.residue(center_residue).atom( "O" ).xyz() );
		numeric::xyzVector<core::Real> CO_plane_vector = center_O_xyz - center_C_xyz;
		numeric::xyzVector<core::Real> anti_vector = cross( CO_plane_vector, parl_vector );

#ifndef NDEBUG
		TR << "Anti axis is from pose residue: " <<jj << " to " << this_end <<"\n"
			<< "The normal to these two residues is: " << anti_vector << "\n"
			<< "The midpoint is : "<< midpoint << std::endl;
#endif

		//define some movers
		Real const min_angle (180.0);
		Real const max_angle(180.0);
		//moves the input pose, not the copy
		pose::Pose anti_pose (pose), parl_pose (pose);
		rigid::RollMoverOP parl_roll_mover( new rigid::RollMover(1 /*start_res*/, n_residues /*stop*/, min_angle, max_angle, parl_vector, center_xyz ) );
		rigid::RollMoverOP anti_roll_mover( new rigid::RollMover(1 /*start_res*/, n_residues /*stop*/, min_angle, max_angle, anti_vector, center_xyz ) );

		//apply to pose and output
		anti_roll_mover->apply( anti_pose );
		parl_roll_mover->apply( parl_pose );
		//std::string anti_string = pdb_file_name.base() + "_anti_complex.pdb" ;
		//std::string parl_string = pdb_file_name.base() + "_parl_complex.pdb" ;

		//now move chains apart
		Real anti_dist (6.0); //distance to move 2 CA atoms apart (largest for antiparallel seen)
		Real parl_dist (5.5); //similar number for parallel sheetsm

		rigid::RigidBodyTransMoverOP push_apart_mover( new rigid::RigidBodyTransMover );
		//figure out which way to push
		//check if the center residue is in bb:bb hbonds (lame way but will do for now)
		TR << "C-O vector on center res is" << CO_plane_vector << ",  ";
		Real scaler (-1); //will define direction
		if ( hbond_set.acc_bbg_in_bb_bb_hbond(center_residue) && hbond_set.don_bbg_in_bb_bb_hbond(center_residue) ) {
			scaler = 1;
		}
		CO_plane_vector *= scaler;
		TR << "Direction of translation is: " << CO_plane_vector << std::endl;
		push_apart_mover->trans_axis(CO_plane_vector);
		//now do push
		push_apart_mover->step_size(anti_dist);
		push_apart_mover->apply(anti_pose);
		push_apart_mover->step_size(parl_dist);
		push_apart_mover->apply(parl_pose);

		//chains are now apart
		//now explore space along the sheet
		//need to figure out how far along to move along the sheet.
		//Figure we need at LEAST two bb hbonding residues (4 total hbonds)

		//first recalculate the vector along the sheet as it has moved
		start_xyz = pose.residue(jj).atom(atom_to_use).xyz() ;
		end_xyz = pose.residue(this_end).atom(atom_to_use).xyz() ;
		parl_vector = end_xyz - start_xyz;

		//figure out how many steps are possible given the lenght of the sheet
		//doubt that sheet definition will ever be greater than 12
		int numsteps (0);
		if ( window <=4 ) {
			numsteps = 0;
		} else if ( window <= 6 ) {
			numsteps = 1;
		} else if ( window <= 8 ) {
			numsteps = 2;
		} else if ( window <= 10 ) {
			numsteps = 3;
		} else if ( window <= 12 ) {
			numsteps = 4;
		} else numsteps = 5;

		TR<<"Number of RB steps to search: "<< numsteps << std::endl;

		//now set up for searching along strand
		rigid::RigidBodyTransMoverOP sheet_trans_mover( new rigid::RigidBodyTransMover );
		//quick bump to line up parallel sheet
		sheet_trans_mover->trans_axis(parl_vector * (-1));
		sheet_trans_mover->step_size(3.6); //aprox only
		sheet_trans_mover->apply(parl_pose);

		//save poses
		pose::Pose saved_anti( anti_pose );
		pose::Pose saved_parl( parl_pose );
		//now do searching
		sheet_trans_mover->trans_axis(parl_vector);
		for ( int ii = (-1)*numsteps ; ii <= numsteps; ++ii ) {
			Real trans_dist (7.0); //aprox dist of CA from i to i+2 in sheet
			sheet_trans_mover->step_size(trans_dist*ii);
			sheet_trans_mover->apply(anti_pose);
			sheet_trans_mover->apply(parl_pose);
			TR<< "Translate step: "<< ii << " distance: " << trans_dist * ii << std::endl;

			//bb score the complex
			Real const anti_bb_score (bb_score(anti_pose,2/*assumtiuon!!!*/, scorefxn_));
			Real const parl_bb_score (bb_score(parl_pose,2/*assumtiuon!!!*/, scorefxn_));

			//naming for output
			std::stringstream anti_ss;
			std::stringstream parl_ss;
			int pdb_res( pose.pdb_info()->number(jj) );
			char pdb_chain( pose.pdb_info()->chain(jj) );
			anti_ss << pdb_file_name.base() <<"_"<< pdb_chain << pdb_res <<"_anti_wind_"<< window_num << "_step_" << ii   << ".pdb";
			parl_ss << pdb_file_name.base() <<"_"<< pdb_chain << pdb_res <<"_parl_wind_" << window_num << "_step_" << ii   << ".pdb";
			std::string anti_filename (anti_ss.str());
			std::string parl_filename (parl_ss.str());

			std::cout << "FILE:  "<< anti_filename << "  bb_score: " << anti_bb_score << std::endl;
			std::cout << "FILE:  "<< parl_filename << "  bb_score: " << parl_bb_score << std::endl;

			//dump pdbs only if score is good enough...
			if ( anti_bb_score < maxE_ ) {
				anti_pose.dump_pdb(anti_filename);
			}
			if ( parl_bb_score < maxE_ ) {
				parl_pose.dump_pdb(parl_filename);
			}

			//reverts to original
			anti_pose = saved_anti;
			parl_pose = saved_parl;
		}

		//   anti_pose.dump_scored_pdb(anti_string, *(scorefxn_));
		//   parl_pose.dump_scored_pdb(parl_string, *(scorefxn_));

		//set up for next itteration
		++window_num;
		pose = saved_pose;
	}//end window loop
}//end apply

//begin main
int
main( int argc, char* argv[] ) {
	try {
		//options
		option.add(sheet_start, "Start of beta sheet (PDBNum) to do rolls and translates about");
		option.add(sheet_stop, "End of beta sheet (PDBNum) to do rolls and translates about");
		option.add(E_cutoff, "Max E for output" );
		option.add( struct_file, "file with info about chains and such");
		option.add( window_size, "file with info about chains and such");

		//init
		devel::init( argc, argv );
		//making our own output
		basic::options::option[ OptionKeys::jd2::no_output ].def(true);

		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new HDmakerMover ) );

		TR<< "Complete." << std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}//main
