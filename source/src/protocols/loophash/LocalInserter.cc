// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LocalInserter.cc
/// @brief
/// @author Mike Tyka

#include <protocols/loophash/LocalInserter.hh>
// AUTO-REMOVED #include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/BackboneDB.hh>

#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/conformation/util.hh>

#include <protocols/relax/cst_util.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/Pose.hh>

using namespace core;
using namespace core::pose;

namespace protocols {
namespace loophash {

static basic::Tracer TR("LocalInserter");

/// @details Auto-generated virtual destructor
LocalInserter::~LocalInserter() {}


///@brief set the score function for the first round of minimization
///during a loophash insert
void
LocalInserter_SimpleMin::scorefxn_rama_cst(
	core::scoring::ScoreFunction const & scorefxn
){
	scorefxn_rama_cst_=scorefxn.clone();
}

///@brief set the score function for the second round of minimization
///during a loophash insert
void
LocalInserter_SimpleMin::scorefxn_cen_cst(
	core::scoring::ScoreFunction const & scorefxn
){
	scorefxn_cen_cst_=scorefxn.clone();
}

core::Real
LocalInserter_SimpleMin::make_local_bb_change(
	core::pose::Pose &newpose,
	const core::pose::Pose &original_pose,
	const protocols::loophash::BackboneSegment &new_bs,
	core::Size res_pos
)
{
	using namespace core;
	using namespace optimization;

	// set newpose
	protocols::loops::Loops exclude_region;
	exclude_region.add_loop( protocols::loops::Loop( res_pos, res_pos + new_bs.length() ) );
	//core::pose::Pose newpose( original_pose );
	transfer_phi_psi( original_pose, newpose );
	add_coordinate_constraints_to_pose( newpose, original_pose, exclude_region );
	new_bs.apply_to_pose( newpose, res_pos );
	pose::set_ss_from_phipsi( newpose );

	//scorefxn_rama_cst->show( TR.Info, *newpose );

	kinematics::MoveMap final_mm;
	final_mm.set_bb(true);
	// setup movemap & minimisation

	// just for comparison with cut!
	//core::pose::PoseOP newpose2( new Pose( original_pose ) );
	//new_bs.apply_to_pose( *newpose2, ir, true );
	//newpose2->dump_pdb("rep_" + utility::to_string( ir ) + "_" + utility::to_string( jr ) + "_" + utility::to_string( int(xyzdist) ) + "_" + utility::to_string( int(angdist) ) + ".cut.pdb" );
	//scorefxn_rama_cst->show( TR.Info, *newpose );
	//newpose->dump_pdb("rep_" + utility::to_string( ir ) + "_" + utility::to_string( jr ) + "_" + utility::to_string( int(xyzdist) ) + "_" + utility::to_string( int(angdist) ) + ".bef.pdb" );
	AtomTreeMinimizer().run( newpose, final_mm, *scorefxn_rama_cst_, options_ );
	//scorefxn_rama_cst->show( TR.Info, *newpose );
	//newpose->dump_pdb("rep_" + utility::to_string( ir ) + "_" + utility::to_string( jr ) + "_" + utility::to_string( int(xyzdist) ) + "_" + utility::to_string( int(angdist) ) + ".aft.pdb" );
	//newpose->dump_pdb("rep_" + utility::to_string( ir ) + "_" + utility::to_string( jr ) + "_" + utility::to_string( int(xyzdist) ) + "_" + utility::to_string( int(angdist) ) + ".pdb" );

	core::Real premin_rms = core::scoring::CA_rmsd( newpose, original_pose );
	//scorefxn_cen_cst->show( TR.Info, *newpose );
	AtomTreeMinimizer().run( newpose, final_mm, *scorefxn_cen_cst_, options2_ );
	//scorefxn_cen_cst->show( TR.Info, *newpose );

	// get final RMS
	core::Real final_rms = core::scoring::CA_rmsd( newpose, original_pose );
	TR.Debug << "Premin RMS: " << premin_rms << "Min Score3 " << "Final RMS: " << final_rms << std::endl;


	//transfer_phi_psi( newpose, start_pose );

	core::Real final_score = scorefxn_cen_cst_->score(newpose);

	TR.Debug << "INSERTRESULT: " << final_rms << "	" << final_score << std::endl;
	core::pose::setPoseExtraScores( newpose, "lh_censcore", final_score );

	//transfer_phi_psi( newpose, start_pose );

	return final_rms;
}

core::Real
LocalInserter_SimpleMin::make_local_bb_change_close_gaps(
	core::pose::Pose &newpose,
	const core::pose::Pose &original_pose,
	const protocols::loophash::BackboneSegment &new_bs,
	core::Size res_pos
)
{
	using namespace core;
	using namespace optimization;

	// set newpose
	protocols::loops::Loops exclude_region;
	utility::vector1 < Size > excluded_res;
	exclude_region.add_loop( protocols::loops::Loop( res_pos, res_pos + new_bs.length() ) );
	//core::pose::Pose newpose( original_pose );
	transfer_phi_psi( original_pose, newpose );
	add_coordinate_constraints_to_pose( newpose, original_pose, exclude_region );


	// fix gaps between ir and jr if it exists by idealizing every position
	for( Size idx = res_pos, end_pos = res_pos + new_bs.length(); idx <= end_pos; idx ++ ) {
		conformation::insert_ideal_bonds_at_polymer_junction( idx, newpose.conformation() );
		excluded_res.push_back(idx);
	}
	new_bs.apply_to_pose( newpose, res_pos );

	pose::set_ss_from_phipsi( newpose );


	kinematics::MoveMap final_mm;
	final_mm.set_bb(true);
	// setup movemap & minimisation

	AtomTreeMinimizer().run( newpose, final_mm, *scorefxn_rama_cst_, options_ );

	core::Real premin_rms = core::scoring::CA_rmsd( newpose, original_pose );
	//scorefxn_cen_cst->show( TR.Info, *newpose );
	AtomTreeMinimizer().run( newpose, final_mm, *scorefxn_cen_cst_, options2_ );
	//scorefxn_cen_cst->show( TR.Info, *newpose );



	// get final RMS
	// since we're closing gaps, we only want to make sure the rmsd of the non loophash region is under some cutoff
	core::Real final_rms = core::scoring::CA_rmsd( newpose, original_pose, 1, original_pose.total_residue(), excluded_res );
	TR.Debug << "Premin RMS: " << premin_rms << "Min Score3 " << "Final RMS: " << final_rms << std::endl;



	core::Real final_score = scorefxn_cen_cst_->score(newpose);

	TR.Debug << "INSERTRESULT: " << final_rms << "	" << final_score << std::endl;
	core::pose::setPoseExtraScores( newpose, "censcore", final_score );

	// remove that extra virtual atom cooordinate constraints adds on
	protocols::relax::delete_virtual_residues( newpose );

	return final_rms;
}

core::Real
LocalInserter_SimpleMin::make_local_bb_change_include_cut(
	core::pose::Pose &newpose,
	const core::pose::Pose &original_pose,
	const protocols::loophash::BackboneSegment &new_bs,
	core::Size res_pos
)
{
	using namespace core;
	using namespace optimization;

	// set newpose
	protocols::loops::Loops exclude_region;
	exclude_region.add_loop( protocols::loops::Loop( res_pos, res_pos + new_bs.length() ) );
	transfer_phi_psi( original_pose, newpose );
	add_coordinate_constraints_to_pose( newpose, original_pose, exclude_region );
	new_bs.apply_to_pose( newpose, res_pos, true );
	pose::set_ss_from_phipsi( newpose );

	// get final RMS
	core::Real final_rms = core::scoring::CA_rmsd( newpose, original_pose );
	TR << "Final RMS: " << final_rms << std::endl;
	TR.Debug << "Final RMS: " << final_rms << std::endl;

	//transfer_phi_psi( newpose, start_pose );

	core::Real final_score = scorefxn_cen_cst_->score(newpose);

	TR.Debug << "INSERTRESULT: " << final_rms << "	" << final_score << std::endl;
	core::pose::setPoseExtraScores( newpose, "censcore", final_score );
	core::pose::setPoseExtraScores( newpose, "rms", final_rms );

	//transfer_phi_psi( newpose, start_pose );

	return final_rms;
}

void
LocalInserter_SimpleMin::set_default_score_functions(){
	using namespace core::scoring;

	scorefxn_rama_cst_->set_weight( coordinate_constraint, 0.5 );
	scorefxn_rama_cst_->set_weight( rama		, 1.0 );

	scorefxn_cen_cst_->set_weight( coordinate_constraint, 0.05 );
	scorefxn_cen_cst_->set_weight( env			, 1.0);
	scorefxn_cen_cst_->set_weight( pair		 , 1.0);
	scorefxn_cen_cst_->set_weight( cbeta		, 1.0);
	scorefxn_cen_cst_->set_weight( vdw			, 1.0);
	scorefxn_cen_cst_->set_weight( rg			 , 3.0);
	scorefxn_cen_cst_->set_weight( cenpack	, 1.0);
	scorefxn_cen_cst_->set_weight( hs_pair	, 1.0);
	scorefxn_cen_cst_->set_weight( ss_pair	, 1.0);
	scorefxn_cen_cst_->set_weight( rsigma	 , 1.0);
	scorefxn_cen_cst_->set_weight( sheet		, 1.0);
}

} // namespace loophash
} // namespace protocols




