// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/util/disulfide_util.cc
/// @brief A collection of procedures for manipulating disulfide bonds
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date 4/30/2009

// Unit Headers
#include <core/util/disulfide_util.hh>

// Project Headers
// Package Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/make_symmetric_task.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ headers

#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>


namespace core {
namespace util {

using namespace core;
using namespace std;

using utility::vector1;
using core::pose::Pose;
using core::pose::PoseOP;
using namespace core::conformation;
using core::pack::task::PackerTaskOP;
using core::scoring::ScoreFunctionOP;
using core::kinematics::MoveMap;
using core::kinematics::MoveMapOP;
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.disulfide_util" );

/// @details A convenience function for calling core::util:: rebuild_disulfide() with
///  only a single disulfide bond.  Supports symmetric poses.
void
rebuild_disulfide( Pose & pose, Size lower_res, Size upper_res,
	PackerTaskOP packer_task, ScoreFunctionOP packer_score,
	MoveMapOP mm, ScoreFunctionOP minimizer_score )
{
	vector1< pair<Size,Size> > disulfides;

	//If this is a symmetric pose, add the equivalent positions:
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformationOP conf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation >( pose.conformation_ptr() ) );
		bool r1_is_indep = conf->Symmetry_Info()->bb_is_independent( lower_res );
		bool r2_is_indep = conf->Symmetry_Info()->bb_is_independent( upper_res );

		Size s1master = lower_res, s2master = upper_res;
		if ( !r1_is_indep && !r2_is_indep ) {
			s1master = conf->Symmetry_Info()->bb_follows(lower_res);
			s2master = conf->Symmetry_Info()->bb_follows(upper_res);
		}
		disulfides.push_back(std::make_pair(s1master,s2master));

		// special logic if this crosses a symm boundary
		if ( r1_is_indep != r2_is_indep ) {
			Size A = r1_is_indep? lower_res : upper_res;
			Size Bprime = r1_is_indep? upper_res : lower_res;
			Size B = conf->Symmetry_Info()->bb_follows(Bprime);

			core::Size nclones = conf->Symmetry_Info()->num_bb_clones( );
			conformation::symmetry::SymmetryInfo::Clones clonesA = conf->Symmetry_Info()->bb_clones( A );
			bool found = false;
			numeric::xyzVector< core::Real > xyzA = conf->residue( A ).xyz(1);
			for ( int i=1; i<=(int)nclones && !found; ++i ) {
				Size Aprime = clonesA[i];
				numeric::xyzVector< core::Real > xyzAstar = conf->apply_transformation( conf->residue(Aprime).xyz(1), A, Bprime );
				core::Real dist = (xyzAstar - xyzA).length();
				if ( dist < 1e-4 ) {
					disulfides.push_back(std::make_pair(Aprime,B));
					found = true;
				}
			}
			if ( !found ) {
				TR << "Error in rebuild_disulfide: unable to find symmetry partner!" << std::endl;
			}
		}
	} else {
		// nonsymmetric, do nothing
		disulfides.push_back(std::make_pair(lower_res,upper_res));
	}

	core::util:: rebuild_disulfide(pose, disulfides, packer_task, packer_score, mm, minimizer_score);
}


/// @details
/// @pre The two residues specified should already be cysteines with a bond
///  defined between the SG atoms. This can be accomplished by first calling
///  protocols::toolbox::form_disulfide().
///
///  To provide the most flexibility, this function does not restrict
///  the degrees of freedom in the repacking and minimization steps. This is
///  needed in some cases where the backbone position is not optimal for making
///  a disulfide bond. However, if this function needs to be reasonably fast
///  you should prevent repacking on all residues except the two targets.
///
/// @param pose[in,out] The pose to modify with a disulfides
/// @param lower_res[in] The first residue of the disulfide
/// @param upper_res[in] The second residue of the disulfide
/// @param packer_task[in] A task to control repacking. Defaults to repacking
///  lower_res \& upper_res if ommitted or NULL.
/// @param packer_score[in] A scoring function to use while repacking. Use default
///  score function if ommitted or NULL.
/// @param mm[in] A MoveMap to control minimization. Defaults to full degrees
///  of freedom if ommitted or NULL.
/// @param minimizer_score[in] The score to use for minimization. Defaults to
///  packer_score if ommitted or NULL.
void
rebuild_disulfide( core::pose::Pose & pose,
	utility::vector1<std::pair<core::Size, core::Size> > disulfides,
	core::pack::task::PackerTaskOP packer_task,
	core::scoring::ScoreFunctionOP packer_score,
	core::kinematics::MoveMapOP mm,
	core::scoring::ScoreFunctionOP minimizer_score )
{
	// Quick lookup of whether a residue is a disulfide or not
	vector1<bool> is_disulf(pose.size(), false);

	for ( vector1<pair<Size, Size> >::const_iterator
			disulf(disulfides.begin()), end_disulf(disulfides.end());
			disulf != end_disulf; ++disulf ) {
		is_disulf[disulf->first] = true;
		is_disulf[disulf->second] = true;

		// Verify precondition
		if ( ! core::conformation::is_disulfide_bond(pose.conformation(), disulf->first, disulf->second) ) {
			TR.Error << "Disulfide bond required between " << disulf->first
				<< " and " << disulf->second << "." << std::endl;
			utility_exit();
		}
	}

	// Set up any NULL parameters
	if ( !packer_task ) {
		packer_task = core::pack::task::TaskFactory::create_packer_task( pose );
		packer_task->initialize_from_command_line().or_include_current( true );
		packer_task->restrict_to_repacking();

		// Restrict repacking to the targets
		for ( Size i(1); i <= pose.size(); ++i ) {
			if ( !is_disulf[i] ) {
				packer_task->nonconst_residue_task(i).prevent_repacking();
			}
		}
	}
	if ( !packer_score ) {
		packer_score = scoring::get_score_function();
	}
	if ( !mm ) {
		mm = MoveMapOP( new MoveMap );
		mm->set_bb( true );
		mm->set_chi( true );
	}
	if ( !minimizer_score ) {
		minimizer_score = packer_score;
	}

	// Extend rotamers for the disulfide
	for ( vector1<pair<Size, Size> >::const_iterator
			disulf(disulfides.begin()), end_disulf(disulfides.end());
			disulf != end_disulf; ++disulf ) {
		packer_task->nonconst_residue_task(disulf->first).and_extrachi_cutoff( 0 );
		packer_task->nonconst_residue_task(disulf->second).and_extrachi_cutoff( 0 );
		packer_task->nonconst_residue_task(disulf->first).or_ex1_sample_level(
			pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
		packer_task->nonconst_residue_task(disulf->second).or_ex1_sample_level(
			pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
	}

	// REPACK
	(*packer_score)(pose); // structure must be scored before rotamer_trials can be called
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		TR << "Performing symmetric packing on disulfides." << std::endl;
	}
	core::pack::pack_rotamers(pose, *packer_score, packer_task);

	using namespace core::optimization;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		TR << "Performing symmetric minimization on disulfides." << std::endl;
		core::kinematics::MoveMapOP symm_mm( mm->clone() );
		core::pose::symmetry::make_symmetric_movemap( pose, *symm_mm );
		core::optimization::symmetry::SymAtomTreeMinimizer symm_minimizer;
		(*minimizer_score)(pose);
		symm_minimizer.run( pose, *symm_mm, *minimizer_score, MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/ ) );
	} else {
		AtomTreeMinimizer().run( pose, *mm, *minimizer_score, MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true/*nblist*/, false/*deriv_check*/ ) );
	}

	// update score
	pose.update_residue_neighbors();
	(*minimizer_score)( pose );

}

} // util
} // core
