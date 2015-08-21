// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MiniRelax.cc
/// @brief
/// @details
/// @author James Thompson

#include <protocols/relax/MiniRelax.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/cst_util.hh>


#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/pose/Pose.hh>

#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.relax.MiniRelax" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {

MiniRelax::MiniRelax(
	core::scoring::ScoreFunctionOP scorefxn_in
) :
	parent(),
	scorefxn_(scorefxn_in)
{}

MiniRelax::MiniRelax( MiniRelax const & other ) :
	//utility::pointer::ReferenceCount(),
	parent( other ),
	scorefxn_( other.scorefxn_ )
{}

MiniRelax::~MiniRelax() {}

protocols::moves::MoverOP
MiniRelax::clone() const {
	return protocols::moves::MoverOP( new MiniRelax(*this) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void MiniRelax::apply( core::pose::Pose & pose ) {
	if ( !pose.is_fullatom() ) {
		core::util::switch_to_residue_type_set(
			pose, core::chemical::FA_STANDARD
		);
		std::cerr << "Fullatom mode .... " << std::endl;
	}

	core::scoring::constraints::ConstraintSetOP orig_cst_set(
		pose.constraint_set()->clone()
	);

	// add coordinate constraints
	using core::Real;
	using core::sequence::alignment_from_pose;
	using namespace basic::options;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	Real const default_coord_sdev( option[ OptionKeys::relax::minirelax_sdev ]() );
	{
		coordinate_constrain_selection(
			pose, alignment_from_pose(pose),
			default_coord_sdev
		);
		scorefxn_->set_weight( coordinate_constraint, 1.0 );
	}

	FastRelax relaxer( scorefxn_ );

	int const n_repeats( option[ OptionKeys::relax::minirelax_repeats ]() );
	Real multiplier( 1.0 );
	Real const growth_factor( 3.0 );

	core::id::SequenceMapping map
		= alignment_from_pose(pose).sequence_mapping(1,2);
	for ( int ii = 1; ii <= n_repeats; ++ii ) {
		Real const current_round_sdev( default_coord_sdev * multiplier );
		// relax for a single round
		apply_disulfides(pose);
		relaxer.apply(pose);

		// generate a new set of coordinate constraints
		using utility::vector1;
		vector1< Real > fa_reps = get_per_residue_scores(
			pose, core::scoring::fa_rep
		);
		Real const clash_cutoff(2); // min fa-rep score for a clash

		vector1< Real > coord_sdevs;
		for ( Size idx = 1; idx <= pose.total_residue(); ++ii ) {
			Real coord_sdev(0);
			if ( map[idx] != 0 ) {
				// residue is aligned, sdev = default_coord_sdev
				coord_sdev = default_coord_sdev;
			}
			if ( fa_reps[idx] > clash_cutoff ) {
				// residue is aligned, sdev = default_coord_sdev * multiplier
				coord_sdev = current_round_sdev;
			}
			++idx;
			coord_sdevs.push_back(coord_sdev);
		}
		// apply coordinate constraints
		ConstraintSetOP cst_set = generate_bb_coordinate_constraints(
			pose, coord_sdevs
		);

		// increase multiplier by growth_factor
		multiplier *= growth_factor;
	} // n_repeats

	// clean up hydrogens
	//{
	// using namespace core::conformation;
	// typedef core::conformation::ResidueOPs::iterator iter;
	// for ( iter it = pose.res_begin(), end = pose.res_end();
	//    it != end; ++it
	// ) {
	//  idealize_hydrogens( *(*it), pose.conformation() );
	// }
	//}

	delete_virtual_residues(pose);

	// one more round with no constraints
	{
		using namespace core::scoring;
		scorefxn_->set_weight( coordinate_constraint, 0.0 );
		FastRelax no_cst_relaxer( scorefxn_ , 1 );
		no_cst_relaxer.apply(pose);
	}

	apply_disulfides(pose);
	pose.constraint_set( orig_cst_set );
} // apply

std::string
MiniRelax::get_name() const {
	return "MiniRelax";
}

} // namespace relax
} // namespace protocols
