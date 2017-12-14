// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    FlexBBDesignMeanField.cc

/// @brief   Method definitions for FlexBBDesignMeanField.
/// @author  arubenstein

// Unit headers
#include <protocols/mean_field/FlexBBDesignMeanField.hh>

// Package headers
#include <protocols/mean_field/AAProb.hh>
#include <protocols/mean_field/AAMatrix.hh>
#include <protocols/mean_field/RotMatrix.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers


// C++ headers
#include <iostream>


// Construct tracer.
static basic::Tracer TR("protocols.mean_field.FlexBBDesignMeanField");


namespace protocols {
namespace mean_field {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
FlexBBDesignMeanField::FlexBBDesignMeanField( Size option,
	pose::PoseOPs & poses,
	utility::vector1 < core::pack::task::PackerTaskOP > tasks,
	scoring::ScoreFunctionOP scfxn
) :
	FlexBBMeanField::FlexBBMeanField( option, poses, tasks, scfxn ),
	exp_aa_matrix_( new AAMatrix() )
{
	aa_matrices_.resize( poses.size() );
}

FlexBBDesignMeanField::FlexBBDesignMeanField( Size option,
	pose::PoseOPs & poses,
	utility::vector1 < core::pack::task::PackerTaskOP > tasks,
	scoring::ScoreFunctionOP scfxn,
	Real lambda_mem,
	Real tolerance,
	Real temp,
	Real threshold
) :
	FlexBBMeanField::FlexBBMeanField( option, poses, tasks, scfxn, lambda_mem, tolerance, temp, threshold ),
	exp_aa_matrix_( new AAMatrix() )
{
	aa_matrices_.resize( poses.size() );
}

// Destructor
FlexBBDesignMeanField::~FlexBBDesignMeanField() = default;


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods

void
FlexBBDesignMeanField::show(std::ostream & output) const
{
	output << "Averaged Conformational Matrix (RotMatrix): " << std::endl;
	rot_matrix()->show( output );
	output << "Averaged Specificity Profile (AAMatrix): " << std::endl;
	aa_matrix()->show( output );
}

/// @details to implement for design adds step to calculate AAMatrixs
void
FlexBBDesignMeanField::process()
{
	//calculates rot_matrices, bb_boltz_probs, and exp_rot_matrix
	FlexBBMeanField::process();

	for ( Size bb = 1; bb <= num_poses(); ++bb ) {
		aa_matrices_[ bb ] = AAMatrix ( rot_matrices()[ bb ], energy_matrices()[ bb ], temperature() );
	}
	calc_exp_value_aa_matrix();
}

/// @details resizes the aa_matrices_ vector to allow for the deleted pose
/// @remarks called before aa_matrices_ has been fully assigned so simply truncating it is okay
void
FlexBBDesignMeanField::delete_pose( Size pose_ind )
{

	FlexBBMeanField::delete_pose( pose_ind );
	aa_matrices_.resize( num_poses() );
}

/// @details erases bb_boltz_probs_ for positions that aren't designed
/// @remarks this doesn't affect the averaging of the RotMatrix (which is done for positions that are packed
/// @remarks and for positions that are designed) because the RotMatrix is averaged using the bb_boltz_weights_
void
FlexBBDesignMeanField::calc_bb_boltz_probs()
{
	FlexBBMeanField::calc_bb_boltz_probs();
	utility::vector1 < bool > is_designed = rot_matrices()[ 1 ].is_designed();

	for ( Size pos = 0; pos < is_designed.size(); ++pos ) {
		if ( ! is_designed[ pos+1 ] ) {
			bb_boltz_probs().erase( bb_boltz_probs().begin() + pos );
			//Added per_aa
			for ( core::Size bb = 1; bb <= bb_boltz_probs_per_aa_.size(); ++bb ) {
				bb_boltz_probs_per_aa_[bb].erase( bb_boltz_probs_per_aa_[bb].begin()+pos);
			}
		}
	}
}

AAMatrixCOP
FlexBBDesignMeanField::aa_matrix() const
{
	return exp_aa_matrix_;
}

// Accessors/Mutators



// Private methods /////////////////////////////////////////////////////////////

/// @details calculate "expected value" or averaged AAMatrix
/// @details weight each position-backbone-aa probability by the Boltzmann probability of that backbone occurring for that position
/// @details and then sum them up across all backbones for that position
void
FlexBBDesignMeanField::calc_exp_value_aa_matrix()
{
	//no need for following line unless bb boltz per aa

	utility::vector1 < Real > totals( aa_matrices_[1].size(), Real( 0.0 ) );
	for ( Size pos = 1; pos <= aa_matrices_[1].size(); ++pos ) {
		exp_aa_matrix_->push_back( utility::vector1 < AAProb> ( aa_matrices_[1][ pos ].size() ) );

		for ( Size bb = 1; bb <= aa_matrices_.size(); ++bb ) {
			for ( Size aa = 1; aa <= aa_matrices_[ bb ][ pos ].size(); ++aa ) {
				//Added per_aa
				//This AAProb wasn't yet initialized so initialize it to value of current RotProb multiplied by the bb probability
				if ( (*exp_aa_matrix_)[ pos ][ aa ].nrot() == 0 ) {
					//      (*exp_aa_matrix_)[ pos ][ aa ] = AAProb( aa_matrices_[ bb ][ pos ][ aa ] ) *
					//       bb_boltz_probs_per_aa()[ bb ][ pos ][ aa ];
					(*exp_aa_matrix_)[ pos ][ aa ] = AAProb( aa_matrices_[ bb ][ pos ][ aa ] ) *
						bb_boltz_probs()[ pos ][ bb ];
				} else {
					//Add second RotProb multiplied by its bb probability
					//     (*exp_aa_matrix_)[ pos ][ aa ] += aa_matrices_[ bb ][ pos ][ aa ] *
					//       bb_boltz_probs_per_aa()[ bb ][ pos ][ aa ];
					(*exp_aa_matrix_)[ pos ][ aa ] += aa_matrices_[ bb ][ pos ][ aa ] *
						bb_boltz_probs()[ pos ][ bb ];
				}
				//no need for following line unless bb boltz per aa

				//    totals [ pos ] += aa_matrices_[ bb ][ pos ][ aa ].probability() *
				//      bb_boltz_probs_per_aa()[ bb ][ pos ][ aa ];
			} // loops through rot
		} //loops through bb
	} //loops through pos
	//no need for following line unless bb boltz per aa
	// (*exp_aa_matrix_) /= totals;


}

// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that FlexBBDesignMeanField can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, FlexBBDesignMeanField const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace mean_field
}  // namespace protocols
