// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    FlexBBMeanField.cc

/// @brief   Method definitions for FlexBBMeanField.
/// @author  arubenstein

// Unit headers
#include <protocols/mean_field/FlexBBMeanField.hh>

// Package headers
#include <protocols/mean_field/RotProb.hh>
#include <protocols/mean_field/RotMatrix.hh>
#include <protocols/mean_field/ResHashMap.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mean_field.OptionKeys.gen.hh>

// Numeric headers

// Boost headers
#include <boost/unordered_map.hpp>

// C++ headers
#include <iostream>


// Construct tracer.
static basic::Tracer TR("protocols.mean_field.FlexBBMeanField");


namespace protocols {
namespace mean_field {

using namespace core;
using namespace utility;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
/// @details standard constructor initializes using base class constructor and resizes all jagged_arrays/vectors
FlexBBMeanField::FlexBBMeanField(core::uint option, core::pose::PoseOPs & poses, utility::vector1 < core::pack::task::PackerTaskOP > tasks, core::scoring::ScoreFunctionOP scfxn ) :
	MeanField::MeanField( option, poses, tasks, scfxn ),
	exp_rot_matrix_( new RotMatrix )
{
	rot_matrices_.resize( poses.size() );
	energy_matrices_.resize( poses.size() );
	bb_boltz_weights_.resize( tasks[ 1 ]->num_to_be_packed(), vector1 < Real > ( poses.size(), Real( 0.0 ) ) );
	bb_boltz_probs_.resize( tasks[ 1 ]->num_to_be_packed(), vector1 < Real > ( poses.size(), Real( 0.0 ) ) );
	nrot_per_pos_.resize( tasks[ 1 ]->num_to_be_packed() );
	bb_boltz_probs_per_aa_.resize( poses.size(), jagged_array < Real > ( tasks[ 1 ]->num_to_be_packed(),
		vector1 < Real > ( chemical::num_canonical_aas, Real( 0.0 ) ) ) );
}

/// @details constructor initializes using base class constructor and resizes all jagged_arrays/vectors
/// @details initializes parameters to input values
FlexBBMeanField::FlexBBMeanField(core::uint option, core::pose::PoseOPs & poses, utility::vector1 < core::pack::task::PackerTaskOP > tasks, core::scoring::ScoreFunctionOP scfxn,
	Real lambda_mem, Real tolerance, Real temp, Real threshold) :
	MeanField::MeanField( option, poses, tasks, scfxn, lambda_mem, tolerance, temp, threshold ),
	exp_rot_matrix_( new RotMatrix )
{
	rot_matrices_.resize( poses.size() );
	energy_matrices_.resize( poses.size() );
	bb_boltz_weights_.resize( tasks[ 1 ]->num_to_be_packed(), vector1 < Real > ( poses.size(), Real( 0.0 ) ) );
	bb_boltz_probs_.resize( tasks[ 1 ]->num_to_be_packed(), vector1 < Real > ( poses.size(), Real( 0.0 ) ) );
	nrot_per_pos_.resize( tasks[ 1 ]->num_to_be_packed() );
	bb_boltz_probs_per_aa_.resize( poses.size(), jagged_array < Real > ( tasks[ 1 ]->num_to_be_packed(),
		vector1 < Real > ( chemical::num_canonical_aas, Real( 0.0 ) ) ) );
}

// Destructor
FlexBBMeanField::~FlexBBMeanField() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods

/// @details prints out the exp_rot_matrix_
void
FlexBBMeanField::show(std::ostream & output) const
{
	output << "Averaged Conformational Matrix (RotMatrix): " << std::endl;
	rot_matrix()->show( output );
}

/// @details for each pose, tries to build a RotMatrix, converge, and save values of RotMatrix and energy_matrix
/// @details to rot_matrices and energy_matrices, then renumbers the rotamers, calculates the Backbone Boltzmann Probabilities
/// @details and calculates the averaged RotMatrix
/// @remarks if there are errors in building a RotMatrix or if the pose did not converge, that pose is deleted
void
FlexBBMeanField::process()
{

	for ( core::Size bb = 1; bb <= num_poses(); ++bb ) {
		try
{
			build_rot_matrix_for_pose( bb );
			converge();
			rot_matrices_[ bb ] = *( MeanField::rot_matrix() );
			energy_matrices_[ bb ] = energies_matrix();
		} catch ( utility::excn::EXCN_Msg_Exception &excn )
{
			TR.Info << "Pose did not converge, deleted." << std::endl;
			delete_pose( bb );
			--bb;
		}
	}

	if ( num_poses() == Size( 0 ) ) {
		std::stringstream error_message;
		error_message
			<< "No poses reached convergence." << std::endl;
		throw utility::excn::EXCN_Msg_Exception( error_message.str() );
	}

	TR.Info << "Total number of final poses " << num_poses() << std::endl;

	//renumber_rotamers();
	calc_bb_boltz_probs();
	//calc_exp_value_rot_matrix();
}


// Private methods /////////////////////////////////////////////////////////////

/// @details overrides MeanField::delete_pose to resize appropriate data members
/// @details used to delete a pose which does not reach convergence
/// @remarks meant to be called before rot_matrices_ and energy_matrices_ have been fully assigned
/// @remarks and before any bb_boltz_weights or nrot_per_pos have been assigned at all
void
FlexBBMeanField::delete_pose( core::Size pose_ind )
{
	MeanField::delete_pose( pose_ind );

	rot_matrices_.resize( num_poses() );
	energy_matrices_.resize( num_poses() );

	bb_boltz_weights_.assign( num_packed(), vector1 < Real > ( num_poses(), Real( 0.0 ) ) );
	bb_boltz_probs_.assign( num_packed(), vector1 < Real > ( num_poses(), Real( 0.0 ) ) );
	bb_boltz_probs_per_aa_.assign( num_poses(), jagged_array < Real > ( num_packed(),
		vector1 < Real > ( chemical::num_canonical_aas, Real( 0.0 ) ) ) );
	nrot_per_pos_.resize( num_packed() );
}

/// @details overrides MeanField::rot_matrix to return the averaged RotMatrix
RotMatrixCOP
FlexBBMeanField::rot_matrix() const
{
	return exp_rot_matrix_;
}

/// @details renumbers the rotamers using a ResHashMap
/// @details inserts rotamers into the ResHashMap using a hash function that rounds chi angles down to nearest 10 degrees
/// @details before attempting an insert, check if a rotamer has already been inserted and if so, what its rotamer index is
void
FlexBBMeanField::renumber_rotamers()
{

	//iterates over the positions in the RotMatrix
	for ( Size pos = 1; pos <= rot_matrices_[1].size(); ++pos ) {

		ResHashMap rothash;

		//iterates over all backbones
		for ( Size bb = 1; bb <= rot_matrices_.size(); ++bb ) {

			//iterates over all rots for pos in bb
			for ( Size rot = 1; rot <= rot_matrices_[ bb ][ pos ].size(); ++rot ) {
				SSize rot_ind = rothash.attempt_insert( rot_matrices_[ bb ][ pos ][ rot ].res() );

				rot_matrices_[ bb ][ pos ][ rot ].rot_ind( rot_ind );

			}

		}
		nrot_per_pos_[ pos ] = rothash.last_ind_assigned();
	}


}

/// @details calculates Boltzmann probabilities of each position in each backbone, relative to the other backbones for that position
/// @details uses the sum of energies across all rotamers for each position at a backbone as the E(i) in P(i) = -e^(E(i)/kT)
/// @remarks bb_boltz_weights_ = -e^(E(i)/kT)
/// @remarks bb_boltz_probs_ = -e^(E(i)/kT) / Z where Z is the sum of -e^(E(i)/kT) for that position at each backbone
void
FlexBBMeanField::calc_bb_boltz_probs()
{

	// can use rot_matrices_[1] to determine # of positions since that is constant across the rot_matrices
	vector1 < Real > total_bb_boltz_prob( rot_matrices_[1].size(), Real( 0.0 ) );

	// necessary for nrot normalization of bb_boltz_probs_per_aa
	jagged_array < Size > nrot_per_aa_per_pos( rot_matrices_[1].size(),
		vector1< Size > ( chemical::num_canonical_aas, Size( 0 ) ) );
	jagged_array < Real > total_bb_boltz_prob_per_aa( rot_matrices_[1].size(),
		vector1< Real > ( chemical::num_canonical_aas, Real( 0.0 ) ) );

	//loops through backbones
	for ( Size bb = 1; bb <= rot_matrices_.size(); ++bb ) {
		//loops through positions
		for ( Size pos = 1; pos <= energy_matrices_[ bb ].size() ; ++pos ) {
			//loops through rotamers
			for ( Size rot = 1; rot <= energy_matrices_[ bb ][ pos ].size(); ++rot ) {

				core::chemical::AA aa = rot_matrices_[bb][pos][rot].aa_ind();

				//    bb_boltz_weights_[ pos ][ bb ] += energy_matrices_[ bb ][ pos ][ rot ];

				bb_boltz_probs_per_aa_[ bb ][ pos ][ aa ] += energy_matrices_[ bb ][ pos ][ rot ] * rot_matrices_[bb][pos][rot].probability();
				nrot_per_aa_per_pos[ pos ][ aa ] += 1;
			}

			for ( core::Size aa = 1; aa <= bb_boltz_probs_per_aa_[bb][pos].size(); ++aa ) {
				//if nrot is still 0 set it to 1 to avoid dividing by 0
				nrot_per_aa_per_pos[ pos ][ aa ] = nrot_per_aa_per_pos[ pos ][ aa ] == 0 ? 1 : nrot_per_aa_per_pos[ pos ][ aa ];
				bb_boltz_weights_[ pos ][ bb ] += bb_boltz_probs_per_aa_[bb][pos][aa] /
					pow( nrot_per_aa_per_pos[ pos ][ aa ], basic::options::option[ basic::options::OptionKeys::mean_field::bb_average_weight ] ) ;
				bb_boltz_probs_per_aa_[ bb ][ pos ][ aa ] = exp( ( bb_boltz_probs_per_aa_[bb][pos][aa] /
					pow( nrot_per_aa_per_pos[ pos ][ aa ], basic::options::option[ basic::options::OptionKeys::mean_field::bb_average_weight ] ) ) *
					core::Real(-1.0) / temperature() );
				total_bb_boltz_prob_per_aa[ pos ][ aa ] += bb_boltz_probs_per_aa_[bb][pos][aa];
			}

			//E/K
			//   bb_boltz_weights_[ pos ][ bb ] = exp( ( bb_boltz_weights_[ pos ][ bb ] / energy_matrices_[ bb ][ pos ].size() ) * Real( -1.0 ) / temperature() );
			//E(truncated)/K^RNW
			//            bb_boltz_weights_[ pos ][ bb ] = ( bb_boltz_weights_[ pos ][ bb ] < threshold() ? bb_boltz_weights_[ pos ][ bb ] : threshold() );
			//   bb_boltz_weights_[ pos ][ bb ] = exp( ( bb_boltz_weights_[ pos ][ bb ] / pow ( energy_matrices_[ bb ][ pos ].size(),
			//     basic::options::option[ basic::options::OptionKeys::mean_field::rot_norm_weight ] ) ) * Real( -1.0 ) / temperature() );
			//E/K^RNW * 100.0
			//   bb_boltz_weights_[ pos ][ bb ] = exp( ( bb_boltz_weights_[ pos ][ bb ] / pow ( energy_matrices_[ bb ][ pos ].size(),
			//     basic::options::option[ basic::options::OptionKeys::mean_field::rot_norm_weight ] ) ) * Real( -1.0 ) / ( Real( 100.0 ) * ( temperature() ) ) );
			bb_boltz_weights_[ pos ][ bb ] = exp( ( bb_boltz_weights_[ pos ][ bb ] ) * Real( -1.0 ) / ( ( temperature() ) ) );
			//GC expression - not using for now
			//   bb_boltz_weights_[ pos ][ bb ] = exp( ( energy_matrices_[ bb ][ pos ].size() *
			//     basic::options::option[ basic::options::OptionKeys::mean_field::rot_norm_weight ] -
			//     bb_boltz_weights_[ pos ][ bb ] ) /
			//     ( temperature() * 100.0 ) );

			//            TR << "bb_boltz post-: " << pos << " " << bb << " " << bb_boltz_weights_[ pos ][ bb ] << std::endl;
			//            TR << "power: " << pos << " " << bb << " " << pow ( energy_matrices_[ bb ][ pos ].size(), basic::options::option[ basic::options::OptionKeys::mean_field::rot_norm_weight ] ) << std::endl;
			//            TR << "Temp: " << pos << " " << bb << " " <<  temperature() << std::endl;
			total_bb_boltz_prob[ pos ] += bb_boltz_weights_[ pos ][ bb ];
			//            TR << "total_bb_boltz prob: " << pos << " " << bb << " " << total_bb_boltz_prob[ pos ] << std::endl;
		}
	}

	bb_boltz_probs_ = bb_boltz_weights_ / total_bb_boltz_prob;


	for ( core::Size bb = 1; bb <= bb_boltz_probs_per_aa_.size(); ++bb ) {

		bb_boltz_probs_per_aa_[bb] /= total_bb_boltz_prob_per_aa;
	}

}

/// @details calculate "expected value" or averaged RotMatrix
/// @details weight each position-backbone-rot probability by the Boltzmann probability of that backbone occurring for that position
/// @details and then sum them up across all backbones for that position
void
FlexBBMeanField::calc_exp_value_rot_matrix()
{
	jagged_array < Real > total_bb_boltz_wt_per_rot;

	//added 6/26/15
	exp_rot_matrix_->is_designed( rot_matrices_[1].is_designed() );

	for ( Size pos = 1; pos <= rot_matrices_[1].size(); ++pos ) {
		exp_rot_matrix_->push_back( vector1 < RotProb> ( nrot_per_pos_[ pos ] ) );
		//added 6/23/15
		exp_energy_matrix_.push_back( vector1 < Real > ( nrot_per_pos_[ pos ], Real( 0.0 ) ) );
		total_bb_boltz_wt_per_rot.push_back( vector1 < Real> ( nrot_per_pos_[ pos ], Real( 0.0 ) ) );

		for ( Size bb = 1; bb <= rot_matrices_.size(); ++bb ) {
			for ( Size rot = 1; rot <= rot_matrices_[ bb ][ pos ].size(); ++rot ) {
				Size rot_ind = rot_matrices_[ bb ][ pos ][ rot ].rot_ind();

				//This RotProb wasn't yet initialized so initialize it to value of current RotProb multiplied by the bb probability
				if ( ! (*exp_rot_matrix_)[ pos ][ rot_ind ].res() ) {
					// uses bb_boltz_weights_ since each rotamer has to be normalized separately
					(*exp_rot_matrix_)[ pos ][ rot_ind ] = RotProb( rot_matrices_[ bb ][ pos ][ rot ] ) * bb_boltz_weights_[ pos ][ bb ];
				} else {
					//Add second RotProb multiplied by its bb probability
					(*exp_rot_matrix_)[ pos ][ rot_ind ] += rot_matrices_[ bb ][ pos ][ rot ] * bb_boltz_weights_[ pos ][ bb ];
				}
				//added 6/23/15
				exp_energy_matrix_[ pos ][ rot_ind ] += energy_matrices_[ bb ][ pos ][ rot ] * bb_boltz_weights_[ pos ][ bb ];

				total_bb_boltz_wt_per_rot[ pos ][ rot_ind ] += bb_boltz_weights_[ pos ][ bb ];
			} // loops through rot
		} //loops through bb
	} //loops through pos

	*exp_rot_matrix_ /= total_bb_boltz_wt_per_rot;
	//added 6/23/15
	exp_energy_matrix_ /= total_bb_boltz_wt_per_rot;

}



// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that FlexBBMeanField can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, FlexBBMeanField const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace mean_field
}  // namespace protocols
