// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    MeanField.cc

/// @brief   Method definitions for MeanField.
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

// Unit headers
#include <protocols/mean_field/MeanField.hh>

// Package headers
#include <protocols/mean_field/RotProb.hh>
#include <protocols/mean_field/RotMatrix.hh>

// Project headers
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunction.hh>


// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <protocols/mean_field/jagged_array.functions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>

// Numeric headers


// C++ headers
#include <iostream>


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR("protocols.mean_field.MeanField");


namespace protocols {
namespace mean_field {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////

MeanField::MeanField( Size option,
	core::pose::PoseOPs & poses,
	utility::vector1 < pack::task::PackerTaskOP > tasks,
	scoring::ScoreFunctionOP scfxn ) :
	utility::pointer::ReferenceCount(),
	rot_matrix_( new RotMatrix() ),
	alt_rot_matrix_( new RotMatrix() ),
	lambda_memory_( Real( 0.5 ) ),
	tolerance_ ( Real( 0.0001 ) ),
	temperature_( Real( 0.6 ) ),
	init_option_( option ),
	threshold_( Real( 10.0 ) ) ,
	poses_( poses ),
	tasks_( tasks ),
	scfxn_( scfxn )
{}

MeanField::MeanField( Size option,
	pose::PoseOPs & poses,
	utility::vector1 < pack::task::PackerTaskOP > tasks,
	scoring::ScoreFunctionOP scfxn,
	Real lambda_mem,
	Real tolerance,
	Real temp,
	Real threshold ) :
	utility::pointer::ReferenceCount(),
	rot_matrix_( new RotMatrix() ),
	alt_rot_matrix_( new RotMatrix() ),
	lambda_memory_( lambda_mem ),
	tolerance_ ( tolerance ),
	temperature_( temp ),
	init_option_( option ),
	threshold_( threshold ),
	poses_( poses ),
	tasks_( tasks ),
	scfxn_( scfxn )

{}

// Destructor
MeanField::~MeanField() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods

/// @details prints out RotMatrix, Energies Matrix, and Alt RotMatrix
void
MeanField::show(std::ostream & output) const
{
	output << "Conformational Matrix (RotMatrix): " << std::endl;
	rot_matrix()->show( output );
	output << "Energies Matrix: " << std::endl;
	energies_matrix_.show( output );
	output << "Alternate Conformational Matrix (RotMatrix): " << std::endl;
	alt_rot_matrix_->show( output );
}

/// @details runs mean-field algorithm until convergence is reached or 1000 rounds of convergence has been run
/// @remarks 1000 is currently arbitrary - should possibly be changed later
void
MeanField::converge ()
{

	calc_alt_rot_matrix();

	Size counter( 1 );

	while ( ! has_converged() && counter <= Size( 1000 ) ) //doesn't allow it to try to converge ad infinitum. may have to raise max rounds later.
			{



		//can technically use utility overloaded operators in jagged_array.functions.hh
		//but they may be less efficient since requires three nested-for-loops
		for ( Size pos = 1 ; pos <= rot_matrix_->size() ; ++pos ) {
			for ( Size rot = 1 ; rot <= (*rot_matrix_)[ pos ].size() ; ++rot ) {
				(*rot_matrix_)[ pos ][ rot ] *= ( Real (1.0) - lambda_memory() );
				(*alt_rot_matrix_)[ pos ][ rot ] *= lambda_memory();
				(*rot_matrix_)[ pos ][ rot ] += (*alt_rot_matrix_)[ pos ][ rot ];
			}
		}

		calc_alt_rot_matrix();
		++counter;

	}

	if ( has_converged() ) {
		TR.Info << "Following Round " << counter << " of convergence, converged!" << std::endl;

	} else {
		//throws exception which can be caught by FlexBBMeanField's
		std::stringstream error_message;
		error_message
			<< "Unable to converge after " << counter << " rounds of convergence." << std::endl;
		throw utility::excn::EXCN_Msg_Exception(error_message.str());

	}

}

/// @details required so that it can be overridden in derived classes which perform many more functions while processing
void
MeanField::process()
{

	build_rot_matrix_for_pose( 1 );

	converge();

}

/// @details builds a RotMatrix for a specific pose by initializing RotamerSets and using pack_rotamers_setup
/// @details also builds an IG for use in calculating energies_matrix
void
MeanField::build_rot_matrix_for_pose( Size pose_ind )
{

	TR.Info << "Processing pose: " << poses_[ pose_ind ]->pdb_info()->name() << std::endl;

	core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets() );
	core::pack::interaction_graph::AnnealableGraphBaseOP ig = NULL;

	core::pack::pack_rotamers_setup( *poses_[ pose_ind ], *scfxn_, tasks_[ pose_ind ], rotsets, ig );

	// error catching - if someone changes types of graphs returned by interaction_graph_factory assert fails
	//TODO: AR: what if someone uses linmem_ig?
	assert( utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph >( ig ) );
	pig_ = utility::pointer::static_pointer_cast< core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph >( ig );

	rot_matrix_->build_rot_matrix( init_option_, rotsets );
	alt_rot_matrix_->build_rot_matrix( init_option_, rotsets );
	energies_matrix_.clear();

	for ( Size pos = 1; pos <= rot_matrix_->size(); ++pos ) {
		energies_matrix_.push_back( utility::vector1 < Real > ( (*rot_matrix_)[ pos ].size() ) );
	}

}

/// @details delete pose from the list of poses
/// @remarks used if a pose fails to converge
void
MeanField::delete_pose( Size pose_ind )
{
	poses_.erase( poses_.begin() + pose_ind - 1);
	tasks_.erase( tasks_.begin() + pose_ind - 1);
}

/// @remarks virtual because can be overridden by FlexBBMeanField
RotMatrixCOP
MeanField::rot_matrix() const
{
	return rot_matrix_;
}

/// @details uses task_ to determine number of residues to be packed
Size
MeanField::num_packed() const
{
	return tasks_[ 1 ]->num_to_be_packed();
}

// Private methods /////////////////////////////////////////////////////////////
/// @details calculates the next conformational matrix (RM)
/// @details i.e. if RotMatrix is RM(i), prints out RM(i+1)
/// @details first calculates energies matrix based on RM, then converts energies to Boltzmann weights (RM)
void
MeanField::calc_alt_rot_matrix ()
{
	calc_energies_matrix();

	convert_energies_to_alt_rm();
}

/// @details calculates energies matrix from RM and IG
/// @details uses equation EM[i][j] = one_body_energy[ii][jj] + sum(two_body_energy[jj][n] where:
/// @details ii = node (position) of residue
/// @details jj = state (rotamer) of residue
/// @details n = other rotamers that interact with [ii][jj]
/// @details sum is over all rotamers of other residues that interact with [ii][jj]
void
MeanField::calc_energies_matrix ()
{
	//initialize energies_matrix to one_body_energies
	//separate for-loop so as to initialize all values before adding two_body_energies
	for ( Size node = 1; node <= rot_matrix_->size(); ++node ) {
		for ( Size state = 1 ; state <= (*rot_matrix_)[ node ].size(); ++state ) {
			energies_matrix_[ node ][ state ] = pig_->get_one_body_energy_for_node_state( node, state );
			energies_matrix_[node][state] = ( energies_matrix_[node][state] < threshold_ ? energies_matrix_[node][state] : threshold_ );
		}
	}

	//two_body_energy iterate through one node->edge at a time and sum all energies associated with both sides of the edge
	for ( Size node = 1; node <= rot_matrix_->size(); ++node ) {

		for (  pig_->reset_edge_list_iterator_for_node( node ) ; !pig_->edge_list_iterator_at_end() ; pig_->increment_edge_list_iterator() ) {
			//no need for assert on dynamic cast as already checked that pig_ is FixedBBInteractionGraph
			core::pack::interaction_graph::FixedBBEdge const & pedge =
				static_cast< core::pack::interaction_graph::FixedBBEdge const & > ( pig_->get_edge() );
			Size node2 = pedge.get_other_ind( node );

			//no double-counting edges
			if ( node2 == node || node2 < node ) continue;

			for ( Size state = 1; state <= (*rot_matrix_)[ node ].size(); ++state ) {

				for ( Size state2 = 1 ; state2 <= (*rot_matrix_)[ node2 ].size(); ++state2 ) {

					Real two_body_energy = pedge.get_two_body_energy( state, state2 );
					two_body_energy = ( two_body_energy < threshold_ ? two_body_energy : threshold_ );

					energies_matrix_[ node ][ state ] += two_body_energy * (*rot_matrix_)[ node2 ][ state2 ].probability();

					energies_matrix_[ node2 ][ state2 ] += two_body_energy * (*rot_matrix_)[ node ][ state ].probability();
				}
			}
		}

	}

}

/// @details calculates Boltzmann weights from energies_matrix_ to form alt_rot_matrix_
void
MeanField::convert_energies_to_alt_rm ()
{
	//doesn't use jagged_array get_totals_columns because this is faster
	utility::vector1 < Real > totals( energies_matrix_.size(), core::Real( 0.0 ) );

	for ( Size pos = 1 ; pos <= energies_matrix_.size() ; ++pos ) {
		for ( Size rot = 1 ; rot <= energies_matrix_[ pos ].size() ; ++rot ) {
			Real energy = energies_matrix_[ pos ][ rot ] < threshold_ ? energies_matrix_[ pos ][ rot ] : threshold_;

			(*alt_rot_matrix_)[ pos ][ rot ].probability( exp( ( energy * Real( -1.0 ) ) / temperature_ ) );

			totals[ pos ] += (*alt_rot_matrix_)[ pos ][ rot ].probability();
		}
	}

	*alt_rot_matrix_ /= totals;
}

/// @details checks if algorithm has converged by checking if RMSmat < tolerance
bool
MeanField::has_converged () const
{

	Real total = 0;

	for ( Size pos = 1 ; pos <= energies_matrix_.size() ; ++pos ) {
		for ( Size rot = 1 ; rot <= energies_matrix_[ pos ].size() ; ++rot ) {
			total += pow ( ( (*alt_rot_matrix_)[ pos ][ rot ].probability() - (*rot_matrix_)[ pos ][ rot ].probability() ), 2.0 );
		}
	}

	Real result = sqrt( total );

	return result < tolerance_;
}

// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that MeanField can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, MeanField const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace mean_field
}  // namespace protocols
