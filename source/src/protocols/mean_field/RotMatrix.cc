// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    RotMatrix.cc

/// @brief   Method definitions for RotMatrix.
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

// Unit headers
#include <protocols/mean_field/RotMatrix.hh>

// Package headers

// Project headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/conformation/Residue.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
RotMatrix::RotMatrix() : jagged_array< RotProb > ()
{
	curr_rot();
}

// Standard methods ////////////////////////////////////////////////////////////
/// @details standard constructor calls private member function
/// @param [in] option - determines which numbers to initialize RotMatrix with
/// @param [in] rs - RotamerSetsOP used to gather the rest of the information used to build the RotMatrix
RotMatrix::RotMatrix( Size option, core::pack::rotamer_set::RotamerSetsOP rs ) : jagged_array < RotProb > ()
{
	init(option, rs);
}

/// @details Copy constructor calls private member function
RotMatrix::RotMatrix(RotMatrix const & object_to_copy) : jagged_array < RotProb > ()
{
	copy_data( *this, object_to_copy );
}

/// @details Assignment operator calls private member function
RotMatrix &
RotMatrix::operator=(RotMatrix const & object_to_copy)
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	copy_data( *this, object_to_copy );
	return *this;
}

// Destructor
RotMatrix::~RotMatrix() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods

/// @details prints representation of RotMatrix using show() function of RotProbs
/// @details in the outputted table, each column represents the probabilities of rotamers at that position
/// @details if current rotamers are set, outputs a vector of current rotamers and their probability
void
RotMatrix::show( std::ostream & output ) const
{
	output << "Conformational Matrix (RotMatrix)---------------------------------" << std::endl;

	//uses derived method of jagged_array and show function of RotProb to print
	jagged_array<RotProb>::show( output );

	//outputs probabilities of current rotamers if current rotamer vector is set (i.e. -use_input_sc is used)
	utility::vector1 < RotProb > curr_rot_p = curr_rot_prob();

	bool printed_header = false;

	for ( utility::vector1<RotProb>::const_iterator iter = curr_rot_p.begin(); iter != curr_rot_p.end(); ++iter ) {
		//checks that RotProb is meaningful before outputting data
		if ( iter->pos() != Size( 0 ) ) {
			if ( ! printed_header ) {
				output << "Curr Rotamer Probability------------------------------" << std::endl;
				printed_header = true;
			}
			output << ( *iter ) << "\t";
		}
	}
	output << std::endl;
}

/// @details counts the number of positions for which design is on
Size
RotMatrix::n_designed() const
{
	Size n_des( 0 );

	for ( Size n = 1; n <= is_designed_.size(); ++n ) {
		if ( is_designed_[ n ] ) ++n_des;
	}

	return n_des;
}

/// @details returns vector with probabilities of the current rotamers based on curr_rot_ vector
utility::vector1 < RotProb >
RotMatrix::curr_rot_prob() const
{
	utility::vector1 < RotProb > rp ( curr_rot_.size() );
	for ( Size pos = 1; pos <= curr_rot_.size(); ++pos ) {
		if ( curr_rot_[ pos ] != 0 ) {
			rp[ pos ] = ( *this )[ pos ][ curr_rot_[ pos ] ];
		}
	}
	return rp;
}

void
RotMatrix::build_rot_matrix( core::Size const option, core::pack::rotamer_set::RotamerSetsOP rs )
{
	init( option, rs );
}

// Private methods /////////////////////////////////////////////////////////////
// Initialize data members.
/// @details private method used to initialize RotMatrix as jagged_array of RotProbs
/// @details called by constructors and build_rot_matrix
/// @details also builds curr_rot_ and is_designed_ vectors
/// @param [in] option - determines which numbers to initialize RotMatrix with
/// @param [in] rs - RotamerSetsOP used to gather the rest of the information used to build the RotMatrix
void
RotMatrix::init( Size const option, core::pack::rotamer_set::RotamerSetsOP rs )
{

	//clears variables in order to reinitialize
	clear();
	curr_rot_.clear();
	is_designed_.clear();

	//resizes vectors to appropriate sizes
	is_designed_.resize( rs->nmoltenres(), false );

	//checks to see if current rotamers are set (i.e. -use_input_sc is true)
	if ( rs->rotamer_set_for_moltenresidue( 1 )->id_for_current_rotamer() != 0 ) {
		curr_rot_.resize( rs->nmoltenres() );
	}

	//iterates through molten residues to create a vector for each residue of RotProbs
	for ( Size pos = 1; pos <= rs->nmoltenres(); ++pos ) {
		//if ( rs->task()->design_residue( rs->moltenres_2_resid( pos ) ) )

		Size nrot = rs->rotamer_set_for_moltenresidue( pos )->num_rotamers();
		utility::vector1 < RotProb > rps ( nrot );

		push_back( rps );

		//setup current rotamer vector
		if ( curr_rot_.size() > 0 ) {
			Size cr = rs->rotamer_set_for_moltenresidue( pos )->id_for_current_rotamer();
			if ( cr == 0 ) {
				curr_rot_.clear();
			} else {
				curr_rot_[ pos ] = cr;
			}
		}

		//sets initial value based on option
		Real init_val;
		switch ( option ) {
		case 1 :
			init_val = Real( 1.0 ) / nrot; //option 1 sets initial value to 1/nrot
			break;
		default :
			init_val = Real( 0.0 );
			break;
		}

		//initializes RotProb and sets is_designed_ accordingly
		core::chemical::AA aa = rs->rotamer_for_moltenres( pos, 1 )->type().aa();

		for ( Size rot = 1; rot <= nrot ; ++rot ) {
			( *this )[ pos ][ rot ] = RotProb( init_val, rot, rs->moltenres_2_resid( pos ), rs->rotamer_for_moltenres( pos, rot ) );
			if ( aa != rs->rotamer_for_moltenres( pos, rot )->type().aa() ) {
				is_designed_[ pos ] = true;
			}
		}
	}

}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RotMatrix::copy_data(
	RotMatrix & object_to_copy_to,
	RotMatrix const & object_to_copy_from)
{
	object_to_copy_to.clear();
	for ( Size i = 1; i <= object_to_copy_from.size(); ++i ) {
		object_to_copy_to.push_back ( object_to_copy_from[i] );
	}
	object_to_copy_to.curr_rot_ = object_to_copy_from.curr_rot_;
	object_to_copy_to.is_designed_ = object_to_copy_from.is_designed_;
}

// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that RotMatrix can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, RotMatrix const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}  // namespace mean_field
}  // namespace protocols
