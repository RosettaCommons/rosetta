// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/RotamerAssigningAnnealer.cc
/// @brief  Annealer class that assigns rotamers implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/annealer/RotamerAssigningAnnealer.hh>

// Package headers
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/Fmath.hh>

/// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

/// Utility headers
#include <utility/exit.hh>

/// C++ headers
#include <iostream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1.hh>


using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace annealer {

RotamerAssigningAnnealer::RotamerAssigningAnnealer(
	int num_rots_to_pack_in,
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current,
	FixbbRotamerSetsCOP rotamer_sets,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq
) :
	SimAnnealerBase(
	num_rots_to_pack_in,
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	),
	rotamer_sets_( rotamer_sets ),
	rot_to_pack_( num_rots_to_pack_in ),
	current_to_pick_( 1 ),
	n_assigned_at_start_( 0 ),
	assign_state_to_all_nodes_immediately_( false )
{
	for ( unsigned int ii(1); ii <= rot_to_pack_.size(); ++ii ) {
		rot_to_pack_[ ii-1 ] = ii;
	}
	setup_rots_for_node( rotamer_sets );
}


RotamerAssigningAnnealer::RotamerAssigningAnnealer(
	utility::vector0< int > & rot_to_pack,
	int num_rots_to_pack_in,
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current,
	FixbbRotamerSetsCOP rotamer_sets,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq
) :
	SimAnnealerBase(
	num_rots_to_pack_in,
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	),
	rotamer_sets_( rotamer_sets ),
	current_to_pick_( 1 ),
	n_assigned_at_start_( 0 ),
	assign_state_to_all_nodes_immediately_( false )
{
	//ja no need for overloading: just resort to generic behavior if rot_to_pack is empty.  Should just have one constructor with an empty default argument for rot_to_pack
	if ( !rot_to_pack.empty() ) rot_to_pack_ = rot_to_pack;
	else {
		// rot_to_pack obtained was empty, use all rotamers
		num_rots_to_pack( rotamer_sets->nrotamers() ); // important
		rot_to_pack_.reserve( num_rots_to_pack() );
		for ( Size rot = 1; rot <= num_rots_to_pack(); ++rot ) {
			rot_to_pack_.push_back( rot );
		}
	}
	setup_rots_for_node( rotamer_sets );
}

/// @details if you want to make a single rotamer substitution for a single node
/// then you'll want to grab a random rotamer for that node -- this function creates
/// a mapping from moltenres id to the set of accessible rotamers for that moltenres.
/// It also makes sure that the user has provided a valid set of rotamers: each molten
/// residue must have at least one pickable rotamer, or the final rotamer assignment
/// will have some molten residues in an unassigned state.  Bad user input is most common
/// when the annealer is constructed using the rot_to_pack vector.
void
RotamerAssigningAnnealer::setup_rots_for_node(
	FixbbRotamerSetsCOP rotamer_sets
)
{
	int const nmoltres = rotamer_sets->nmoltenres();
	utility::vector1< int > npickable_rots_for_node( nmoltres, 0 );
	rots_for_nodes_.resize( nmoltres );
	for ( Size ii = 0; ii < rot_to_pack_.size(); ++ii ) {
		++npickable_rots_for_node[ rotamer_sets->moltenres_for_rotamer( rot_to_pack_[ ii ] ) ];
	}
	for ( int ii = 1; ii <= nmoltres; ++ii ) {
		/// Error detection -- make sure each molten residue has at least one accessible rotamer
		/// before SA begins.  No point in waiting until after packing has ended to realize that
		/// the final state assignment has unassigned states.
		if ( npickable_rots_for_node[ ii ] == 0 ) {
			std::cerr << "ERROR: No pickable rotamers for molten residue " << ii << " ( seqpos= " <<
				rotamer_sets->moltenres_2_resid( ii ) << " ) found in rots_to_pack_." << std::endl;
			std::cerr << "ERROR: Each molten residue must have at least one pickable rotamer." << std::endl;
			std::cerr << "Provided rotamers: (rotid, moltenresid, moltenres_rotid ) " << std::endl;
			for ( Size jj = 0; jj < rot_to_pack_.size(); ++jj ) {
				std::cerr << "( " << jj <<", "
					<< rotamer_sets->moltenres_for_rotamer( rot_to_pack_[ ii ] ) << ", "
					<< rotamer_sets->rotid_on_moltenresidue( rot_to_pack_[ ii ] ) << ") ";
				if ( jj % 4 == 3 ) std::cerr << "\n";
			}
			std::cerr << "\n Exiting." << std::endl;
			utility_exit();
		}
		rots_for_nodes_.reserve( npickable_rots_for_node[ ii ] );
	}
	for ( Size ii = 0; ii < rot_to_pack_.size(); ++ii ) {
		rots_for_nodes_[ rotamer_sets->moltenres_for_rotamer( rot_to_pack_[ ii ] ) ].
			push_back( rot_to_pack_[ ii ] );
	}
}

RotamerAssigningAnnealer::~RotamerAssigningAnnealer()
{}

/// @brief pick a rotamer from a list
///
/// @details if no rotamer is available, return a nonsense number, like -1
int RotamerAssigningAnnealer::pick_a_rotamer( int cycle )
{
	bool start_with_current = get_start_with_current();
	int ranrotamer = -1;
	if ( assign_state_to_all_nodes_immediately_ && core::Size(n_assigned_at_start_) < rots_for_nodes_.size() ) {
		++n_assigned_at_start_;
		return pick_a_rotamer_for_node( n_assigned_at_start_ );
	}

	//bk if quench cycle, pass through all rotamers before
	//bk repeating a rotamer
	if ( quench() ) {
		int num =  mod( cycle - 1, (int) num_rots_to_pack());
		if ( num == 0 ) {
			numeric::random::random_permutation( rot_to_pack_, numeric::random::rg() );
		}
		ranrotamer = rot_to_pack_.at(num);
		//bk if start of run and start_with_current is true then first nres
		//bk iterations will be used to place the current rotamers
	} else if ( start_with_current && (unsigned int) current_to_pick_ <= current_rot_index().size() ) {

		if ( current_to_pick_ == 1  ) {
			std::cout << "RotamerAssigningAnnealer: recovering current rotamers for start." << std::endl;
		}
		if ( current_rot_index()(current_to_pick_) != -1 ) {
			ranrotamer = current_rot_index()(current_to_pick_);
		} else {
			std::cout << "RotamerAssigningAnnealer: unable to recover current rotamer for moltenres: " << current_to_pick_ << std::endl;
			ranrotamer = -1;
		}
		++current_to_pick_;
	} else {
		ranrotamer = rot_to_pack_.at(static_cast<int>( num_rots_to_pack() * numeric::random::rg().uniform() ));
	}
	//std::cerr << "ranrotamer: " << ranrotamer << " total: " << current_rot_index_.size() << std::endl;
	return ranrotamer;
}

int RotamerAssigningAnnealer::pick_a_rotamer_for_node( int node ) const
{
	return rots_for_nodes_[node][ static_cast<int> ( rots_for_nodes_[ node ].size() * numeric::random::rg().uniform() + 1 ) ];
}

void RotamerAssigningAnnealer::set_assign_state_to_all_nodes_immediately( bool setting )
{
	assign_state_to_all_nodes_immediately_ = setting;
}

utility::vector0< int > const &
RotamerAssigningAnnealer::rot_to_pack() const
{
	return rot_to_pack_;
}

RotamerAssigningAnnealer::FixbbRotamerSetsCOP
RotamerAssigningAnnealer::rotamer_sets() const { return rotamer_sets_; }


} //end namespace annealer
} //end namespace pack
} //end namespace core
