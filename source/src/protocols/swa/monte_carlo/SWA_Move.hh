// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/monte_carlo/SWA_Move.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_monte_carlo_SWA_Move_HH
#define INCLUDED_protocols_swa_monte_carlo_SWA_Move_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/monte_carlo/SWA_Move.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <string>
#include <ostream>

using namespace core;

namespace protocols {
namespace swa {
namespace monte_carlo {

	typedef  utility::vector1< Size >  Chunk; // move to fwd file?

	enum MovingResidueCase { NO_CASE = 0, CHAIN_TERMINUS_5PRIME, CHAIN_TERMINUS_3PRIME, INTERNAL, FLOATING_BASE, LAST_MOVING_RESIDUE_CASE };

	enum AddOrDeleteChoice{ NO_ADD_OR_DELETE = 0, ADD, DELETE, LAST_ADD_OR_DELETE_CHOICE };

	std::string to_string( MovingResidueCase const & moving_residue_case );

	std::string to_string( AddOrDeleteChoice const & add_or_delete_choice_name );

	class SWA_Move: public utility::pointer::ReferenceCount {

	public:

		//constructor
		SWA_Move();

		SWA_Move( Chunk const & chunk,
									MovingResidueCase const & moving_residue_case,
									AddOrDeleteChoice const & add_or_delete_choice );

		SWA_Move( SWA_Move const & src );

		SWA_Move &
		operator=( SWA_Move const & src );

		//destructor
		~SWA_Move();

	public:

		void set_chunk( Chunk const & setting ){ chunk_ = setting; }
		Chunk chunk() const{ return chunk_; }

		void set_moving_residue_case( MovingResidueCase const & setting ){ moving_residue_case_ = setting; }
		MovingResidueCase moving_residue_case() const{ return moving_residue_case_; }

		void set_add_or_delete_choice( AddOrDeleteChoice const & setting ){ add_or_delete_choice_ = setting; }
		AddOrDeleteChoice add_or_delete_choice() const{ return add_or_delete_choice_; }

	private:

		Chunk chunk_;
		// later may want to add information on residues that are attached at 3' or 5' ends -- following wouldn't be a single Case
		// but instead some kind of a map or list.
		MovingResidueCase moving_residue_case_;
		AddOrDeleteChoice add_or_delete_choice_;

	};

	std::ostream &
	operator <<( std::ostream & os, SWA_Move const swa_move );

} //monte_carlo
} //swa
} //protocols

#endif
