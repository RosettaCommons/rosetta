// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /protocols/simple_moves/chiral/ChiralMover.hh
/// @brief
/// @author Kevin Drew, kdrew@nyu.edu
#ifndef INCLUDED_protocols_simple_moves_chiral_ChiralMover_hh
#define INCLUDED_protocols_simple_moves_chiral_ChiralMover_hh
// Unit Headers
#include <protocols/simple_moves/chiral/ChiralMover.fwd.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Patch.hh>

#include <map>

// Utility Headers
//#include <core/types.hh>
//#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace chiral {

enum Chirality {
	L_CHIRALITY=1,
	D_CHIRALITY,
	FLIP_CHIRALITY //flip to the other chirality, ie L->D or D->L
};


bool is_d_chiral( core::chemical::ResidueType restype );
bool is_l_chiral( core::chemical::ResidueType restype );
core::chemical::ResidueType const & get_chiral_residue_type( core::chemical::ResidueType const & , Chirality );

///@details
class ChiralMover : public protocols::moves::Mover {

public:

	///@brief
	ChiralMover( core::Size chiral_seq_position );
	ChiralMover( core::Size chiral_seq_position, Chirality chirality );
	ChiralMover( core::Size chiral_seq_position, bool orient_functional_group );
	ChiralMover( core::Size chiral_seq_position, Chirality chirality, bool orient_functional_group );

	virtual ~ChiralMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	
private:

	core::Size const chiral_seq_pos_;
	Chirality const chirality_;
	bool orient_functional_group_;

};//end ChiralMover


}//namespace chiral
}//namespace simple_moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_chiral_ChiralMover_hh

