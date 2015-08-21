// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/methods/FreeDOF_Options.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_methods_FreeDOF_Options_HH
#define INCLUDED_core_scoring_methods_FreeDOF_Options_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/methods/FreeDOF_Options.fwd.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace methods {

class FreeDOF_Options: public utility::pointer::ReferenceCount {

public:

	//constructor
	FreeDOF_Options();

	//destructor
	~FreeDOF_Options();

	void initialize_from_options();

public:

	core::Real free_suite_bonus() const { return free_suite_bonus_; }
	void free_suite_bonus( core::Real setting ){ free_suite_bonus_ = setting; }

	core::Real free_2HOprime_bonus() const { return free_2HOprime_bonus_; }
	void free_2HOprime_bonus( core::Real setting ){ free_2HOprime_bonus_ = setting; }

	core::Real free_sugar_bonus() const { return free_sugar_bonus_; }
	void free_sugar_bonus( core::Real setting ){ free_sugar_bonus_ = setting; }

	core::Real pack_phosphate_penalty() const { return pack_phosphate_penalty_; }
	void pack_phosphate_penalty( core::Real setting ){ pack_phosphate_penalty_ = setting; }

	core::Real free_side_chain_bonus() const { return free_side_chain_bonus_; }
	void free_side_chain_bonus( core::Real setting ){ free_suite_bonus_ = setting; }

	friend
	bool
	operator==( FreeDOF_Options const & a, FreeDOF_Options const & b );

	friend
	std::ostream &
	operator<< ( std::ostream & out, const FreeDOF_Options & options );

	void
	show( std::ostream & out ) const;

private:

	core::Real free_suite_bonus_;
	core::Real free_2HOprime_bonus_;
	core::Real free_sugar_bonus_;
	core::Real pack_phosphate_penalty_;
	core::Real free_side_chain_bonus_;

};

} //methods
} //scoring
} //core

#endif
