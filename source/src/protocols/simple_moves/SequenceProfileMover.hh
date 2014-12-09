// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/SequenceProfileMover.hh
/// @brief  BS mover to get around a stupid "mover" that was embedded in the parser
/// @author Brian Weitzner brian.weitzner@gmail.com, Steven Lewis smlewi@gmail.com
/// @date   Rebased to next year.
///

#ifndef INCLUDED_protocols_simple_moves_SequenceProfileMover_HH
#define INCLUDED_protocols_simple_moves_SequenceProfileMover_HH

// Unit Headers
#include <protocols/simple_moves/SequenceProfileMover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <string>

namespace protocols {
namespace simple_moves {

class SequenceProfileMover : public moves::Mover {
public:
	~SequenceProfileMover();
	SequenceProfileMover();
	virtual void apply( core::pose::Pose& pose );
	virtual std::string get_name() const;

	// function for the parser with lots of accessors
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	void set_cst_file_name( std::string const& cst_file_name ) { cst_file_name_ = cst_file_name; }
	void set_profile_wgt( core::Real const profile_wgt ) { profile_wgt_ = profile_wgt; }

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	core::Real get_profile_wgt() const { return profile_wgt_; }
	std::string const& get_cst_file_name() const { return cst_file_name_; }

private:
	core::Real profile_wgt_;
	std::string cst_file_name_;

};


} //simple_moves
} //protocols
#endif //protocols_simple_moves_SequenceProfileMover_HH
