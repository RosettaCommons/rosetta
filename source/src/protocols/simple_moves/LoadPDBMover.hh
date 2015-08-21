// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoadPDBMover.hh
/// @brief switch the chain order

#ifndef INCLUDED_protocols_simple_moves_LoadPDBMover_hh
#define INCLUDED_protocols_simple_moves_LoadPDBMover_hh

#include <protocols/simple_moves/LoadPDBMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMapObj.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class LoadPDBMover : public moves::Mover {
public:
	LoadPDBMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	std::string filename() const;
	void filename( std::string const f );

	bool append() const { return append_; };
	void append( bool setting ){ append_ = setting; }
private:
	std::string filename_;
	bool append_;
};


} // simple_moves
} // protocols

#endif
