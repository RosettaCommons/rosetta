// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#ifndef INCLUDED_protocols_simple_moves_MembraneTopology_hh
#define INCLUDED_protocols_simple_moves_MembraneTopology_hh

#include <protocols/simple_moves/MembraneTopology.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/AA.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

#include <utility/vector1.hh>


// Utility headers

// C++ headers

// Unit headers

// @brief Simple wrapper to MembraneTopology which lives in core/scoring

namespace protocols {
namespace simple_moves {

/// @brief A mover to change one torsion angle
class MembraneTopology : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	///@brief default ctor
	MembraneTopology();
	virtual ~MembraneTopology();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return (protocols::moves::MoverOP( new protocols::simple_moves::MembraneTopology( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new MembraneTopology );
	}

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	std::string span_file() const { return span_file_; }
  void span_file( std::string const s ){ span_file_ = s; }

private:
	std::string span_file_;
};

} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_MembraneTopology_HH_
