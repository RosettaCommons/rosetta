// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/NamedMover.hh
/// @brief NamedMover functions for building structures from components
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_NamedMover_hh
#define INCLUDED_protocols_denovo_design_components_NamedMover_hh

// Unit headers
#include <protocols/denovo_design/components/NamedMover.fwd.hh>

// Project headers
#include <protocols/denovo_design/types.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>

// Core headers

// Basic/Numeric/Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers

namespace protocols {
namespace denovo_design {
namespace components {

/// @brief manages information about segments of residues
class NamedMover : public protocols::moves::Mover {
public:
	NamedMover();
	NamedMover( std::string const & id, std::string const & parent_id );

	virtual ~NamedMover() {}

	/// @brief setup the parameters via an xml tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	/// @brief returns an identifier for this mover
	std::string const & id() const;

	/// @brief returns the name of the parent mover
	std::string const & parent_id() const;

	/// @brief sets an identifier for this mover
	void set_id( std::string const & idval );

	/// @brief sets an identifier for the parent of this mover, for nested moves
	void set_parent_id( std::string const & parent );

protected:
	/// @brief adds prefix if necessary, returns result
	std::string add_parent_prefix( std::string const & s ) const;

private:
	std::string id_;
	std::string parent_id_;
};

} // components
} // denovo_design
} // protocols

#endif
