// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/matdes/InterfacePackingFilter.hh
/// @brief Calculates the RosettaHoles score for the atoms within a user-defined cutoff distance
/// (default = 9.0) of the subpose interface.  Filters based on lower and upper thresholds,
/// which are set to -/+5 by default.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_matdes_InterfacePackingFilter_hh
#define INCLUDED_protocols_matdes_InterfacePackingFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/matdes/InterfacePackingFilter.fwd.hh>

#include <utility/vector1.hh>

// Unit headers

namespace protocols {
namespace matdes {

class InterfacePackingFilter : public protocols::filters::Filter

{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	InterfacePackingFilter();
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~InterfacePackingFilter();
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	core::Real distance_cutoff() const;
	core::Real contact_dist() const;
	core::Real lower_threshold() const;
	core::Real upper_threshold() const;
	bool multicomp() const;
	std::string sym_dof_names() const;

	void distance_cutoff( core::Real const d );
	void contact_dist( core::Real const c );
	void lower_threshold( core::Real const l );
	void upper_threshold( core::Real const u );
	void sym_dof_names( std::string const s );
	void multicomp( bool const multicomp );

private:

	core::Real distance_cutoff_, contact_dist_, lower_threshold_, upper_threshold_; // distance within which atoms must be across interface in order to be scored. Lower and upper thresholds for the RosettaHoles score.

	bool multicomp_;
	std::string sym_dof_names_;
};

} // matdes
} // protocols

#endif

