// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/filters/CoreResiduesPerElementFilter.hh
/// @brief Tom's Denovo Protocol. This is freely mutable and used for playing around with stuff
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_devel_denovo_design_filters_CoreResiduesPerElementFilter_hh
#define INCLUDED_devel_denovo_design_filters_CoreResiduesPerElementFilter_hh

// Unit headers
#include <devel/denovo_design/filters/CoreResiduesPerElementFilter.fwd.hh>

// Project headers
#include <devel/denovo_design/scoring/SideChainNeighborsEnergy.fwd.hh>

// Protocol headers
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// C++ headers

namespace devel {
namespace denovo_design {
namespace filters {

class CoreResiduesPerElementFilter : public protocols::filters::Filter {
public:

	/// @brief Initialize CoreResiduesPerElementFilter
	CoreResiduesPerElementFilter();

	/// @brief virtual constructor to allow derivation
	virtual ~CoreResiduesPerElementFilter();

	/// @brief Parses the CoreResiduesPerElementFilter tags
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	/// @brief Return the name of this mover.
	virtual std::string get_name() const;

	/// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::filters::FilterOP clone() const;

	/// @brief Apply the CoreResiduesPerElementFilter. Overloaded apply function from filter base class.
	virtual protocols::filters::FilterOP fresh_instance() const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual bool apply( core::pose::Pose const & pose ) const;

	core::Real compute( core::pose::Pose const & pose ) const;

	// mutators
public:
	/// @brief sets the residue selector used to decide which residues to evaluate
	void set_selector( core::select::residue_selector::ResidueSelectorCOP rs );
	/// @brief sets the cutoff used to define a residue as core
	void set_core_cutoff( core::Real const core_cutoff );

	// protected functions
protected:
	core::Size evaluate_element(
		core::pose::Pose const & pose,
		core::scoring::methods::SideChainNeighborsEnergy const & scn,
		core::Size const start,
		core::Size const stop ) const;

	core::Size evaluate_element(
		core::pose::Pose const & pose,
		core::scoring::methods::SideChainNeighborsEnergy const & scn,
		core::select::residue_selector::ResidueSubset const & subset ) const;

	core::Size evaluate_element(
		core::pose::Pose const & pose,
		core::scoring::methods::SideChainNeighborsEnergy const & scn,
		utility::vector1< core::Size > const & residues ) const;

private:   // options

private:   // other data
	core::Real core_cutoff_;
	/// @brief residue selector to identify positions to rebuild
	core::select::residue_selector::ResidueSelectorCOP selector_;
};


} // filters
} // denovo_design
} // devel

#endif
