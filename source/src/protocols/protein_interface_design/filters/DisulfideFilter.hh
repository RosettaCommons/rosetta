// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/protein_interface_design/filters/DisulfideFilter.hh
/// @brief Filters for interfaces which could form a disulfide bond between
/// docking partners.
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date Created 4/30/2009

#ifndef INCLUDED_protocols_protein_interface_design_filters_DisulfideFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_DisulfideFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers

#include <core/scoring/disulfides/CentroidDisulfidePotential.fwd.hh>
#include <utility/vector1.hh>


// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace filters {

/**
 * @brief Filters for structures which could form a disulfide bond across the
 * docking interface.
 *
 * @details Use this filter when you are trying to design one docking member
 * so that it forms a disulfide bond to one or more target residues of the other
 * docking partner. The filter does not consider the indentities of the residues
 * involved, only their Cb position.
 *
 * This filter only applies to centroid poses. Calling it with a full atom pose
 * will result in everything failing.
 *
 * @author Spencer Bliven <blivens@u.washington.edu>
 */
class DisulfideFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	DisulfideFilter();
	/// @brief copy ctor
	DisulfideFilter(DisulfideFilter const& df);
	///@brief Constructor with a single target residue
	DisulfideFilter( core::Size targetResidue );
	///@brief Constructor with multiple target residues
	/// @details targets may come from either binding partner. If no targets
	///   are specified for one target, all residues on the interface will be
	///   concidered.
	DisulfideFilter( utility::vector1<core::Size> const& targetResidues );
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const {
		return new DisulfideFilter( *this );
	}
	virtual protocols::filters::FilterOP fresh_instance() const{
		return new DisulfideFilter();
	}

	virtual ~DisulfideFilter();
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
private:
	/// @brief a list of residues which may participate in the disulfide.
	/// @details If either docking partner has no target residues specified,
	///   all interface residues will be allowed to disulfide bond.
	utility::vector1<core::Size> targets_;
	/// @brief for calculating centroid disulfide energies
	static const core::scoring::disulfides::CentroidDisulfidePotential potential_;
	/// @brief Which jump defines the interface where the targets lie?
	Size rb_jump_;
};

} // filters
} // protein_interface_design
} // devel

#endif //INCLUDED_protocols_protein_interface_design_filters_DisulfideFilter_HH_
