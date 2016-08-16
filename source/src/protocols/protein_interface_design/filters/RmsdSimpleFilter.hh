// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/RmsdSimpleFilter.hh
/// @brief Simply calculate RMSD, doesn't assume that poses are statick unmutable objects
/// @author Daniel Silva (dadriano@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_RmsdSimpleFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_RmsdSimpleFilter_hh

//Include Rosetta Core Stuff
#include <core/pose/Pose.hh>
#include <core/types.hh>

//Include Rosetta protocols
#include <protocols/filters/Filter.hh>

//Include ObjexxFCL
#include <ObjexxFCL/FArray.all.hh>

//Include Rosetta utilities
#include <utility/vector1.hh>

//Include Rosetta XML tag reader
#include <utility/tag/Tag.fwd.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

class RmsdSimpleFilter : public protocols::filters::Filter
{
public:
	RmsdSimpleFilter();

	RmsdSimpleFilter(
		core::Real const threshold,
		core::pose::PoseOP reference_pose
	);

	bool apply( core::pose::Pose const & pose ) const;

	protocols::filters::FilterOP clone() const;
	protocols::filters::FilterOP fresh_instance() const;

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~RmsdSimpleFilter();

	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & reference_pose
	);

private:

	void
	superposition_transform(
		core::Size natoms,
		ObjexxFCL::FArray1_double const& weights,
		ObjexxFCL::FArray2_double& ref_coords,
		ObjexxFCL::FArray2_double& coords,
		numeric::xyzMatrix< core::Real > &RotM,
		numeric::xyzVector< core::Real > &TvecA,
		numeric::xyzVector< core::Real > &TvecB
	) const;

	core::Real
	rmsd_bb(
		core::pose::Pose const & poseA,
		utility::vector1< core::Size > const & positions_to_alignA,
		core::pose::Pose const & poseB,
		utility::vector1< core::Size > const & positions_to_alignB
	) const;

	core::Real
	dist_bb(
		core::pose::Pose const & poseA,
		utility::vector1< core::Size > const & positions_to_alignA,
		core::pose::Pose const & poseB,
		utility::vector1< core::Size > const & positions_to_alignB
	) const;

private:

	core::Real threshold_;
	core::pose::PoseOP reference_pose_;
	core::Size target_chain_;
	core::Size do_align_;
	bool b_target_chain_;
};

} // filters
} // protein_interface_design
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_filters_RmsdSimpleFilter_HH_
