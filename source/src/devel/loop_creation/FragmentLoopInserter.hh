// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FragmentLoopInserter.hh
///
/// @brief LoopInserter that identifies fragments with ends matching RMSD
/// @author Tim Jacobs

#ifndef INCLUDED_devel_loop_creation_FragmentLoopInserter_HH
#define INCLUDED_devel_loop_creation_FragmentLoopInserter_HH

//Unit
#include <devel/loop_creation/LoopInserter.hh>

//Core
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>

//Numeric
#include <numeric/xyzVector.fwd.hh>

//Protocols
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loops/Loop.fwd.hh>

namespace devel {
namespace loop_creation {

class FragmentLoopInserter : public LoopInserter
{
public:

	FragmentLoopInserter();

	protocols::moves::MoverOP
	clone() const override;

	protocols::moves::MoverOP
	fresh_instance() const override;


	void
	apply(
		core::pose::Pose & pose
	) override;

	void
	find_loop_fragments(
		core::pose::Pose & pose
	);

	void
	build_fragment_loop(
		core::pose::Pose & pose,
		core::fragment::FragDataCOP fragment
	);

	utility::vector1< numeric::xyzVector<core::Real> >
	get_pose_coords_to_match(
		core::pose::Pose const & pose
	);

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	utility::vector1<core::fragment::FragSetOP> frag_sets_;

	core::Real max_rms_;

	utility::vector1<core::Size> loop_sizes_;
	bool pack_;
	bool design_;

	//a map from loop anchor size to fragments that pass RMS requirements. loop anchor size
	//is needed so that if the loop anchor is changed, invalid fragments aren't returned.
	std::map<core::Size, utility::vector1<core::fragment::FragDataCOP> > anchor_frags_;

	//number of flanking residues for each segment that must match the
	//torsions of the input pose (within loophash min and max rms)
	core::Size num_flanking_residues_to_match_;

	bool modify_flanking_regions_;
};

} //loop creation
} //devel

#endif
