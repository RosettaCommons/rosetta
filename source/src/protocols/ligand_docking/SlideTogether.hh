// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_SlideTogether_hh
#define INCLUDED_protocols_ligand_docking_SlideTogether_hh

#include <protocols/moves/Mover.hh>

//// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>


// Unit Headers

//STL headers
#include <string>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {


/// @brief
class SlideTogether : public protocols::moves::Mover
{
public:
	SlideTogether();
	SlideTogether(std::string chain);

	~SlideTogether() override;
	SlideTogether(SlideTogether const & that);

	void apply( core::pose::Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
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
	std::string chain_;
	utility::vector1<core::Size> jumps_; // these guys tag along, such as waters and metals
};

} //namespace ligand_docking
} //namespace protocols

#endif
