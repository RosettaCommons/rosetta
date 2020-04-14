// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/VLB.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_VLB_hh
#define INCLUDED_protocols_protein_interface_design_movers_VLB_hh

//#include <protocols/protein_interface_design/movers/DesignRepackMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/forge/build/BuildManager.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/ResidueIndexDescription.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

enum class VLBInstructionType {
	Bridge,
	ConnectRight,
	GrowLeft,
	GrowRight,
	SegmentInsert,
	SegmentRebuild,
	SegmentSwap
};

struct VLBInstruction {
	VLBInstructionType type;

	// Not all parameters are used by all types.
	core::pose::ResidueIndexDescriptionCOP res1;
	core::pose::ResidueIndexDescriptionCOP res2;

	std::string aa;
	std::string ss;
	std::string filename;
	std::string side;

	bool keep_bb = false;
};

/// @brief user interface for YAB's Variable Length Build.
class VLB : public protocols::moves::Mover
{
public:
	VLB(); // default ctor//design_ = true;

	void apply( core::pose::Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	protocols::forge::build::BuildManagerOP make_manager(core::pose::Pose const & pose) const;

private:
	utility::vector1< VLBInstruction > instruction_list_;

	core::scoring::ScoreFunctionOP scorefxn_;
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_VLB_HH*/
