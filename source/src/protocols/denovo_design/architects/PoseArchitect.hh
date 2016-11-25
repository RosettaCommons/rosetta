// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/PoseArchitect.hh
/// @brief Design segments based on a pose
/// @author Tom Linsky (tlinsky@uw.edu)
#ifndef INCLUDED_protocols_denovo_design_architects_PoseArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_PoseArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/PoseArchitect.fwd.hh>
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace denovo_design {
namespace architects {

///@brief Design segments based on a pose
class PoseArchitect : public DeNovoArchitect {
public:
	PoseArchitect( std::string const & id_value );

	virtual ~PoseArchitect();

	virtual DeNovoArchitectOP
	clone() const;

	static std::string
	architect_name() { return "PoseArchitect"; }

	virtual std::string
	type() const;

	virtual components::StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
protected:
	/// @brief Configuration by XML
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

private:
	/// @brief If set, this secondary structure will be used to override what is in the input pose
	std::string secstruct_;

	/// @brief if true, a one-residue "pad" will be added to the beginning and end of each chain
	///        if false, the pad will not be added
	bool add_padding_;

};

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_PoseArchitect_hh

