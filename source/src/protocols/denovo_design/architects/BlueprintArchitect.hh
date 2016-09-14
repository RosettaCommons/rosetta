// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/BlueprintArchitect.hh
/// @brief Designs a structure using a Blueprint file
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_architects_BlueprintArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_BlueprintArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/BlueprintArchitect.fwd.hh>
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace denovo_design {
namespace architects {

///@brief Designs a structure using a Blueprint file
class BlueprintArchitect : public DeNovoArchitect {
public:
	typedef components::StructureData StructureData;

public:
	BlueprintArchitect( std::string const & id_value );

	virtual ~BlueprintArchitect();

	static std::string
	class_name() { return "BlueprintArchitect"; }

	virtual std::string
	type() const;

	DeNovoArchitectOP
	clone() const;

	virtual StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

public:
	void
	set_blueprint( protocols::jd2::parser::BluePrint const & bp );

private:
	/// @brief add templated segments from pose. These are designated by a residue number
	/// in the blueprint instead of 0
	/// @param[in,out] sd    StructureData to be modified
	/// @param[in]     pose  Pose containing template residues
	void
	set_template_segments( StructureData & sd, core::pose::Pose const & pose ) const;

	/// @brief gets names of templated segments from pose. These are designated by a residue number
	/// in the blueprint instead of 0.  All residues in the segment must be sequential and not built
	/// denovo
	/// @param[in] sd    StructureData to be scanned
	/// @returns SegmentNameSet of segments with templates
	SegmentNames
	get_template_segments( StructureData const & sd ) const;

	/// @brief Adds helix pairings to the given SD using the HHPAIR line of the blueprint
	/// @param[in,out] sd StructureData object to be modified
	void
	set_helix_pairings( StructureData & sd ) const;

	/// @brief Converts HHPAIR data in the blueprint into usable information
	/// @param[in]  hhpair_str  Blueprint HHPAIR string
	/// @param[out] h1          Helix number (N-->C ordering) for first paired helix
	/// @param[out] h2          Helix number (N-->C ordering) for first paired helix
	/// @param[out] orientation Orientation of the helices (P or A)
	void
	get_helix_pairings(
		std::string const & hhpair_str,
		core::Size & h1,
		core::Size & h2,
		char & orientation ) const;

private:
	protocols::jd2::parser::BluePrintCOP blueprint_;

};

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_BlueprintArchitect_hh

