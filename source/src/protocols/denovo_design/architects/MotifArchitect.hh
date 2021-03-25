// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/denovo_design/architects/MotifArchitect.hh
/// @brief Designs de novo motifs
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_architects_MotifArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_MotifArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/MotifArchitect.fwd.hh>
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// Boost headers
#include <boost/lexical_cast.hpp>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace architects {

/// @brief for planning arbitrary motifs
class MotifArchitect : public DeNovoArchitect {
public:
	typedef components::Segment Motif;
	typedef components::SegmentOP MotifOP;
	typedef components::SegmentCOP MotifCOP;
	typedef utility::vector1< MotifCOP > MotifCOPs;

public:
	MotifArchitect( std::string const & id );

	~MotifArchitect() override;

	DeNovoArchitectOP
	clone() const override;

	static std::string
	architect_name() { return "DeNovoMotif"; }

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	std::string
	type() const override;

	components::StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const override;

protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

public:
	/// @brief Returns iterator to the start of the list of possible motifs that could be used during design
	MotifCOPs::const_iterator
	motifs_begin() const;

	/// @brief Returns end iterator to the start of the list of possible motifs that could be used during design
	MotifCOPs::const_iterator
	motifs_end() const;

	/// @brief Returns the list of possible motifs that could be used during design
	MotifCOPs const &
	motifs() const;

	/// @brief Sets the list of allowed motifs using a string
	void
	set_motifs( std::string const & motif_str );

	/// @brief Sets the list of allowed motifs from a list of Motif objects
	void
	set_motifs( MotifCOPs const & motifs );

private:
	MotifCOPs motifs_;
};

/// @brief Creates a "len"-residue long SecStructInfo object with SS and ABEGO as specified here
SecStructInfo
generate_secstruct_for_length(
	char const ss_char,
	std::string const & abego,
	core::Size const len );

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_MotifArchitect_hh
