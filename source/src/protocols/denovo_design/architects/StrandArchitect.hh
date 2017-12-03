// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/StrandArchitect.hh
/// @brief Architect that creates a beta strand
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_architects_StrandArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_StrandArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/StrandArchitect.fwd.hh>
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>
#include <protocols/denovo_design/components/SegmentPairing.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace architects {

typedef utility::vector1< core::Size > StrandExtensions;

typedef components::SegmentResid SegmentResid;
typedef utility::vector1< SegmentResid > StrandBulges;
typedef utility::vector1< StrandBulges > AllowedStrandBulges;
typedef utility::vector1< SegmentResid > StrandExtended;
typedef utility::vector1< StrandExtended > AllowedStrandExtended;

///@brief Architect that creates a beta strand
class StrandArchitect : public protocols::denovo_design::architects::DeNovoArchitect {
public:
	typedef std::set< core::Size > LengthSet;
	typedef protocols::denovo_design::architects::DeNovoArchitect DeNovoArchitect;
	typedef protocols::denovo_design::architects::DeNovoArchitectOP DeNovoArchitectOP;
	typedef protocols::denovo_design::components::StructureData StructureData;
	typedef protocols::denovo_design::components::StructureDataOP StructureDataOP;
	typedef components::RegisterShift RegisterShift;
	typedef components::RegisterShifts RegisterShifts;
	typedef components::StrandOrientation StrandOrientation;
	typedef components::StrandOrientations StrandOrientations;

public:
	StrandArchitect( std::string const & id_value );

	virtual ~StrandArchitect();

	static std::string
	class_name() { return "StrandArchitect"; }

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
	components::StructureDataCOPs::const_iterator
	motifs_begin() const;

	components::StructureDataCOPs::const_iterator
	motifs_end() const;

	void
	set_length( std::string const & length_str );

	void
	set_length( Lengths const & lengths_val );

	void
	set_bulges( std::string const & bulges_str );

	void
	set_bulges( AllowedStrandBulges const & bulges );

	void
	set_extended( std::string const & extended_str );

	void
	set_extended( AllowedStrandExtended const & extended );

	void
	enumerate_permutations();

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public:
	static std::string const
	bulge_keyname();

	static StrandOrientation
	int_to_orientation( int const integer );

	static StrandBulges
	retrieve_bulges_s( StructureData const & sd, std::string const & segment_id ); // Fix for PyRosetta: renaming static version to retrieve_bulges_s to avoid name clashed with member function

public:
	StrandBulges
	retrieve_bulges( StructureData const & sd ) const;

private:
	void
	store_bulges( StructureData & sd, StrandBulges const & bulge ) const;

	/// @brief Given a list of bulges, build a motif
	/// @param[in] bulges    List of bulge positions to place onto the strand
	/// @param[in] secstruct Secondary structure for the new segment
	/// @param[in] abego     Abego for the new segment, without bulges placed
	StructureDataOP
	create_motif(
		StrandBulges const & bulges,
		StrandExtended const & extended,
		std::string const & secstruct,
		std::string const & abego ) const;

private:
	void
	needs_update();

	/// @brief If architect is updated, simply exit without error.  If not, exit with an error message
	void
	check_updated() const;

	components::StructureDataCOPs
	compute_permutations() const;

private:
	components::StructureDataCOPs motifs_;
	Lengths lengths_;
	AllowedStrandBulges bulges_;
	AllowedStrandExtended extended_;
	bool updated_;

};

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_StrandArchitect_hh
