// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/SegmentPairing.hh
/// @brief Handles user-specified pairing between/among segments
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_denovo_design_components_SegmentPairing_hh
#define INCLUDED_protocols_denovo_design_components_SegmentPairing_hh

// Unit headers
#include <protocols/denovo_design/components/SegmentPairing.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace denovo_design {
namespace components {

class SegmentPairing : public utility::pointer::ReferenceCount {
public:
	// types
	enum PairingType {
		HELIX = 1,
		STRAND = 2,
		HELIX_SHEET = 3,
		UNKNOWN
	};

public:
	// constants
	static std::string
		TAG_NAME;

public:
	SegmentPairing();

	SegmentPairing( SegmentNames const & paired_segments );

	virtual
	~SegmentPairing() {};

	virtual SegmentPairingOP
	clone() const = 0;

	virtual PairingType
	type() const = 0;

	static std::string
	complex_type_name_for_pairing( std::string const & pairing_name );

protected:
	static void
	add_common_xml_elements(
		utility::tag::XMLSchemaDefinition & xsd,
		std::string const & class_name,
		std::string const & description,
		utility::tag::AttributeList & attlist );

	virtual void
	parse_tag( utility::tag::Tag const & tag ) = 0;

	virtual void
	to_xml( utility::tag::Tag & tag ) const = 0;

public:
	virtual std::string
	pairing_string( StructureData const & sd ) const = 0;

	void
	parse_my_tag( utility::tag::Tag const & tag );

	friend std::ostream &
	operator<<( std::ostream & os, SegmentPairing const & pairing );

	bool
	has_segment( std::string const & segment ) const;

	SegmentNames const &
	segments() const;

	void
	set_segments( std::string const & segments_str );

	void
	set_segments( SegmentNames const & segments );

public:
	// Static Functions

	static SegmentPairingOP
	create( std::string const & type_name );

	/// @brief Gets string for all strand pairings from a StructureData
	static std::string
	get_strand_pairings( StructureData const & sd );

	/// @brief Gets string for all helix pairings from a StructureData
	static std::string
	get_helix_pairings( StructureData const & sd );

	/// @brief Gets string for all helix-strand-strand triplets in a StructureData
	static std::string
	get_hss_triplets( StructureData const & sd );

	/// @brief Gets pairs of residues involved in strand pairing
	static ResiduePairs
	get_strand_residue_pairs( StructureData const & sd );

private:
	// static functions

	static std::string
	type_to_str( PairingType const & type );

	static SegmentPairingCOPs
	get_pairings( StructureData const & sd, PairingType const & type );

	static std::string
	get_pairing_str( StructureData const & sd, PairingType const & type );

private:
	SegmentNames segments_;

};

class HelixPairing : public SegmentPairing {
public:
	HelixPairing();

	HelixPairing( SegmentName const & h1, SegmentName const & h2, bool const is_parallel );

	virtual
	~HelixPairing() {};

	virtual SegmentPairingOP
	clone() const;

	static std::string
	class_name() { return "HelixPairing"; }

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	virtual PairingType
	type() const { return HELIX; }

protected:
	virtual void
	parse_tag( utility::tag::Tag const & tag );

	virtual void
	to_xml( utility::tag::Tag & tag ) const;

public:
	virtual std::string
	pairing_string( StructureData const & sd ) const;

private:
	bool parallel_;
};

class StrandPairing : public SegmentPairing {
public:
	StrandPairing();

	StrandPairing(
		SegmentName const & s1,
		SegmentName const & s2,
		StrandOrientation const & orient1,
		StrandOrientation const & orient2,
		RegisterShift const & shift );

	virtual
	~StrandPairing() {};

	virtual SegmentPairingOP
	clone() const;

	static std::string
	class_name() { return "StrandPairing"; }

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	virtual PairingType
	type() const { return STRAND; }

protected:
	virtual void
	parse_tag( utility::tag::Tag const & tag );

	virtual void
	to_xml( utility::tag::Tag & tag ) const;

public:
	virtual std::string
	pairing_string( StructureData const & sd ) const;

	bool
	parallel() const;

	StrandOrientation
	orient1() const;

	StrandOrientation
	orient2() const;

	RegisterShift
	shift() const;

private:
	RegisterShift
	nobu_register_shift(
		StructureData const & sd,
		protocols::fldsgn::topology::Strand const & s1,
		protocols::fldsgn::topology::Strand const & s2,
		RegisterShift const nc_order_shift,
		StrandOrientation const & nc_order_orient ) const;

private:
	StrandOrientation orient1_;
	StrandOrientation orient2_;
	RegisterShift shift_;
};

class HelixSheetPairing : public SegmentPairing {
public:
	HelixSheetPairing();

	HelixSheetPairing(
		SegmentName const & helix,
		SegmentName const & s1,
		SegmentName const & s2 );

	virtual
	~HelixSheetPairing() {};

	virtual SegmentPairingOP
	clone() const;

	static std::string
	class_name() { return "HelixSheetPairing"; }

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	virtual PairingType
	type() const { return HELIX_SHEET; }

protected:
	virtual void
	parse_tag( utility::tag::Tag const & tag );

	virtual void
	to_xml( utility::tag::Tag & tag ) const;

public:
	virtual std::string
	pairing_string( StructureData const & sd ) const;

};


} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_SegmentPairing_hh
