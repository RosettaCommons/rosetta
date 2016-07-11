// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/SegmentPairing.hh
/// @brief Handles user-specified pairing between/among segments
/// @author Tom Linsky (tlinsky@uw.edu)


#ifndef INCLUDED_protocols_denovo_design_components_SegmentPairing_hh
#define INCLUDED_protocols_denovo_design_components_SegmentPairing_hh

// Unit headers
#include <protocols/denovo_design/components/SegmentPairing.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace denovo_design {
namespace components {

enum PairingType {
	HELIX = 1,
	STRAND = 2,
	UNKNOWN = 3
};

SegmentPairingOP
create_segment_pairing( std::string const & type_name );

class SegmentPairing : public utility::pointer::ReferenceCount {
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

protected:
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

	SegmentNames const &
	segments() const;

	void
	set_segments( std::string const & segments_str );

	void
	set_segments( SegmentNames const & segments );

private:
	// static functions
	static std::string
	type_to_str( PairingType const & type );

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
		architects::StrandOrientation const & orient1,
		architects::StrandOrientation const & orient2,
		architects::RegisterShift const & shift );

	virtual
	~StrandPairing() {};

	virtual SegmentPairingOP
	clone() const;

	static std::string
	class_name() { return "StrandPairing"; }

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

	//architects::RegisterShift
	//nobu_register_shift() const;

private:
	architects::StrandOrientation orient1_;
	architects::StrandOrientation orient2_;
	architects::RegisterShift shift_;

};

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_SegmentPairing_hh
