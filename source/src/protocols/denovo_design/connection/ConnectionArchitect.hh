// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/connection/ConnectionArchitect.hh
/// @brief Architect for covalently joining two segments of a pose
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_connection_ConnectionArchitect_hh
#define INCLUDED_protocols_denovo_design_connection_ConnectionArchitect_hh

// Unit headers
#include <protocols/denovo_design/connection/ConnectionArchitect.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/architects/StructureArchitect.hh>
#include <protocols/denovo_design/components/IdealAbegoGenerator.fwd.hh>
#include <protocols/denovo_design/components/Segment.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/excn/EXCN_Base.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace connection {

typedef components::Segment Motif;
typedef components::SegmentOP MotifOP;
typedef components::SegmentCOP MotifCOP;
typedef utility::vector1< MotifOP > MotifOPs;
typedef utility::vector1< MotifCOP > MotifCOPs;

///@brief Architect for covalently joining two segments of a pose
class ConnectionArchitect : public architects::StructureArchitect {
public:
	typedef std::pair< std::string, std::string > SegmentPair;
	typedef utility::vector1< SegmentPair > SegmentPairs;
	typedef std::set< core::Size > LengthSet;
	typedef architects::Lengths Lengths;

public:
	//constants
	static std::string const
	DATA_MAP_NAME;

public:
	ConnectionArchitect( std::string const & id_value );

	virtual ~ConnectionArchitect();

	static std::string
	class_name() { return "Connection"; }

	virtual std::string
	type() const;

	virtual ConnectionArchitectOP
	clone() const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

public:
	/// @brief does the work of modifying the StructureData
	void
	apply( components::StructureData & sd ) const;

	/// @brief applies with a specific random number
	void
	apply( components::StructureData & sd, core::Real & random ) const;

	/// @brief computes a list of possible motifs
	MotifOPs
	compute_connection_candidates( components::StructureData const & sd ) const;

	/// @brief returns list of allowed segment 1 ids
	SegmentNames const &
	segment1_ids() const;

	/// @brief sets list of segment1 ids from string
	void
	set_segment1_ids( std::string const & segment1_str );

	/// @brief sets list of segment1 ids from list
	void
	set_segment1_ids( SegmentNames const & segments );

	/// @brief returns list of allowed segment2 ids
	SegmentNames const &
	segment2_ids() const;

	/// @brief sets list of segment2 ids from string
	void
	set_segment2_ids( std::string const & segment2_str );

	/// @brief sets list of segment1 ids from list
	void
	set_segment2_ids( SegmentNames const & segments );

	/// @brief gets user-specified chain number for the lower chain to be connected. 0 if not specified
	core::Size
	user_chain1() const;

	/// @brief sets user-specified chain number for the lower chain to be connected. 0 if unspecified
	void
	set_user_chain1( core::Size const chain );

	/// @brief gets user-specified chain number for the upper chain to be connected. 0 if not specified
	core::Size
	user_chain2() const;

	/// @brief sets user-specified chain number for the upper chain to be connected. 0 if unspecified
	void
	set_user_chain2( core::Size const chain );

	/// @brief sets motifs using a motif string and a string of cutpoints
	void
	set_motifs( std::string const & motif_str, std::string const & cutpoints_str );

	/// @brief sets motifs via a vector
	void
	set_motifs( MotifCOPs const & motifs );

	/// @brief sets whether to use "ideal abego" loops according to Koga papers
	void
	set_ideal_abego( bool const use_ideal_abego, bool const extend_ss );

	/// @brief sets whether to always try to bridge.  If true, a random cutpoint will be selected in the connection
	///        if the chains to be connected have different movable groups
	void
	set_bridge( bool const bridge );

private:
	/// @brief Get list of segment pairs allowed to be connected
	MotifOPs
	compute_connection_candidates(
		components::StructureData const & sd,
		AreConnectablePredicate const & connectable ) const;

	SegmentPairs
	segment_pairs( components::StructureData const & sd ) const;

	SegmentPairs
	combine_segment_names(
		SegmentNames const & seg1s,
		SegmentNames const & seg2s ) const;

	MotifOP
	choose_motif( MotifOPs const & motifs, core::Real & random ) const;

	/// @brief creates StructureData from given Motif
	void
	connect( components::StructureData & sd, Motif & motif ) const;

	LengthSet
	lengths() const;

	/// @brief returns a set of valid loop index cutpoints
	LengthSet
	cutpoints() const;

	MotifOPs
	motifs_for_pair(
		SegmentPair const & pair,
		components::StructureData const & sd,
		LengthSet const & length_set,
		LengthSet const & cutpoint_set	) const;

	MotifCOPs
	parse_motif_string( std::string const & motifs ) const;

	SegmentNames
	available_upper_termini( components::StructureData const & sd ) const;

	SegmentNames
	available_lower_termini( components::StructureData const & sd ) const;

private:
	bool bridge_;
	components::IdealAbegoGeneratorCOP ideal_abego_;
	MotifCOPs motifs_;
	SegmentNames segment1_ids_;
	SegmentNames segment2_ids_;
	core::Size chain1_;
	core::Size chain2_;
};

class AreConnectablePredicate {
public:
	typedef std::set< core::Size > MovableGroupSet;

public:
	AreConnectablePredicate( bool const allow_cyclic );

	virtual bool
	operator()(
		components::StructureData const & sd,
		Motif const & motif ) const;

private:
	bool
	check_distance( components::StructureData const & sd, Motif const & motif ) const;

	bool
	check_movable_groups( components::StructureData const & sd, Motif const & motif ) const;

	MovableGroupSet
	connected_movable_groups( components::StructureData const & sd, std::string const & seg_name ) const;

private:
	bool allow_cyclic_;
	AreConnectablePredicate();
};
typedef utility::pointer::shared_ptr< AreConnectablePredicate > AreConnectablePredicateOP;
typedef utility::pointer::shared_ptr< AreConnectablePredicate const > AreConnectablePredicateCOP;

SegmentNames
parse_segment_names( std::string const & segment_name_str );

class EXCN_ConnectionSetupFailed : public utility::excn::EXCN_Base {
public:
	EXCN_ConnectionSetupFailed( std::string const & msg ):
		utility::excn::EXCN_Base(),
		message_( msg ) {}
	std::string const & message() const { return message_; }
	virtual void show( std::ostream & os ) const { os << message_; }
private:
	std::string const message_;
};


} //connection
} //denovo_design
} //protocols

#endif //INCLUDED_protocols_denovo_design_connection_ConnectionArchitect_hh

