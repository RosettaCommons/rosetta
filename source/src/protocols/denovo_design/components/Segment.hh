// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/components/Segment.hh
/// @brief Segment functions for building structures from components
/// @details
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_Segment_hh
#define INCLUDED_protocols_denovo_design_components_Segment_hh

// Unit headers
#include <protocols/denovo_design/components/Segment.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/types.hh>

// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Basic/Numeric/Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

// C++ Headers

namespace protocols {
namespace denovo_design {
namespace components {

class ResidueDihedrals : public utility::pointer::ReferenceCount {
public:
	ResidueDihedrals();

	ResidueDihedrals( core::pose::Pose const & input, core::Size const lower_resid );

	core::Real
	phi() const;

	core::Real
	psi() const;

	core::Real
	omega() const;

	core::Real
	lower_phi() const;

	core::Real
	upper_psi() const;

	core::Real
	upper_omega() const;

	void
	set_in_pose( core::pose::Pose & pose, core::Size const lower_resid ) const;

private:
	core::Real lower_phi_;
	core::Real psi_;
	core::Real omega_;
	core::Real phi_;
	core::Real upper_psi_;
	core::Real upper_omega_;

#ifdef    SERIALIZATION
public:
	template< class Archive >
	void
	save( Archive & arc ) const;

	template< class Archive >
	void
	load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief manages information about segments of residues
class Segment : public utility::pointer::ReferenceCount {
public:
	Segment();

	Segment(
		std::string const & ss_val,
		std::string const & abego_val,
		bool const start_inc,
		bool const stop_inc );

	virtual ~Segment() {}

	virtual SegmentOP
	clone() const;

	/// @brief construct from xml tag
	void
	parse_tag( utility::tag::TagCOP tag );

	/// @brief construct from motif string (i.e. '1LX-10HA-1LA-1LB')
	void
	parse_motif( std::string const & motif_str );

	inline bool nterm_included() const { return nterm_included_; }
	inline bool cterm_included() const { return cterm_included_; }

public:
	//////////////////////////////////////////
	// pose relationship
	//////////////////////////////////////////
	void
	set_pose_start( core::Size const pose_resid );

	core::Size
	lower() const;

	core::Size
	upper() const;

	core::Size
	start() const;

	core::Size
	stop() const;

	core::Size
	safe() const;

	core::Size
	cutpoint() const;

	bool
	contains( core::Size const pose_resid ) const;

	/// @brief converts a internal segment resid to a pose resid
	///        1 --> segment_start - 1 + 1
	///        2 --> segment_start - 1 + 2
	///        N --> segment_start - 1 + N
	core::Size
	segment_to_pose( SegmentResid const segment_resid ) const;

	/// @brief converts a internal segment resid to a pose resid
	///        segment_start - 1 + 1 --> 1
	///        segment_start - 1 + 2 --> 2
	///        segment_start - 1 + N --> N
	SegmentResid
	pose_to_segment( core::Size const pose_resid ) const;


	//////////////////////////////////////////
	// properties
	//////////////////////////////////////////
	std::string const &
	ss() const;

	void
	set_ss( SegmentResid const segment_resid, char ss_type );

	std::string const &
	abego() const;

	char
	abego( SegmentResid const segment_resid ) const;

	void
	set_abego( std::string const & abego_str );

	void
	set_abego( utility::vector1< std::string > const & abego );

	core::Size
	length() const;

	core::Size
	elem_length() const;

	/// @brief returns template pose
	core::pose::PoseCOP
	template_pose() const;

	/// @brief sets template pose to be the given residues from template_pose
	void
	set_template_pose(
		core::pose::Pose const & template_pose,
		core::Size const start_resid,
		core::Size const stop_resid );

	/// @brief appends segment by the given secstruct and abego
	void
	extend( std::string const & secstruct, std::string const & abego );

	core::Size
	movable_group() const;

	void
	set_movable_group( core::Size const mg );

	/// @brief segment residue number for "safe" residue
	SegmentResid
	safe_segment() const;

	/// @brief sets safe residue for this segment to be the ith residue in the segment
	///        safe = segment_start - 1 + cut_res
	void
	set_safe( SegmentResid const segment_resid );

	/// @brief segment residue number for cutpoint
	SegmentResid
	cutpoint_segment() const;

	/// @brief sets cutpoint for this segment to be the ith residue in the segment
	///        cut = segment_start - 1 + cut_res
	void
	set_cutpoint( SegmentResid const segment_resid );

	/// @brief number of residues before the cutpoint, 0 if cutpoint not set
	core::Size
	n_residues_before_cutpoint() const;

	/// @brief number of residues after the cutpoint, length() if cutpoint not set
	core::Size
	n_residues_after_cutpoint() const;

	//////////////////////////////////////////
	// connectivity
	//////////////////////////////////////////

	std::string const &
	upper_segment() const { return upper_segment_; }

	std::string const &
	lower_segment() const { return lower_segment_; }

	bool
	has_free_lower_terminus() const { return ( lower_segment_ == "" ); }

	bool
	has_free_upper_terminus() const { return ( upper_segment_ == "" ); }

	void
	set_upper_segment( std::string const & comp );

	void
	set_lower_segment( std::string const & comp );

	//////////////////////////////////////////
	// padding
	//////////////////////////////////////////

	ResidueDihedrals const &
	lower_dihedrals() const;

	ResidueDihedrals const &
	upper_dihedrals() const;

	core::conformation::Residue const &
	lower_residue() const;

	core::conformation::Residue const &
	upper_residue() const;

	core::Size
	lower_padding() const;

	core::Size
	upper_padding() const;

	/// @brief adds "padding" residue before the fixed portion of this segment
	void
	add_lower_padding();

	/// @brief adds "padding" residue after the fixed portion of this segment
	void
	add_upper_padding();

	/// @brief deletes dummy residues before the fixed portion of this segment
	void
	delete_lower_padding();

	/// @brief deletes dummy residues after the fixed portion of this segment
	void
	delete_upper_padding();

	/// @brief given a segment residue number, delete that residue. Resid for start_local() == 1
	void
	delete_residue( SegmentResid const segment_resid );

	/// on the chopping block
public:
	// local resids
	core::Size
	lower_local() const;

	core::Size
	upper_local() const;

	core::Size
	start_local() const;

	core::Size
	stop_local() const;

private:
	/// @brief given a residue number range local to this 1=start, length=end, delete the residue
	void
	delete_residues( core::Size const local_resnum_start, core::Size const local_resnum_stop );

	// I/O
public:
	/// output residueinfo
	friend std::ostream &
	operator<<( std::ostream & os, Segment const & res );

	/// @brief converts a segment resid to a "local" resid
	core::Size
	segment_to_local( SegmentResid const segment_resid ) const;

private:
	/// @brief converts a pose resid to a "local" resid
	core::Size
	pose_to_local( core::Size const pose_resid ) const;

	/// @brief converts a local resid to a segment resid
	SegmentResid
	local_to_segment( core::Size const local_resid ) const;

	/// @brief converts a internal segment "local" resid to a pose resid
	///        1 --> posestart - 1 + 1
	///        2 --> posestart - 1 + 2
	///        N --> posestart - 1 + N
	core::Size
	local_to_pose( core::Size const local_resid ) const;

	void
	clear();

	core::Size
	template_resid( SegmentResid const segment_resid ) const;

private:
	core::Size posestart_;
	core::Size movable_group_;
	core::Size saferes_;
	core::Size cutpoint_;
	std::string ss_;
	std::string abego_;
	bool nterm_included_;
	bool cterm_included_;
	std::string lower_segment_;
	std::string upper_segment_;
	core::pose::PoseOP template_pose_;
	ResidueDihedrals lower_dihedrals_;
	ResidueDihedrals upper_dihedrals_;
	core::conformation::ResidueCOP lower_residue_;
	core::conformation::ResidueCOP upper_residue_;

#ifdef    SERIALIZATION
public:
	template< class Archive >
	void
	save( Archive & arc ) const;

	template< class Archive >
	void
	load( Archive & arc );
#endif // SERIALIZATION

};

class NamedSegment : public std::pair< std::string, Segment > {
public:
	NamedSegment( std::string const & name, Segment const & res ):
		std::pair< std::string, Segment >( name, res ) {}
	NamedSegment( std::pair< std::string, Segment > const & resp ):
		std::pair< std::string, Segment >( resp ) {}
	friend std::ostream & operator<<( std::ostream & os, NamedSegment const & resis );
};

} // components
} // denovo_design
} // protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_denovo_design_components_Segment )
CEREAL_FORCE_DYNAMIC_INIT( protocols_denovo_design_components_ResidueDihedrals )
#endif // SERIALIZATION

#endif
