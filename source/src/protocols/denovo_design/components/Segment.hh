// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/Segment.hh
/// @brief Segment functions for building structures from components
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_Segment_hh
#define INCLUDED_protocols_denovo_design_components_Segment_hh

// Unit headers
#include <protocols/denovo_design/components/Segment.fwd.hh>

// Project headers
#include <protocols/denovo_design/types.hh>

// Protocol headers

// Core headers

// Basic/Numeric/Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers

namespace protocols {
namespace denovo_design {
namespace components {

/// @brief manages information about segments of residues
class Segment : public utility::pointer::ReferenceCount {
public:
	Segment();

	Segment(
		core::Size const pose_length_val,
		core::Size const local_saferes,
		core::Size const local_cutpoint,
		core::Size const movable_group_val,
		bool const is_loop_val,
		bool const start_inc,
		bool const stop_inc,
		std::string const & lower,
		std::string const & upper,
		std::string const & ss_val,
		utility::vector1< std::string > const & abego_val );

	virtual ~Segment() {}

	/// @brief construct from xml tag
	void parse_tag( utility::tag::TagCOP tag );

	core::Size resid( core::Size const local_resnum ) const;
	core::Size start() const;
	core::Size stop() const;
	core::Size safe() const;
	core::Size cutpoint() const;

	inline std::string const & ss() const { return ss_; }
	inline utility::vector1< std::string > const & abego() const { return abego_; }
	inline utility::vector1< std::string > & abego_nonconst() { return abego_; }

	inline bool nterm_included() const { return nterm_included_; }
	inline bool cterm_included() const { return cterm_included_; }

	void set_nterm_included( bool const ntermval );
	void set_cterm_included( bool const ctermval );

	inline bool has_free_lower_terminus() const { return ( lower_segment_ == "" ); }
	inline bool has_free_upper_terminus() const { return ( upper_segment_ == "" ); }

	inline std::string const & upper_segment() const { return upper_segment_; }
	inline std::string const & lower_segment() const { return lower_segment_; }


	bool contains( core::Size const res ) const
	{
		return ( ( nterm_resi() <= res ) && ( res <= cterm_resi() ) );
	}

	core::Size nterm_resi() const;
	core::Size cterm_resi() const;

	core::Size nterm_pad() const
	{
		return start_ - 1;
	}

	core::Size cterm_pad() const
	{
		return cterm_resi() - stop();
	}

	core::Size length() const
	{
		return cterm_resi() - nterm_resi() + 1;
	}

	core::Size elem_length() const
	{
		assert( stop_ >= start_ );
		return stop_ - start_ + 1;
	}

	inline void set_pose_start( core::Size const res )
	{
		posestart_ = res;
	}

	void set_upper_segment( std::string const & comp ) { upper_segment_ = comp; }
	void set_lower_segment( std::string const & comp ) { lower_segment_ = comp; }

	/// @brief sets cutpoint for this segment to be the ith residue
	void set_cutpoint( core::Size const cut_res ) { cutpoint_ = cut_res; }

	/// @brief deletes dummy residues before the fixed portion of this segment
	void delete_leading_residues();

	/// @brief deletes dummy residues after the fixed portion of this segment
	void delete_trailing_residues();

	/// @brief expands this residue set to include the dummy trailing residues
	void engulf_leading_residues();

	/// @brief expands this residue set to include the dummy trailing residues
	void engulf_trailing_residues();

	/// @brief given a residue number range local to this 1=start, length=end, delete the residue
	void delete_residues( core::Size const local_resnum_start, core::Size const local_resnum_stop );

	// I/O
public:
	std::string serialize() const;

	/// output residueinfo
	friend std::ostream &
	operator<<( std::ostream & os, Segment const & res );

	core::Size movable_group;
	bool is_loop;

private:
	core::Size posestart_;
	core::Size start_;
	core::Size stop_;
	core::Size saferes_;
	core::Size cutpoint_;
	std::string ss_;
	utility::vector1< std::string > abego_;
	bool nterm_included_;
	bool cterm_included_;
	std::string lower_segment_;
	std::string upper_segment_;
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

#endif
