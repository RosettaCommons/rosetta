// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/WriteUpstreamCoordinateKinemage.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_WriteUpstreamCoordinateKineamge_hh
#define INCLUDED_protocols_match_output_WriteUpstreamCoordinateKineamge_hh

// Unit headers
#include <protocols/match/output/WriteUpstreamCoordinateKineamge.fwd.hh>

// Package headers
#include <protocols/match/upstream/UpstreamBuilder.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/ClassicMatchAlgorithm.fwd.hh>


// Project headers
#if defined(WIN32) || defined(PYROSETTA)
#include <core/id/AtomID.hh>
#endif

// C++ headers
#include <fstream>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace match {
namespace output {


class ResidueKinemageWriter
{
public:
	typedef core::Size   Size;

public:
	ResidueKinemageWriter();

	/// @brief Write out the coordinates for a particular residue; the kinemage tag
	/// is assumed to have been writen already.
	void
	write_rsd_coords(
		std::ostream & ostr,
		Size const scaffold_build_point_id,
		Size const upstream_conf_id,
		core::conformation::Residue const & rsd,
		bool is_instance = false
	) const;

	void dominant( bool setting );
	void animate(  bool setting );
	void group(    bool setting );
	void write_virtual_atoms( bool setting );

	void master( std::string const & setting );

private:
	std::string master_;
	bool dominant_;
	bool animate_;
	bool group_; // false for subgroup
	bool write_virtual_atoms_;
};

class WriteUpstreamCoordinateKinemage : public downstream::DownstreamAlgorithm
{
public:
	typedef core::Size   Size;
	typedef core::Vector Vector;

public:
	WriteUpstreamCoordinateKinemage();
	WriteUpstreamCoordinateKinemage( std::string const & fname );
	WriteUpstreamCoordinateKinemage( std::ostream & ostr );

	virtual ~WriteUpstreamCoordinateKinemage();

	virtual
	downstream::DownstreamAlgorithmOP
	clone() const;

	virtual
	std::list< Hit >
	build(
		Size const scaffold_build_point_id,
		Size const upstream_conf_id,
		core::conformation::Residue const & upstream_residue
	) const;

	/// @brief This method returns 'true' whether or not it's ClassicMatchAlgorithm is set
	/// as it should not have its hits_to_include_with_partial_match method invoked.
	virtual
	bool
	upstream_only() const;

	/// @brief This method returns 'true' since when it does return hits, it's those generated
	/// by the ClassicMatchAlgorithm
	virtual
	bool
	generates_primary_hits() const;


	/// @brief This method should not be invoked on this class,
	/// since it returns "true" in its upstream_only
	/// method.
	virtual
	HitPtrListCOP
	hits_to_include_with_partial_match( match_dspos1 const & m ) const;

	virtual
	Size
	n_possible_hits_per_upstream_conformation() const;

	void
	set_kinemage_file_name( std::string const & filename );

	void
	set_match_algorithm( downstream::ClassicMatchAlgorithmCOP algorithm );

	void
	set_downstream_writer( DownstreamCoordinateKinemageWriterCOP dswriter );

	void
	set_n_downstream_to_output( Size n_downstream_to_output );

	bool
	return_pseudo_hits() const {
		return return_pseudo_hits_;
	}

	void
	return_pseudo_hits( bool setting ) {
		return_pseudo_hits_ = setting;
	}

private:
	std::string kinemage_file_name_;
	std::ofstream file_out_;
	std::ostream & ostr_;
	mutable Size last_scaffold_build_point_;
	mutable Size nkins_;

	downstream::ClassicMatchAlgorithmCOP match_algorithm_;
	DownstreamCoordinateKinemageWriterCOP dswriter_;
	Size n_downstream_to_output_;
	mutable Size n_output_so_far_;

	bool return_pseudo_hits_;
};

class WriteUpstreamHitKinemage : public upstream::UpstreamResidueProcessor
{
public:
	WriteUpstreamHitKinemage();
	WriteUpstreamHitKinemage( std::string const & fname );
	WriteUpstreamHitKinemage( std::ostream & ostr );

	virtual ~WriteUpstreamHitKinemage();

	virtual
	void
	process_hit(
		Hit const & hit,
		core::conformation::Residue const & upstream_conformation
	);

	/// @brief Non-virtual method to write out a kinemage for an upstream residue
	void
	output_hit(
		Hit const & hit,
		core::conformation::Residue const & upstream_conformation
	);

	/// @brief non-virtual method to write out a kinemage for an upstream residue;
	/// without invoking dswriter_->write_downstream_coordinates().
	/// Only information about the upstream portion of the hit is needed.
	void
	output_upstream_coordinates(
		upstream_hit const & hit,
		core::conformation::Residue const & upstream_conformation
	);

	void
	start_new_match();

	void
	set_kinemage_file( std::string const & fname );

	void
	set_dswriter( DownstreamCoordinateKinemageWriterOP dswriter );

	void geom_id( Size setting );

	/// @brief Set the kinemage master for the upstream residue, overriding
	/// the default master, which is "geom#"
	void
	set_master( std::string const & master );

	/// @brief Set whether the default master should be used, or whether the
	/// user-defined master should be used.  If the user defined master is the
	/// empty string, no master will be writen.
	void
	default_master( bool setting );

	/// @brief Returns whether or not the default master is being used.
	bool
	default_master() const;

	void animate( bool setting );
	void dominant( bool setting );
	void group( bool setting );
	void write_virtual_atoms( bool setting );

private:
	Size matches_output_count_;

	bool use_default_master_;
	std::string master_;

	bool animate_;
	bool dominant_;
	bool group_;
	bool write_virtual_atoms_;

	Size geom_id_;
	std::string kinemage_file_name_;
	std::ofstream file_out_;
	std::ostream & ostr_;

	DownstreamCoordinateKinemageWriterOP dswriter_;
};

class DownstreamCoordinateKinemageWriter : public utility::pointer::ReferenceCount
{
public:
	DownstreamCoordinateKinemageWriter();
	virtual ~DownstreamCoordinateKinemageWriter();

	virtual
	void
	write_downstream_coordinates(
		Hit const & hit,
		std::ostream & ostr
	) const = 0;

	virtual
	void
	set_downstream_master( std::string const & str ) = 0;

};

/// @brief Class for writing conformations of the downstream partner in a kinemage
/// description.
class SingleDownstreamResidueWriter : public DownstreamCoordinateKinemageWriter
{
public:
	SingleDownstreamResidueWriter();
	virtual ~SingleDownstreamResidueWriter();

	virtual
	void
	write_downstream_coordinates(
		Hit const & hit,
		std::ostream & ostr
	) const;

	void
	set_restype( core::chemical::ResidueTypeCOP );

	void
	set_downstream_builder( downstream::DownstreamBuilderCOP dsbuilder );

	virtual
	void
	set_downstream_master( std::string const & master );

private:

	std::string master_;
	core::chemical::ResidueTypeCOP restype_;
	utility::vector1< core::id::AtomID >       all_atom_inds_;
	downstream::DownstreamBuilderCOP           dsbuilder_;

};

}
}
}

#endif
