// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/Hit.hh
/// @brief  Hit typedef
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_Hit_hh
#define INCLUDED_protocols_match_Hit_hh

// Unit headers
#include <protocols/match/Hit.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace match {

class Hit
{
public:
	core::Size
	scaffold_build_id() const {
		return first_[ 1 ];
	}

	core::Size
	upstream_conf_id() const {
		return first_[ 2 ];
	}

	core::Size
	external_geom_id() const {
		return first_[ 3 ];
	}


	core::Size
	downstream_conf_id() const {
		return first_[ 4 ];
	}

	Size4 const & first() const {
		return first_;
	}

	Size4 & first() {
		return first_;
	}

	Real6 const & second() const {
		return second_;
	}

	Real6 & second() {
		return second_;
	}

private:

	Size4 first_;
	Real6 second_;

};

class upstream_hit
{
public:
	upstream_hit() {}

	upstream_hit( Hit const & source ) {
		first_[ 1 ] = source.first()[ 1 ];
		first_[ 2 ] = source.first()[ 2 ];
		first_[ 3 ] = source.first()[ 3 ];
	}

	void
	copy_hit( Hit const & source ) {
		first_[ 1 ] = source.first()[ 1 ];
		first_[ 2 ] = source.first()[ 2 ];
		first_[ 3 ] = source.first()[ 3 ];
	}

	core::Size
	scaffold_build_id() const {
		return first_[ 1 ];
	}

	core::Size
	upstream_conf_id() const {
		return first_[ 2 ];
	}

	core::Size
	external_geom_id() const {
		return first_[ 3 ];
	}

	void
	scaffold_build_id( core::Size setting ) {
		first_[ 1 ] = setting;
	}

	void
	upstream_conf_id( core::Size setting ) {
		first_[ 2 ] = setting;
	}

	void
	external_geom_id( core::Size setting ) {
		first_[ 3 ] = setting;
	}

	bool
	operator < ( upstream_hit const & rhs ) const;

	bool
	operator == ( upstream_hit const & rhs ) const;


private:
	Size3 first_;
};


class downstream_hit
{

public:

	downstream_hit(){}

	downstream_hit( Hit const & source ) {
		downstream_conf_id_ = source.first()[ 4 ];
		second_ = source.second();
	}

	void
	copy_hit( Hit const & source ) {
		downstream_conf_id_ = source.first()[ 4 ];
		second_ = source.second();
	}

	core::Size
	downstream_conf_id() const {
		return downstream_conf_id_;
	}

	Real6 const & second() const {
		return second_;
	}

	Real6 & second() {
		return second_;
	}

	bool
	operator < ( downstream_hit const & rhs ) const;

	bool
	operator == ( downstream_hit const & rhs ) const;


private:
	core::Size downstream_conf_id_;
	Real6 second_;
};


/// @brief Describe a match as n_geometric_constraint upstream residue conformations and
/// one positioning of the downstream partner ( "dspos1" = 1 downstrem position)
struct match_dspos1
{
public:
	match_dspos1();
	match_dspos1( core::Size n_geometric_constraints );
	match_dspos1( match const & m, core::Size geomcst_specifying_dspos );

	utility::vector1< upstream_hit > upstream_hits;
	core::Size originating_geom_cst_for_dspos;
	core::Size downstream_conf_id;
	Real6 dspos;
};

/// @brief Create a fake hit from an upstream_hit where hit.first()[4] and hit.second() are 0's.
Hit fake_hit( upstream_hit const & );

/// @brief Create a fake hit from a downstream_hit where hit.first()[1-3] are 0's.
Hit fake_hit( downstream_hit const & );

/// @brief Create a hit with the full data from a given match_dspos1 representing
/// the upstream conformation from the originating_geom_cst and its
/// description of the downstream position.
Hit full_hit( match_dspos1 const & m );

}
}

#endif
