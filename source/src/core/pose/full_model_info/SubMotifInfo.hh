// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/full_model_info/SubMotifInfo.hh
/// @brief  Stores information about submotifs in a pose.
/// @author Caleb Geniesse

#ifndef INCLUDED_core_pose_full_model_info_SubMotifInfo_hh
#define INCLUDED_core_pose_full_model_info_SubMotifInfo_hh


// Project headers
#include <core/pose/full_model_info/SubMotifInfo.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>


namespace core {
namespace pose {
namespace full_model_info {

//////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of information about a submotif added to a pose.
class SubMotifInfo : public utility::pointer::ReferenceCount {

public:

	SubMotifInfo();

	SubMotifInfo(
		utility::vector1< Size > const & res_list,
		std::string const & tag,
		bool const & seed = false
	);

	SubMotifInfo( SubMotifInfo const & src );

	~SubMotifInfo();

	SubMotifInfoOP
	clone() const
	{
		return SubMotifInfoOP( new SubMotifInfo( *this ) );
	}

	// setters/getters
	utility::vector1< Size > const & res_list() const { return res_list_; }
	Size const & res_list( Size const & index ) const { return res_list().at( index ); }
	void res_list( utility::vector1< core::Size > const & res_list ) { res_list_ = res_list; }

	utility::vector1< Size > sorted_res_list() const;
	Size sorted_res_list( Size const & index ) { return sorted_res_list().at( index ); }

	std::string const & tag() const { return tag_; }
	void tag( std::string const & tag ) { tag_ = tag; }

	bool const & seed() const { return seed_; }
	void seed( bool const & seed ) { seed_ = seed; }

public:
	/// @brief Equality comparator
	friend bool operator ==( SubMotifInfoOP lhs, SubMotifInfoOP rhs );

	/// @brief Equality comparator
	friend bool operator !=( SubMotifInfoOP lhs, SubMotifInfoOP rhs );

	/// @brief << operator
	friend std::ostream &
	operator <<( std::ostream & os, SubMotifInfoOP submotif_info );

	/// @brief << operator
	friend std::istream &
	operator >>( std::istream & is, SubMotifInfoOP submotif_info );


private:

	utility::vector1< Size > res_list_;
	std::string tag_;
	bool seed_;

};


} //full_model_info
} //pose
} //core
#endif
