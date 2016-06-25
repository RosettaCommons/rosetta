// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/libraries/BasePairStepLibrary.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rna_BasePairStepLibrary_HH
#define INCLUDED_protocols_rna_BasePairStepLibrary_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/farna/base_pairs/BasePairStep.hh>
#include <protocols/farna/libraries/BasePairStepLibrary.fwd.hh>
#include <protocols/farna/util.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <iostream>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace protocols::farna::base_pairs;

namespace protocols {
namespace farna {
namespace libraries {

int const MAX_BULGE_LENGTH( 3 );

////////////////////////////////////////////////////////////////////////
// BasePairStepSequence stores sequence of first dinucleotide strand and second dinucleotide strand --
// used as a 'key' for BasePairStepLibrary below.
class BasePairStepSequence: public utility::pointer::ReferenceCount {

public:

	//constructor
	BasePairStepSequence( char const nt_i, char const nt_i_next, char const nt_j, char const nt_j_next, Size const num_bulge = 0 );

	BasePairStepSequence( std::string const & sequence, Size const i, Size const i_next, Size const j, Size const j_next );

	BasePairStepSequence( std::string const & sequence, BasePairStep const & base_pair_step );

	//destructor
	~BasePairStepSequence(){}

	friend
	bool operator < (BasePairStepSequence const & lhs, BasePairStepSequence const & rhs )
	{
		//There must be a more elegant way to do this...
		if ( lhs.base_pair_step_sequence_.first.first < rhs.base_pair_step_sequence_.first.first ) {
			return true;
		} else if ( lhs.base_pair_step_sequence_.first.first == rhs.base_pair_step_sequence_.first.first ) {
			if ( lhs.base_pair_step_sequence_.first.second < rhs.base_pair_step_sequence_.first.second ) {
				return true;
			} else if ( lhs.base_pair_step_sequence_.first.second == rhs.base_pair_step_sequence_.first.second ) {
				if ( lhs.base_pair_step_sequence_.second.first < rhs.base_pair_step_sequence_.second.first ) {
					return true;
				} else if ( lhs.base_pair_step_sequence_.second.first == rhs.base_pair_step_sequence_.second.first ) {
					if ( lhs.base_pair_step_sequence_.second.second < rhs.base_pair_step_sequence_.second.second ) {
						return true;
					} else if  ( lhs.base_pair_step_sequence_.second.second == rhs.base_pair_step_sequence_.second.second ) {
						if ( lhs.num_bulge_ < rhs.num_bulge_ ) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}

	friend
	std::ostream &
	operator <<( std::ostream & os, BasePairStepSequence const & bps ){
		os << bps.base_pair_step_sequence_.first.first << "-" << bps.base_pair_step_sequence_.first.second;
		os << " ";
		os << bps.base_pair_step_sequence_.second.first;
		for ( Size n = 1; n <= Size(bps.num_bulge_); n++ ) os << "-n";
		os << "-" << bps.base_pair_step_sequence_.second.second;
		return os;
	}

	std::string
	tag() const {
		std::string bps_seq;
		bps_seq += base_pair_step_sequence_.first.first;
		bps_seq += base_pair_step_sequence_.first.second;
		bps_seq += "_";
		bps_seq += base_pair_step_sequence_.second.first;
		for ( Size n = 1; n <= Size(num_bulge_); n++ ) bps_seq += 'n';
		bps_seq += base_pair_step_sequence_.second.second;
		return bps_seq;
	}

	std::string
	subdir() const;

	int const &
	num_bulge() const { return num_bulge_; }

private:

	typedef std::pair< char, char > DinucleotideStrandSequence;
	std::pair< DinucleotideStrandSequence, DinucleotideStrandSequence > base_pair_step_sequence_;
	int num_bulge_;

};


////////////////////////////////////////////////////////////////////////
class BasePairStepLibrary: public utility::pointer::ReferenceCount {

public:

	//constructor
	BasePairStepLibrary( bool const canonical = true );

	//destructor
	~BasePairStepLibrary();

public:

	bool has_value( BasePairStepSequence const & base_pair_step_sequence ) const;

	// List of poses for each base pair step -- converted to
	// 'mini-pose' format with just coordinates.
	utility::vector1< core::pose::MiniPoseOP > const &
	mini_pose_list( BasePairStepSequence const & base_pair_step_sequence ) const;

	// example of a full pose for each base pair step
	pose::PoseCOP const &
	scratch_pose( BasePairStepSequence const & base_pair_step_sequence ) const;

	bool
	canonical() const { return canonical_; }

	std::string
	database_dir() const;

private:

	void
	initialize();

	void
	initialize_data( BasePairStepSequence const & base_pair_step_sequence,
		bool const load_in_poses = true ) const;

private:

	bool const canonical_;
	mutable std::map< BasePairStepSequence, utility::vector1< core::pose::MiniPoseOP > > mini_pose_lists_;
	mutable std::map< BasePairStepSequence, pose::PoseCOP > scratch_poses_;

};

} //libraries
} //farna
} //protocols

#endif
