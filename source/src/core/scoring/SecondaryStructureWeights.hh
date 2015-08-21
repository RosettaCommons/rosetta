// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_SecondaryStructureWeights_hh
#define INCLUDED_core_scoring_SecondaryStructureWeights_hh

#include <core/types.hh>


#include <utility/pointer/ReferenceCount.hh>


// ObjexxFCL Headers

// C++ headers

namespace core {
namespace scoring {


/// @brief  Holds weights and flags for configuring a SecondaryStructureEnergy evaluation

class SecondaryStructureWeights : public utility::pointer::ReferenceCount {

public: // construct/destruct

	/// @brief default constructor
	inline
	SecondaryStructureWeights()
	{
		initialize();
	}

	/// @brief copy constructor
	inline
	SecondaryStructureWeights(
		SecondaryStructureWeights const & s
	) : ReferenceCount(),
		parallel_weight_( s.parallel_weight_ ),
		antiparallel_weight_( s.antiparallel_weight_ ),

		ss_lowstrand_( s.ss_lowstrand_ ),
		ss_cutoff_( s.ss_cutoff_ ),
		localstrandpair_penalty_( s.localstrandpair_penalty_ ),
		seq_sep_scale_( s.seq_sep_scale_ ),
		max_strand_dist_cutoff_( s.max_strand_dist_cutoff_ ),
		strand_dist_cutoff_( s.strand_dist_cutoff_ ),
		stretch_strand_dist_cutoff_( s.stretch_strand_dist_cutoff_ ),
		handedness_score_flag_( s.handedness_score_flag_ )
	{}

	/// @brief default destructor
	virtual
	~SecondaryStructureWeights();


public: // assignment

	/// @brief copy assignment
	inline
	SecondaryStructureWeights &
	operator =( SecondaryStructureWeights const & s )
	{
		if ( this != &s ) {
			parallel_weight_ = s.parallel_weight_;
			antiparallel_weight_ = s.antiparallel_weight_;

			ss_lowstrand_ = s.ss_lowstrand_;
			ss_cutoff_ = s.ss_cutoff_;
			localstrandpair_penalty_ = s.localstrandpair_penalty_;
			seq_sep_scale_ = s.seq_sep_scale_;
			max_strand_dist_cutoff_ = s.max_strand_dist_cutoff_;
			strand_dist_cutoff_ = s.strand_dist_cutoff_;
			stretch_strand_dist_cutoff_ = s.stretch_strand_dist_cutoff_;
			handedness_score_flag_ = s.handedness_score_flag_;
		}
		return *this;
	}


public: // bulk setup

	/// @brief setup parallel/antiparallel weights according to scheme
	void
	setup_parallel_antiparallel_weights(
		bool const & randomize_weights = false
	);

public: // get accessors

	/// @brief operator==
	friend
	inline
	bool
	operator ==( SecondaryStructureWeights const & a, SecondaryStructureWeights const & b )
	{
		return ( ( a.parallel_weight_            == b.parallel_weight_ ) &&
			( a.antiparallel_weight_        == b.antiparallel_weight_ ) &&
			( a.ss_lowstrand_               == b.ss_lowstrand_ ) &&
			( a.ss_cutoff_                  == b.ss_cutoff_ ) &&
			( a.localstrandpair_penalty_    == b.localstrandpair_penalty_ ) &&
			( a.seq_sep_scale_              == b.seq_sep_scale_ ) &&
			( a.max_strand_dist_cutoff_     == b.max_strand_dist_cutoff_ ) &&
			( a.strand_dist_cutoff_         == b.strand_dist_cutoff_ ) &&
			( a.stretch_strand_dist_cutoff_ == b.stretch_strand_dist_cutoff_ ) &&
			( a.handedness_score_flag_      == b.handedness_score_flag_ ) );
	}


	/// @brief parallel strand weight
	inline
	Real const &
	get_parallel_weight() const
	{
		return parallel_weight_;
	}

	/// @brief antiparallel strand weight
	inline
	Real const &
	get_antiparallel_weight() const
	{
		return antiparallel_weight_;
	}


	inline
	int
	get_ss_lowstrand() const
	{
		return ss_lowstrand_;
	}


	inline
	int
	get_ss_cutoff() const
	{
		return ss_cutoff_;
	}

	/// @brief local strand pair penalty
	inline
	Real const &
	get_localstrandpair_penalty() const
	{
		return localstrandpair_penalty_;
	}

	/// @brief sequence separation scale
	inline
	Real const &
	get_seq_sep_scale() const
	{
		return seq_sep_scale_;
	}

	/// @brief max strand distance cutoff
	inline
	Real const &
	get_max_strand_dist_cutoff() const
	{
		return max_strand_dist_cutoff_;
	}

	/// @brief strand distance cutoff
	inline
	Real const &
	get_strand_dist_cutoff() const
	{
		return strand_dist_cutoff_;
	}

	/// @brief stretching the strand distance cutoff?
	inline
	bool const &
	get_stretch_strand_dist_cutoff() const
	{
		return stretch_strand_dist_cutoff_;
	}

	/// @brief using handedness score?
	inline
	bool const &
	get_handedness_score_flag() const
	{
		return handedness_score_flag_;
	}


public: // set accessors

	/// @brief set parallel strand weight
	inline
	void
	set_parallel_weight(
		Real const & weight
	)
	{
		parallel_weight_ = weight;
	}

	/// @brief set antiparallel strand weight
	inline
	void
	set_antiparallel_weight(
		Real const & weight
	)
	{
		antiparallel_weight_ = weight;
	}


	void
	set_ss_lowstrand( int const setting )
	{
		ss_lowstrand_ = setting;
	}


	void
	set_ss_cutoff( int const setting )
	{
		ss_cutoff_ = setting;
	}


	/// @brief set local strand pair penalty
	inline
	void
	set_localstrandpair_penalty(
		Real const & penalty
	)
	{
		localstrandpair_penalty_ = penalty;
	}

	/// @brief set sequence separation scale
	inline
	void
	set_seq_sep_scale(
		Real const & scale
	)
	{
		seq_sep_scale_ = scale;
	}

	/// @brief set max strand distance cutoff
	inline
	void
	set_max_strand_dist_cutoff(
		Real const & cutoff
	)
	{
		max_strand_dist_cutoff_ = cutoff;
	}

	/// @brief set strand distance cutoff
	inline
	void
	set_strand_dist_cutoff(
		Real const & cutoff
	)
	{
		strand_dist_cutoff_ = cutoff;
	}

	/// @brief stretch the strand distance cutoff
	inline
	void
	set_stretch_strand_dist_cutoff(
		bool const & flag
	)
	{
		stretch_strand_dist_cutoff_ = flag;
	}

	/// @brief use the handedness score
	inline
	void
	set_handedness_score_flag(
		bool const & flag
	)
	{
		handedness_score_flag_ = flag;
	}


private: // initialization

	/// @brief initialize to default values, also load proper data bins from external files
	void
	initialize();

private: // weights

	Real parallel_weight_;
	Real antiparallel_weight_;


private: // additional settings

	int ss_lowstrand_;
	int ss_cutoff_;

	Real localstrandpair_penalty_;
	Real seq_sep_scale_;
	Real max_strand_dist_cutoff_;
	Real strand_dist_cutoff_;
	bool stretch_strand_dist_cutoff_;
	bool handedness_score_flag_;


};


} // ns scoring
} // ns core

#endif
