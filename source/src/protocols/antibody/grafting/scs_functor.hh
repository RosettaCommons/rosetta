
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/scs_functor.hh
/// @brief Structural Component Selector (SCS): implemetation of predicates for filtering and sorting
/// @author Sergey Lyskov


#ifndef INCLUDED_protocols_antibody_grafting_scs_functor_hh
#define INCLUDED_protocols_antibody_grafting_scs_functor_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/scs_functor.fwd.hh>

#include <protocols/antibody/grafting/antibody_sequence.hh>
#include <protocols/antibody/grafting/scs_blast.hh>

#include <string>

namespace protocols {
namespace antibody {
namespace grafting {


class SCS_Functor
{
public:
	virtual ~SCS_Functor() {}

	virtual void apply(AntibodySequence const &antibody_sequence, SCS_ResultsOP results) const = 0;
};


class SCS_Comparator : public SCS_Functor
{
public:
	virtual void apply(AntibodySequence const &antibody_sequence, SCS_ResultsOP results) const override;
	virtual bool compare(AntibodySequence const &antibody_sequence, SCS_Result const &a, SCS_Result const &b) const = 0;
};


class SCS_BlastComparator : public SCS_Comparator
{
public:
    virtual bool compare(AntibodySequence const &antibody_sequence, SCS_Result const &a, SCS_Result const &b) const override;
	virtual bool compare(AntibodySequence const &antibody_sequence, SCS_BlastResult const &a, SCS_BlastResult const &b) const = 0;
};



class SCS_BlastComparator_BitScore_Resolution : public SCS_BlastComparator
{
public:
	bool compare(AntibodySequence const &, SCS_BlastResult const &a, SCS_BlastResult const &b) const override {
		if( a.bit_score > b.bit_score ) return true;
		if( a.bit_score == b.bit_score ) return a.resolution < b.resolution;
		return false;
	}
};


class SCS_BlastFilter_by_sequence_length : public SCS_Functor
{
public:
	void apply(AntibodySequence const &antibody_sequence, SCS_ResultsOP results) const override;
};


class SCS_BlastFilter_by_alignment_length : public SCS_Functor
{
 public:
	void apply(AntibodySequence const &antibody_sequence, SCS_ResultsOP results) const override;
};


/// @details filter helper function: generate string with results sizes
std::string result_sizes(SCS_ResultsOP r, int width=4);


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_scs_functor_hh
