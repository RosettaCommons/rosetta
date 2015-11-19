
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



/* Old desging with clear separation between filters and comparators
/// @brief Base class for comparison and filter operators
class SCS_Predicate
{
	mutable SCS_ResultsCOP results_; /// reference to whole result set in case operator need it
	mutable std::string region_;     /// name of currently compared region

	friend class SCS_Base;

public:
	//SCS_BlastOperator() {}
	virtual ~SCS_Predicate() {}


	// @brief return name of currently filterred region. Could be: h1, h2, h3, l1, l2, l3, frh, frl
	std::string region() const { return region_; }


	/// @brief Return reference to Blast+ results set containing results for all regions
	///        Note: results (length and order of elements) could change between different compare invocation
	/// @throw _AE_scs_failed_ if no results set was yet set (this should not occur under normal circumstances
	SCS_Results const & results() const;
};


class SCS_Filter : public SCS_Predicate
{
public:
	virtual bool operator()(SCS_ResultCOP &) const = 0; /// return true if results need to be filterd-out (removed) from results set
};



class SCS_Comparator : public SCS_Predicate
{
public:
	/// @brief comparison operator, should return true if first element is less then second (orderd before) and false otherwise
	virtual bool compare(SCS_ResultCOP &a, SCS_ResultCOP &b) const = 0;
};





/// Various Predicated for Blast+ SCS ----------------------------------------------------------
class SCS_BlastComparator : public SCS_Comparator
{
public:
	bool compare(SCS_ResultCOP &a, SCS_ResultCOP &b) const override {
		SCS_BlastResult const *aa = dynamic_cast< SCS_BlastResult const *>( a.get() );
		SCS_BlastResult const *bb = dynamic_cast< SCS_BlastResult const *>( b.get() );
		if( aa and bb ) return compare(*aa, *bb);
		else throw _AE_scs_failed_("SCS_BlastComparator::compare: Error! Could not cast SCS_Results to SCS_BlastResult!");
	}

	/// @brief comparison operator, should return true if first element is less then second (orderd before) and false otherwise
	virtual bool compare(SCS_BlastResult const &a, SCS_BlastResult const &b) const = 0;
};




class SCS_BlastFilter : public SCS_Filter
{
public:
	/// @brief return true if results need to be filterd-out (removed) from results set
    bool operator()(SCS_ResultCOP &o) const override {
		SCS_BlastResult const *oo = dynamic_cast< SCS_BlastResult const *>( o.get() );
		if(oo) return filter(*oo);
		else throw _AE_scs_failed_("SCS_BlastFilter::filter: Error! Could not cast SCS_Results to SCS_BlastResult!");
	}

	virtual bool filter(SCS_BlastResult const &) const = 0;
};




/// @brief filter results by sequence length
class SCS_BlastFilter_by_sequence_length : public SCS_BlastFilter
{
 public:
	bool filter(SCS_BlastResult const &r) const override;
};




/// @brief filter by alignment length
class SCS_BlastFilter_by_alignment_length : public SCS_BlastFilter
{
 public:
	bool filter(SCS_BlastResult const &r) const override;
};

*/

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_scs_functor_hh
