// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/scs_blast.hh
/// @brief Structural Component Selector (SCS) implementation with NCBI-BLAST+
/// @author Sergey Lyskov


#ifndef INCLUDED_protocols_antibody_grafting_scs_blast_hh
#define INCLUDED_protocols_antibody_grafting_scs_blast_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/scs_blast.fwd.hh>

#include <protocols/antibody/grafting/antibody_sequence.hh>
#include <protocols/antibody/grafting/chothia_numberer.hh>
#include <protocols/antibody/grafting/scs_functor.fwd.hh>
#include <protocols/antibody/grafting/exception.hh>

#include <core/types.hh>

#include <map>

namespace protocols {
namespace antibody {
namespace grafting {

/// @brief Base struct for any SCS single result. Hold information of selected template
struct SCS_Result
{
	virtual ~SCS_Result() {}  // This datatype is pure struct but we need to add at least one virtual method so we can use dynamic_cast later

	std::string pdb;
};


struct SCS_BlastResult : public SCS_Result
{
	int alignment_length;
	core::Real resolution, identity, bit_score;
	std::string bio_type, light_type, struct_source;

	/// sequences of selected template (do not confuse with querry sequences! 'sequence' will hold value corresponding to region results column)
	std::string /*sequence, */h1, h2, h3, frh, l1, l2, l3, frl;
};



struct SCS_ResultsVector : public utility::vector0< SCS_ResultOP >
{
	//std::string name, sequence;  // query info: name of the region and it sequence
};


struct SCS_Results
{
	//AntibodySequence const antibody_sequence; // our initial query for reference

	SCS_ResultsVector h1, h2, h3, l1, l2, l3, frh, frl, orientation;

	SCS_Results& operator=(SCS_Results&&) = default;
};


class SCS_Base {
	std::vector<SCS_FunctorCOP> filters_;
	SCS_FunctorCOP sorter_;

public:
	virtual ~SCS_Base() {}

	/// Results post-processing options: use filter to filter out resutls and sorting for chnaging results priority
	//void clear_filters();
	void add_filter(SCS_FunctorCOP filter);
	void set_sorter(SCS_FunctorCOP sorter);

	/// @brief Select CDR's template without filtering or sorting. In general you probably need to call select(...) instead
	/// @throw _AE_scs_failed_ on failure
	virtual SCS_ResultsOP raw_select(AntibodySequence const &) = 0;

	/// @brief Select CDR's template, filter it and sort
	/// @throw _AE_scs_failed_ on failure
	virtual SCS_ResultsOP select(AntibodySequence const &);
};




class SCS_BlastPlus : public SCS_Base {
public:

	/// @brief Select CDR's template
	/// @throw _AE_scs_failed_ on failure
	SCS_ResultsOP raw_select(AntibodySequence const &);  //, AntibodyNumbering const &);
};

struct FRH_FRL
{
	std::string frh1, frh2, frh3, frh4,
				frl1, frl2, frl3, frl4;

	// struct Residue {
	// 	char aa;
	// 	std::string number;
	// };
	// utility::vector0< Residue > heavy, light;
};

/// @brief Calculate antibody framework regions (frh1...frl4) and Trim it by removing less-preserved elements
void trim_framework(AntibodySequence const &A, AntibodyFramework &heavy_fr, AntibodyFramework &light_fr);


/// @brief Calculate frh+frl pair by trimming fr* and removing 'less-preserved' regions
FRH_FRL calculate_frh_frl(AntibodySequence const &);


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__


#endif // INCLUDED_protocols_antibody_grafting_scs_blast_hh
