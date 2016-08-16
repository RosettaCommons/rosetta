// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/scs_subst_mat.hh
/// @brief Structural Component Selector (SCS) implementation using substition matrices
/// @author Brian D. Weitzner


#ifndef INCLUDED_protocols_antibody_grafting_scs_subst_mat_hh
#define INCLUDED_protocols_antibody_grafting_scs_subst_mat_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__


#include <protocols/antibody/grafting/scs_subst_mat.fwd.hh>
#include <protocols/antibody/grafting/scs_blast.hh>


namespace protocols {
namespace antibody {
namespace grafting {

struct SCS_SubstitutionMatrixResult : public SCS_Antibody_Database_Result
{
	core::Real sid;
	core::Real score;

	/// sequences of selected template (do not confuse with query sequences! 'sequence' will hold value corresponding to region results column)
	std::string h1, h2, h3, frh, l1, l2, l3, frl;
};

class SCS_SubstitutionMatrix : public SCS_LoopOverSCs
{

public:

	/// @brief set working dir/output-prefix for intermediate files based on command-line options
	void init_from_options();

private:
	void select_template( Result & j,
	                      std::string const & db_to_query, std::map< std::string,
												std::map< std::string, std::string> > const & ab_db ) const override;
};



} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__


#endif // INCLUDED_protocols_antibody_grafting_scs_blast_hh
