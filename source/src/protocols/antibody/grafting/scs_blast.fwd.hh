// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/scs_blast.fwd.hh
/// @brief Structural Component Selector (SCS) implementation with NCBI-BLAST+
/// @author Sergey Lyskov


#ifdef CXX11

#ifndef INCLUDED_protocols_antibody_grafting_scs_blast_fwd_hh
#define INCLUDED_protocols_antibody_grafting_scs_blast_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace antibody {
namespace grafting {


struct SCS_Result;
typedef utility::pointer::shared_ptr< SCS_Result > SCS_ResultOP;
typedef utility::pointer::shared_ptr< SCS_Result const > SCS_ResultCOP;


struct SCS_BlastResult;
typedef utility::pointer::shared_ptr< SCS_BlastResult > SCS_BlastResultOP;
typedef utility::pointer::shared_ptr< SCS_BlastResult const > SCS_BlastResultCOP;


struct SCS_Results;
typedef utility::pointer::shared_ptr< SCS_Results > SCS_ResultsOP;
typedef utility::pointer::shared_ptr< SCS_Results const > SCS_ResultsCOP;


struct SCS_ResultSet;
typedef utility::pointer::shared_ptr< SCS_ResultSet > SCS_ResultSetOP;
typedef utility::pointer::shared_ptr< SCS_ResultSet const > SCS_ResultSetCOP;

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // INCLUDED_protocols_antibody_grafting_scs_blast_fwd_hh

#endif // CXX11
