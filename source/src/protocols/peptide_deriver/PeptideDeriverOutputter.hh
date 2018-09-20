// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverOutputter.hh
/// @brief an abstract base class for handling calculation outputs from PeptideDeriverFilter
/// @author orlypolo (orlymarcu@gmail.com)


#ifndef INCLUDED_protocols_peptide_deriver_PeptideDeriverOutputter_hh
#define INCLUDED_protocols_peptide_deriver_PeptideDeriverOutputter_hh

#include <protocols/peptide_deriver/PeptideDeriverOutputter.fwd.hh>
#include <protocols/peptide_deriver/PeptideDeriverFilter.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace protocols {
namespace peptide_deriver {

/// @brief an abstract base class for handling calculation outputs from
///        PeptideDeriverFilter
/// Since PeptideDeriverFilter might have a set of outputs for each residue,
/// for each chain in each chain-pair, outputting is quite elaborate. This is
/// an attempt to decouple the calculation from the representation of results.
/// Representation of results is delegated to implementors of this class.
class PeptideDeriverOutputter {
public:

	/// @brief default constructor
	PeptideDeriverOutputter();

	/// @brief destructor
	virtual ~PeptideDeriverOutputter();

	/// @brief called by PeptideDeriverFilter when processing of a strucuture (possibly multi-chain) starts
	virtual void begin_structure(core::pose::Pose const &, std::string const &) = 0;

	/// @brief called by PeptideDeriverFilter when a chain-pair pose is prepared (minimized and disulfide-bridge resolved)
	virtual void chain_pair_pose_prepared(core::pose::Pose const & pose) = 0;

	/// @brief called by PeptideDeriverFilter when calculation commences on a receptor-partner pair
	virtual void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const total_isc,
		std::string const & options_string) = 0;

	/// @brief called by PeptideDeriverFilter when calculation commences for
	///        the specified peptide length
	/// this is called for every peptide length specified to
	/// PeptideDeriverFilter, per each call to
	/// PeptideDeriverOutputter::begin_receptor_partner_pair()
	virtual void peptide_length(core::Size const pep_length) = 0;

	/// @brief called by PeptideDeriverFilter for each 'peptide' entry that should be output.
	/// besides the peptides for each sliding window, this is being called for the 'best' peptides (@see PeptideDeriverEntryType)
	/*
	* virtual void peptide_entry(core::pose::Pose const & pose, PeptideDeriverEntryType const entry_type, core::Size const pep_start,
	core::Real const linear_isc, core::Real const binding_contribution_fraction, std::string const & disulfide_info, bool const was_cyclic_pep_modeled,
	core::pose::Pose const & cyclic_pose, core::Real const cyclic_isc) = 0;
	*/
	virtual void peptide_entry(PeptideDeriverEntryType entry_type, core::Real const total_isc, DerivedPeptideEntryCOP entry) = 0;

	/// @brief called by PeptideDeriverFilter when calculation concludes for a receptor-partner pair (for all the different peptide lengths)
	virtual void end_receptor_partner_pair() = 0;

	/// @brief called by PeptideDeriverFilter when processing of a strucuture (all chain-pairs considered) ends
	virtual void end_structure() = 0;

	core::Real avoid_negative_zero(core::Real const value, core::Real const threshold);

};

} //protocols
} //peptide_deriver



#endif //INCLUDED_protocols_peptide_deriver_PeptideDeriverOutputter_hh
