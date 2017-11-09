// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverOutputterContainer.hh
/// @brief a container class which holds a list of PeptideDeriverOutputter instances and delegates calls to those
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orly.marcu@mail.huji.ac.il)


#ifndef INCLUDED_protocols_peptide_deriver_PeptideDeriverOutputterContainer_hh
#define INCLUDED_protocols_peptide_deriver_PeptideDeriverOutputterContainer_hh

#include <protocols/peptide_deriver/PeptideDeriverOutputterContainer.fwd.hh>
#include <protocols/peptide_deriver/PeptideDeriverOutputter.hh>


// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <string>

namespace protocols {
namespace peptide_deriver {

///@brief a container class which holds a list of PeptideDeriverOutputter instances and delegates calls to those
class PeptideDeriverOutputterContainer : public PeptideDeriverOutputter {
public:

	PeptideDeriverOutputterContainer();
	PeptideDeriverOutputterContainer(PeptideDeriverOutputterContainer const & src);

	virtual ~PeptideDeriverOutputterContainer();

	PeptideDeriverOutputterContainerOP
	clone() const;


	/// @brief add an outputter to the list.
	void push_back(PeptideDeriverOutputterOP const item);

	/// @brief clear all outputters in the list.
	void clear();

	void begin_structure(core::pose::Pose const & pose, std::string const &) override;

	void chain_pair_pose_prepared(core::pose::Pose const & pose) override;

	void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const total_isc,
		std::string const & options_string) override;

	void peptide_length(core::Size const pep_length) override;

	void peptide_entry(PeptideDeriverEntryType entry_type,  core::Real const total_isc, DerivedPeptideEntryCOP entry) override;


	void end_receptor_partner_pair() override;

	void end_structure() override;

private:
	utility::vector1<PeptideDeriverOutputterOP> list_;

}; // PeptideDeriverOutputterContainer

} //protocols
} //peptide_deriver


#endif //INCLUDED_protocols_peptide_deriver_PeptideDeriverOutputterContainer_hh

