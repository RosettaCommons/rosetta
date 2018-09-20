// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverBasicStreamOutputter.hh
/// @brief outputs a Peptiderive report to a stream in a basic, easily parsable format
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orly.marcu@mail.huji.ac.il)


#ifndef INCLUDED_protocols_peptide_deriver_PeptideDeriverBasicStreamOutputter_hh
#define INCLUDED_protocols_peptide_deriver_PeptideDeriverBasicStreamOutputter_hh

#include <protocols/peptide_deriver/PeptideDeriverBasicStreamOutputter.fwd.hh>
#include <protocols/peptide_deriver/PeptideDeriverOutputter.hh>
#include <protocols/peptide_deriver/PeptideDeriverFilter.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace protocols {
namespace peptide_deriver {

///@brief outputs a Peptiderive report to a stream in a basic, easily parsable format
class PeptideDeriverBasicStreamOutputter : public PeptideDeriverOutputter {

public:
	/// @param out    an output stream to which output will be diverted
	///               should belong to the thread on which the outputter
	///               methods are called in, and outlive the outputter object.
	/// @param prefix a string prefix to prepend to each output line
	PeptideDeriverBasicStreamOutputter(utility::io::orstream & out, std::string prefix);

	PeptideDeriverBasicStreamOutputter(PeptideDeriverBasicStreamOutputter const & src);

	~PeptideDeriverBasicStreamOutputter() override;

	PeptideDeriverBasicStreamOutputterOP clone() const;

	void begin_structure(core::pose::Pose const &, std::string const &) override;

	void chain_pair_pose_prepared(core::pose::Pose const & ) override {
		// do nothing
	}
	void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const total_isc,
		std::string const & options_string) override;

	void peptide_length(core::Size const pep_length) override;

	// TODO : consider using a PeptideDeriverEntry struct to make this signature less verbose and more tolerant to changes in the protocol
	//        e.g.,  virtual void peptide_entry(core::pose::Pose const & , PeptideDeriverEntry const & peptide_deriver_entry) --yuvals
	void peptide_entry(PeptideDeriverEntryType entry_type, core::Real const total_isc, DerivedPeptideEntryCOP entry) override;


	void end_receptor_partner_pair() override;

	void end_structure() override;

private:
	utility::io::orstream * out_p_;
	std::string prefix_;

}; // PeptideDeriverBasicStreamOutputter

} //protocols
} //peptide_deriver



#endif //INCLUDED_protocols_peptide_deriver_PeptideDeriverBasicStreamOutputter_hh
