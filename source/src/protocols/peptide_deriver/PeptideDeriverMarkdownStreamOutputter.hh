// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverMarkdownStreamOutputter.hh
/// @brief outputs a Markdown formatted Peptiderive report to a stream
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orlymarcu@gmail.com)


#ifndef INCLUDED_protocols_peptide_deriver_PeptideDeriverMarkdownStreamOutputter_hh
#define INCLUDED_protocols_peptide_deriver_PeptideDeriverMarkdownStreamOutputter_hh

#include <protocols/peptide_deriver/PeptideDeriverMarkdownStreamOutputter.fwd.hh>
#include <protocols/peptide_deriver/PeptideDeriverOutputter.hh>

// Utility headers
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace protocols {
namespace peptide_deriver {

// TODO : allow easy switching with reST format
// TODO : conside creating a separate class(es) for simple Markdown/reST/HTML constructs (headers, tables, lists)
/// @brief outputs a Markdown formatted Peptiderive report to a stream
class PeptideDeriverMarkdownStreamOutputter : public PeptideDeriverOutputter {
public:

	/// @param out    an output stream to which output will be diverted
	///               should belong to the thread on which the outputter
	///               methods are called in, and outlive the outputter object.
	/// @param prefix a string prefix to prepend to each output line
	PeptideDeriverMarkdownStreamOutputter(utility::io::orstream & out, std::string prefix);

	PeptideDeriverMarkdownStreamOutputter(PeptideDeriverMarkdownStreamOutputter const & src);

	~PeptideDeriverMarkdownStreamOutputter() override;

	PeptideDeriverMarkdownStreamOutputterOP
	clone() const;

	void begin_structure(core::pose::Pose const & pose, std::string const &) override;

	// TODO : do we want to output anything here? Probably not.
	void chain_pair_pose_prepared(core::pose::Pose const &) override {
		// do nothing
	}

	void end_structure() override;

	void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const total_isc,
		std::string const & options_string) override;

	void peptide_length(core::Size const pep_length) override;

	void peptide_entry(PeptideDeriverEntryType entry_type, core::Real const total_isc, DerivedPeptideEntryCOP entry) override;

	void end_receptor_partner_pair() override;



private:

	class CyclizedReportInfo {
	public:

		core::Size best_cyclic_pep_start;
		bool cyclic_peptide_encountered_for_current_pep_length;
		std::ostringstream best_cyclic_peptides_;
		std::ostringstream cyclic_peptides_;
	};

	typedef utility::pointer::shared_ptr<CyclizedReportInfo> CyclizedReportInfoOP;
	typedef utility::pointer::shared_ptr<CyclizedReportInfo const> CyclizedReportInfoCOP;

	CyclizedReportInfoOP cyc_report_info_set_[NUM_CYCLIZATION_METHODS];

	/// clear accumulating member variables
	void clear_buffers();

	utility::io::orstream * out_p_;
	std::string prefix_;

	core::Size current_pep_length_;
	core::Real current_total_isc_;
	char current_receptor_chain_letter_;
	char current_partner_chain_letter_;
	std::string current_sequence_;

	std::ostringstream header_;
	std::ostringstream best_linear_peptides_;
	std::ostringstream all_peptides_;
	std::ostringstream footer_;

	/// store location of best cyclic peptide for comparison reasons

}; // PeptideDeriverMarkdownStreamOutputter


} //protocols
} //peptide_deriver



#endif //INCLUDED_protocols_peptide_deriver_PeptideDeriverMarkdownStreamOutputter_hh
