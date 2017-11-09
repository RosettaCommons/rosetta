// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverPoseOutputter.hh
/// @brief outputs poses at different points of the Peptiderive protocol, according to the given set of options
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orlymarcu@gmail.com)


#ifndef INCLUDED_protocols_peptide_deriver_PeptideDeriverPoseOutputter_hh
#define INCLUDED_protocols_peptide_deriver_PeptideDeriverPoseOutputter_hh

#include <protocols/peptide_deriver/PeptideDeriverPoseOutputter.fwd.hh>
#include <protocols/peptide_deriver/PeptideDeriverOutputter.hh>

// Project headers
#include <core/scoring/ScoreFunction.hh>
#include <protocols/filters/Filter.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iosfwd>
#include <string>

namespace protocols {
namespace peptide_deriver {


/// @brief outputs poses at different points of the Peptiderive protocol, according to the given set of options
class PeptideDeriverPoseOutputter : public PeptideDeriverOutputter {
public:
	/// @param is_dump_best_peptide_pose flag to output the best (in terms of linea/cyclic isc) poses. Max. 4 poses per chain per.
	/// @param is_dump_prepared_pose     flag to output each chain-pair pose as Peptiderive sees it
	/// @param is_dump_cyclic_poses      flag to output each modeled cyclic peptide (@see basic::options::OptionKeys::peptide_deriver::optimize_cyclic_threshold)
	/// @param scorefxn                  score function to score the pose with before outputting
	PeptideDeriverPoseOutputter( bool const is_dump_best_peptide_pose,
		bool const is_dump_prepared_pose, bool const is_dump_cyclic_poses,
		core::scoring::ScoreFunctionCOP scorefxn);

	PeptideDeriverPoseOutputter(PeptideDeriverPoseOutputter const & src);

	virtual ~PeptideDeriverPoseOutputter();

	PeptideDeriverPoseOutputterOP
	clone() const;

	void begin_structure(core::pose::Pose const &, std::string const &) override {
		// do nothing
	}

	void chain_pair_pose_prepared(core::pose::Pose const & pose) override;

	void begin_receptor_partner_pair(char const receptor_chain_letter,
		char const partner_chain_letter, core::Real const,
		std::string const &) override;

	void peptide_length(core::Size const pep_length) override;

	void end_receptor_partner_pair() override {
		// do nothing
	}

	void end_structure() override {
		// do nothing
	}

	void peptide_entry(PeptideDeriverEntryType entry_type, core::Real const total_isc, DerivedPeptideEntryCOP entry) override;

private:
	/// output a post to the given file
	void output_pose( core::pose::Pose & pose, std::string const & pose_name );

	bool is_dump_best_peptide_pose_;
	bool is_dump_prepared_pose_;
	bool is_dump_cyclic_poses_;

	char current_receptor_chain_letter_;
	char current_partner_chain_letter_;
	core::Size current_peptide_length_;
	core::pose::PoseOP current_chain_pair_pose_;
	bool is_chain_pair_new_;
	core::scoring::ScoreFunctionCOP scorefxn_;

}; // PeptideDeriverPoseOutputter

} //protocols
} //peptide_deriver



#endif //INCLUDED_protocols_peptide_deriver_PeptideDeriverPoseOutputter_hh





