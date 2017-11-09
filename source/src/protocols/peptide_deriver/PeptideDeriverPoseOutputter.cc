// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/peptide_deriver/PeptideDeriverPoseOutputter.cc
/// @brief outputs poses at different points of the Peptiderive protocol, according to the given set of options
/// @author orlypolo (orlymarcu@gmail.com)

#include <protocols/peptide_deriver/PeptideDeriverPoseOutputter.hh>
#include <protocols/peptide_deriver/PeptideDeriverFilter.hh>
#include <basic/Tracer.hh>

// Project headers
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <protocols/jd2/util.hh>

// Utility headers
#include <utility/io/ozstream.hh>
#include <utility/io/ocstream.hh>

// C++ headers
#include <string>

// External headers
#include <boost/format.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.peptide_deriver.PeptideDeriverPoseOutputter" );


namespace protocols {
namespace peptide_deriver {



PeptideDeriverPoseOutputter::PeptideDeriverPoseOutputter( bool const is_dump_best_peptide_pose,
	bool const is_dump_prepared_pose, bool const is_dump_cyclic_poses,
	core::scoring::ScoreFunctionCOP scorefxn) :
	is_dump_best_peptide_pose_(is_dump_best_peptide_pose),
	is_dump_prepared_pose_(is_dump_prepared_pose),
	is_dump_cyclic_poses_(is_dump_cyclic_poses),
	scorefxn_(std::move(scorefxn)) { }

PeptideDeriverPoseOutputter::~PeptideDeriverPoseOutputter(){}

PeptideDeriverPoseOutputter::PeptideDeriverPoseOutputter( PeptideDeriverPoseOutputter const & src) {

	is_dump_best_peptide_pose_=src.is_dump_best_peptide_pose_;
	is_dump_prepared_pose_=src.is_dump_prepared_pose_;
	is_dump_cyclic_poses_=src.is_dump_cyclic_poses_;
	scorefxn_=src.scorefxn_;

}

PeptideDeriverPoseOutputterOP
PeptideDeriverPoseOutputter::clone() const {
	return PeptideDeriverPoseOutputterOP( new PeptideDeriverPoseOutputter( *this ) );
}

void PeptideDeriverPoseOutputter::output_pose( core::pose::Pose & pose, std::string const & pose_name ) {
	(*scorefxn_)(pose);
	if ( protocols::jd2::jd2_used() ) {
		protocols::jd2::output_intermediate_pose(pose, pose_name );
	} else {
		std::string const file_name(pose_name + ".pdb");
		utility::io::ozstream out_stream( file_name );
		core::io::pdb::dump_pdb( pose, out_stream );
	}
}


void PeptideDeriverPoseOutputter::chain_pair_pose_prepared(core::pose::Pose const & pose) {
	// note: we use the term 'prepared' and 'minimized' interchangeably
	//       anticipating that this step might involve other things
	//       in the meantime, it makes sense to talk about the complex just before
	//       peptiderive goes over it

	// save a copy of the prepared pose aside, but we output only on begin_receptor_partner_pair()
	// since we want to include chain letter in output, and we don't get them from here
	current_chain_pair_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
	is_chain_pair_new_ = true;
}

void PeptideDeriverPoseOutputter::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const,
	std::string const &) {
	current_receptor_chain_letter_ = receptor_chain_letter;
	current_partner_chain_letter_ = partner_chain_letter;

	// NOTE : this assumes the PeptideDeriver uses the same prepared
	//        structure when switching roles (receptor and partner)
	//        between partners in a chain pair, in case both roles are
	//        indeed evaluated.
	if ( is_dump_prepared_pose_ && is_chain_pair_new_ ) {
		std::string pose_name = ( boost::format("%1%%2%")
			% ( (current_receptor_chain_letter_ < current_partner_chain_letter_ )? current_receptor_chain_letter_ : current_partner_chain_letter_  )
			% ( (current_receptor_chain_letter_ < current_partner_chain_letter_ )? current_partner_chain_letter_  : current_receptor_chain_letter_ ) ).str();

		output_pose( *current_chain_pair_pose_, pose_name );
		is_chain_pair_new_ = false;
	}
}

void PeptideDeriverPoseOutputter::peptide_length(core::Size const pep_length) {
	current_peptide_length_ = pep_length;
}

void PeptideDeriverPoseOutputter::peptide_entry(PeptideDeriverEntryType const entry_type, core::Real const total_isc, DerivedPeptideEntryCOP entry) {

	TR.Debug << "structure total_isc: " << total_isc << std::endl;

	// handles best linear peptide (first if statement) and all possible cyclic peptides generated for this peptide
	if ( is_dump_best_peptide_pose_ && entry_type == ET_BEST_LINEAR ) {
		std::string const lin_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_linear_linear_peptide")
			% current_receptor_chain_letter_
			% current_partner_chain_letter_
			% current_peptide_length_).str();
		core::pose::Pose temp_pose(*(entry->lin_pose));
		output_pose(temp_pose, lin_pose_name);

		for ( core::Size method = 0; method<NUM_CYCLIZATION_METHODS; ++method ) {
			if ( entry->cyc_info_set[method]->was_cyclic_model_created ) {
				std::string const cyc_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_linear_cyclic_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_).str();
				core::pose::Pose temp_cyclic_pose(*(entry->cyc_info_set[method]->cyc_pose));
				output_pose(temp_cyclic_pose, cyc_pose_name);
			}
		}
	} else {
		for ( core::Size method = 0; method<NUM_CYCLIZATION_METHODS; ++method ) {
			// handles best cyclic peptide (by any one of the cyclization methods) and its linear / other possible cyclization variants
			if ( is_dump_best_peptide_pose_ && entry_type == ET_BEST_CYCLIC && entry->cyc_info_set[method]->was_cyclic_model_created ) {
				std::string const lin_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_cyclic_linear_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_).str();
				core::pose::Pose temp_pose(*(entry->cyc_info_set[method]->pre_cyc_pose));
				output_pose(temp_pose, lin_pose_name);

				std::string const cyc_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_cyclic_cyclic_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_).str();
				core::pose::Pose temp_cyclic_pose(*(entry->cyc_info_set[method]->cyc_pose));
				output_pose(temp_cyclic_pose, cyc_pose_name);
			} else if ( is_dump_cyclic_poses_ && entry_type == ET_GENERAL && entry->cyc_info_set[method]->was_cyclic_model_created ) {
				// handles all non-best cyclic peptides that did meet the criteria to be generated
				std::string const pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_%4%_cyclic_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_
					% entry->cyc_info_set[method]->cyc_comment ).str();
				core::pose::Pose temp_pose(*(entry->cyc_info_set[method]->cyc_pose));
				output_pose(temp_pose, pose_name);
			}

		}
	}

}



} //protocols
} //peptide_deriver






