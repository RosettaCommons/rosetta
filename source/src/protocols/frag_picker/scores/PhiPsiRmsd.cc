// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/FragmentCrmsd.cc
/// @brief  Object that scores a fragment by root mean square deviation of Phi && Psi dihedrals
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/PhiPsiRmsd.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/PhiPsiTalosIO.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>


// utils
#include <ObjexxFCL/FArray1D.hh>
#include <basic/prof.hh>

// C++

// Boost
#include <boost/tuple/tuple.hpp>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/Pose.hh>
//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
//Auto using namespaces end


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer trPhiPsiRmsd(
	"protocols.frag_picker.scores.PhiPsiRmsd");

PhiPsiRmsd::PhiPsiRmsd(Size priority, Real lowest_acceptable_value, bool use_lowest,
	core::pose::PoseOP reference_pose) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "PhiPsiRmsd") {

	n_atoms_ = reference_pose->total_residue();
	query_phi_.redimension(n_atoms_);
	query_psi_.redimension(n_atoms_);
	existing_data_.resize(n_atoms_);
	for ( Size i = 1; i <= n_atoms_; ++i ) {
		query_phi_(i) = reference_pose->phi(i);
		query_psi_(i) = reference_pose->psi(i);
		existing_data_[i] = true;
	}
}

PhiPsiRmsd::PhiPsiRmsd(Size priority, Real lowest_acceptable_value, bool use_lowest,
	PhiPsiTalosIO& reader) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "PhiPsiRmsd") {

	n_atoms_ = reader.get_sequence().length();
	query_phi_.redimension(n_atoms_);
	query_psi_.redimension(n_atoms_);
	existing_data_.resize(n_atoms_);
	for ( Size i = 1; i <= n_atoms_; ++i ) {
		if ( !reader.has_entry(i) ) {
			existing_data_[i] = false;
			trPhiPsiRmsd.Warning << "Lack of data for position " << i
				<< std::endl;
			continue;
		}
		boost::tuple<Size, char, Real, Real, Real, Real, Real, Real, Size,
			std::string> entry = reader.get_entry(i);
		if ( (entry.get<2> () > 181) || (entry.get<2> () < -181)
				|| (entry.get<3> () > 181) || (entry.get<3> () < -181) ) {
			existing_data_[i] = false;
			trPhiPsiRmsd.Warning
				<< "Unphysical Phi/Psi observation at position " << i
				<< std::endl;
			continue;
		}
		existing_data_[i] = true;
		query_phi_(i) = entry.get<2> ();
		query_psi_(i) = entry.get<3> ();
	}
}

void PhiPsiRmsd::do_caching(VallChunkOP current_chunk) {

	chunk_phi_.redimension(current_chunk->size());
	chunk_psi_.redimension(current_chunk->size());
	for ( Size i = 1; i <= current_chunk->size(); ++i ) {
		VallResidueOP r = current_chunk->at(i);
		chunk_phi_(i) = r->phi();
		chunk_psi_(i) = r->psi();
	}
}

bool PhiPsiRmsd::score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {
	return cached_score( fragment, scores );
}

bool PhiPsiRmsd::cached_score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {

	std::string tmp = fragment->get_chunk()->chunk_key();
	if ( tmp.compare(cached_scores_id_) != 0 ) {
		do_caching(fragment->get_chunk());
		cached_scores_id_ = tmp;
	}

	PROF_START( basic::FRAGMENTPICKING_PHIPSI_SCORE );
	Size offset_q = fragment->get_first_index_in_query() - 1;
	Size offset_v = fragment->get_first_index_in_vall() - 1;
	Real score = 0.0;
	Real stmp = 0.0;
	for ( Size i = 1; i <= fragment->get_length(); ++i ) {
		if ( !existing_data_[i + offset_q] ) {
			continue;
		}
		stmp = std::abs(chunk_phi_(i + offset_v) - query_phi_(i + offset_q));
		if ( stmp > 180.0 ) stmp = std::abs(360.0 - stmp);
		score += stmp * stmp;
		stmp = std::abs(chunk_psi_(i + offset_v) - query_psi_(i + offset_q));
		if ( stmp > 180.0 ) stmp = std::abs(360.0 - stmp);
		score += stmp * stmp;
	}

	score = sqrt(score) / ((Real) fragment->get_length());
	score /= (Real) fragment->get_length();
	scores->set_score_component(score, id_);
	PROF_STOP( basic::FRAGMENTPICKING_PHIPSI_SCORE );
	if ( (score > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}

bool PhiPsiRmsd::describe_score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map, std::ostream& out) {

	Real totalScore = 0;

	out << f->get_chunk()->get_pdb_id() << " " << I(5,
		f->get_first_index_in_vall()) << " ";
	for ( Size i = 1; i <= f->get_length(); i++ ) {
		out << "   " << f->get_residue(i)->aa() << "   ";
	}
	out << std::endl << "            ";
	for ( Size i = 1; i <= f->get_length(); i++ ) {
		out << F(6, 1, chunk_phi_(f->get_first_index_in_vall() + i - 1)) << " ";
	}
	out << std::endl << "query " << I(5, f->get_first_index_in_query()) << " ";
	for ( Size i = 1; i <= f->get_length(); i++ ) {
		out << F(6, 1, query_phi_(f->get_first_index_in_query() + i - 1))
			<< " ";
	}

	Size offset_q = f->get_first_index_in_query() - 1;
	Size offset_v = f->get_first_index_in_vall() - 1;
	Real score = 0.0;
	Real stmp = 0.0;
	for ( Size i = 1; i <= f->get_length(); ++i ) {
		stmp = std::abs(chunk_phi_(i + offset_v) - query_phi_(i + offset_q));
		score += stmp * stmp;
		if ( stmp > 180.0 ) stmp = std::abs(360.0 - stmp);
		stmp = std::abs(chunk_psi_(i + offset_v) - query_psi_(i + offset_q));
		if ( stmp > 180.0 ) stmp = std::abs(360.0 - stmp);
		score += stmp * stmp;
	}

	score = sqrt(score) / ((Real) f->get_length());

	if ( (score > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		out << "\nTotal score " << F(5, 3, score) << " REJECTED" << std::endl;
	} else {
		out << "\nTotal score " << F(5, 3, score) << " ACCEPTED" << std::endl;
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (score > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

void PhiPsiRmsd::clean_up() {
}

/// @brief Creates a PhiPsiRmsd scoring method
/// @param priority - priority of the scoring method. The higher value the earlier the score
///  will be evaluated
/// @param lowest_acceptable_value - if a calculated score is higher than this value,
///  fragment will be neglected
/// @param FragmentPickerOP object - not used
/// @param line - the relevant line extracted from the scoring configuration file that defines this scoring method
///   It could look like: "PhiPsiRmsd                140     -5.0     100.0 additional_string"
///  where 140, -5.0 && 100.0 are priority, weight && treshold, respectively.
///  The additional string may be:
///  - empty: then the maker tries to create a scoring object from a TALOS file
///   trying in::file::talos_phi_psi flag. If fails, will try to use a pose from in::file::s
///  - a pdb file, pdb extension is necessary. This will create a pose && steal Phi && Psi
///  - a TALOS file with Phi/Psi prediction (tab extension is necessary)
FragmentScoringMethodOP MakePhiPsiRmsd::make(Size priority,
	Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
	, std::string input_file) {

	//std::istringstream line_stream(config_line);
	//std::string score_name;
	//Size p;
	//Real weight;
	//Real lowest;
	//std::string input_file;
	//line_stream >> score_name >> p >> weight >> lowest;
	//if (!line_stream.eof()) {
	if ( input_file != "" ) {
		//line_stream >> input_file;
		Size pos = input_file.find(".pdb");
		if ( pos != std::string::npos ) {
			core::pose::PoseOP nativePose( new core::pose::Pose );
			core::import_pose::pose_from_pdb(*nativePose, input_file);
			trPhiPsiRmsd
				<< "Reference structure for Phi,Psi scoring loaded from a PDB file: "
				<< input_file << std::endl;
			trPhiPsiRmsd.Debug << "its sequence is:\n"
				<< nativePose->sequence() << std::endl;

			return (FragmentScoringMethodOP) FragmentScoringMethodOP( new PhiPsiRmsd(priority,
				lowest_acceptable_value, use_lowest, nativePose) );
		}
		pos = input_file.find(".tab");
		if ( pos != std::string::npos ) {
			trPhiPsiRmsd
				<< "Reference file for Phi,Psi scoring loaded from a TALOS file: "
				<< input_file << std::endl;
			PhiPsiTalosIO in(input_file);
			in.write(trPhiPsiRmsd.Debug);
			return (FragmentScoringMethodOP) FragmentScoringMethodOP( new PhiPsiRmsd(priority,
				lowest_acceptable_value, use_lowest, in) );
		}

	}

	if ( option[in::file::talos_phi_psi].user() ) {
		trPhiPsiRmsd
			<< "Reference file for Phi,Psi scoring loaded from a TALOS file: "
			<< input_file << std::endl;
		PhiPsiTalosIO in(option[in::file::talos_phi_psi]());
		in.write(trPhiPsiRmsd.Debug);
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new PhiPsiRmsd(priority,
			lowest_acceptable_value, use_lowest, in) );
	}
	if ( option[in::file::native].user() ) {
		core::pose::PoseOP nativePose( new core::pose::Pose );
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::native]());
		trPhiPsiRmsd << "Reference file for Phi,Psi scoring loaded from "
			<< option[in::file::native]() << std::endl;
		trPhiPsiRmsd.Debug << "its sequence is:\n" << nativePose->sequence()
			<< std::endl;

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new PhiPsiRmsd(priority,
			lowest_acceptable_value, use_lowest, nativePose) );
	}
	if ( option[in::file::s].user() ) {
		core::pose::PoseOP nativePose( new core::pose::Pose );
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::s]()[1]);
		trPhiPsiRmsd << "Reference file for Phi,Psi scoring loaded from "
			<< option[in::file::s]()[1] << std::endl;
		trPhiPsiRmsd.Debug << "its sequence is:\n" << nativePose->sequence()
			<< std::endl;

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new PhiPsiRmsd(priority,
			lowest_acceptable_value, use_lowest, nativePose) );
	}

	utility_exit_with_message(
		"Can't read a reference Phi,Psi data. Provide a structure with in::file::s flag\n\t or a Phi-Psi prediction in TALOS format.");

	return NULL;
}

} // scores
} // frag_picker
} // protocols
