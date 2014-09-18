// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite && is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions && developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/frag_picker/scores/FragmentCrmsd.cc
/// @brief  Object that scores a fragment by root mean square deviation of Phi && Psi dihedrals
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/PhiPsiSquareWell.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
// AUTO-REMOVED #include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#ifdef WIN32
#include <protocols/frag_picker/FragmentPicker.hh>
#endif

#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// utils
#include <ObjexxFCL/FArray1D.hh>
#include <basic/prof.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/Pose.hh>


// C++

// Boost

namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer tr(
		"protocols.frag_picker.scores.PhiPsiSquareWell");

	PhiPsiSquareWell::PhiPsiSquareWell(Size priority, Real lowest_acceptable_value, bool use_lowest,
		core::pose::PoseOP reference_pose) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "PhiPsiSquareWell") {

	Size n_atoms_ = reference_pose->total_residue();
	query_phi_.redimension(n_atoms_);
	query_psi_.redimension(n_atoms_);
	query_d_phi_.redimension(n_atoms_);
	query_d_psi_.redimension(n_atoms_);
	existing_data_.resize(n_atoms_);
	for (Size i = 1; i <= n_atoms_; ++i) {
		query_phi_(i) = reference_pose->phi(i);
		query_psi_(i) = reference_pose->psi(i);
		query_d_phi_(i) = 0.0;
		query_d_psi_(i) = 0.0;
		existing_data_[i] = true;
	}
}

	PhiPsiSquareWell::PhiPsiSquareWell(Size priority, Real lowest_acceptable_value, bool use_lowest,
		PhiPsiTalosIO& reader) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "PhiPsiSquareWell") {

	Size n_atoms_ = reader.get_sequence().length();
	query_phi_.redimension(n_atoms_);
	query_psi_.redimension(n_atoms_);
	query_d_phi_.redimension(n_atoms_);
	query_d_psi_.redimension(n_atoms_);
	query_dist_.redimension(n_atoms_);
	query_s2_.redimension(n_atoms_);
	query_cnt_.resize(n_atoms_);
	query_class_.resize(n_atoms_);
	existing_data_.resize(n_atoms_);
	for (Size i = 1; i <= n_atoms_; ++i) {

		if (!reader.has_entry(i)) {
			existing_data_[i] = false;
			tr.Warning << "Lack of data for position " << i
					<< std::endl;
			continue;
		}

		if ( (reader.quality(i) == "Warn")
				 || (reader.quality(i) == "Dyn")
				 || (reader.quality(i) == "None") ) {
			existing_data_[i] = false;
			tr.Info	<< "Untrustworthy data at position " << i << " due to prediction class " << reader.quality(i)
							<< std::endl;
			continue;
		}

		//boost::tuple<Size, char, Real, Real, Real, Real, Real, Real, Size,
		//		std::string> entry = reader.get_entry(i);
		if ((reader.phi(i) > 181) || (reader.phi(i) < -181)
				|| (reader.psi(i) > 181) || (reader.psi(i) < -181)) {
			existing_data_[i] = false;
			tr.Warning
					<< "Unphysical Phi/Psi observation at position " << i
					<< std::endl;
			continue;
		}

		existing_data_[i] = true;
		query_phi_(i) = reader.phi(i);
		query_psi_(i) = reader.psi(i);
		query_d_phi_(i) = reader.d_phi(i);
		query_d_psi_(i) = reader.d_psi(i);
	}
}

void PhiPsiSquareWell::do_caching(VallChunkOP current_chunk) {

	chunk_phi_.redimension(current_chunk->size());
	chunk_psi_.redimension(current_chunk->size());
	for (Size i = 1; i <= current_chunk->size(); ++i) {
		VallResidueOP r = current_chunk->at(i);
		chunk_phi_(i) = r->phi();
		chunk_psi_(i) = r->psi();
	}
}

bool PhiPsiSquareWell::cached_score(FragmentCandidateOP fragment,
		FragmentScoreMapOP scores) {

	return score(fragment, scores);

// 	std::string ctmp = fragment->get_chunk()->chunk_key();
// 	if (ctmp.compare(cached_scores_id_) != 0) {
// 		do_caching(fragment->get_chunk());
// 		cached_scores_id_ = ctmp;
// 	}

// 	PROF_START( basic::FRAGMENTPICKING_PHIPSI_SCORE );
// 	Size offset_q = fragment->get_first_index_in_query() - 1;
// 	Size offset_v = fragment->get_first_index_in_vall() - 1;
// 	Real score = 0.0;
// 	Real tmp = 0.0;
// 	for (Size i = 1; i < fragment->get_length(); ++i) {
// 		if (!existing_data_[i + offset_q])
// 			continue;
// 		Real d = 0.0;

// 		tmp = std::abs(chunk_phi_(i + offset_v) - query_phi_(i + offset_q));

// 		if ( tmp > 180.0 ) tmp = std::abs(tmp - 360.0);

// 		if ( tmp > query_d_phi_(i + offset_q) ) {
// 			tmp = tmp - query_d_phi_(i + offset_q);
// 		} else {
// 			tmp = 0.0;
// 		}

// 		d += tmp * tmp;

// 		tmp = std::abs(chunk_psi_(i + offset_v) - query_psi_(i + offset_q));

// 		if ( tmp > 180.0 ) tmp = std::abs(tmp - 360.0);

// 		if ( tmp > query_d_psi_(i + offset_q) ) {
// 			tmp = tmp - query_d_psi_(i + offset_q);
// 		} else {
// 			tmp = 0.0;
// 		}

// 		d += tmp * tmp;

// 		score += std::sqrt(d);
// 	}

// 	score = score / ((Real) fragment->get_length());
// 	PROF_STOP( basic::FRAGMENTPICKING_PHIPSI_SCORE );

// 	scores->set_score_component(score, id_);
// 	if ((score > lowest_acceptable_value_) && (use_lowest_ == true))
// 		return false;

// 	return true;
}


bool PhiPsiSquareWell::score(FragmentCandidateOP fragment,
		FragmentScoreMapOP scores) {

	PROF_START( basic::FRAGMENTPICKING_PHIPSI_SCORE );
	Size offset_q = fragment->get_first_index_in_query() - 1;
	Size offset_v = fragment->get_first_index_in_vall() - 1;
	VallChunkOP chunk = fragment->get_chunk();
	Real score = 0.0;
	Real tmp = 0.0;

	utility::vector1< Real > values;
	values.resize(fragment->get_length());

	for (Size i = 1; i < fragment->get_length(); ++i) {
		Real d = 0.0;
		VallResidueOP r = chunk->at(i + offset_v);

		values[i] = 0;
		if (!existing_data_[i + offset_q])
 			continue;

		tmp = std::abs(r->phi() - query_phi_(i + offset_q));

		if ( tmp > 180.0 ) tmp = std::abs(tmp - 360.0);

		if ( tmp > query_d_phi_(i + offset_q) ) {
			tmp = tmp - query_d_phi_(i + offset_q);
		} else {
			tmp = 0.0;
		}

		if ( tmp > 0.0 ) {
			d += 1 / ( 1 + exp( -2.5*( tmp / query_d_phi_(i + offset_q) ) + 5) );
		}

		tmp = std::abs(r->psi() - query_psi_(i + offset_q));

		if ( tmp > 180.0 ) tmp = std::abs(tmp - 360.0);

		if ( tmp > query_d_psi_(i + offset_q) ) {
			tmp = tmp - query_d_psi_(i + offset_q);
		} else {
			tmp = 0.0;
		}

		if ( tmp > 0.0 ) {
			d += 1 / ( 1 + exp( -2.5*( tmp / query_d_psi_(i + offset_q) ) + 5) );
		}

		values[i] = std::sqrt(d);
	}


	std::sort( values.begin(), values.end() );

	score = 0.0;
	for (Size i = 1; i <= fragment->get_length(); i++) {
		//~1 at 1, ~0.05 at f->get_length, 0.5 at 0.7*f->get_length()
		Real sigmoid_weight( 1 / ( 1 + exp( (10*( (Real) i ) / fragment->get_length()) - 7 ) ) );

		score += sigmoid_weight*values[i];
	}


	score = score / ((Real) fragment->get_length());
	PROF_STOP( basic::FRAGMENTPICKING_PHIPSI_SCORE );

	scores->set_score_component(score, id_);
	if ((score > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;

	return true;
}

void PhiPsiSquareWell::clean_up() {
}

/// @brief Creates a PhiPsiSquareWell scoring method
/// @param priority - priority of the scoring method. The higher value the earlier the score
///		will be evaluated
/// @param lowest_acceptable_value - if a calculated score is higher than this value,
///		fragment will be neglected
/// @param FragmentPickerOP object - not used
/// @param extras - additional parameters to create a new object. Allowed values are:
///		- empty: then the maker tries to create a scoring object from a TALOS file
///			trying in::file::talos_phi_psi flag. If fails, will try to use a pose from in::file::s
///		- a pdb file, pdb extension is necessary. This will create a pose && steal Phi && Psi
///		- a TALOS file with Phi/Psi prediction (tab extension is necessary)
FragmentScoringMethodOP MakePhiPsiSquareWell::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
		, std::string input_file) {

	if (input_file != "") {
		Size pos = input_file.find(".pdb");
		if (pos != std::string::npos) {
			core::pose::PoseOP nativePose = new core::pose::Pose;
			core::import_pose::pose_from_pdb(*nativePose, input_file);
			tr.Info
					<< "Reference file for Phi,Psi scoring loaded from "
					<< input_file << std::endl;
			tr.Debug << "its sequence is:\n"
					<< nativePose->sequence() << std::endl;

			return (FragmentScoringMethodOP) new PhiPsiSquareWell(priority,
					lowest_acceptable_value, use_lowest, nativePose);
		}
		pos = input_file.find(".tab");
		if (pos != std::string::npos) {
			tr.Info
					<< "Reference file for Phi,Psi scoring loaded from a TALOS file: "
					<< input_file << std::endl;
			PhiPsiTalosIO in(input_file);
			in.write(tr.Debug);
			return (FragmentScoringMethodOP) new PhiPsiSquareWell(priority,
					lowest_acceptable_value, use_lowest, in);
		}
	}

	if (option[in::file::talos_phi_psi].user()) {
		tr.Info
				<< "Reference file for Phi,Psi scoring loaded from a TALOS file: "
				<< input_file << std::endl;
		PhiPsiTalosIO in(option[in::file::talos_phi_psi]());
		in.write(tr.Debug);
		return (FragmentScoringMethodOP) new PhiPsiSquareWell(priority,
				lowest_acceptable_value, use_lowest, in);
	}
	if (option[in::file::s].user()) {
		core::pose::PoseOP nativePose = new core::pose::Pose;
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::s]()[1]);
		tr.Info << "Reference file for Phi,Psi scoring loaded from "
				<< option[in::file::s]()[1] << std::endl;
		tr.Debug << "its sequence is:\n"
				<< nativePose->sequence() << std::endl;

		return (FragmentScoringMethodOP) new PhiPsiSquareWell(priority,
				lowest_acceptable_value, use_lowest, nativePose);
	}
	if (option[in::file::native].user()) {
		core::pose::PoseOP nativePose = new core::pose::Pose;
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::native]());
		tr.Info << "Reference file for Phi,Psi scoring loaded from "
				<< option[in::file::native]() << std::endl;
		tr.Debug << "its sequence is:\n"
				<< nativePose->sequence() << std::endl;

		return (FragmentScoringMethodOP) new PhiPsiSquareWell(priority,
				lowest_acceptable_value, use_lowest, nativePose);
	}

	utility_exit_with_message(
			"Can't read a reference Phi,Psi data. Provide a structure with in::file::s flag\n\t or a Phi-Psi prediction in TALOS format.");

	return NULL;
}

} // scores
} // frag_picker
} // protocols
