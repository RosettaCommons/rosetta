// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite && is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/AtomPairConstraintsScore.cc
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/AtomPairConstraintsScore.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

// mini headers
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/io/izstream.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

// @brief Auto-generated virtual destructor
AtomPairConstraintsData::~AtomPairConstraintsData() {}

static thread_local basic::Tracer trAtomPairConstraintsScore(
		"fragment.picking.scores.AtomPairConstraintsScore");

/// @brief Prepare an AtomPair constraint - based score that utilizes some user-defined atoms.
/// @param priority - the priority for this scoring method. The lower the priority, the later the score will be evaluated
/// Because a fragment may be discarded when a score is too low, the most accurate && meaningful scores should have the highest priority
/// @param lowest_acceptable_value - a fragment for which this score is below a certain threshold will be discarded
/// @param constraints_file_name - from this file AtomPair constraints will be obtained
/// @param query_size - the number of residues in the query sequence
/// @param constrainable_atoms - a vector of strings providing names of constrained atoms.
/// On every do_cahing() event these && only these atoms will be cached from a chunk's pose
AtomPairConstraintsScore::AtomPairConstraintsScore(Size priority,
		Real lowest_acceptable_value, bool use_lowest, std::string constraints_file_name,
		Size query_size, utility::vector1<std::string> constrainable_atoms) :
	AtomBasedConstraintsScore(priority, lowest_acceptable_value, use_lowest, query_size,
			constrainable_atoms, "AtomPairConstraintsScore") {

	data_.resize(get_query_size());
	read_constraints(constraints_file_name);
}

/// @brief Prepare an AtomPair constraint - based score that utilizes N, C, CA, O && CB atoms
/// @detailed These atoms that will be cached when a new
/// chunk is considered (i.e. at every do_caching() call)
/// @param priority - the priority for this scoring method. The lower the priority, the later the score will be evaluated
/// Because a fragment may be discarded when a score is too low, the most accurate && meaningful scores should have the highest priority
/// @param lowest_acceptable_value - a fragment for which this score is below a certain threshold will be discarded
/// @param constraints_file_name - from this file AtomPair constraints will be obtained
/// @param query_size - the number of residues in the query sequence
AtomPairConstraintsScore::AtomPairConstraintsScore(Size priority,
		Real lowest_acceptable_value, bool use_lowest, std::string constraints_file_name,
		Size query_size) :
	AtomBasedConstraintsScore(priority, lowest_acceptable_value, use_lowest, query_size,
			"AtomPairConstraintsScore") {

	data_.resize(get_query_size());
	read_constraints(constraints_file_name);
}

/// @brief Calculates a score for a given fragment.
/// @detailed Resulting value is written into a given score map.
/// @param fragment - fragment to be scored
/// @param scores - resulting score will be stored here
bool AtomPairConstraintsScore::cached_score(FragmentCandidateOP fragment,
		FragmentScoreMapOP scores) {

	PROF_START( basic::FRAGMENTPICKING_ATOMPAIR_SCORE );

	Size frag_len = fragment->get_length();
	Size vi = fragment->get_first_index_in_vall();
	Size qi = fragment->get_first_index_in_query();

	Real total_score = 0;
	for (Size i = 0; i < frag_len; ++i) {
		for (Size c = 1; c <= data_[i + qi].size(); ++c) {
			Size firstQueryResidueIndex = qi + i;
			AtomPairConstraintsDataOP r = data_[firstQueryResidueIndex][c];
			if (r->get_offset() >= frag_len - i)
				continue;
			Size firstVallResidueIndex = vi + i;
			Size secondVallResidueIndex = firstVallResidueIndex
					+ r->get_offset();
			if (!has_atom(firstVallResidueIndex, r->get_first_atom()))
				continue;
			if (!has_atom(secondVallResidueIndex, r->get_second_atom()))
				continue;

			assert (fragment->get_chunk()->size() >=secondVallResidueIndex);

			numeric::xyzVector<Real> v1 = get_atom_coordinates(
					firstVallResidueIndex, r->get_first_atom());
			numeric::xyzVector<Real> v2 = get_atom_coordinates(
					secondVallResidueIndex, r->get_second_atom());
			Real d = v1.distance(v2);
			total_score += r->get_function()->func(d);
		}
	}
	total_score /= (Real) frag_len;
	scores->set_score_component(total_score, id_);
	PROF_STOP( basic::FRAGMENTPICKING_ATOMPAIR_SCORE );
	if ((total_score > lowest_acceptable_value_) && (use_lowest_ == true)) {
		trAtomPairConstraintsScore.Debug << "Trashing a fragment: "
				<< *fragment << " because its score is: " << total_score
				<< std::endl;
		return false;
	}

	return true;
}

void AtomPairConstraintsScore::read_constraints(
		std::string constraints_file_name) {
	utility::io::izstream data(constraints_file_name.c_str());
	trAtomPairConstraintsScore.Info << "reading constraints from "
			<< constraints_file_name << std::endl;
	if (!data) {
		utility_exit_with_message("[ERROR] Unable to open constraints file: "
				+ constraints_file_name);
	}

	std::string line;
	getline(data, line); // header line
	std::string tag;
	Size n_constr = 0;
	while (!data.fail()) {
		char c = data.peek();
		if (c == '#' || c == '\n') {
			getline(data, line); //comment
			continue;
		}
		data >> tag;
		if (data.fail()) {
			trAtomPairConstraintsScore.Debug << constraints_file_name
					<< " end of file reached" << std::endl;
			break;
		}
		if (tag == "AtomPair") {
			std::string name1, name2, func_type;
			Size id1, id2;
			data >> name1 >> id1 >> name2 >> id2 >> func_type;
			trAtomPairConstraintsScore.Debug << "read: " << name1 << " " << id1
					<< " " << name2 << " " << id2 << " func: " << func_type
					<< std::endl;
			core::scoring::func::FuncOP func = factory_.new_func(
					func_type);
			func->read_data(data);
			std::map<std::string, Size> constr_atoms =
					get_constrainable_atoms_map();
			std::map<std::string, Size>::iterator it = constr_atoms.find(name1);
			if (it == constr_atoms.end()) {
				trAtomPairConstraintsScore.Warning << "Unknown backbone atom: "
						<< name1
						<< "\nThe following constraint will NOT be used:\n"
						<< line << std::endl;
				continue;
			}
			Size a1 = it->second;
			it = constr_atoms.find(name2);
			if (it == constr_atoms.end()) {
				trAtomPairConstraintsScore.Warning << "Unknown backbone atom: "
						<< name2
						<< "\nThe following constraint will NOT be used:\n"
						<< line << std::endl;
				continue;
			}
			Size a2 = it->second;
			AtomPairConstraintsDataOP dat;
			if (id2 > id1) {
				if (id2 > get_query_size()) {
					trAtomPairConstraintsScore.Warning
							<< "Constrained atom id exceeds the length of a query sequence. The following constraint will NOT be used:\n"
							<< line << std::endl;
					continue;
				}
				dat = AtomPairConstraintsDataOP( new AtomPairConstraintsData(func, id2 - id1, a1, a2) );
				data_[id1].push_back(dat);
				n_constr++;
			} else {
				if (id1 > get_query_size()) {
					trAtomPairConstraintsScore.Warning
							<< "Constrained atom id exceeds the length of a query sequence. The following constraint will NOT be used:\n"
							<< line << std::endl;
					continue;
				}
				dat = AtomPairConstraintsDataOP( new AtomPairConstraintsData(func, id1 - id2, a2, a1) );
				data_[id2].push_back(dat);
				n_constr++;
			}
		}
	}
	trAtomPairConstraintsScore << n_constr << " constraints loaded from a file"
			<< std::endl;
}

FragmentScoringMethodOP MakeAtomPairConstraintsScore::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if (option[constraints::cst_file].user()) {
		trAtomPairConstraintsScore << "Constraints loaded from: "
				<< option[constraints::cst_file]()[1] << std::endl;

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new AtomPairConstraintsScore(priority,
				lowest_acceptable_value, use_lowest, option[constraints::cst_file]()[1],
				picker->size_of_query()) );
	}
	utility_exit_with_message(
			"Can't read a constraints file. Provide it with constraints::cst_file flag");

	return NULL;
}

}
} // frag_picker
} // protocols
