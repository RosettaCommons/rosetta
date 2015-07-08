// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/InterbondAngleScore.cc
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/InterbondAngleScore.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

// mini headers
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/func/MinMultiHarmonicFunc.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/io/izstream.hh>
#include <numeric/NumericTraits.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <basic/prof.hh>

#include <protocols/frag_picker/scores/FourAtomsConstraintData.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

static thread_local basic::Tracer trInterbondAngleScore(
		"fragment.picking.scores.InterbondAngleScore");

InterbondAngleScore::InterbondAngleScore(Size priority,
		Real lowest_acceptable_value, bool use_lowest, std::string constraints_file_name,
		Size query_size, utility::vector1<std::string> constrainable_atoms) :
	AtomBasedConstraintsScore(priority, lowest_acceptable_value, use_lowest, query_size,
			constrainable_atoms, "InterbondAngleScore") {

	utility::vector1<Real> x0(2);
	utility::vector1<Real> sd(2);
	using core::scoring::func::FuncOP;
	factory_.add_type("MINMULTIHARMONIC", FuncOP( new core::scoring::func::MinMultiHarmonicFunc(
			x0,sd) ));

	data_.resize(get_query_size());
	read_constraints(constraints_file_name);
}

InterbondAngleScore::InterbondAngleScore(Size priority,
		Real lowest_acceptable_value, bool use_lowest, std::string constraints_file_name,
		Size query_size) :
	AtomBasedConstraintsScore(priority, lowest_acceptable_value, use_lowest, query_size,
			"InterbondAngleScore") {

	utility::vector1<Real> x0(2);
	utility::vector1<Real> sd(2);
	using core::scoring::func::FuncOP;
	factory_.add_type("MINMULTIHARMONIC", FuncOP( new core::scoring::func::MinMultiHarmonicFunc(
			x0,sd) ));

	data_.resize(get_query_size());
	read_constraints(constraints_file_name);
}

bool InterbondAngleScore::cached_score(FragmentCandidateOP fragment,
		FragmentScoreMapOP scores) {

	PROF_START( basic::FRAGMENTPICKING_DIHEDRALCONSTR_SCORE );
	Size frag_len = fragment->get_length();
	Size vi = fragment->get_first_index_in_vall();
	Size qi = fragment->get_first_index_in_query();

	Real total_score = 0;
	for (Size i = 0; i < frag_len; ++i) {
		for (Size c = 1; c <= data_[i + qi].size(); ++c) {
			Size firstQueryResidueIndex = qi + i;
			Size firstVallResidueIndex = vi + i;
			FourAtomsConstraintDataOP r = data_[firstQueryResidueIndex][c];
			if (r->get_second_offset() >= frag_len - i)
				continue;
			if (r->get_third_offset() >= frag_len - i)
				continue;
			if (r->get_fourth_offset() >= frag_len - i)
				continue;
			Size secondVallResidueIndex = firstVallResidueIndex
					+ r->get_second_offset();
			Size thirdVallResidueIndex = firstVallResidueIndex
					+ r->get_third_offset();
			Size fourthVallResidueIndex = firstVallResidueIndex
					+ r->get_fourth_offset();
			if (!has_atom(firstVallResidueIndex, r->get_first_atom()))
				continue;
			if (!has_atom(secondVallResidueIndex, r->get_second_atom()))
				continue;
			if (!has_atom(thirdVallResidueIndex, r->get_third_atom()))
				continue;
			if (!has_atom(fourthVallResidueIndex, r->get_fourth_atom()))
				continue;
			numeric::xyzVector<Real> v1(
				get_atom_coordinates(secondVallResidueIndex, r->get_second_atom()) -
				get_atom_coordinates(firstVallResidueIndex, r->get_first_atom())
			);
			numeric::xyzVector<Real> v2(
				get_atom_coordinates(fourthVallResidueIndex, r->get_fourth_atom()) -
				get_atom_coordinates(thirdVallResidueIndex, r->get_third_atom())
			);

			double angle = angle_of(v1, v2)*180.0/numeric::NumericTraits<Real>::pi();
			total_score += r->get_function()->func(angle);
		}
	}
	total_score /= (Real) frag_len;
	scores->set_score_component(total_score, id_);
	PROF_STOP( basic::FRAGMENTPICKING_DIHEDRALCONSTR_SCORE );

	if ((total_score > lowest_acceptable_value_) && (use_lowest_ == true)) {
		trInterbondAngleScore.Debug << "Trashing a fragment: "
				<< *fragment << " because its score is: " << total_score
				<< std::endl;
		return false;
	}

	return true;
}

void InterbondAngleScore::read_constraints(
		std::string constraints_file_name) {
	utility::io::izstream data(constraints_file_name.c_str());
	trInterbondAngleScore.Info << "read constraints from "
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
			trInterbondAngleScore.Debug << constraints_file_name
					<< " end of file reached" << std::endl;
			break;
		}
		if (tag == "InterbondAngleScore") {
			Size res1, res2, res3, res4;
			std::string name1, name2, name3, name4;
			std::string func_type;
			//std::string type;

			data >> name1 >> res1 >> name2 >> res2 >> name3 >> res3 >> name4
					>> res4 >> func_type;

			trInterbondAngleScore.Debug << "read: " << name1 << " "
					<< name2 << " " << name3 << " " << name4 << " " << res1
					<< " " << res2 << " " << res3 << " " << res4 << " func: "
					<< func_type << std::endl;
			core::scoring::func::FuncOP func = factory_.new_func(
					func_type);
			func->read_data(data);
			std::map<std::string, Size> constr_atoms =
					get_constrainable_atoms_map();
			std::map<std::string, Size>::iterator it = constr_atoms.find(name1);
			if (it == constr_atoms.end()) {
				trInterbondAngleScore.Warning << "Unknown atom: " << name1
						<< "\nThe following constraint will NOT be used:\n"
						<< line << std::endl;
				continue;
			}
			Size a1 = it->second;
			it = constr_atoms.find(name2);
			if (it == constr_atoms.end()) {
				trInterbondAngleScore.Warning << "Unknown  atom: "
						<< name2
						<< "\nThe following constraint will NOT be used:\n"
						<< line << std::endl;
				continue;
			}
			Size a2 = it->second;
			it = constr_atoms.find(name3);
			if (it == constr_atoms.end()) {
				trInterbondAngleScore.Warning << "Unknown atom: " << name3
						<< "\nThe following constraint will NOT be used:\n"
						<< line << std::endl;
				continue;
			}
			Size a3 = it->second;
			it = constr_atoms.find(name4);
			if (it == constr_atoms.end()) {
				trInterbondAngleScore.Warning << "Unknown  atom: "
						<< name4
						<< "\nThe following constraint will NOT be used:\n"
						<< line << std::endl;
				continue;
			}
			Size a4 = it->second;
			Size o2 = res2 - res1;
			Size o3 = res3 - res1;
			Size o4 = res4 - res1;
			if ((res2 < res1) || (res3 < res1) || (res4 < res1)) {
				trInterbondAngleScore.Warning
						<< "The residue of the first constrained atoms must precede all the other three.\n\t\t"
						<< "Check residue indexes && redefine the constraint if necessary.\n\t\t"
						<< "The following constraint will NOT be used:\n"
						<< line << std::endl;
				continue;
			}

			FourAtomsConstraintDataOP dat;
			dat = FourAtomsConstraintDataOP( new FourAtomsConstraintData(func, a1, o2, a2, o3, a3, o4, a4) );
			if (res1 > data_.size()) {
				trInterbondAngleScore.Warning
						<< "Skipping a constraint that involves residue "
						<< res1 << " that does not exist in a query"
						<< std::endl;
				continue;
			}
			data_[res1].push_back(dat);
			n_constr++;
		}
	}
	trInterbondAngleScore << n_constr << " constraints loaded from a file"
			<< std::endl;
}

FragmentScoringMethodOP MakeInterbondAngleScore::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

        // Hack for VEAN experiment scoring. In the future atom names should come from command-line
        utility::vector1<std::string> constrainable_atoms;
	constrainable_atoms.push_back("N");
	constrainable_atoms.push_back("CA");
	constrainable_atoms.push_back("H");
	constrainable_atoms.push_back("HA");

	if (option[constraints::cst_file].user()) {
		trInterbondAngleScore << "Constraints loaded from: "
				<< option[constraints::cst_file]()[1] << std::endl;

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new InterbondAngleScore(priority,
				lowest_acceptable_value, use_lowest, option[constraints::cst_file]()[1],
				picker->size_of_query(),constrainable_atoms) );
	}
	utility_exit_with_message(
			"Can't read a constraints file. Provide it with constraints::cst_file flag");

	return NULL;
}

}
} // frag_picker
} // protocols
