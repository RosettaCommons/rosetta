// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FragmentAllAtomCrmsd.cc
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/FragmentAllAtomCrmsd.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#ifdef WIN32
#include <protocols/frag_picker/FragmentPicker.hh>
#endif

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>


#include <numeric/model_quality/rms.hh>

// utils
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <fstream>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL;
using namespace core;

static basic::Tracer trRmsScore(
	"protocols.frag_picker.scores.FragmentAllAtomCrmsd");

FragmentAllAtomCrmsd::~FragmentAllAtomCrmsd() {}

FragmentAllAtomCrmsd::FragmentAllAtomCrmsd(core::Size priority, core::Real lowest_acceptable_value,
	bool use_lowest, std::string query_sequence, core::pose::PoseOP reference_pose) :
	FragmentScoringMethod(priority, lowest_acceptable_value, use_lowest, "FragmentAllAtomCrmsd"), query_sequence_(query_sequence) {

	reference_pose_ = reference_pose;
	fragment_pose_ = pose::PoseOP( new core::pose::Pose );
	total_ref_resids_ = reference_pose_->size();
	reference_coordinates_.redimension(3, 4 * total_ref_resids_, 0.0);
	fill_coords(*reference_pose_, reference_coordinates_, total_ref_resids_, reference_pose_->sequence());
	weights_.redimension(reference_pose_->size() * 4, 1.0);

	core::pose::make_pose_from_sequence(*fragment_pose_, reference_pose_->sequence(), *(chemical::ChemicalManager::get_instance()->residue_type_set("centroid")));
}

FragmentAllAtomCrmsd::FragmentAllAtomCrmsd(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
	std::string query_sequence, utility::vector1< utility::vector1<core::Real> > xyz)  :
	FragmentScoringMethod(priority, lowest_acceptable_value, use_lowest, "FragmentAllAtomCrmsd"),
	query_sequence_(query_sequence) {

	total_ref_resids_ = xyz.size() / 5;
	fragment_pose_ = pose::PoseOP( new core::pose::Pose );
	reference_coordinates_.redimension(3, 4 * total_ref_resids_, 0.0);
	weights_.redimension(total_ref_resids_ * 4, 1.0);
	for ( core::Size i = 1; i <= total_ref_resids_ * 4; i++ ) {
		trRmsScore.Trace << i<<" ";
		for ( core::Size d = 1; d <= 3; ++d ) {
			reference_coordinates_(d, i) = xyz[i][d];
			trRmsScore.Trace <<xyz[i][d]<<" ";
		}
		trRmsScore.Trace <<std::endl;
	}
	core::pose::make_pose_from_sequence(*fragment_pose_, reference_pose_->sequence(), *(chemical::ChemicalManager::get_instance()->residue_type_set("centroid")));
}

void FragmentAllAtomCrmsd::fill_coords(core::pose::Pose const& pose,
	FArray2_double& coords, core::Size n_res, std::string ) {

	trRmsScore.Debug << "Copying coordinates from ... The first residues are: "
		<< pose.residue(1).name3() << " " << pose.residue(2).name3() << " "
		<< pose.residue(3).name3() << std::endl;

	core::Size n_at = 1;
	for ( core::Size i = 1; i <= n_res; i++ ) {
		id::NamedAtomID idN("N", i);
		PointPosition const& xyzN = pose.xyz(idN);
		for ( core::Size d = 1; d <= 3; ++d ) {
			coords(d, n_at) = xyzN[d - 1];
		}
		n_at++;

		id::NamedAtomID idCA("CA", i);
		PointPosition const& xyzCA = pose.xyz(idCA);
		for ( core::Size d = 1; d <= 3; ++d ) {
			coords(d, i) = xyzCA[d - 1];
		}
		n_at++;

		id::NamedAtomID idC("C", i);
		PointPosition const& xyzC = pose.xyz(idC);
		for ( core::Size d = 1; d <= 3; ++d ) {
			coords(d, n_at) = xyzC[d - 1];
		}
		n_at++;

		id::NamedAtomID idO("O", i);
		PointPosition const& xyzO = pose.xyz(idO);
		for ( core::Size d = 1; d <= 3; ++d ) {
			coords(d, i) = xyzO[d - 1];
		}
		n_at++;
	}
}

bool FragmentAllAtomCrmsd::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {


	if ( (core::Size) fragment_coordinates_.size2() < f->get_length() * 4 ) {
		fragment_coordinates_.redimension(3, f->get_length()*4, 0.0);
		frag_pos_coordinates_.redimension(3, f->get_length()*4, 0.0);
	}
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		fragment_pose_->set_phi( i, f->get_residue(i)->phi() );
		fragment_pose_->set_psi( i, f->get_residue(i)->psi() );
		fragment_pose_->set_omega( i, f->get_residue(i)->omega() );

		for ( core::Size d = 1; d <= 3; ++d ) {
			//         std::cerr << i << " " << reference_coordinates_.size2() << " " << frag_pos_coordinates_.size2() << " " << f->get_first_index_in_query()<< std::endl;
			frag_pos_coordinates_(d,i) = reference_coordinates_(d,i + f->get_first_index_in_query() - 1);
		}
	}

	fill_coords(*fragment_pose_, fragment_coordinates_, f->get_length(), f->sequence());

	core::Real rms = numeric::model_quality::rms_wrapper(f->get_length(),
		fragment_coordinates_, frag_pos_coordinates_);

	empty_map->set_score_component(rms, id_);

	if ( (rms > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}


FragmentScoringMethodOP MakeFragmentAllAtomCrmsd::make(core::Size priority,
	core::Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
	, std::string) {

	if ( option[in::file::native].user() ) {
		trRmsScore
			<< "Reference structure to score fragments by crmsd loaded from: "
			<< option[in::file::native]() << std::endl;
		core::pose::PoseOP nativePose( new core::pose::Pose );
		core::import_pose::pose_from_file(*nativePose, option[in::file::native](), core::import_pose::PDB_file);

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new FragmentAllAtomCrmsd(priority,
			lowest_acceptable_value, use_lowest, nativePose->sequence(), nativePose) );
	}
	if ( option[in::file::s].user() ) {
		trRmsScore
			<< "Reference structure to score fragments by crmsd loaded from: "
			<< option[in::file::s]()[1] << std::endl;
		core::pose::PoseOP nativePose( new core::pose::Pose );
		core::import_pose::pose_from_file(*nativePose, option[in::file::s]()[1], core::import_pose::PDB_file);

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new FragmentAllAtomCrmsd(priority,
			lowest_acceptable_value, use_lowest, nativePose->sequence(), nativePose) );
	}
	if ( option[in::file::xyz].user() ) {

		trRmsScore
			<< "Reference structure to score fragments by crmsd loaded from: "
			<< option[in::file::xyz]().c_str() << std::endl;
		std::string line;
		utility::vector1< utility::vector1< core::Real > > xyz;
		// std::istream & input = utility::io::izstream( option[ in::file::xyz ]().c_str() );
		std::ifstream input( option[ in::file::xyz ]().c_str() );

		std::string seq;
		while ( getline( input, line ) ) {
			trRmsScore.Warning << line<<std::endl;
			if ( line.substr(0,1) == "#" ) continue;
			std::istringstream line_stream( line );
			utility::vector1<core::Real> row;
			core::Real x,y,z;
			char c;
			line_stream >> c >> x >> y >> z;
			seq += c;
			row.push_back(x);
			row.push_back(y);
			row.push_back(z);
			xyz.push_back( row );
		}
		trRmsScore <<  xyz.size() << " atoms found in the reference" << std::endl;

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new FragmentAllAtomCrmsd(priority,
			lowest_acceptable_value, use_lowest, seq, xyz) );
	}
	utility_exit_with_message(
		"Can't read a reference structure. Provide it with in::file::s flag");

	return NULL;
}

} // scores
} // frag_picker
} // protocols
