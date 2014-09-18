// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite && is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions && developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/frag_picker/scores/FragmentCrmsdResDepth.cc
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author David E Kim

#include <protocols/frag_picker/scores/FragmentCrmsdResDepth.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>

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
#include <basic/prof.hh>

//Auto Headers
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <string>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer trRmsScore(
		"protocols.frag_picker.scores.FragmentCrmsdResDepth");

FragmentCrmsdResDepth::~FragmentCrmsdResDepth() {}

FragmentCrmsdResDepth::FragmentCrmsdResDepth(Size priority, Real lowest_acceptable_value,
		bool use_lowest, core::pose::PoseOP reference_pose, utility::vector1<core::Real> & query_residue_depth) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "FragmentCrmsdResDepth") {
	reference_pose_ = reference_pose;
	n_atoms_ = reference_pose_->total_residue();
	reference_coordinates_.redimension(3, n_atoms_, 0.0);
	fill_CA_coords(*reference_pose_, reference_coordinates_, n_atoms_);
	weights_.redimension(reference_pose_->total_residue(), 1.0);
	query_residue_depth_ = query_residue_depth;
}

FragmentCrmsdResDepth::FragmentCrmsdResDepth(Size priority, Real lowest_acceptable_value, bool use_lowest,
	utility::vector1< utility::vector1<Real> > xyz, utility::vector1<core::Real> & query_residue_depth)  :
        CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "FragmentCrmsdResDepth") {

	n_atoms_ = xyz.size();
	reference_coordinates_.redimension(3, n_atoms_, 0.0);
	weights_.redimension(n_atoms_, 1.0);
	for (core::Size i = 1; i <= n_atoms_; i++) {
		trRmsScore.Debug << i<<" ";
		for (core::Size d = 1; d <= 3; ++d) {
			reference_coordinates_(d, i) = xyz[i][d];
			trRmsScore.Debug <<xyz[i][d]<<" ";
		}
		trRmsScore.Debug <<std::endl;
	}
	query_residue_depth_ = query_residue_depth;
}

void FragmentCrmsdResDepth::fill_CA_coords(core::pose::Pose const& pose,
		FArray2_double& coords, Size n_atoms) {

	trRmsScore.Debug << "Copying coordinates from ... The first residues are: "
			<< pose.residue(1).name3() << " " << pose.residue(2).name3() << " "
			<< pose.residue(3).name3() << std::endl;

	for (core::Size i = 1; i <= n_atoms; i++) {
		id::NamedAtomID idCA("CA", i);
		PointPosition const& xyz = pose.xyz(idCA);
		for (core::Size d = 1; d <= 3; ++d) {
			coords(d, i) = xyz[d - 1];
		}
	}
}

bool FragmentCrmsdResDepth::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	PROF_START( basic::CA_RMSD_EVALUATION );

	Real depth_score = 0;

	if ((Size) fragment_coordinates_.size2() < f->get_length()) {
		fragment_coordinates_.redimension(3, f->get_length(), 0.0);
		frag_pos_coordinates_.redimension(3, f->get_length(), 0.0);
	}
	for (Size i = 1; i <= f->get_length(); i++) {
		Size qindex = i + f->get_first_index_in_query() - 1;
		VallChunkOP chunk = f->get_chunk();
		VallResidueOP r = chunk->at( f->get_first_index_in_vall() + i - 1 );
		fragment_coordinates_(1, i) = r->x();
		frag_pos_coordinates_(1, i) = reference_coordinates_(1, qindex);
		fragment_coordinates_(2, i) = r->y();
		frag_pos_coordinates_(2, i) = reference_coordinates_(2, qindex);
		fragment_coordinates_(3, i) = r->z();
		frag_pos_coordinates_(3, i) = reference_coordinates_(3, qindex);
		// residue depth similarity score from Zhou and Zhou, Proteins 2005
		// (exp( - Dquery / 2.8) - exp( - Dtemplate / 2.8))**2
		Real depthdiff = exp(-1*query_residue_depth_[qindex]/2.8) - exp(-1*r->depth()/2.8);
		depth_score += depthdiff*depthdiff;
	}
	Real rms = numeric::model_quality::rms_wrapper(f->get_length(),
			fragment_coordinates_, frag_pos_coordinates_);

	// similarity score from Zhou and Zhou, Proteins 2005
	// rmsd**2 + weight*depth_score  with weight = 10
	Real totalScore = rms*rms + 10*depth_score;

	empty_map->set_score_component(totalScore, id_);
	PROF_STOP( basic::CA_RMSD_EVALUATION );
	if ((rms > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;

	return true;
}


void FragmentCrmsdResDepth::do_caching(VallChunkOP current_chunk) {

	chunk_coordinates_.redimension(3, current_chunk->size());
	for (core::Size i = 1; i <= current_chunk->size(); i++) {
		VallResidueOP r = current_chunk->at(i);
		chunk_coordinates_(1, i) = r->x();
		chunk_coordinates_(2, i) = r->y();
		chunk_coordinates_(3, i) = r->z();
	}
}

bool FragmentCrmsdResDepth::cached_score(FragmentCandidateOP fragment,
		FragmentScoreMapOP scores) {

	std::string tmp = fragment->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0)
		do_caching(fragment->get_chunk());

	Real depth_score = 0;

	PROF_START( basic::CA_RMSD_EVALUATION );
	Size frag_len = fragment->get_length();
	if ((Size) fragment_coordinates_.size2() < frag_len)
		fragment_coordinates_.redimension(3, frag_len, 0.0);
	if ((Size) frag_pos_coordinates_.size2() < frag_len)
		frag_pos_coordinates_.redimension(3, frag_len, 0.0);

	for (core::Size i = 1; i <= frag_len; ++i) {
		Size qindex = i + fragment->get_first_index_in_query() - 1;
		Size vindex = i + fragment->get_first_index_in_vall() - 1;
		for (core::Size d = 1; d <= 3; ++d) {
			fragment_coordinates_(d, i) = chunk_coordinates_(d, vindex);
			frag_pos_coordinates_(d, i) = reference_coordinates_(d, qindex);
		}
		// residue depth similarity score from Zhou and Zhou, Proteins 2005
		// (exp( - Dquery / 2.8) - exp( - Dtemplate / 2.8))**2
		Real depthdiff = exp(-1*query_residue_depth_[qindex]/2.8) -
				exp(-1*fragment->get_chunk()->at( vindex )->depth()/2.8);
		depth_score += depthdiff*depthdiff;
	}

	Real rms = numeric::model_quality::rms_wrapper(frag_len,
			fragment_coordinates_, frag_pos_coordinates_);

	// similarity score from Zhou and Zhou, Proteins 2005
	// rmsd**2 + weight*depth_score  with weight = 10
	Real totalScore = rms*rms + 10*depth_score;

	scores->set_score_component(totalScore, id_);
	PROF_STOP( basic::CA_RMSD_EVALUATION );
	if ((rms > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;

	return true;
}

void FragmentCrmsdResDepth::clean_up() {
}

FragmentScoringMethodOP MakeFragmentCrmsdResDepth::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string) {

	if (option[in::file::native].user()) {
		trRmsScore
				<< "Reference structure to score fragments by crmsd loaded from: "
				<< option[in::file::native]() << std::endl;
		core::pose::PoseOP nativePose = new core::pose::Pose;
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::native]());

		if (nativePose->total_residue() != picker->get_query_residue_depth().size())
			utility_exit_with_message("MakeFragmentCrmsdResDepth native total residue != query residue depth length");

		return (FragmentScoringMethodOP) new FragmentCrmsdResDepth(priority,
				lowest_acceptable_value, use_lowest, nativePose, picker->get_query_residue_depth());
	}
	if (option[in::file::s].user()) {
		trRmsScore
				<< "Reference structure to score fragments by crmsd loaded from: "
				<< option[in::file::s]()[1] << std::endl;
		core::pose::PoseOP nativePose = new core::pose::Pose;
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::s]()[1]);

		if (nativePose->total_residue() != picker->get_query_residue_depth().size())
			utility_exit_with_message("MakeFragmentCrmsdResDepth native total residue != query residue depth length");

		return (FragmentScoringMethodOP) new FragmentCrmsdResDepth(priority,
				lowest_acceptable_value, use_lowest, nativePose, picker->get_query_residue_depth());
	}
	if (option[in::file::xyz].user()) {

		trRmsScore
				<< "Reference structure to score fragments by crmsd loaded from: "
				<< option[in::file::xyz]().c_str() << std::endl;
                std::string line;
                utility::vector1< utility::vector1< Real > > xyz;
                // std::istream & input = utility::io::izstream( option[ in::file::xyz ]().c_str() );
		std::ifstream input( option[ in::file::xyz ]().c_str() );

		while( getline( input, line ) ) {
			trRmsScore.Warning << line<<std::endl;
			if ( line.substr(0,1) == "#" ) continue;
			std::istringstream line_stream( line );
			utility::vector1<Real> row;
			Real x,y,z;
			line_stream >> x >> y >> z;
			row.push_back(x);
			row.push_back(y);
			row.push_back(z);
			xyz.push_back( row );
		}
		trRmsScore <<  xyz.size() << " atoms found in the reference" << std::endl;

		if (xyz.size() != picker->get_query_residue_depth().size())
			utility_exit_with_message("MakeFragmentCrmsdResDepth xyz size != query residue depth length");

		return (FragmentScoringMethodOP) new FragmentCrmsdResDepth(priority,
				lowest_acceptable_value, use_lowest, xyz, picker->get_query_residue_depth());
	}
	utility_exit_with_message(
			"Can't read a reference structure. Provide it with in::file::s flag");

	return NULL;
}

} // scores
} // frag_picker
} // protocols
