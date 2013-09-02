// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite && is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions && developers.
// For more information, see http://www.rosettacommons.org/.

///// @file   protocols/frag_picker/scores/FragmentChunkCrms.cc
///// @brief  Object that scores a fragment by its chunkrms to the native
///// @author Lei Shi (shilei@uw.edu)

#include <protocols/frag_picker/scores/FragmentChunkCrms.hh>
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
#include <basic/prof.hh>

//Auto Headers
#include <core/id/NamedAtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <string>

//Sequence alignment
#include <basic/database/open.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/NWAligner.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/MatrixScoringScheme.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer trTmScore(
		"protocols.frag_picker.scores.FragmentChunkCrms");

FragmentChunkCrms::~FragmentChunkCrms() {}

FragmentChunkCrms::FragmentChunkCrms(Size priority, Real lowest_acceptable_value,
		bool use_lowest, std::string query_sequence, core::pose::PoseOP reference_pose, FArray1D_int& seqmapping) :
	FragmentScoringMethod(priority, lowest_acceptable_value, use_lowest, "FragmentChunkCrms") {
	reference_pose_ = reference_pose;
  fragment_pose_ = new core::pose::Pose;
	n_atoms_ = reference_pose_->total_residue();
	reference_coordinates_.redimension(3, 4*n_atoms_, 0.0);
	fill_bb_coords(*reference_pose_, reference_coordinates_, n_atoms_);
	weights_.redimension(reference_pose_->total_residue()*4, 1.0);
	seqmapping_=seqmapping;
  core::pose::make_pose_from_sequence(*fragment_pose_, query_sequence, *(chemical::ChemicalManager::get_instance()->residue_type_set("centroid")));

  //seqmapping_.redimension(query_sequence.length(), 0);
}

void FragmentChunkCrms::fill_bb_coords(core::pose::Pose const& pose, FArray2_double& coords, Size n_atoms) {

	trTmScore.Debug << "Copying coordinates from ... The first residues are: "
			<< pose.residue(1).name3() << " " << pose.residue(2).name3() << " "
			<< pose.residue(3).name3() << std::endl;

  core::Size n_at = 1;
	for (core::Size i = 1; i <= n_atoms; i++) {
    id::NamedAtomID idN("N", i);
    PointPosition const& xyzN = pose.xyz(idN);
    for (core::Size d = 1; d <= 3; ++d) {
      coords(d, n_at) = xyzN[d - 1];
    }
    n_at++;

    id::NamedAtomID idCA("CA", i);
    PointPosition const& xyzCA = pose.xyz(idCA);
    for (core::Size d = 1; d <= 3; ++d) {
      coords(d, n_at) = xyzCA[d - 1];
    }
    n_at++;

    id::NamedAtomID idC("C", i);
    PointPosition const& xyzC = pose.xyz(idC);
    for (core::Size d = 1; d <= 3; ++d) {
      coords(d, n_at) = xyzC[d - 1];
    }
    n_at++;

    id::NamedAtomID idO("O", i);
    PointPosition const& xyzO = pose.xyz(idO);
    for (core::Size d = 1; d <= 3; ++d) {
      coords(d, n_at) = xyzO[d - 1];
    }
    n_at++;
	}
}

void FragmentChunkCrms::fill_bb_coords(core::pose::Pose const& pose, FArray2_double& coords, FArray1D_int& seqmapping) {

  trTmScore.Debug << "Copying coordinates from according to seqmapping " << std::endl;

  core::Size n_at = 1;
  for (core::Size i = 1; i <= pose.total_residue(); i++) {
		if (seqmapping(i)==1) {
    id::NamedAtomID idN("N", i);
    PointPosition const& xyzN = pose.xyz(idN);
    for (core::Size d = 1; d <= 3; ++d) {
      coords(d, n_at) = xyzN[d - 1];
    }
    n_at++;

    id::NamedAtomID idCA("CA", i);
    PointPosition const& xyzCA = pose.xyz(idCA);
    for (core::Size d = 1; d <= 3; ++d) {
      coords(d, n_at) = xyzCA[d - 1];
    }
    n_at++;

    id::NamedAtomID idC("C", i);
    PointPosition const& xyzC = pose.xyz(idC);
    for (core::Size d = 1; d <= 3; ++d) {
      coords(d, n_at) = xyzC[d - 1];
    }
    n_at++;

    id::NamedAtomID idO("O", i);
    PointPosition const& xyzO = pose.xyz(idO);
    for (core::Size d = 1; d <= 3; ++d) {
      coords(d, n_at) = xyzO[d - 1];
    }
    n_at++;
  }
	}
}

bool FragmentChunkCrms::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	if ((Size) fragment_coordinates_.size2() < n_atoms_*4) {
		fragment_coordinates_.redimension(3, n_atoms_*4, 0.0);
	}

  for (Size i = 1; i <= f->get_length(); i++) {
      fragment_pose_->set_phi( i, f->get_residue(i)->phi() );
      fragment_pose_->set_psi( i, f->get_residue(i)->psi() );
      fragment_pose_->set_omega( i, f->get_residue(i)->omega() );
  }

  fill_bb_coords(*fragment_pose_, fragment_coordinates_, seqmapping_);

	Real chunkrms = numeric::model_quality::rms_wrapper(4*n_atoms_, fragment_coordinates_, reference_coordinates_);

	empty_map->set_score_component(chunkrms, id_);
	if ((chunkrms > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;

	return true;
}

void sequencealign(core::sequence::SequenceOP seq1, core::sequence::SequenceOP seq2, FArray1D_int& seqmapping) {

  seqmapping.redimension(std::max(seq1->length(),seq2->length()), 0);
	//std::cout << "map size: " << std::min(seq1->length(),seq2->length()) << std::endl;
  utility::file::FileName blosum62( basic::database::full_name("sequence/substitution_matrix/BLOSUM62"));
  core::sequence::ScoringSchemeOP blosum_score( new core::sequence::MatrixScoringScheme( -11, -1, blosum62 ) );

	//alignment
	core::sequence::NWAligner nw_aligner;
	core::sequence::SequenceAlignment global_align = nw_aligner.align( seq1, seq2, blosum_score ) ;
	core::id::SequenceMapping mapping = global_align.sequence_mapping(1, 2);
  //std::cout << mapping;
  //std::cout << std::endl;

  for (core::Size i = 1; i <= std::max(seq1->length(),seq2->length()); ++i) {
           if (mapping[i]) {
							seqmapping(i)=1;
           } else {
							seqmapping(i)=0;
					}
   }

}

FragmentScoringMethodOP MakeFragmentChunkCrms::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker
		, std::string) {

  //trTmScore << "QUERY_SEQUENCE " << picker->get_query_seq_string() << std::endl;
  //Initialized an array and align to each other
  core::sequence::SequenceOP seq1 ( new core::sequence::Sequence(picker->get_query_seq_string(),"blank",1));
	core::sequence::SequenceOP seq2;
	FArray1D_int seqmapping;

	if (option[in::file::native].user()) {
		trTmScore
				<< "Reference structure to score fragments by chunkrms loaded from: "
				<< option[in::file::native]() << std::endl;
		core::pose::PoseOP nativePose = new core::pose::Pose;
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::native]());
    seq2 = new core::sequence::Sequence(*nativePose);
		seqmapping.redimension(std::max(seq1->length(),seq2->length()), 0);
		sequencealign(seq1,seq2,seqmapping);

		return (FragmentScoringMethodOP) new FragmentChunkCrms(priority, lowest_acceptable_value, use_lowest, picker->get_query_seq_string(), nativePose, seqmapping);

	} else if (option[in::file::s].user()) {
		trTmScore
				<< "Reference structure to score fragments by chunkrms loaded from: "
				<< option[in::file::s]()[1] << std::endl;
		core::pose::PoseOP nativePose = new core::pose::Pose;
		core::import_pose::pose_from_pdb(*nativePose, option[in::file::s]()[1]);
    seq2 = new core::sequence::Sequence(*nativePose);
		seqmapping.redimension(std::max(seq1->length(),seq2->length()), 0);
		sequencealign(seq1,seq2,seqmapping);

		return (FragmentScoringMethodOP) new FragmentChunkCrms(priority, lowest_acceptable_value, use_lowest, picker->get_query_seq_string(), nativePose, seqmapping);

	} else {
		utility_exit_with_message("Can't read a reference structure. Provide it with in::file::s or in:file:native flag");

	return NULL;
	}
}

} // scores
} // frag_picker
} // protocols
