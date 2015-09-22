// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/RDCScore.cc
/// @brief  Object that scores a fragment by its RDC
/// @author Ray Wang (wangyr@u.washington.edu)

#include <protocols/frag_picker/scores/RDCScore.hh>
#include <core/scoring/ResidualDipolarCoupling.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#ifdef WIN32
#include <protocols/frag_picker/FragmentPicker.hh>
#endif

// option key includes
#include <basic/options/keys/OptionKeys.hh>

#include <core/chemical/ChemicalManager.hh>

#include <basic/Tracer.hh>

// utils
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer trRDCScore(
	"protocols.frag_picker.scores.RDCScore");


RDCScore::RDCScore(
	Size priority,
	Real lowest_acceptable_value,
	bool use_lowest
) :
	CachingScoringMethod(
	priority,
	lowest_acceptable_value,
	use_lowest,
	"RDCScore"
	) {
	rdc_file_ = core::scoring::ResidualDipolarCouplingOP( new core::scoring::ResidualDipolarCoupling );
	rdc_raw_data_ = rdc_file_->get_RDC_data();
}

void RDCScore::do_caching(VallChunkOP current_chunk) {
	std::string ctmp = current_chunk->chunk_key();
	if ( ctmp.compare("change to 'cached_scores_id_' when ready") != 0 ) {
		return; // CACHING NOT BUILT IN YET
	}
}

bool RDCScore::cached_score(
	FragmentCandidateOP fragment,
	FragmentScoreMapOP empty_map
) {
	return score( fragment, empty_map );
}

void RDCScore::clean_up() {
}

bool RDCScore::score(
	FragmentCandidateOP fragment,
	FragmentScoreMapOP empty_map
) {
	////////////////////////////////////////////////////////////////
	//
	// make a frag_pose
	//
	Size offset_query = fragment->get_first_index_in_query(); //index of the 1st rsd in a query sequence covered by this fragment
	Size offset_vall  = fragment->get_first_index_in_vall();  //index of the 1st rsd in a vall chunk covered by this fragment

	core::pose::Pose frag_pose;
	core::chemical::ResidueTypeSetCOP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	//get the whole vall chunk where the fragment from
	VallChunkOP chunk = fragment->get_chunk();

	std::string fragment_seq ( fragment->get_length(), 'A' );
	core::pose::make_pose_from_sequence( frag_pose, fragment_seq, *rsd_set);

	for ( Size i = 1; i <= fragment->get_length(); ++i ) {
		Size vall_rsd_num  = offset_vall  + i - 1;

		VallResidueOP vall_rsd = chunk->at( vall_rsd_num );

		frag_pose.set_phi  ( i, vall_rsd->phi()   );
		frag_pose.set_psi  ( i, vall_rsd->psi()   );
		frag_pose.set_omega( i, vall_rsd->omega() );
	}

	//////////////////////////////////////////////////////////////
	//
	// store rdc for a given query_sequence


	// rdc data handling
	core::scoring::ResidualDipolarCoupling::RDC_lines rdc_data_given_fragment_lines;
	for ( core::scoring::ResidualDipolarCoupling::RDC_lines::const_iterator it = rdc_raw_data_.begin(); it != rdc_raw_data_.end(); ++it ) {
		//try to put RDC of a fragment into a place
		if ( it->res1() >= offset_query && it->res1() <= offset_query + fragment->get_length() - 1 ) {
			if ( it->res2() >= offset_query && it->res2() <= offset_query + fragment->get_length() - 1 ) {
				core::scoring::RDC selected ( it->res1() - offset_query + 1, it->atom1(), it->res2() - offset_query + 1, it->atom2(), it->Jdipolar());
				rdc_data_given_fragment_lines.push_back( selected );
			}
		}
	}

	//read rdc data as lines
	core::scoring::ResidualDipolarCouplingOP rdc_data_given_fragment( new core::scoring::ResidualDipolarCoupling ( rdc_data_given_fragment_lines ) );
	core::scoring::store_RDC_in_pose( rdc_data_given_fragment, frag_pose );

	///////////////////////////////////////////////////////////////
	//
	// scoring the frag_pose using RDC in range
	//
	//Real rdc_score = 0.0;
	Real rdc_score = rdc_data_given_fragment->compute_dipscore( frag_pose );

	empty_map->set_score_component( rdc_score, id_ );
	if ( ( rdc_score > lowest_acceptable_value_ ) && ( use_lowest_ == true ) ) {
		return false;
	}
	return true;
}


FragmentScoringMethodOP MakeRDCScore::make(
	Size priority,
	Real lowest_acceptable_value,
	bool use_lowest,
	FragmentPickerOP, //picker,
	std::string
) {
	return ( FragmentScoringMethodOP )
		FragmentScoringMethodOP( new RDCScore(
		priority,
		lowest_acceptable_value,
		use_lowest
		) );
}

} // scores
} // frag_picker
} // protocols
