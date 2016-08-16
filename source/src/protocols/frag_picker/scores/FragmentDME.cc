// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FragmentDME.cc
/// @brief  Object that scores a fragment by its Distance Matrix Error (DME) to the native
/// @author David E Kim

#include <protocols/frag_picker/scores/FragmentDME.hh>
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

#include <core/pose/Pose.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <string>

using namespace ObjexxFCL;

namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer trDMEScore(
	"protocols.frag_picker.scores.FragmentDME");

FragmentDME::~FragmentDME() {}

FragmentDME::FragmentDME(Size priority, Real lowest_acceptable_value,
	bool use_lowest, core::pose::PoseOP reference_pose) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "FragmentDME") {
	reference_pose_ = reference_pose;
	n_atoms_ = reference_pose_->total_residue();
	reference_coordinates_.redimension(3, n_atoms_, 0.0);
	fill_CA_coords(*reference_pose_, reference_coordinates_, n_atoms_);
	weights_.redimension(reference_pose_->total_residue(), 1.0);
}

FragmentDME::FragmentDME(Size priority, Real lowest_acceptable_value, bool use_lowest,
	utility::vector1< utility::vector1<Real> > xyz)  :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "FragmentDME") {

	n_atoms_ = xyz.size();
	reference_coordinates_.redimension(3, n_atoms_, 0.0);
	weights_.redimension(n_atoms_, 1.0);
	for ( core::Size i = 1; i <= n_atoms_; i++ ) {
		trDMEScore.Debug << i<<" ";
		for ( core::Size d = 1; d <= 3; ++d ) {
			reference_coordinates_(d, i) = xyz[i][d];
			trDMEScore.Debug <<xyz[i][d]<<" ";
		}
		trDMEScore.Debug <<std::endl;
	}
}

void FragmentDME::fill_CA_coords(core::pose::Pose const& pose,
	FArray2_double& coords, Size n_atoms) {

	trDMEScore.Debug << "Copying coordinates from ... The first residues are: "
		<< pose.residue(1).name3() << " " << pose.residue(2).name3() << " "
		<< pose.residue(3).name3() << std::endl;

	for ( core::Size i = 1; i <= n_atoms; i++ ) {
		id::NamedAtomID idCA("CA", i);
		PointPosition const& xyz = pose.xyz(idCA);
		for ( core::Size d = 1; d <= 3; ++d ) {
			coords(d, i) = xyz[d - 1];
		}
	}
}

bool FragmentDME::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	PROF_START( basic::CA_DME_EVALUATION );
	Real dme=0.0;
	Real dme_native =99.0;
	Size tot_pair=0;
	for ( Size i = 1; i<=f->get_length(); ++i ) {
		tot_pair = tot_pair+i;
	}
	if ( tot_pair > 0 ) {
		VallChunkOP chunk = f->get_chunk();
		for ( Size i = 1; i<=f->get_length(); ++i ) {
			VallResidueOP ri = chunk->at( f->get_first_index_in_vall() + i - 1 );
			Size qindexi = i + f->get_first_index_in_query() - 1;
			for ( Size j = i+1; j<=f->get_length(); ++j ) {
				VallResidueOP rj = chunk->at( f->get_first_index_in_vall() + j - 1 );
				Real xdiff = ri->x()-rj->x();
				Real ydiff = ri->y()-rj->y();
				Real zdiff = ri->z()-rj->z();
				Real dist1 = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;
				Size qindexj = j + f->get_first_index_in_query() - 1;
				Real xrdiff = reference_coordinates_(1, qindexi)-reference_coordinates_(1, qindexj);
				Real yrdiff = reference_coordinates_(2, qindexi)-reference_coordinates_(2, qindexj);
				Real zrdiff = reference_coordinates_(3, qindexi)-reference_coordinates_(3, qindexj);
				Real dist2 = xrdiff*xrdiff + yrdiff*yrdiff + zrdiff*zrdiff;
				dme += (sqrt(dist1)-sqrt(dist2))*(sqrt(dist1)-sqrt(dist2));
			}
		}
		dme_native=sqrt(dme/tot_pair);
	}
	empty_map->set_score_component(dme_native, id_);
	PROF_STOP( basic::CA_DME_EVALUATION );
	if ( (dme_native > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}


void FragmentDME::do_caching(VallChunkOP current_chunk) {

	chunk_coordinates_.redimension(3, current_chunk->size());
	for ( core::Size i = 1; i <= current_chunk->size(); i++ ) {
		VallResidueOP r = current_chunk->at(i);
		chunk_coordinates_(1, i) = r->x();
		chunk_coordinates_(2, i) = r->y();
		chunk_coordinates_(3, i) = r->z();
	}
}

bool FragmentDME::cached_score(FragmentCandidateOP fragment,
	FragmentScoreMapOP scores) {

	std::string tmp = fragment->get_chunk()->chunk_key();
	if ( tmp.compare(cached_scores_id_) != 0 ) {
		do_caching(fragment->get_chunk());
	}

	PROF_START( basic::CA_DME_EVALUATION );

	Real dme=0.0;
	Real dme_native =99.0;
	Size tot_pair=0;
	for ( Size i = 1; i<=fragment->get_length(); ++i ) {
		tot_pair = tot_pair+i;
	}
	if ( tot_pair > 0 ) {
		for ( Size i = 1; i<=fragment->get_length(); ++i ) {
			Size qindexi = i + fragment->get_first_index_in_query() - 1;
			Size vindexi = i + fragment->get_first_index_in_vall() - 1;
			for ( Size j = i+1; j<=fragment->get_length(); ++j ) {
				Size qindexj = j + fragment->get_first_index_in_query() - 1;
				Size vindexj = j + fragment->get_first_index_in_vall() - 1;
				Real xdiff = chunk_coordinates_(1, vindexi)-chunk_coordinates_(1, vindexj);
				Real ydiff = chunk_coordinates_(2, vindexi)-chunk_coordinates_(2, vindexj);
				Real zdiff = chunk_coordinates_(3, vindexi)-chunk_coordinates_(3, vindexj);
				Real dist1 = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;
				Real xrdiff = reference_coordinates_(1, qindexi)-reference_coordinates_(1, qindexj);
				Real yrdiff = reference_coordinates_(2, qindexi)-reference_coordinates_(2, qindexj);
				Real zrdiff = reference_coordinates_(3, qindexi)-reference_coordinates_(3, qindexj);
				Real dist2 = xrdiff*xrdiff + yrdiff*yrdiff + zrdiff*zrdiff;
				dme += (sqrt(dist1)-sqrt(dist2))*(sqrt(dist1)-sqrt(dist2));
			}
		}
		dme_native=sqrt(dme/tot_pair);
	}

	scores->set_score_component(dme_native, id_);
	PROF_STOP( basic::CA_DME_EVALUATION );
	if ( (dme_native > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}

void FragmentDME::clean_up() {
}

FragmentScoringMethodOP MakeFragmentDME::make(Size priority,
	Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
	, std::string) {

	if ( option[in::file::native].user() ) {
		trDMEScore
			<< "Reference structure to score fragments by DME loaded from: "
			<< option[in::file::native]() << std::endl;
		core::pose::PoseOP nativePose( new core::pose::Pose );
		core::import_pose::pose_from_file(*nativePose, option[in::file::native](), core::import_pose::PDB_file);

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new FragmentDME(priority,
			lowest_acceptable_value, use_lowest, nativePose) );
	}
	if ( option[in::file::s].user() ) {
		trDMEScore
			<< "Reference structure to score fragments by DME loaded from: "
			<< option[in::file::s]()[1] << std::endl;
		core::pose::PoseOP nativePose( new core::pose::Pose );
		core::import_pose::pose_from_file(*nativePose, option[in::file::s]()[1], core::import_pose::PDB_file);

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new FragmentDME(priority,
			lowest_acceptable_value, use_lowest, nativePose) );
	}
	if ( option[in::file::xyz].user() ) {

		trDMEScore
			<< "Reference structure to score fragments by DME loaded from: "
			<< option[in::file::xyz]().c_str() << std::endl;
		std::string line;
		utility::vector1< utility::vector1< Real > > xyz;
		// std::istream & input = utility::io::izstream( option[ in::file::xyz ]().c_str() );
		std::ifstream input( option[ in::file::xyz ]().c_str() );

		while ( getline( input, line ) ) {
			trDMEScore.Warning << line<<std::endl;
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
		trDMEScore <<  xyz.size() << " atoms found in the reference" << std::endl;

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new FragmentDME(priority,
			lowest_acceptable_value, use_lowest, xyz) );
	}
	utility_exit_with_message(
		"Can't read a reference structure. Provide it with in::file::s flag");

	return NULL;
}

} // scores
} // frag_picker
} // protocols
