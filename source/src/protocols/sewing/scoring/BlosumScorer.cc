// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file BlosumScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/scoring/BlosumScorer.hh>

//Core headers
#include <core/types.hh>

//Utility headers
#include <basic/Tracer.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

static THREAD_LOCAL basic::Tracer TR("protocols.sewing.scoring.BlosumScorer");

BlosumScorer::BlosumScorer() {}

core::Real
BlosumScorer::score( AssemblyCOP /*assembly*/ ){
	return 0;
}

//bool
//Assembly::check_blosum(
// AtomMap const & atom_map,
// std::set<SewSegment> const & reference_segments,
// Model const & mobile_model
//){
// core::Real blosum_score = 0;
//
// core::sequence::MatrixScoringScheme blosum_reader;
// utility::file::FileName fn(basic::database::full_name("sequence/substitution_matrix/BLOSUM62_aa"));
// blosum_reader.read_from_file(fn);
// utility::vector1< utility::vector1<core::Real> > blosum_matrix = blosum_reader.scoring_matrix();
//
// AtomMap::const_iterator it = atom_map.begin();
// AtomMap::const_iterator it_end = atom_map.end();
// core::Size mobile_first_res = it->second.rsd();
//
// SewSegment matched_mobile_seg;
// for(core::Size i=1; i <= mobile_model.segments_.size(); ++i){
//  for(core::Size j=1; j <= mobile_model.segments_[i].residues_.size(); ++j){
//   if(mobile_model.segments_[i].residues_[j].resnum_ == mobile_first_res){
//    matched_mobile_seg = mobile_model.segments_[i];
//   }
//  }
// }
//
// SewSegment reference_seg = *(reference_segments.begin());
//
// core::Size k1_old = 0;
// core::Size k2_old = 0;
// core::Size res_counter = 0;
// for(; it != it_end; ++it){
//  std::string res1;
//  std::string res2;
//  core::Size k1 = it->first.rsd();
//  core::Size k2 = it->second.rsd();
//  if (k1_old == k1 && k2_old == k2){
//   continue;
//  }
//  k1_old = k1;
//  k2_old = k2;
//  for(core::Size j=1; j <= reference_seg.residues_.size(); ++j){
//   if(reference_seg.residues_[j].resnum_ == k1){
//    res1 = reference_seg.residues_[j].residue_type_;
//
//    TR.Debug << "Residue number " << reference_seg.residues_[j].resnum_ << " of type "
//     << reference_seg.residues_[j].residue_type_ << " from model "
//     << reference_seg.model_id_ << std::endl;
//
//    if(res1 == "HIS_D"){
//     res1 = "HIS";
//    }
//    if(res1 == "CYD"){
//     res1 = "CYS";
//    }
//    if(res1 == "SER_p:phosphorylated"){
//     res1 = "SER";
//    }
//   }
//  }
//  for(core::Size k=1; k <= matched_mobile_seg.residues_.size(); ++k){
//   if(matched_mobile_seg.residues_[k].resnum_ == k2){
//    res2 = matched_mobile_seg.residues_[k].residue_type_;
//
//    TR.Debug << "Residue number " << matched_mobile_seg.residues_[k].resnum_ << " of type "
//    << matched_mobile_seg.residues_[k].residue_type_ << " from model "
//    << matched_mobile_seg.model_id_ << std::endl;
//
//    if(res2 == "HIS_D"){
//     res2 = "HIS";
//    }
//    if(res2 == "CYD"){
//     res2 = "CYS";
//    }
//    if(res2 == "SER_p:phosphorylated"){
//     res2 = "SER";
//    }
//   }
//  }
//  core::Size r1 = core::chemical::aa_from_name(res1);
//  core::Size r2 = core::chemical::aa_from_name(res2);
//  core::Real blosum_value = blosum_matrix[r1][r2];
//  blosum_score += blosum_value;
//  res_counter++;
// }
//
// //Tracer for testing purposes
// TR.Debug << "Blosum score: " << blosum_score << " Number of Residues: " << res_counter << std::endl;
//
// blosum_history_.push_back(std::make_pair(blosum_score, reference_seg.residues_.size()));
//
// //BLOSUM score is being used as a weighting factor for the current edge.
// if(blosum_score < -8){
//  TR << "Edge failed to pass check_blosum" << std::endl;
//  return false;
// }
// return true;
//}


///////RIPPED OUT OF ASSEMBLY


} //scoring namespace
} //sewing namespace
} //protocols namespace
