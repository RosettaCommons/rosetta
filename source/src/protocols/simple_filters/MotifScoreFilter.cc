// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MotifScoreFilter.cc
/// @brief Scores poses using will's motifScore
/// @author TJ Brunette

// Unit headers
#include <protocols/simple_filters/MotifScoreFilter.hh>
#include <protocols/simple_filters/MotifScoreFilterCreator.hh>


// Project Headers
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/motif/reference_frames.hh>
#include <numeric/xyzTransform.hh>
#include <core/pose/util.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <string>

using namespace core;
using namespace std;
using utility::vector1;

namespace protocols {
namespace simple_filters {
typedef numeric::xyzTransform<double> Xform;
static basic::Tracer TR( "protocols.simple_moves.MotifScoreFilter" );

protocols::filters::FilterOP
MotifScoreFilterCreator::create_filter() const { return protocols::filters::FilterOP( new MotifScoreFilter ); }

std::string
MotifScoreFilterCreator::keyname() const { return "MotifScore"; }

MotifScoreFilter::MotifScoreFilter(){
	mman_ = core::scoring::motif::MotifHashManager::get_instance();
}

bool MotifScoreFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	TR << "Motif score is " << score << ". ";
	if ( score >= score_threshold_ ) {
		TR <<"passing." << std::endl;
		return true;
	} else {
		TR <<"failing." << std::endl;
		return false;
	}
}

void MotifScoreFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	out<<"motif score: "<< score <<'\n';
}

core::Real MotifScoreFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	return( score );
}

core::Real MotifScoreFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::scoring;
	using namespace core::scoring::motif;
	double score = 0.0;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	for ( size_t ir = 1; ir <= pose.size(); ++ir ) {
		Xform const ibb_stub = core::pose::motif::get_backbone_reference_frame(pose,ir);
		char ss1 = dssp.get_dssp_secstruct( ir );
		char aa1 = pose.residue(ir).name1();
		for ( size_t jr = ir+1; jr <= pose.size(); ++jr ) {
			Real dist = pose.residue(ir).xyz("CA").distance(pose.residue(jr).xyz("CA"));
			if ( dist < 12 ) {
				char ss2 = dssp.get_dssp_secstruct( jr );
				char aa2 = pose.residue(jr).name1();
				//std::cout << ss1 << ss2 << " " << aa1 << aa2 << "dist" << dist <<  std::endl;
				Xform const jbb_stub = core::pose::motif::get_backbone_reference_frame(pose,jr);
				Xform const Xbb = ibb_stub.inverse() * jbb_stub;
				core::scoring::motif::XformScoreCOP xs_bb_fxn1(mman_->get_xform_score_BB_BB(ss1,ss2,aa1,aa2));
				core::scoring::motif::XformScoreCOP xs_bb_fxn2(mman_->get_xform_score_BB_BB(ss2,ss1,aa2,aa1));
				if ( xs_bb_fxn1 != nullptr ) {
					score += xs_bb_fxn1->score_of_bin(Xbb);
					score += xs_bb_fxn2->score_of_bin(Xbb.inverse());
				}
				//std::cout << "score pos" << ir << "," << jr << "," << score << std::endl;
			}
		}
	}
	return(score);

}

void MotifScoreFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & //pose
){
	if ( ! tag->hasOption( "threshold" ) ) {
		score_threshold_ = -999;
	} else {
		score_threshold_ = tag->getOption<core::Real>( "threshold" );
	}
}

} // filters
} // protocols
