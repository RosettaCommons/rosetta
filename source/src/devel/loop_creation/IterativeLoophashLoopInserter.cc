// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IterativeLoophashLoopInserter.cc
///
/// @brief Use loophash to find a fragment from loop_anchor to loop_anchor+1. Build ideal residues
/// using the loophash sequence.
/// @author Tim Jacobs

//Unit
//#include <protocols/loophash/IterativeLoophashLoopInserter.hh>
//#include <protocols/loophash/IterativeLoophashLoopInserterCreator.hh>
#include <devel/loop_creation/IterativeLoophashLoopInserter.hh>
#include <devel/loop_creation/IterativeLoophashLoopInserterCreator.hh>

//Numeric
#include <numeric/random/random.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

//protocols
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/ccd_closure.hh>

//utility
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


//namespace protocols {
//namespace loophash {
namespace devel {
namespace loop_creation {

static basic::Tracer TR( "protocols.loophash.IterativeLoophashLoopInserter" );

//****CREATOR METHODS****//
// XRW TEMP std::string
// XRW TEMP IterativeLoophashLoopInserterCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return IterativeLoophashLoopInserter::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP IterativeLoophashLoopInserterCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new IterativeLoophashLoopInserter );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP IterativeLoophashLoopInserter::mover_name()
// XRW TEMP {
// XRW TEMP  return "IterativeLoophashLoopInserter";
// XRW TEMP }
//****END CREATOR METHODS****//

IterativeLoophashLoopInserter::IterativeLoophashLoopInserter():
	LoophashLoopInserter(),
	max_closure_deviation_(5.0),
	max_insertions_(10)
{}

protocols::moves::MoverOP
IterativeLoophashLoopInserter::clone() const {
	return( protocols::moves::MoverOP( new IterativeLoophashLoopInserter( *this ) ) );
}
protocols::moves::MoverOP
IterativeLoophashLoopInserter::fresh_instance() const {
	return protocols::moves::MoverOP( new IterativeLoophashLoopInserter );
}

// XRW TEMP std::string
// XRW TEMP IterativeLoophashLoopInserter::get_name() const {
// XRW TEMP  return "IterativeLoophashLoopInserter";
// XRW TEMP }

void
IterativeLoophashLoopInserter::apply(
	core::pose::Pose & pose
){
	using namespace core;
	using namespace protocols::loophash;

	init(pose);

	//Build the first fragment
	core::Size lh_fragment_begin = loop_anchor()-num_flanking_residues_to_match_+1;
	core::Size lh_fragment_end = loop_anchor()+num_flanking_residues_to_match_;

	HashBuckets hash_buckets = find_fragments(pose, lh_fragment_begin, lh_fragment_end);

	std::pair<core::Size, core::Size> random_fragment =
		get_random_fragment(hash_buckets);

	std::pair<core::Real,core::Real> deviations =
		build_loop(pose, lh_fragment_begin, lh_fragment_end, random_fragment.first, random_fragment.second);

	TR << "Deviations after initial loop insert: " << deviations.first << " " << deviations.second << std::endl;

	lh_fragment_end+=created_loop_.size();
	core::Size loop_cut = created_loop_.cut();

	//core::Size max_lh_fragment_size = created_loop_.size() + (2*num_flanking_residues_to_match_);
	core::Size min_lh_fragment_size=100000;
	for ( core::Size const & i : lh_library_->hash_sizes() ) {
		min_lh_fragment_size = std::min(min_lh_fragment_size, i);
	}

	//restrict future loophash lookups to be maximally the size we just created
	//Pick two random loophash fragment ends
	for ( core::Size i=2; i<=max_insertions_; ++i ) {
		if ( deviations.first <= max_closure_deviation_ && deviations.second <= max_closure_deviation_ ) {
			TR.Debug << "Insertion deviations satisfied after " << i-1 << " insertions." << std::endl;
			break;
		}

		core::Size new_begin = numeric::random::random_range(lh_fragment_begin, loop_cut);
		core::Size new_end = numeric::random::random_range(loop_cut+1, lh_fragment_end);
		core::Size temp_size = new_end-new_begin+1;
		while ( temp_size < min_lh_fragment_size )
				{
			new_begin = numeric::random::random_range(lh_fragment_begin, loop_cut);
			new_end = numeric::random::random_range(loop_cut+1, lh_fragment_end);
			temp_size = new_end-new_begin+1;
		}
		TR.Debug << "new begin, new end, size = " << new_begin << ", " << new_end << ", " << temp_size << std::endl;

		hash_buckets = find_fragments(pose, new_begin, new_end, temp_size, temp_size);

		//If there are no loops for this particular break, then skip
		if ( hash_buckets.size()==0 ) {
			TR << "No loophash fragments found for transform between residues " << new_begin << " and " << new_end << std::endl;
			continue;
		}

		std::pair<core::Size, core::Size> random_fragment =
			get_random_fragment(hash_buckets);

		//apply random fragment and get deviations
		LoopHashMap & hashmap = lh_library_->gethash( random_fragment.first );
		LeapIndex cp = hashmap.get_peptide( random_fragment.second );

		BackboneSegment lh_fragment_bs;
		lh_library_->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , lh_fragment_bs );

		protocols::loops::Loop temp_loop(new_begin, new_end, loop_cut);
		protocols::loops::set_single_loop_fold_tree(pose, temp_loop);

		TR.Debug << "Applying " << lh_fragment_bs.length() << " lh fragment torsions starting at residue " << new_begin << std::endl;
		lh_fragment_bs.apply_to_pose(pose, new_begin, false);
		deviations = protocols::loops::loop_closure::ccd::get_deviation(pose, loop_cut);
		TR.Debug << "Deviations after insertion " << i << ": " << deviations.first << " " << deviations.second << std::endl;

		//  core::Size fewer_n_term = new_begin - lh_fragment_begin;
		//  core::Size fewer_c_term = lh_fragment_end - new_end;
		//  int cur_max_size = max_lh_fragment_size - fewer_n_term - fewer_c_term;
		//  TR.Debug << "min_lh_fragment_size, max_lh_fragment_size, fewer_n_term, fewer_c_term -- " <<
		//   min_lh_fragment_size << ", " << max_lh_fragment_size << ", " << fewer_n_term << ", " << fewer_c_term << std::endl;
		//  TR.Debug << "cut, max_size, new_begin, new_end -- " << loop_cut << ", " << cur_max_size << ", " << new_begin << ", " << new_end << std::endl;
		//  while(cur_max_size < (int)min_lh_fragment_size)
		//  {
		//   new_begin = numeric::random::random_range(lh_fragment_begin, loop_cut);
		//   new_end = numeric::random::random_range(loop_cut+1, lh_fragment_end);
		//
		//   fewer_n_term = new_begin - lh_fragment_begin;
		//   fewer_c_term = lh_fragment_end - new_end;
		//   cur_max_size = max_lh_fragment_size - fewer_n_term - fewer_c_term;
		//
		//   TR.Debug << "min_lh_fragment_size, max_lh_fragment_size, fewer_n_term, fewer_c_term -- " <<
		//    min_lh_fragment_size << ", " << max_lh_fragment_size << ", " << fewer_n_term << ", " << fewer_c_term << std::endl;
		//   TR.Debug << "cut, max_size, new_begin, new_end -- " << loop_cut << ", " << cur_max_size << ", " << new_begin << ", " << new_end << std::endl;
		//  }
	}
	TR << "Final loop deviations " << deviations.first << " " << deviations.second << std::endl;
}

void
IterativeLoophashLoopInserter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){
	LoophashLoopInserter::parse_my_tag(tag, data, filters, movers, pose);

	if ( tag->hasOption("max_closure_deviation") ) {
		max_closure_deviation_ = tag->getOption<core::Real>("max_closure_deviation");
	}
	if ( tag->hasOption("max_insertions") ) {
		max_insertions_ = tag->getOption<core::Real>("max_insertions");
	}
}

std::string IterativeLoophashLoopInserter::get_name() const {
	return mover_name();
}

std::string IterativeLoophashLoopInserter::mover_name() {
	return "IterativeLoophashLoopInserter";
}

void IterativeLoophashLoopInserter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	LoophashLoopInserter::attributes_for_tag( attlist );
	attlist
		//+ XMLSchemaAttribute( "max_closure_deviation", xsct_real, "Why is this separately specified in the derived class?")
		+ XMLSchemaAttribute( "max_insertions", xsct_real, "Probably-integral description of how many insertions are appropriate over the course of simulation" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string IterativeLoophashLoopInserterCreator::keyname() const {
	return IterativeLoophashLoopInserter::mover_name();
}

protocols::moves::MoverOP
IterativeLoophashLoopInserterCreator::create_mover() const {
	return protocols::moves::MoverOP( new IterativeLoophashLoopInserter );
}

void IterativeLoophashLoopInserterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	IterativeLoophashLoopInserter::provide_xml_schema( xsd );
}


} //loop creation
} //devel
//} //protocols
//} //loophash
