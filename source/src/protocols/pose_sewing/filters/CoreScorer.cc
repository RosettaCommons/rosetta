// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/filters/CoreScorer.cc
/// @brief a filter that evaluates pairwise MotifScores
/// @author frankdt (frankdt@email.unc.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/pose_sewing/filters/CoreScorer.hh>
#include <protocols/pose_sewing/filters/CoreScorerCreator.hh>
#include <core/select/residue_selector/BlockSelector.hh>
#include <protocols/pose_sewing/util.hh>
#include <protocols/pose_sewing/data_storage/DsspShiftArray.hh>
#include <protocols/sewing/scoring/MotifScorer.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/motif/util.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/selection.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <utility/string_util.hh>






static basic::Tracer TR( "protocols.pose_sewing.filters.CoreScorer" );

namespace protocols {
namespace pose_sewing {
namespace filters {

CoreScorer::CoreScorer():
	protocols::filters::Filter( "CoreScorer" )
{
}

CoreScorer::~CoreScorer()
{}

void
CoreScorer::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{

	if ( tag->hasOption("score_cutoff") ) {
		score_cutoff_ = tag->getOption<core::Real>("score_cutoff");
	}
	if ( tag->hasOption("distance_cutoff") ) {
		distance_cutoff_ = tag->getOption<core::Real>("distance_cutoff");
	}
	if ( tag->hasOption("link_cutoff") ) {
		link_cutoff_ = tag->getOption<core::Real>("link_cutoff");
	}

	if ( tag->hasOption("window_width") ) {
		window_width_ = tag->getOption<core::Size>("window_width");
	}

	if ( tag->hasOption("sum") ) {
		sum_ = tag->getOption<bool>("sum");
	}
	if ( tag->hasOption("distance_mode") ) {
		distance_mode_ = tag->getOption<bool>("distance_mode");
	}
	if ( tag->hasOption( "selector" ) ) {
		selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "selector" ), datamap );
	}
	if ( ! selector_ ) {
		utility_exit_with_message("CoreScorer requires a selector for blockwise calculations!");
	}
}

protocols::filters::FilterOP
CoreScorer::clone() const
{
	return protocols::filters::FilterOP( new CoreScorer( *this ) );
}


protocols::filters::FilterOP
CoreScorer::fresh_instance() const
{
	return protocols::filters::FilterOP( new CoreScorer );
}

bool
CoreScorer::apply( core::pose::Pose const & pose ) const
{
	core::Real sm = this->report_sm(pose);
	return(sm >= link_cutoff_);
}

/*
* The idea is to ensure that the linkages between blocks -- which is really to say the packable cross-block residue pairs -- are colocated. We want A's links to B to be in the same place, or nearly so, to A's links to C.
* Therefore, once we assign blocks, we build a map of residues to sets of links, where each residue either does or does not link to a given block. If a residue links to multiple residues in the same block, that link is counted only once.
* We then iterate through a sliding window in that map, summing the number of links to each block found in a given window.
* Then we examine the second-highest number of links for that window, and pass iff that is greater than the listed cutoff.
* The idea is that, in a hypothetical triangular arrangement of long helices, A's links to B will be distant from A's links to C, so the second-highest link count for both windows will be 0.
*/

core::Real
CoreScorer::report_sm( core::pose::Pose const & pose) const
{


	protocols::sewing::scoring::MotifScorer motif_scorer = protocols::sewing::scoring::MotifScorer();

	// make blockselectors here
	std::map<core::Size,core::Size> block_assignments;
	core::Size block_count = 1;

	core::select::residue_selector::ResidueSubset selection = selector_->apply( pose );
	calculate_blocks_from_subset(block_assignments, selection );

	for ( auto current_pair : block_assignments ) {
		if ( current_pair.second > block_count ) {
			block_count = current_pair.second;
		}
	}

	if ( block_count <= 2 ) {
		TR << "Blockwise invalid with <2 blocks!" << std::endl;
		return 0.0;
	}


	core::Size reverse_it = 0;
	if ( block_count >= 2 ) {
		reverse_it = 1;
	}


	core::Real running_score = 0.0;
	std::map<core::Size,std::set<core::Size>> block_links; //residue-wise list of links by block
	std::set<core::Size> empty_set;
	for ( core::Size N_resnum = 1; N_resnum <= pose.size(); ++N_resnum ) {
		if ( !pose.residue(N_resnum).is_virtual_residue() && (block_assignments[N_resnum] > 0) ) {

			block_links[N_resnum] = empty_set;
			core::conformation::Residue N_res = pose.residue(N_resnum);
			utility::vector1< core::conformation::Atom > const & N_atoms = N_res.atoms();
			numeric::xyzTransform< core::Real > stub1 = core::pose::motif::get_backbone_reference_frame(N_atoms[ 1 ].xyz(), N_atoms[ 2 ].xyz(), N_atoms[ 3 ].xyz() );
			char aa1 = N_res.name1();
			char ss1 = pose.secstruct(N_resnum);

			//Look at all other residues in the pose
			//  Make sure the residue is assigned a block
			//  Make sure N and C resnum are not in the same block
			//  If the motif score is less than the score cutoff, count it as a link.
			//  Store the block number in block_links for each N_resnum
			//Strand Blocks:
			// Only count a link of a strand to paired strands or other elements

			for ( core::Size C_resnum = 1; C_resnum <= pose.size(); ++C_resnum ) {
				if ( !pose.residue(C_resnum).is_virtual_residue() && (block_assignments[C_resnum] > 0)  && (block_assignments[C_resnum] != block_assignments[N_resnum]) ) {
					char ss2 = pose.secstruct(C_resnum);

					core::conformation::Residue C_res = pose.residue(C_resnum);
					utility::vector1< core::conformation::Atom > const & C_atoms = C_res.atoms();
					numeric::xyzTransform< core::Real > stub2 = core::pose::motif::get_backbone_reference_frame(C_atoms[ 1 ].xyz(), C_atoms[ 2 ].xyz(), C_atoms[ 3 ].xyz() );
					char aa2 = C_res.name1();
					if ( distance_mode_ ) {
						if ( N_atoms[2].xyz().distance(C_atoms[2].xyz()) <= distance_cutoff_ ) {
							block_links[N_resnum].insert(block_assignments[C_resnum]);
						}
					} else {
						running_score = -1 * motif_scorer.get_score(stub1, ss1, aa1, stub2, ss2, aa2);
						if ( running_score <= score_cutoff_ ) { // only count if motifscore is good enough, but only count once
							block_links[N_resnum].insert(block_assignments[C_resnum]);
						}
					}
				}
			}
		}
	}
	// now we have a map of all the links for each residue

	std::list<core::Size> residues_to_consider;
	std::map<core::Size,core::Size> block_windows;

	for ( core::Size count_resnum = 1; count_resnum <= pose.size(); ++count_resnum ) {
		if ( block_assignments[count_resnum] != 0 && block_links.count(count_resnum) != 0 ) {

			residues_to_consider.push_back(count_resnum);
			if ( residues_to_consider.size() > window_width_ ) {
				residues_to_consider.pop_front();
			}
			if ( residues_to_consider.size() == window_width_ || ((block_assignments[count_resnum] == 0 || count_resnum == pose.size() )&& !residues_to_consider.empty()) ) {
				// score here
				utility::vector1<core::Size> link_counts;
				for ( core::Size work_block = 1; work_block <= block_count; ++work_block ) {
					link_counts.push_back(0);
				}
				for ( auto list_it = residues_to_consider.begin(); list_it != residues_to_consider.end(); ++list_it ) {
					for ( auto const & link : block_links[*list_it] ) {
						++link_counts[link]; // we end up with the sum total of all the links from a particular window, sorted by block

					}
				}
				link_counts[block_assignments[*residues_to_consider.begin()]] = 0; //Null self-block counting
				std::sort(link_counts.begin(),link_counts.end());

				if ( block_windows.count(block_assignments[*residues_to_consider.begin()]) == 0 ) {
					block_windows[block_assignments[*residues_to_consider.begin()]] = 0;
				}

				//Reverse IT is one - helps to figure out what the highest and less highest is.
				if ( block_windows.count(block_assignments[*residues_to_consider.begin()]) < link_counts[link_counts.size()-reverse_it] ) {
					block_windows[block_assignments[*residues_to_consider.begin()]] = link_counts[link_counts.size()-reverse_it];
				}

				// end scoring
			}
		} else {
			residues_to_consider.clear();
		}
	}

	//Now we have a list of the second highest links (blockwise)
	if ( sum_ ) {
		core::Size total = 0;
		for ( auto current_pair : block_windows ) {
			TR << "block score: "<< current_pair.first<<" " << current_pair.second << std::endl;
			total = total + current_pair.second;
		}
		return total;
	} else {
		core::Size min_links = 9999;
		for ( auto current_pair : block_windows ) {
			TR << "block score: "<< current_pair.first<<" " << current_pair.second << std::endl;
			TR << current_pair.first << ":" << current_pair.second << std::endl;
			if ( min_links > current_pair.second ) {
				min_links = current_pair.second;
			}
		}

		return min_links;
	}
	return 0;
}

void
CoreScorer::report( std::ostream &, core::pose::Pose const & pose) const
{
	TR << this->report_sm(pose) << std::endl;
}

std::string CoreScorer::name() const {
	return class_name();
}

std::string CoreScorer::class_name() {
	return "CoreScorer";
}

void
CoreScorer::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}
void
CoreScorer::set_score_cutoff(core::Real score_cutoff){
	score_cutoff_ = score_cutoff;
}
void
CoreScorer::set_distance_cutoff(core::Real distance_cutoff){
	distance_cutoff_ = distance_cutoff;
}
void
CoreScorer::set_window_width(core::Size window_width){
	window_width_ = window_width;
}
void
CoreScorer::set_link_cutoff(core::Size link_cutoff){
	link_cutoff_ = link_cutoff;
}
void
CoreScorer::set_sum(bool sum){
	sum_ = sum;
}
void
CoreScorer::set_distance_mode(bool distance_mode){
	distance_mode_ = distance_mode;
}

void CoreScorer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaAttribute score_cutoff = XMLSchemaAttribute::attribute_w_default( "score_cutoff", xsct_real, "score per residue cutoff","-0.2" );

	XMLSchemaAttribute distance_cutoff = XMLSchemaAttribute::attribute_w_default( "distance_cutoff", xsct_real, "cutoff for distance mode","4" );

	XMLSchemaAttribute link_cutoff = XMLSchemaAttribute( "link_cutoff", xsct_non_negative_integer, "The number of links required of a window to other blocks. Now defaults to per-ss settings unless this is overriden");

	XMLSchemaAttribute window_width = XMLSchemaAttribute( "window_width", xsct_non_negative_integer, "The number of residues that constitute a window in a particular block. Now defaults to per-ss settings unlesss this is overridden");

	AttributeList attlist;
	attlist
		+ score_cutoff
		+ distance_cutoff
		+ link_cutoff
		+ window_width
		+ XMLSchemaAttribute::attribute_w_default( "sum", xs_boolean, "Are we taking the sum across a window, as opposed to the maximum?" , "false")
		+ XMLSchemaAttribute::attribute_w_default( "distance_mode", xs_boolean, "do we care about motifscore or raw distance" , "false");

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "selector",  "Selector to create SS blocks from. Requred for blockwise, used as a mask for elementwise." );

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.


	std::string docs = "Details:\n The idea is to ensure that the linkages between blocks -- which is really to say the packable cross-block residue pairs -- are colocated. We want A's links to B to be in the same place, or nearly so, to A's links to C.\n"
		"Therefore, once we assign blocks, we build a map of residues to sets of links, where each residue either does or does not link to a given block. If a residue links to multiple residues in the same block, that link is counted only once.\n"
		"We then iterate through a sliding window in that map, summing the number of links to each block found for each residue in a given window.\n"
		"What we then have is a map of blocks to counts.  \n"
		"For a window of 3, and two possible blocks where res 1 is bound to 1 block, 2 bound to 2, and 3 bound to 1, we have a list of [2, 2].  We then use the second-highest as number to use for filtering for that window (Block count in subset of 2, we use the highest). If we only have connections to a single block, then the total is 0.\n"
		"We then use this for pass/fail of the window.  If any of the windows fail, the filter fails.\n"
		"The idea is that, in a hypothetical triangular arrangement of long helices, A's links to B will be distant from A's links to C, so the second-highest link count for both windows will be 0.\n";

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "A filter for the core of a denovo design. Only works for blocks or elements larger than 1. For n blocks of 2, it basically becomes WindowPoseCompatMotif.\n\n." + docs, attlist );



}

/////////////// Creator ///////////////

protocols::filters::FilterOP
CoreScorerCreator::create_filter() const
{
	return protocols::filters::FilterOP( new CoreScorer );
}

std::string
CoreScorerCreator::keyname() const
{
	return CoreScorer::class_name();
}

void CoreScorerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CoreScorer::provide_xml_schema( xsd );
}

} //protocols
} //pose_sewing
} //filters
