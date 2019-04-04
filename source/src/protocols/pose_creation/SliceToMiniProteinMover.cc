// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SliceToMiniProteinMover.cc
/// @brief This mover chops a big protein into miniProteins
/// @author TJ Brunette (tjbrunette@gmail.com)
///

// Unit headers
#include <protocols/pose_creation/SliceToMiniProteinMoverCreator.hh>
#include <protocols/pose_creation/SliceToMiniProteinMover.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.pose_creation.SliceToMiniProteinMover" );

#include <utility/tag/Tag.hh>

#include <core/pose/PDBInfo.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//for parsing residue_selectors
#include <core/conformation/ResidueFactory.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>
#include <protocols/rosetta_scripts/util.hh>

//filters
#include <protocols/simple_ddg/SSElementBisectddGFilter.hh>

#include <core/types.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/selection.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace pose_creation {

using std::map;
using namespace core;
using core::pose::Pose;
using utility::vector1;

SliceToMiniProteinMover::SliceToMiniProteinMover()
: moves::Mover("SliceToMiniProteinMover"),
	selector_()
{
}

moves::MoverOP
SliceToMiniProteinMover::clone() const
{
	return moves::MoverOP( utility::pointer::make_shared<SliceToMiniProteinMover>( *this ) );
}

moves::MoverOP
SliceToMiniProteinMover::fresh_instance() const
{
	return moves::MoverOP( utility::pointer::make_shared<SliceToMiniProteinMover>() );
}


utility::vector1<SliceToMiniProteinMover::SSElement> SliceToMiniProteinMover::parse_ss(core::pose::Pose const & pose) const{
	utility::vector1<SliceToMiniProteinMover::SSElement> ss_elements;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	core::select::residue_selector::ResidueSelectorCOP labeled_res_selector(utility::pointer::make_shared<core::select::residue_selector::ResiduePDBInfoHasLabelSelector>("miniCore"));
	core::select::residue_selector::ResidueSubset residues( pose.total_residue(), false );
	residues = labeled_res_selector->apply( pose );
	Size reassign_short_terminal_loop_=2;
	Size first_residue = 1;
	Size last_residue = residues.size();
	bool first_found =false;
	bool last_found=false;

	if ( residues[1]==true ) {
		first_found=true;
	}
	for ( Size ii=1; ii<=residues.size(); ++ii ) {
		if ( residues[ii]==true && first_found==false ) {
			first_residue=ii;
			first_found=true;
		}
		if ( residues[ii]==false && first_found==true && last_found==false ) {
			last_residue=ii-1;
			last_found=true;
		}
	}
	//pass 1----------------------------------------
	char lastSecStruct = dssp.get_dssp_secstruct( first_residue );
	Size startSS = first_residue;
	Size endSS = 0;
	for ( core::Size ii = first_residue+1; ii <= last_residue; ++ii ) {
		if ( dssp.get_dssp_secstruct(ii)!=lastSecStruct ) {
			endSS = ii-1;
			TR.Debug << "type: " << lastSecStruct << " startSS: " << startSS << " endSS: " << endSS << std::endl;
			SliceToMiniProteinMover::SSElement ss_element_tmp(startSS,endSS,lastSecStruct);
			ss_elements.push_back(ss_element_tmp);
			startSS=ii;
		}
		lastSecStruct = dssp.get_dssp_secstruct(ii);
	}
	endSS = last_residue;
	SliceToMiniProteinMover::SSElement ss_element_tmp(startSS,endSS,lastSecStruct);
	ss_elements.push_back(ss_element_tmp);
	//pass 2----------------------------------------Remove first and last loop. Could be done faster in above loop, but for future modifications I'm doing it below.
	if ( !reassign_short_terminal_loop_ ) {
		return(ss_elements);
	}

	bool delete_n_term_loop=false;
	bool delete_c_term_loop=false;
	Size ss_length = ss_elements[1].end_res-ss_elements[1].start_res;
	if ( (ss_length<reassign_short_terminal_loop_)&&(ss_elements[1].type=="L") ) {
		ss_elements[2].start_res =ss_elements[1].start_res;
		delete_n_term_loop=true;
	}
	ss_length = ss_elements[ss_elements.size()].end_res-ss_elements[ss_elements.size()].start_res;
	if ( (ss_length<reassign_short_terminal_loop_)&&(ss_elements[ss_elements.size()].type=="L") ) {
		ss_elements[ss_elements.size()-1].end_res =ss_elements[ss_elements.size()].end_res;
		delete_c_term_loop=true;
	}
	Size start_pos = 1;
	if ( delete_n_term_loop ) {
		start_pos=reassign_short_terminal_loop_;
	}
	Size end_pos = ss_elements.size();
	if ( delete_c_term_loop ) {
		end_pos = ss_elements.size()-1;
	}
	utility::vector1<SliceToMiniProteinMover::SSElement> ss_elements_v2;
	for ( Size ii=start_pos; ii<=end_pos; ++ii ) {
		ss_elements_v2.push_back(ss_elements[ii]);
	}
	return(ss_elements_v2);
}


utility::vector1<SliceToMiniProteinMover::Chunk> SliceToMiniProteinMover::parse_chunks(core::pose::Pose const & pose,utility::vector1<SliceToMiniProteinMover::SSElement> ss_elements){
	utility::vector1<SliceToMiniProteinMover::Chunk> chunks;
	for ( Size ii=1; ii<=ss_elements.size(); ++ii ) {
		if ( ss_elements[ii].type=="H" ) {
			bool too_long=false;
			for ( Size jj=ii+1; jj<=ss_elements.size(); ++jj ) {
				if ( ss_elements[jj].type=="H" ) {
					Size start_res=ss_elements[ii].start_res;
					Size end_res=ss_elements[jj].end_res;
					Size numb_elements = (jj-ii+2)/2; //assume loops between secondary structure elements
					Size length =end_res-start_res+1;
					if ( (numb_elements>=min_sse_count_) && (too_long==false) ) { //allow slightly longer chunks. Will be trimmed later.
						if ( length>=max_length_ ) {
							too_long=true;
						}
						utility::vector1<Size> pose_positions;
						for ( Size ll=start_res; ll<=end_res; ++ll ) {
							pose_positions.push_back(ll);
						}
						pose::Pose sliced_pose;
						pdbslice(sliced_pose,pose,pose_positions);
						SliceToMiniProteinMover::Chunk chunk(sliced_pose,start_res,end_res);
						chunks.push_back(chunk);
					}
				}
			}
		}
	}
	return(chunks);
}


utility::vector1<SliceToMiniProteinMover::Chunk> SliceToMiniProteinMover::filter_chunks(utility::vector1<SliceToMiniProteinMover::Chunk> chunks){
	utility::vector1<SliceToMiniProteinMover::Chunk> filtered_chunks;
	bool threshold = ddg_ala_slice_score_;
	bool report_avg = false;
	Size ignore_terminal_SS=1;
	bool only_n_term= false;
	bool only_c_term= false;
	bool skip_ss_element=true;
	bool report_sasa_instead=false;
	bool convert_charged_res_to_ala=false;
	protocols::simple_ddg::SSElementBisectddGFilterOP ddg_bisect_filter( utility::pointer::make_shared<protocols::simple_ddg::SSElementBisectddGFilter>(scorefxn_,threshold,report_avg,ignore_terminal_SS,only_n_term,only_c_term,skip_ss_element,report_sasa_instead,convert_charged_res_to_ala,relax_mover_) );
	//relax_mover_(relax_mover)
	for ( Size ii=1; ii<=chunks.size(); ++ii ) {
		Real ddg = ddg_bisect_filter->compute(chunks[ii].pose);
		chunks[ii].ddg=ddg;
	}
	//for(Size ii=1; ii<=chunks.size(); ++ii){
	// TR <<"pre-filter chunk_options" << chunks[ii].ddg << "," << chunks[ii].start_res << "," << chunks[ii].end_res << std::endl;
	//}
	for ( Size ii=1; ii<=chunks.size(); ++ii ) {
		if ( chunks[ii].ddg<ddg_ala_slice_score_ ) {
			filtered_chunks.push_back(chunks[ii]);
		}
	}
	for ( Size ii=1; ii<=filtered_chunks.size(); ++ii ) {
		TR <<"chunk_options" << filtered_chunks[ii].ddg << "," << filtered_chunks[ii].start_res << "," << filtered_chunks[ii].end_res << std::endl;
	}
	return(filtered_chunks);
}


utility::vector1<SliceToMiniProteinMover::Chunk> SliceToMiniProteinMover::trim_chunks(utility::vector1<SliceToMiniProteinMover::Chunk> chunks){
	utility::vector1<SliceToMiniProteinMover::Chunk> trimmed_chunks;
	//Trim both tail helices equally.
	//If either tail_helices<min_helix_length get rid of this option completely. Trim of a single helix exists in other options.
	//step2: trim the size to insure that it's shorter then max_length
	for ( Size ii=1; ii<=chunks.size(); ++ii ) {
		Size length = chunks[ii].pose.total_residue();
		if ( length<=max_length_ ) {
			trimmed_chunks.push_back(chunks[ii]);
		} else {
			Size residue_to_delete = length-max_length_;
			utility::vector1<SSElement> ss_elements = parse_ss(chunks[ii].pose);
			Size n_helix_length = ss_elements[1].end_res-ss_elements[1].start_res+1;
			Size last_helix = ss_elements.size();
			Size c_helix_length = ss_elements[last_helix].end_res-ss_elements[last_helix].start_res+1;
			Size n_term_trim =0;
			Size c_term_trim =0;
			while ( residue_to_delete>(n_term_trim+c_term_trim) ) {
				if ( (n_helix_length-n_term_trim)>(c_helix_length-c_term_trim) ) {
					n_term_trim++;
				} else {
					c_term_trim++;
				}
			}
			if ( ((n_helix_length-n_term_trim)>=min_sse_length_) && ((c_helix_length-c_term_trim)>=min_sse_length_) ) {
				utility::vector1<Size> pose_positions;
				for ( Size ll=1+n_term_trim; ll<=chunks[ii].pose.total_residue()-c_term_trim; ++ll ) {
					pose_positions.push_back(ll);
				}
				pose::Pose sliced_pose;
				pdbslice(sliced_pose,chunks[ii].pose,pose_positions);
				chunks[ii].pose=sliced_pose;
				chunks[ii].start_res=chunks[ii].start_res-n_term_trim;
				chunks[ii].end_res=chunks[ii].end_res-c_term_trim;
				trimmed_chunks.push_back(chunks[ii]);
			}
		}
	}
	return(trimmed_chunks);
}


void SliceToMiniProteinMover::setup_chunks_for_output(utility::vector1<SliceToMiniProteinMover::Chunk> & chunks){
	if ( output_mode_=="longest" ) {
		//length ties go to the longest ddg
		Size longest=0;
		Size position=0;
		for ( Size ii=1; ii<=chunks.size(); ++ii ) {
			Size length = chunks[ii].pose.total_residue();
			if ( length>longest ) {
				longest=length;
				position=ii;
			}
		}
		for ( Size ii=1; ii<=chunks.size(); ++ii ) {
			if ( ii != position ) {
				chunks[ii].outputed=true;
			}
		}
	}
	if ( output_mode_=="best_ddg" ) {
		Real ddg=99;
		Size position=0;
		for ( Size ii=1; ii<=chunks.size(); ++ii ) {
			if ( ddg>chunks[ii].ddg ) {
				ddg=chunks[ii].ddg;
				position=ii;
			}
		}
		for ( Size ii=1; ii<=chunks.size(); ++ii ) {
			if ( ii != position ) {
				chunks[ii].outputed=true;
			}
		}
	}
}


core::pose::PoseOP SliceToMiniProteinMover::get_additional_output(){
	for ( Size ii=1; ii<=final_chunks_.size(); ++ii ) {
		if ( final_chunks_[ii].outputed == false ) {
			final_chunks_[ii].outputed=true;
			set_last_move_status(protocols::moves::MS_SUCCESS);
			pose::PoseOP poseOP( utility::pointer::make_shared<pose::Pose>( final_chunks_[ii].pose ) );
			TR<<"output_chunk" <<  final_chunks_[ii].ddg << "," << final_chunks_[ii].start_res << "," << final_chunks_[ii].end_res << std::endl;
			return(poseOP);
		}
	}
	set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
	return nullptr;
}

void SliceToMiniProteinMover::apply( Pose & pose )
{
	core::select::residue_selector::ResidueSubset residues( pose.total_residue(), false );
	residues = selector_->apply( pose );
	for ( core::Size ii=1; ii<=pose.size(); ++ii ) {
		if ( residues[ii]==1 ) {
			pose.pdb_info()->add_reslabel(ii, "miniCore");
		}
	}
	utility::vector1<SSElement> ss_elements = parse_ss(pose);
	utility::vector1<Chunk> chunks = parse_chunks(pose,ss_elements);
	utility::vector1<Chunk> trimmed_chunks = trim_chunks(chunks);
	utility::vector1<Chunk> filtered_chunks = filter_chunks(trimmed_chunks);
	setup_chunks_for_output(filtered_chunks);
	final_chunks_ = filtered_chunks;
	core::pose::PoseOP tmpPoseOP=get_additional_output();
	if ( tmpPoseOP!=nullptr ) {
		pose=*tmpPoseOP;
	}
}

void SliceToMiniProteinMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & )
{
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	max_length_ = tag->getOption< core::Size >( "max_length", 119);
	min_sse_count_ = tag->getOption< core::Size >( "min_sse_count", 4);
	min_sse_length_ = tag->getOption< core::Size >( "min_sse_length", 13);
	ddg_ala_slice_score_ = tag->getOption< Real >( "ddg_ala_slice_score", -5.0);
	output_mode_ = tag->getOption< std::string >( "output_mode", "all"); //other options longest,best_ddg

	if ( tag->hasOption("residue_selector") ) {
		selector_ =  protocols::rosetta_scripts::parse_residue_selector( tag, data ) ;
	}

	if ( tag->hasOption( "relax_mover" ) ) {
		relax_mover_ = protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "relax_mover"), movers );
	} else {
		relax_mover_=NULL;
	}


}

void SliceToMiniProteinMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	XMLSchemaRestriction attachment_termini_type;

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist
		+ XMLSchemaAttribute::attribute_w_default("max_length", xsct_non_negative_integer, "max length of miniprotein","119")
		+ XMLSchemaAttribute::attribute_w_default("min_sse_count", xsct_non_negative_integer, "min secondary structure element count","4")
		+ XMLSchemaAttribute::attribute_w_default("min_sse_length", xsct_non_negative_integer, "min secondary structure element length","13")
		+ XMLSchemaAttribute::attribute_w_default("ddg_ala_slice_score", xsct_real, "energy when slicing up the miniprotein","-5.0")
		+ XMLSchemaAttribute( "relax_mover" , xs_string , "Optionally define a mover which will be applied prior to computing the system energy in the unbound state." )
		+ XMLSchemaAttribute::attribute_w_default( "output_mode" , xs_string , "output mode, all,longest,best_ddg","all" )
		+ utility::tag::optional_name_attribute();
	core::select::residue_selector::attributes_for_parse_residue_selector(
		attlist, "residue_selector",
		"name of a residue selector that specifies the subset allowed to be part of the miniprotein. Note. This segment must be contiguous" );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.element_name( SliceToMiniProteinMover::mover_name() )
		.description( "slices up repeat proteins to miniproteins"  )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );

}


std::string SliceToMiniProteinMover::get_name() const {
	return mover_name();
}

std::string SliceToMiniProteinMover::mover_name() {
	return "SliceToMiniProteinMover";
}


std::string SliceToMiniProteinMoverCreator::keyname() const {
	return SliceToMiniProteinMover::mover_name();
}

protocols::moves::MoverOP
SliceToMiniProteinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( utility::pointer::make_shared<SliceToMiniProteinMover>() );
}

void SliceToMiniProteinMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SliceToMiniProteinMover::provide_xml_schema( xsd );
}

} // pose_creation
} // protocols
