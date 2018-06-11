// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/SpliceOut.cc
/// @brief
/// @author Gideon Lapidoth (glapidoth@gmail.com)


// Unit headers
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/splice/SpliceSegment.hh>
#include <protocols/splice/SpliceOut.hh>
#include <protocols/splice/SpliceOutTail.hh>
#include <protocols/splice/SpliceOutAntibody.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/splice/SpliceOutCreator.hh>
#include <protocols/splice/SpliceOutAntibodyCreator.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AA.hh>
#include <protocols/protein_interface_design/filters/TorsionFilter.hh>
#include <protocols/protein_interface_design/util.hh>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/predicate.hpp>//for comparing string case insensitive
#include <protocols/task_operations/RestrictChainToRepackingOperation.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/splice/util.hh>
#include <protocols/moves/Mover.hh>
#include <core/id/SequenceMapping.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/datacache/DataMapObj.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <boost/algorithm/string.hpp>
//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/task_operations/DesignAroundOperation.hh>
#include <protocols/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/task_operations/ThreadSequenceOperation.hh>
#include <protocols/task_operations/SeqprofConsensusOperation.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <numeric/xyzVector.hh>
#include <protocols/loops/FoldTreeFromLoopsWrapper.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <numeric/constants.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/scoring/Energies.hh>
#include <numeric/xyz.functions.hh>
#include <protocols/simple_moves/CutChainMover.hh>
//////////////////////////////////////////////////
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
///////////////////////////////////////////////////
#include <fstream>
#include <ctime>
#include <protocols/splice/RBInMover.hh>
#include <protocols/splice/RBOutMover.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace splice {

using namespace core::conformation;

static  basic::Tracer TR( "protocols.splice.SpliceOutAntibody" );
static  basic::Tracer TR_ccd( "protocols.splice.Splice_ccd" );
static  basic::Tracer TR_constraints( "protocols.splice.Splice_constraints" );
static  basic::Tracer TR_pssm( "protocols.splice.Splice_pssm" );
static  basic::Tracer TR_min( "protocols.splice.Splice_min" );

std::string SpliceOutAntibodyCreator::keyname() const {
	return SpliceOutAntibodyCreator::mover_name();
}

protocols::moves::MoverOP SpliceOutAntibodyCreator::create_mover() const {
	return protocols::moves::MoverOP( new SpliceOutAntibody );
}

std::string SpliceOutAntibodyCreator::mover_name() {
	return "SpliceOutAntibody";
}


SpliceOutAntibody::SpliceOutAntibody() :  Mover(SpliceOutAntibodyCreator::mover_name()),antibody_DB_("additional_protocol_data/splice/antibodies/"),vl_vh_cut_(0)
{
	splicemanager.profile_weight_away_from_interface(1.0);
	design_shell_ = 6.0;
	repack_shell_ = 8.0;
	basic::options::option[basic::options::OptionKeys::out::file::pdb_comments].value(true);
	splicemanager.rb_sensitive(false);
	tolerance_=0.23;
	delete_hairpin( false );
	delete_hairpin_c( 0 );
	delete_hairpin_n( 0 );
	mover_type_={"MinMover","LoopMover_Refine_CCD","TailSegmentMover"};
	call_mover=nullptr;
	protein_family_to_database_["antibodies"] = "additional_protocol_data/splice/antibodies/";
	//Hard coding the correct order of segments to be spliced
	splicemanager.segment_names_ordered().push_back("L1_L2");
	splicemanager.segment_names_ordered().push_back("L3");
	splicemanager.segment_names_ordered().push_back("H1_H2");
	splicemanager.segment_names_ordered().push_back("H3");
	splicemanager.tail_seg("");
}

SpliceOutAntibody::~SpliceOutAntibody() {}

void SpliceOutAntibody::apply(core::pose::Pose & pose) {
	using namespace protocols::rosetta_scripts;
	using core::chemical::DISULFIDE;
	set_last_move_status(protocols::moves::MS_SUCCESS);


	find_disulfide_postions(*splicemanager.template_pose(),template_cys_pos_);
	find_disulfide_postions(pose,pose_cys_pos_);
	find_vl_vh_cut(pose);
	adjust_n_ter_tail_length(pose);
	//Because the pssms are built here we have to make sure that the pose length matches the pssm length. So I am changing the tail length at this point. Gideon, Oct17

	find_disulfide_postions(pose,pose_cys_pos_);//update pose disulfie postions after n-ter change

	assign_from_res_to_res(pose);
	splicemanager.adjust_stem_positions_by_template_=false;

	if ( (splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2") ) {
		splicemanager.tail_seg(""); //Every time I go through another itteration of Splice mover I have to reset the segment type.
		source_to_res(find_nearest_res(*source_pose_, pose, splicemanager.pose_to_res(), 0/*chain*/));
		if ( !(source_pose_->residue(source_to_res()+1).name3()=="CYS") ) {
			TR<<"source_to_res():"<<source_to_res()<<std::endl;
			utility_exit_with_message(" to_res on source pdb is not a cysteine, check alignment between pose and source pdb\n");
		}
	} else if ( (splicemanager.segment_type()=="L3")||(splicemanager.segment_type()=="H3") ) {
		source_to_res(source_pose_->total_residue());
		TR<<"source_to_res:"<<source_to_res()<<std::endl;
	} else {
		utility_exit_with_message("Illegal segment type passed for antibodies\n");
	}
	//From_res for all segments (L1_L2, H1_H2,L3, H3 must be a cysteine)
	source_from_res(find_nearest_res(*source_pose_, pose, splicemanager.pose_from_res(), 0/*chain*/));
	TR<<"source_from_res"<<source_from_res()<<std::endl;
	if ( !(source_pose_->residue(source_from_res()-1).name3()=="CYS") ) {
		utility_exit_with_message(" from_res on source pdb is not a cysteine, check alignmnet between pose and source pdb\n");
	}



	if ( (splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2") ) {
		write_to_database_=false;// Writing to database happens after tail apply
		SpliceOut::apply(pose);
		if ( get_last_move_status() ) { //failed job is true
			return;
		}
		//setup tail segment mover
		splicemanager.tail_seg("n");
		core::Size const prev_template_from_res(splicemanager.template_from_res());
		core::Size const prev_template_to_res(splicemanager.template_to_res());
		core::Size const prev_pose_from_res(splicemanager.pose_from_res());
		core::Size const prev_pose_to_res(splicemanager.pose_to_res());
		core::Size const cut_site(splicemanager.cut_site());

		splicemanager.template_to_res(splicemanager.template_from_res()-2);//remember that from res in the loop segment is one after the cysteine, so the from res for the tail should be one before the cys
		splicemanager.template_from_res(splicemanager.template_from_res()-2);
		if ( splicemanager.segment_type()=="L1_L2" ) {
			splicemanager.adjust_stem_positions_by_template_=true;
		} else {
			splicemanager.adjust_stem_positions_by_template_=false;
			splicemanager.pose_to_res(protocols::rosetta_scripts::find_nearest_res(pose, *splicemanager.template_pose(),splicemanager.template_from_res(), 0/*chain*/));
			splicemanager.pose_from_res(vl_vh_cut_+1);
		}

		protocols::moves::MoverOP tmp_submover_ = submover_;
		call_mover_tmp = call_mover;
		call_mover = call_mover_tail_;
		submover_ = tail_submover_;
		superimposed(true);//no need to superimpose again for L1_L2 tail segment
		SpliceOut::apply(pose);
		submover_ = tmp_submover_;
		call_mover  = call_mover_tmp;
		if ( get_last_move_status() ) { //failed job is true
			return;
		}
		write_to_database_=true;
		SpliceOutTail::write_database_to_file(pose);
		splicemanager.template_to_res(prev_template_to_res);
		splicemanager.template_from_res(prev_template_from_res);
		splicemanager.pose_to_res(prev_pose_to_res);
		splicemanager.pose_from_res(prev_pose_from_res);
		splicemanager.cut_site(cut_site);
		SpliceOut::write_database_to_file(pose);
	} else if ( (splicemanager.segment_type()=="L3")||(splicemanager.segment_type()=="H3") ) {
		splicemanager.tail_seg("c");
		call_mover = call_mover_tail_;
		submover_ = tail_submover_;
		superimposed(true);//no need to superimpose again for L1_L2 tail segment
		write_to_database_=true;
		SpliceOut::apply(pose);
	} else {
		utility_exit_with_message("Illegal segment type passed for antibodies\n");
	}

}//apply


std::string SpliceOutAntibody::get_name() const {
	return SpliceOutAntibodyCreator::mover_name();
}



void SpliceOutAntibody::parse_my_tag(TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const & filters, protocols::moves::Movers_map const & movers, core::pose::Pose const & /*pose*/) {
	if ( tag->hasOption("segment") ) {
		splicemanager.segment_type(tag->getOption<std::string>("segment"));
	} else {
		TR.Error<<TR.Red<<"Must specify antibody segment (L1_L2,H1_H2,H3,L3)"<<TR.Reset<<std::endl;
		runtime_assert((tag->hasOption("segment")));
	}

	utility::vector1<TagCOP> const sub_tags(tag->getTags());

	parse_SpliceOut_tags(tag,movers,filters);
	splicemanager.parse_tags(tag,data);
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		SpliceOut::abstract_parse_tag(tag);
	}

	scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));


	handle_tail_mover_tag(tag,movers);

	if ( tag->hasOption("use_sequence_profile") ) {
		splicemanager.use_sequence_profiles(tag->getOption<bool>("use_sequence_profile"));
	}
	if ( splicemanager.use_sequence_profiles() ) {
		TR << " Using Splice antibody pssms" << std::endl;
		splicemanager.segment_type(tag->getOption<std::string>("segment"));
		auto find_segment = std::find (splicemanager.segment_names_ordered().begin(), splicemanager.segment_names_ordered().end(), splicemanager.segment_type());
		if ( find_segment == splicemanager.segment_names_ordered().end() ) {
			TR.Error << "ERROR !! segment '"<<splicemanager.segment_type()<<"is not a recognized antibody segmetn \n" << std::endl;
			runtime_assert( find_segment != splicemanager.segment_names_ordered().end() );
		}
		for ( auto const segment_type : splicemanager.segment_names_ordered() ) {
			TR<<segment_type<<std::endl;
			SpliceSegmentOP splice_segment( new SpliceSegment );
			splice_segment->read_pdb_profile_file( antibody_DB(), segment_type );
			splice_segment->read_many (antibody_DB(),segment_type);
			TR<<splice_segment->get_cluster_name_by_PDB("1AHWD")<<std::endl;
			splicemanager.splice_segments().insert( std::pair< std::string, SpliceSegmentOP >( segment_type, splice_segment ) );
		}
		splicemanager.use_sequence_profiles(true);
		splicemanager.profile_weight_away_from_interface( tag->getOption< core::Real >( "profile_weight_away_from_interface", 1.0 ) );


		pdb_to_H3_seq_map_=read_pdb_seq(antibody_DB()+"pssm/H3/H3_seq" );

	} else {
		return;
	}

}

protocols::moves::MoverOP SpliceOutAntibody::clone() const {
	return (protocols::moves::MoverOP(new SpliceOutAntibody(*this)));
}

void SpliceOutAntibody::find_disulfide_postions(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos) {
	cys_pos.clear();
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue(i).has_variant_type(core::chemical::DISULFIDE) ) {
			cys_pos.push_back(i);
		}
	}
}

void SpliceOutAntibody::update_vl_vh_cut() {
	if ( splicemanager.segment_type().find_first_of("L") !=std::string::npos ) {
		vl_vh_cut_+=splicemanager.residue_diff();
	}
}

void SpliceOutAntibody::assign_from_res_to_res(core::pose::Pose const pose){
	core::conformation::Conformation const & conf(pose.conformation());
	if ( splicemanager.segment_type()=="L1_L2" ) {
		splicemanager.template_from_res(template_cys_pos_[1]+1);
		splicemanager.template_to_res(template_cys_pos_[2]-1);
		splicemanager.pose_from_res(pose_cys_pos_[1]+1);
		splicemanager.pose_to_res(pose_cys_pos_[2]-1);
	} else if ( splicemanager.segment_type()=="H1_H2" ) {
		splicemanager.template_from_res(template_cys_pos_[3]+1);
		splicemanager.template_to_res(template_cys_pos_[4]-1);
		splicemanager.pose_from_res(pose_cys_pos_[3]+1);
		splicemanager.pose_to_res(pose_cys_pos_[4]-1);
	} else if ( splicemanager.segment_type()=="L3" ) {
		splicemanager.template_from_res(template_cys_pos_[2]+1);
		splicemanager.pose_from_res(pose_cys_pos_[2]+1);
		splicemanager.pose_to_res(vl_vh_cut_);
	} else if ( splicemanager.segment_type()=="H3" ) {
		splicemanager.template_from_res(template_cys_pos_[4]+1);
		splicemanager.pose_from_res(pose_cys_pos_[4]+1);
		splicemanager.pose_to_res(conf.num_chains()==1?pose.total_residue():conf.chain_end(1));

	}

}


std::string SpliceOutAntibody_complex_type_name_for_subsubtag( std::string const & foo ) {
	return "subsubtag_spliceoutAntibody_" + foo + "_type";
}

std::string SpliceOutAntibody_complex_type_name_for_subtag( std::string const & foo ) {
	return "subtag_spliceoutSntibody_" + foo + "_type";
}

void SpliceOutAntibody::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;


	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "tolerance", xsct_real, "XRW TO DO", "0.23" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_chain_break", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "CG_const", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rb_sensitive", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "chain_num", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "cut_site", xsct_non_negative_integer, "residue number of where to place cut", "1" )

		+ XMLSchemaAttribute( "segment", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "superimposed", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_n", xsct_non_negative_integer, "XRW TO DO", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_c", xsct_non_negative_integer, "XRW TO DO", "13" );

	// The "Segments" subtag
	AttributeList segments_subtag_attlist;
	segments_subtag_attlist + XMLSchemaAttribute( "current_segment", xs_string, "XRW TO DO" );

	// The "segment" sub-subtag"
	AttributeList subtag_segments_subtag_attlist;
	subtag_segments_subtag_attlist + XMLSchemaAttribute( "pdb_profile_match", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "profiles", xs_string, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator segment_subsubtag_gen;
	segment_subsubtag_gen.complex_type_naming_func( & SpliceOutAntibody_complex_type_name_for_subsubtag )
		.element_name( "segment" )
		.description( "individual segment tag" )
		.add_attributes( subtag_segments_subtag_attlist )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_already_defined_subelement( "segment", SpliceOutAntibody_complex_type_name_for_subsubtag/*, 0*/ );

	XMLSchemaComplexTypeGenerator segments_subtag_gen;
	segments_subtag_gen.complex_type_naming_func( & SpliceOutAntibody_complex_type_name_for_subtag )
		.element_name( "Segments" )
		.description( "Wrapper for multiple segments tags" )
		.add_attributes( segments_subtag_attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );


	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( "Segments", SpliceOutAntibody_complex_type_name_for_subtag/*, 0*/ );


	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "add_sequence_constraints_only", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "template_file", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "source_pdb", xs_string, "XRW TO DO");
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "from_res", xsct_refpose_enabled_residue_number, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "to_res", xsct_refpose_enabled_residue_number, "XRW TO DO", "0" )
		+ XMLSchemaAttribute( "design_task_operations", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "residue_numbers_setter", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "torsion_database", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "design_shell", xsct_real, "XRW TO DO", "6.0" )
		+ XMLSchemaAttribute::attribute_w_default( "repack_shell", xsct_real, "XRW TO DO", "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "rms_cutoff", xsct_real, "XRW TO DO", "999999" )
		+ XMLSchemaAttribute::attribute_w_default( "rms_cutoff_loop", xsct_real, "XRW TO DO", "999999" )
		+ XMLSchemaAttribute::attribute_w_default( "res_move", xsct_non_negative_integer, "XRW TO DO", "1000" )
		+ XMLSchemaAttribute::attribute_w_default( "randomize_cut", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "cut_secondarystruc", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_ala", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "design", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_original_sequence", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rtmin", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "allow_all_aa", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "locked_residue", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "checkpointing_file", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "splice_filter", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "mover", xs_string, "Which mover to use to close the segment" )
		+ XMLSchemaAttribute( "tail_mover", xs_string, "Which mover to use to close the segment" )
		+ XMLSchemaAttribute::attribute_w_default( "restrict_to_repacking_chain2", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "use_sequence_profile", xsct_rosetta_bool, "XRW TO DO", "true" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, subelements );

}

void SpliceOutAntibody::handle_tail_mover_tag(TagCOP const tag,protocols::moves::Movers_map const & movers){
	std::string mover_name("");
	if ( tag->hasOption("tail_mover") ) {
		mover_name = tag->getOption< std::string >( "tail_mover" );
	} else {
		utility_exit_with_message("Must specify mover to copy source conformation onto template segment\n");
	}
	protocols::moves::Movers_map::const_iterator  find_mover ( movers.find( mover_name ));
	tail_submover_ = find_mover->second;
	TR<<tail_submover_->get_name()<<std::endl;
	if ( find_mover == movers.end() ) {
		TR.Error << "ERROR !! mover '"<<mover_name<<"' not found in map: \n" << tag << std::endl;
		runtime_assert( find_mover != movers.end() );
	}

	if ( std::find(mover_type_.begin(), mover_type_.end(), tail_submover_->get_name()) == mover_type_.end() ) {
		utility_exit_with_message("Please choose only \"Minmover\", \"LoopMover_Refine_CCD\", or \"TailSegmentMover\" for \"tail_mover=\" tag \n");
	} else {
		if ( tail_submover_->get_name()=="MinMover" ) {
			call_mover_tail_ =&SpliceOut::min_mover;
		} else if ( tail_submover_->get_name()=="LoopMover_Refine_CCD" ) {
			call_mover_tail_ = &SpliceOut::ccd_mover;
		} else if ( tail_submover_->get_name()=="TailSegmentMover" ) {
			call_mover_tail_ = &SpliceOut::tail_mover;
		}
	}
}

void SpliceOutAntibodyCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SpliceOutAntibody::provide_xml_schema( xsd );
}

std::string
SpliceOutAntibody::mover_name()  {
	return "SpliceOutAntibody";
}

void
SpliceOutAntibody::find_vl_vh_cut(core::pose::Pose pose)  {
	protocols::simple_moves::CutChainMover ccm;
	vl_vh_cut_ = ccm.chain_cut(pose); //
	TR<<"vl vh cut is: "<<vl_vh_cut_<<std::endl;
}

void SpliceOutAntibody::set_source_from_to_res(){
	using namespace protocols::rosetta_scripts;

	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		SpliceOut::set_source_from_to_res();
	}
	if ( splicemanager.tail_seg()=="n" ) {
		SpliceOutTail::set_source_from_to_res();
	}
	return;
}

void SpliceOutAntibody::place_cut_site_in_segment(core::pose::Pose const & pose){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		SpliceOut::place_cut_site_in_segment(pose);
	} else {
		SpliceOutTail::place_cut_site_in_segment(pose);
	}
}

void SpliceOutAntibody::set_fold_tree_nodes(core::pose::Pose const & pose){
	fold_tree_nodes_.clear();

	/*Build the following Fold tree for the Vl/Vh fragment
	___________________________
	|                         |
	|                         |
	<------*<--------*------>//<-----*<-------*------->

	*/  update_vl_vh_cut();
	find_disulfide_postions(pose,pose_cys_pos_);
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[1],1,-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[2],(int) pose_cys_pos_[1],-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[2],(int) vl_vh_cut_,-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[3],(int) vl_vh_cut_+1,-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[4],(int) pose_cys_pos_[3],-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[4],(int) pose.conformation().chain_end(1),-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[4],(int) pose_cys_pos_[2],1}});

	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		SpliceOut::set_fold_tree_nodes(pose);
		return;
	}

	if ( pose.conformation().num_chains()>1 ) { //if ligand is present we need to add edge between receptor and ligand
		core::Size CoM = (core::Size ) core::pose::residue_center_of_mass( pose, pose.conformation().chain_begin(2), pose.conformation().chain_end(2) );
		fold_tree_nodes_.push_back({{(int)pose_cys_pos_[4],(int) CoM,2}});
		fold_tree_nodes_.push_back({{(int) CoM, (int)pose.conformation().chain_end(2),2}});
		fold_tree_nodes_.push_back({{(int) CoM, (int)pose.conformation().chain_begin(2),2}});
	}
}

core::Size SpliceOutAntibody::set_anchor_res(){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceOut::set_anchor_res();
	}
	return SpliceOutTail::set_anchor_res();
}

void SpliceOutAntibody::set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceOut::set_loop_length_change(llc);
	}
	return SpliceOutTail::set_loop_length_change(llc);
}

void SpliceOutAntibody::superimpose_source_on_pose( core::pose::Pose const & pose, core::pose::Pose & source_pose ){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceOut::superimpose_source_on_pose(pose,source_pose);
	}
	return SpliceOutTail::superimpose_source_on_pose(pose,source_pose);
}


void SpliceOutAntibody::build_ideal_segment(core::pose::Pose & pose){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceOut::build_ideal_segment(pose);
	}
	return SpliceOutTail::build_ideal_segment(pose);
}

void
SpliceOutAntibody::write_database_to_file(core::pose::Pose const & pose){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceOut::write_database_to_file(pose);
	}
	return SpliceOutTail::write_database_to_file(pose);
}

//Change the length of the N-ter tail in the H1_H2/L1_L2 segment to match that of the source pose
void
SpliceOutAntibody::adjust_n_ter_tail_length(core::pose::Pose & pose){
	using namespace protocols::rosetta_scripts;
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	int residue_diff = 0;

	if ( splicemanager.segment_type()=="L1_L2" ) {
		source_from_res(find_nearest_res(*source_pose_, pose, pose_cys_pos_[1]+1, 0/*chain*/));
		if ( !(source_pose_->residue(source_from_res()-1).name3()=="CYS") ) {
			TR<<"source_from_res():"<<source_from_res()<<std::endl;
			utility_exit_with_message(" to_res on source pdb is not a cysteine, check alignment between pose and source pdb\n");
		}
		residue_diff =(source_from_res()-2)-(pose_cys_pos_[1]-1);
		runtime_assert(residue_diff<1000);//sanity check to make sure we are not getting weird residue diffs, gideonla apr18

		TR<<"Tail residue diff: "<<residue_diff<<std::endl;
		llc.loop_start(1);
		llc.loop_end(pose_cys_pos_[1]-1);
	}
	if ( splicemanager.segment_type()=="H1_H2" ) {
		source_from_res(find_nearest_res(*source_pose_, pose, pose_cys_pos_[3]+1, 0/*chain*/));
		if ( !(source_pose_->residue(source_from_res()-1).name3()=="CYS") ) {
			TR<<"source_from_res():"<<source_from_res()<<std::endl;
			utility_exit_with_message(" to_res on source pdb is not a cysteine, check alignment between pose and source pdb\n");
		}
		residue_diff =(source_from_res()-2)-(pose_cys_pos_[3]-(vl_vh_cut_+1));
		llc.loop_start(vl_vh_cut_+1);
		llc.loop_end(pose_cys_pos_[3]-1);
	}
	TR<<"Source from res before tail change length  = "<<source_from_res()<<std::endl;
	llc.delta(residue_diff);
	llc.tail(1);
	llc.direction(1);//N-ter direction
	llc.apply(pose);

}
std::string SpliceOutAntibody::name_for_filter(){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceOut::name_for_filter();
	}
	return SpliceOutTail::name_for_filter();
}


} //splice
} //protocols
