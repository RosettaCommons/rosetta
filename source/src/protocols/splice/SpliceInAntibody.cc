// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/SpliceInAntibody.cc
/// @brief
/// @author Gideon Lapidoth (glapidoth@gmail.com)

// Unit headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/splice/SpliceSegment.hh>
#include <protocols/splice/TailSegmentMover.hh>
#include <protocols/splice/SpliceInAntibody.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/splice/SpliceInAntibodyCreator.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AA.hh>
#include <protocols/protein_interface_design/filters/TorsionFilter.hh>
#include <protocols/protein_interface_design/util.hh>
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
#include <protocols/toolbox/superimpose.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/protein_interface_design/movers/SetAtomTree.hh>


namespace protocols {
namespace splice {

using namespace core::conformation;
static  basic::Tracer TR_constraints("protocols.splice.Splice_constraints");
static  basic::Tracer TR("protocols.splice.SpliceInAntibody");
static  basic::Tracer TR_pssm("protocols.splice.Splice_pssm");
std::string SpliceInAntibodyCreator::keyname() const {
	return SpliceInAntibodyCreator::mover_name();
}

protocols::moves::MoverOP SpliceInAntibodyCreator::create_mover() const {
	return protocols::moves::MoverOP(new SpliceInAntibody);
}

std::string SpliceInAntibodyCreator::mover_name() {
	return "SpliceInAntibody";
}

SpliceInAntibody::SpliceInAntibody() : Mover(SpliceInAntibodyCreator::mover_name()),antibody_DB_("additional_protocol_data/splice/antibodies/")
{
	splicemanager.profile_weight_away_from_interface(1.0);
	splicemanager.design_shell(6.0);
	splicemanager.repack_shell(8.0);
	basic::options::option[basic::options::OptionKeys::out::file::pdb_comments].value(
		true);
	splicemanager.rb_sensitive(false);
	tolerance_ = 0.23;
	allowed_cuts_ = 1;
	dbase_subset_.clear();
	end_dbase_subset_ = DataccacheBoolDataOP( new basic::datacache::DataMapObj<bool> );
	protein_family_to_database_["antibodies"] = "additional_protocol_data/splice/antibodies/";
	//Hard coding the correct order of segments to be spliced
	splicemanager.segment_names_ordered().push_back("L1_L2");
	splicemanager.segment_names_ordered().push_back("L3");
	splicemanager.segment_names_ordered().push_back("H1_H2");
	splicemanager.segment_names_ordered().push_back("H3");

}
SpliceInAntibody::~SpliceInAntibody() = default;


void SpliceInAntibody::apply(core::pose::Pose & pose) {
	using namespace protocols::rosetta_scripts;
	using core::chemical::DISULFIDE;
	set_last_move_status(protocols::moves::MS_SUCCESS);
	splicemanager.scorefxn(scorefxn());
	//Get entry from dbase file
	find_disulfide_postions(pose,pose_cys_pos_);
	find_vl_vh_cut(pose);
	TR<<"segment type:"<<splicemanager.segment_type()<<std::endl;
	TR<<"tail type:"<<splicemanager.tail_seg()<<std::endl;
	core::Size const dbase_entry(find_dbase_entry(pose));
	TR<<"dbase_entry:"<<dbase_entry<<std::endl;
	//if dbase_entry came back with value of zero it mteans that no pdb entry was found. issue warning nad continue
	if ( !dbase_entry ) {
		TR << TR.Red << "ERROR! Could not entry \""<< database_pdb_entry()<<"\" in database file: \""<<splicemanager.dbase_file_name()<<"\""<<TR.Reset << std::endl;
		utility_exit();
	}
	if ( (splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2") ) {
		database_pdb_entry(torsion_database_[dbase_entry].source_pdb());
		splicemanager.dofs(tail_torsion_database_[dbase_entry]);
		find_vl_vh_cut(pose);
		adjust_n_ter_tail_length(pose);
	}
	//assign the template from_res and to_res according to current antibody segment (L1_L2,L3, H1_H2,H3)
	find_vl_vh_cut(pose);
	// find_disulfide_postions(pose,pose_cys_pos_);
	if ( (splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2") ) {
		splicemanager.tail_seg("");
		assign_from_res_to_res(pose);

		SpliceIn::apply(pose);
		database_pdb_entry(splicemanager.source_pdb_name());
		dbase_iterate(false);//using entry found in the previous splicein
		first_pass_=true;

		splicemanager.tail_seg("n");
		SpliceInTail::apply(pose);
		database_pdb_entry("");//need to re-initialize - otherwise will keep looking for entries
		splicemanager.tail_seg("");//re-initilaize tail segment type
		return;
	} else if ( (splicemanager.segment_type()=="L3")||(splicemanager.segment_type()=="H3") ) {
		if ( splicemanager.segment_type()=="L3" ) {
			splicemanager.pose_to_res((vl_vh_cut_));
			SpliceInTail::apply(pose);
			update_vl_vh_cut();
		} else {
			SpliceInTail::apply(pose);
			splicemanager.tail_seg("");//re-initilaize tail segment type
			return;
		}
	} else {
		utility_exit_with_message("Illegal segment type passed for antibodies\n");
	}

} //apply

std::string SpliceInAntibody::get_name() const {
	return SpliceInAntibodyCreator::mover_name();
}

void SpliceInAntibody::parse_my_tag(TagCOP const tag, basic::datacache::DataMap &data,protocols::filters::Filters_map const & /*filters*/,protocols::moves::Movers_map const &/* movers*/,core::pose::Pose const & /*pose*/) {
	if ( tag->hasOption("segment") ) {
		splicemanager.segment_type(tag->getOption<std::string>("segment"));
	} else {
		TR.Error<<TR.Red<<"Must specify antibody segment (L1_L2,H1_H2,H3,L3)"<<TR.Reset<<std::endl;
		runtime_assert((tag->hasOption("segment")));
	}
	typedef utility::vector1<std::string> StringVec;
	utility::vector1<TagCOP> const sub_tags(tag->getTags());
	splicemanager.parse_tags(tag,data);
	tolerance_ = tag->getOption < core::Real > ("tolerance", 0.23); //for debugging purposes
	splicemanager.allowed_cuts(tag->getOption < core::Size > ("allowed_cuts", 1));
	dbase_iterate(tag->getOption<bool>("dbase_iterate", false));
	min_seg(tag->getOption<bool>("min_seg", true));
	//TR<<"READING TORSION DB "<<splicemanager.dbase_file_name()<<std::endl;

	std::string delta;
	if ( tag->hasOption("delta_lengths") ) {
		delta = tag->getOption<std::string>("delta_lengths");
		//check in put is valid, either : "1,2,3,..." or "1:10", Gideon 19may15
		StringVec const lengths_keys(utility::string_split(delta, ','));
		for ( std::string const delta: lengths_keys ) {
			if ( delta == "" ) continue;
			int const delta_i( 1 * atoi( delta.c_str() ) );
			delta_lengths_.push_back( delta_i );
		}//for
	} else { //fi  (tag->hasOption("delta_lengths"))
		delta_lengths_.push_back(0);
	}
	std::sort(delta_lengths_.begin(), delta_lengths_.end());
	std::unique(delta_lengths_.begin(), delta_lengths_.end());


	//Read map between segments and PSSMs
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
			SpliceSegmentOP splice_segment( new SpliceSegment );
			TR<<segment_type<<std::endl;
			if ( data.has("pssms", segment_type) ) { //If I already read pssms from disk, no need to read them again, GDL nov17
				splice_segment = data.get_ptr<SpliceSegment>("pssms", segment_type);
				splicemanager.splice_segments().insert( std::pair< std::string, SpliceSegmentOP >( segment_type, splice_segment ) );
				continue;
			}



			splice_segment->read_pdb_profile_file( antibody_DB(), segment_type );
			splice_segment->read_many (antibody_DB(),segment_type);
			TR<<splice_segment->get_cluster_name_by_PDB("1AHWD")<<std::endl;
			splicemanager.splice_segments().insert( std::pair< std::string, SpliceSegmentOP >( segment_type, splice_segment ) );
			data.add("pssms", segment_type, splice_segment);
		}
		splicemanager.use_sequence_profiles(true);
		splicemanager.profile_weight_away_from_interface( tag->getOption< core::Real >( "profile_weight_away_from_interface", 1.0 ) );

	}

	scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));



	database_entry(tag->getOption<core::Size>("database_entry", 0));
	database_pdb_entry(tag->getOption<std::string>("database_pdb_entry", ""));
	runtime_assert(!(tag->hasOption("database_entry") && tag->hasOption("database_pdb_entry")));
	read_torsion_database();
	TR << "torsion_database: " << splicemanager.dbase_file_name() << " ";
	if ( database_entry() == 0 ) {
		if ( database_pdb_entry_ == "" ) {
			TR << " database entry will be randomly picked at run time. ";
		} else {
			TR << " picking database entry " << database_pdb_entry() << std::endl;
		}
	} else {
		TR << " database_entry: " << database_entry() << " ";
		runtime_assert(database_entry() <= torsion_database_.size());
	}
}

protocols::moves::MoverOP SpliceInAntibody::clone() const {
	return (protocols::moves::MoverOP(new SpliceInAntibody(*this)));
}



std::string SpliceInAntibody::mover_name()  {
	return "SpliceInAntibody";
}

std::string SpliceInAntibody_complex_type_name_for_subsubtag( std::string const & foo ) {
	return "subsubtag_SpliceInAntibody_" + foo + "_type";
}

std::string SpliceInAntibody_complex_type_name_for_subtag( std::string const & foo ) {
	return "subtag_SpliceInAntibody_" + foo + "_type";
}


void SpliceInAntibody::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;


	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "tolerance", xsct_real, "XRW TO DO", "0.23" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_chain_break", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "min_seg", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "CG_const", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rb_sensitive", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "chain_num", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute( "segment", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "superimposed", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_n", xsct_non_negative_integer, "XRW TO DO", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_c", xsct_non_negative_integer, "XRW TO DO", "13" )
		+ XMLSchemaAttribute( "delta_lengths", xsct_int_cslist, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "dbase_iterate", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "database_entry", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "database_pdb_entry", xs_string, "XRW TO DO" );

	// The "Segments" subtag
	AttributeList segments_subtag_attlist;
	segments_subtag_attlist + XMLSchemaAttribute( "current_segment", xs_string, "XRW TO DO" );

	// The "segment" sub-subtag"
	AttributeList subtag_segments_subtag_attlist;
	subtag_segments_subtag_attlist + XMLSchemaAttribute( "pdb_profile_match", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "profiles", xs_string, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator segment_subsubtag_gen;
	segment_subsubtag_gen.complex_type_naming_func( & SpliceInAntibody_complex_type_name_for_subsubtag )
		.element_name( "segment" )
		.description( "individual segment tag" )
		.add_attributes( subtag_segments_subtag_attlist )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_already_defined_subelement( "segment", SpliceInAntibody_complex_type_name_for_subsubtag/*, 0*/ );

	XMLSchemaComplexTypeGenerator segments_subtag_gen;
	segments_subtag_gen.complex_type_naming_func( & SpliceInAntibody_complex_type_name_for_subtag )
		.element_name( "Segments" )
		.description( "Wrapper for multiple segments tags" )
		.add_attributes( segments_subtag_attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );


	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( "Segments", SpliceInAntibody_complex_type_name_for_subtag/*, 0*/ );


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

void SpliceInAntibodyCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SpliceInAntibody::provide_xml_schema( xsd );
}
void SpliceInAntibody::set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceIn::set_loop_length_change(llc);
	}
	return SpliceInTail::set_loop_length_change(llc);
}

void SpliceInAntibody::set_fold_tree_nodes(core::pose::Pose const & pose){
	fold_tree_nodes_.clear();
	update_vl_vh_cut();
	/*Build the following Fold tree for the Vl/Vh fragment
	___________________________
	|                         |
	|                         |
	<------*<--------*------>//<-----*<-------*------->
	*/
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		SpliceIn::set_fold_tree_nodes(pose);
		return;
	}
	find_disulfide_postions(pose,pose_cys_pos_);
	//if ((splicemanager.segment_type()=="L3"))
	// pose.conformation().delete_chain_ending(pose.conformation().chain_end(1));//Need to revert the scfv back to one chain

	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[1],1,-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[2],(int) pose_cys_pos_[1],-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[2],(int) vl_vh_cut_,-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[3],(int) vl_vh_cut_+1,-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[4],(int) pose_cys_pos_[3],-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[4],(int) pose.conformation().chain_end(1),-1}});
	fold_tree_nodes_.push_back({{(int) pose_cys_pos_[4],(int) pose_cys_pos_[2],1}});



	if ( pose.conformation().num_chains()>1 ) { //if ligand is present we need to add edge between receptor and ligand
		core::Size CoM = (core::Size ) core::pose::residue_center_of_mass( pose, pose.conformation().chain_begin(2), pose.conformation().chain_end(2) );
		fold_tree_nodes_.push_back({{(int)pose_cys_pos_[4],(int) CoM,2}});
		fold_tree_nodes_.push_back({{(int) CoM, (int)pose.conformation().chain_end(2),2}});
		fold_tree_nodes_.push_back({{(int) CoM, (int)pose.conformation().chain_begin(2),2}});
	}
}

void
SpliceInAntibody::find_vl_vh_cut(core::pose::Pose pose)  {
	protocols::simple_moves::CutChainMover ccm;
	vl_vh_cut_ = ccm.chain_cut(pose); //
	TR<<"vl vh cut is: "<<vl_vh_cut_<<std::endl;
}

core::Size SpliceInAntibody::set_anchor_res(){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceIn::set_anchor_res();
	}
	return SpliceInTail::set_anchor_res();
}

core::Size SpliceInAntibody::find_dbase_entry(core::pose::Pose const & pose) {
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceIn::find_dbase_entry(pose);
	}
	return SpliceInTail::find_dbase_entry(pose);
}

void SpliceInAntibody::update_vl_vh_cut() {
	if ( splicemanager.segment_type().find_first_of("L") !=std::string::npos ) {
		vl_vh_cut_+=splicemanager.residue_diff();
	}
}

void SpliceInAntibody::find_disulfide_postions(core::pose::Pose const & pose,utility::vector1<core::Size> & cys_pos) {
	cys_pos.clear();
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue(i).has_variant_type(core::chemical::DISULFIDE) ) {
			cys_pos.push_back(i);
		}
	}
	runtime_assert(cys_pos.size() >3 );
}

void SpliceInAntibody::assign_from_res_to_res(core::pose::Pose const & pose){

	find_disulfide_postions(pose,pose_cys_pos_);
	// find_disulfide_postions(*splicemanager.template_pose(),template_cys_pos_);
	core::conformation::Conformation const & conf(pose.conformation());
	if ( splicemanager.segment_type()=="L1_L2" ) {
		if ( splicemanager.tail_seg()=="n" ) {
			splicemanager.pose_from_res(1);
			splicemanager.pose_to_res(pose_cys_pos_[1]-1);
			return;
		}

		splicemanager.pose_from_res(pose_cys_pos_[1]+1);
		splicemanager.pose_to_res(pose_cys_pos_[2]-1);

		//  splicemanager.template_from_res(template_cys_pos_[1]+1);
		//  splicemanager.template_to_res(template_cys_pos_[2]-1);
	} else if ( splicemanager.segment_type()=="H1_H2" ) {
		if ( splicemanager.tail_seg()=="n" ) {
			splicemanager.pose_from_res(vl_vh_cut_+1);
			splicemanager.pose_to_res(pose_cys_pos_[3]-1);
			return;
		}
		//splicemanager.template_from_res(template_cys_pos_[3]+1);
		//splicemanager.template_to_res(template_cys_pos_[4]-1);
		splicemanager.pose_from_res(pose_cys_pos_[3]+1);
		splicemanager.pose_to_res(pose_cys_pos_[4]-1);
		// splicemanager.template_from_res(template_cys_pos_[3]+1);
		// splicemanager.template_to_res(template_cys_pos_[4]-1);
	} else if ( splicemanager.segment_type()=="L3" ) {
		//splicemanager.template_from_res(template_cys_pos_[2]+1);
		splicemanager.pose_from_res(pose_cys_pos_[2]+1);
		splicemanager.pose_to_res(vl_vh_cut_);
		return;
	} else if ( splicemanager.segment_type()=="H3" ) {
		//splicemanager.template_from_res(template_cys_pos_[4]+1);
		splicemanager.pose_from_res(pose_cys_pos_[4]+1);
		splicemanager.pose_to_res(conf.num_chains()==1?pose.total_residue():conf.chain_end(1));
		return;
	}
	splicemanager.cut_site (splicemanager.dofs().cut_site() - splicemanager.dofs().start_loop() + splicemanager.pose_from_res()-tail_diff_);
	splicemanager.residue_diff(splicemanager.dofs().size() - (splicemanager.pose_to_res() - splicemanager.pose_from_res()) -1);

}

//
void SpliceInAntibody::adjust_template_jump(core::pose::Pose & pose) {
	protocols::protein_interface_design::movers::SetAtomTree AB_fold_tree;
	AB_fold_tree.ab_fold_tree(true);
	AB_fold_tree.chain(1);
	AB_fold_tree.apply(pose);
	TR<<"Pose fold tree before adjust jump: "<<pose.fold_tree();
	AB_fold_tree.apply(*(splicemanager.template_pose()));
	TR<<"template fold tree before adjust jump: "<<pose.fold_tree();
	splicemanager.template_pose()->set_jump( 1,pose.jump( 1 ) );
	splicemanager.template_pose()->dump_pdb("template_after_jump_adjust.pdb");
}

void SpliceInAntibody::build_ideal_segment(core::pose::Pose & pose){
	if ( ((splicemanager.segment_type()=="L1_L2")||(splicemanager.segment_type()=="H1_H2"))&&(splicemanager.tail_seg()=="") ) {
		return SpliceIn::build_ideal_segment(pose);
	}
	// if (splicemanager.segment_type()=="L3"){ //need to remove chain ending before building segment
	//  pose.conformation().delete_chain_ending(pose.conformation().chain_end(1));
	// }
	// SpliceInTail::build_ideal_segment(pose);
	// if (splicemanager.segment_type()=="L3"){
	//  pose.conformation().insert_chain_ending(vl_vh_cut_);
	// }
	return SpliceInTail::build_ideal_segment(pose);
}
void
SpliceInAntibody::adjust_n_ter_tail_length(core::pose::Pose & pose){
	using namespace protocols::rosetta_scripts;
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	int residue_diff = 0;

	if ( splicemanager.segment_type()=="L1_L2" ) {
		residue_diff = splicemanager.dofs().size()-(pose_cys_pos_[1]-1);
		runtime_assert(residue_diff<1000);//sanity check to make sure we are not getting weird residue diffs, gideonla apr18
		llc.loop_start(1);
		llc.loop_end(pose_cys_pos_[1]-1);
	}
	if ( splicemanager.segment_type()=="H1_H2" ) {

		residue_diff =splicemanager.dofs().size()-(pose_cys_pos_[3]-(vl_vh_cut_+1));
		llc.loop_start(vl_vh_cut_+1);
		llc.loop_end(pose_cys_pos_[3]-1);
		TR<<"loop_start:"<<llc.loop_start()<<std::endl;
		TR<<"loop_end:"<<llc.loop_end()<<std::endl;
	}
	TR<<"Tail residue diff: "<<residue_diff<<std::endl;
	tail_diff_=residue_diff;
	llc.delta(residue_diff);
	llc.tail(1);
	llc.direction(1);//N-ter direction
	llc.apply(pose);

}



} //splice
} //protocols
