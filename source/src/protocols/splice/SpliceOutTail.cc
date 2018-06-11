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
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/splice/SpliceOutCreator.hh>
#include <protocols/splice/SpliceOutTailCreator.hh>
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
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>

namespace protocols {
namespace splice {

using namespace core::conformation;

static  basic::Tracer TR( "protocols.splice.SpliceOutTail" );
static  basic::Tracer TR_ccd( "protocols.splice.Splice_ccd" );
static  basic::Tracer TR_constraints( "protocols.splice.Splice_constraints" );
static  basic::Tracer TR_pssm( "protocols.splice.Splice_pssm" );
static  basic::Tracer TR_min( "protocols.splice.Splice_min" );

std::string SpliceOutTailCreator::keyname() const {
	return SpliceOutTailCreator::mover_name();
}

protocols::moves::MoverOP SpliceOutTailCreator::create_mover() const {
	return protocols::moves::MoverOP( new SpliceOutTail );
}

std::string SpliceOutTailCreator::mover_name() {
	return "SpliceOutTail";
}


SpliceOutTail::SpliceOutTail() :  Mover(SpliceOutTailCreator::mover_name())
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
	splicemanager.adjust_stem_positions_by_template_=true;
	ignore_chain_break_  = false;
}


SpliceOutTail::~SpliceOutTail() {}



std::string SpliceOutTail::get_name() const {
	return SpliceOutTailCreator::mover_name();
}

void SpliceOutTail::abstract_parse_tag(TagCOP const tag){
	splicemanager.tail_seg(tag->getOption<std::string>("tail_segment", ""));
	if ( (!(boost::iequals(splicemanager.tail_seg(), "c")))&&(!(boost::iequals(splicemanager.tail_seg(), "n"))) ) {
		utility_exit_with_message("\"tail_segment\" tag must get either \"c\" or \"n\" \n");
	}

	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		splicemanager.template_to_res(core::pose::parse_resnum(tag->getOption<std::string>("from_res", "0"), *(splicemanager.template_pose())));
		splicemanager.template_from_res(1);
		TR<<splicemanager.template_to_res()<<std::endl;
	} else /*assming this is a "c" ter tail */{
		core::conformation::Conformation const & conf(splicemanager.template_pose()->conformation());
		splicemanager.template_from_res(core::pose::parse_resnum(tag->getOption<std::string>("from_res", "0"), *(splicemanager.template_pose())));
		splicemanager.template_to_res(conf.chain_end(1));
	}
}





protocols::moves::MoverOP SpliceOutTail::clone() const {
	return (protocols::moves::MoverOP(new SpliceOutTail(*this)));
}



void SpliceOutTail::set_source_from_to_res(){

	using namespace protocols::rosetta_scripts;
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		source_to_res(find_nearest_res(*source_pose_, *splicemanager.template_pose(), splicemanager.template_to_res(), 0/*chain*/));
		if ( source_to_res()==0 ) {
			utility_exit_with_message("Could not find source to_res for n-ter tail segment, check alignmnet between template pose and source pose\n");
		}
		source_from_res(1);
	} else if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		core::conformation::Conformation const & conf(source_pose_->conformation());
		source_from_res(find_nearest_res(*source_pose_, *splicemanager.template_pose(), splicemanager.template_from_res(), 0/*chain*/));
		source_to_res(conf.chain_end(1));
	} else {
		utility_exit_with_message("Tail segment should be noted as \"c\" on \"n\", check XML\n");
	}
}


core::Size SpliceOutTail::set_anchor_res(){
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		return splicemanager.pose_to_res()+1;
	} else if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		return splicemanager.pose_from_res()-1;
	}
	return 0;
}

void SpliceOutTail::set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc){
	llc.loop_start(splicemanager.pose_from_res());
	llc.loop_end(splicemanager.pose_to_res());
	llc.tail(1);
	llc.delta(splicemanager.residue_diff());
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		llc.direction(1);
	} else if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		llc.direction(0);
	}
}

void SpliceOutTail::set_fold_tree_nodes(core::pose::Pose const & pose){
	core::conformation::Conformation const & conf(pose.conformation());
	TR<<conf.num_chains()<<std::endl;
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		fold_tree_nodes_.push_back({{(int) splicemanager.pose_to_res()+1,(int) splicemanager.pose_from_res(),-1}});
		fold_tree_nodes_.push_back({{(int) splicemanager.pose_to_res()+1,(int)conf.chain_end(1),-1}});
	} else if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		fold_tree_nodes_.push_back({{(int) splicemanager.pose_from_res()-1,(int) splicemanager.pose_to_res(),-1}});
		fold_tree_nodes_.push_back({{(int) splicemanager.pose_from_res()-1,1,-1}});
	}
	if ( conf.num_chains()>1 ) { //if ligand is present we need to add edge between receptor and ligand
		core::Size CoM = (core::Size ) core::pose::residue_center_of_mass( pose, conf.chain_begin(2), conf.chain_end(2) );
		fold_tree_nodes_.push_back({{(int) conf.chain_end(1),(int) CoM,1}});
	}
}
void SpliceOutTail::place_cut_site_in_segment(core::pose::Pose const & pose){
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		splicemanager.cut_site(pose.total_residue()+splicemanager.residue_diff()-1);
	} else if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		splicemanager.cut_site(1);
	}
}
void SpliceOutTailCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SpliceOutTail::provide_xml_schema( xsd );
}

std::string
SpliceOutTail::mover_name()  {
	return "SpliceOutTail";
}

std::string SpliceOutTail_complex_type_name_for_subsubtag( std::string const & foo ) {
	return "subsubtag_spliceoutTail_" + foo + "_type";
}

std::string SpliceOutTail_complex_type_name_for_subtag( std::string const & foo ) {
	return "subtag_spliceoutTail_" + foo + "_type";
}

void SpliceOutTail::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction n_or_c;
	n_or_c.name( "n_or_c" );
	n_or_c.base_type( xs_string );
	n_or_c.add_restriction( xsr_enumeration, "n" );
	n_or_c.add_restriction( xsr_enumeration, "c" );
	xsd.add_top_level_element( n_or_c );

	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "CG_const", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rb_sensitive", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "chain_num", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute( "segment", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "superimposed", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_n", xsct_non_negative_integer, "XRW TO DO", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_c", xsct_non_negative_integer, "XRW TO DO", "13" )
		+ XMLSchemaAttribute( "tail_segment", "n_or_c", "XRW TO DO" );

	// The "Segments" subtag
	AttributeList segments_subtag_attlist;
	segments_subtag_attlist + XMLSchemaAttribute( "current_segment", xs_string, "XRW TO DO" );

	// The "segment" sub-subtag"
	AttributeList subtag_segments_subtag_attlist;
	subtag_segments_subtag_attlist + XMLSchemaAttribute( "pdb_profile_match", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "profiles", xs_string, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator segment_subsubtag_gen;
	segment_subsubtag_gen.complex_type_naming_func( & SpliceOutTail_complex_type_name_for_subsubtag )
		.element_name( "segment" )
		.description( "individual segment tag" )
		.add_attributes( subtag_segments_subtag_attlist )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_already_defined_subelement( "segment", SpliceOutTail_complex_type_name_for_subsubtag/*, 0*/ );

	XMLSchemaComplexTypeGenerator segments_subtag_gen;
	segments_subtag_gen.complex_type_naming_func( & SpliceOutTail_complex_type_name_for_subtag )
		.element_name( "Segments" )
		.description( "Wrapper for multiple segments tags" )
		.add_attributes( segments_subtag_attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );


	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( "Segments", SpliceOutTail_complex_type_name_for_subtag/*, 0*/ );

	attlist + XMLSchemaAttribute( "use_sequence_profile", xsct_rosetta_bool, "XRW TO DO" );
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "add_sequence_constraints_only", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "template_file", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "set_fold_tree_only", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "source_pdb", xs_string, "XRW TO DO");
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "from_res", xsct_refpose_enabled_residue_number, "XRW TO DO", "0" )
		+ XMLSchemaAttribute( "design_task_operations", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "residue_numbers_setter", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "torsion_database", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "design_shell", xsct_real, "XRW TO DO", "6.0" )
		+ XMLSchemaAttribute::attribute_w_default( "repack_shell", xsct_real, "XRW TO DO", "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "rms_cutoff", xsct_real, "XRW TO DO", "999999" )
		+ XMLSchemaAttribute::attribute_w_default( "res_move", xsct_non_negative_integer, "XRW TO DO", "1000" )
		+ XMLSchemaAttribute::attribute_w_default( "cut_secondarystruc", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_ala", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "design", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_original_sequence", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rtmin", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "allow_all_aa", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "checkpointing_file", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "splice_filter", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "mover", xs_string, "Which mover to use to close the segment" )
		+ XMLSchemaAttribute::attribute_w_default( "restrict_to_repacking_chain2", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "use_sequence_profiles", xsct_rosetta_bool, "XRW TO DO", "true" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, subelements );
}

void
SpliceOutTail::write_database_to_file(core::pose::Pose const & pose){
	if ( !(write_to_database_) ) {
		return;
	}
	TR<<"Writing to tail database"<<std::endl;
	TR<<"From res: "<<splicemanager.pose_from_res()<<std::endl;
	TR<<"To res: "<<splicemanager.pose_to_res()<<std::endl;
	std::ofstream dbase_file;
	std::stringstream torsion_stringstream; // We cannot write directly into the dbase_file, as the assert statement might stop the run, and thus leave the db poluted, Chris Sep23
	dbase_file.open(splicemanager.dbase_file_name().c_str(), std::ios::app);
	for ( core::Size i = splicemanager.pose_from_res(); i <= splicemanager.pose_to_res(); ++i ) {
		if ( i != splicemanager.pose_to_res()&&i != splicemanager.pose_from_res() ) { // don't do the check for the last residue, which would have psi and omega equals 0 for tail segments
			//TR<<"Checking residue:"<<i<<std::endl;
			runtime_assert( pose.phi(i) && pose.psi(i) && pose.omega(i) );// Make sure that all dihedral angles have non zero values,
		}
		//0 value dihedral angle is probably an error, Gideon Sep14
		torsion_stringstream << pose.phi(i) << ' ' << pose.psi(i) << ' ' << pose.omega(i) << ' ' << pose.residue(i).name3() << ' ';
	}
	core::Size tail_int=0;

	if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		tail_int=1;
	}
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		tail_int=0;
	}
	torsion_stringstream << splicemanager.template_from_res() << ' '<< "0 " << tail_int << ' ';
	torsion_stringstream << splicemanager.source_pdb_name()<<"\n";
	dbase_file << torsion_stringstream.str();
	dbase_file.close();
}


void SpliceOutTail::superimpose_source_on_pose( core::pose::Pose const & pose, core::pose::Pose & source_pose ){
	using namespace protocols::toolbox;
	using namespace protocols::rosetta_scripts;

	core::Size source_anchor(find_nearest_res(source_pose, pose, splicemanager.anchor_res(), 0/*chain*/));
	core::Size template_anchor(find_nearest_res(*splicemanager.template_pose(), pose, splicemanager.anchor_res(), 0/*chain*/));
	if ( superimposed() ) { //Structures are superimposed. Nothing to do
		return;
	}
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		source_anchor = source_anchor -2;
		template_anchor = template_anchor-2;
	}
	utility::vector1< core::Size > pose_positions, template_positions;
	pose_positions.clear(); template_positions.clear();
	pose_positions.push_back( source_anchor );
	pose_positions.push_back( source_anchor +1  );
	pose_positions.push_back(  source_anchor +2 );
	template_positions.push_back( template_anchor );
	template_positions.push_back( template_anchor +1 );
	template_positions.push_back( template_anchor +2 );
	TR<<" template scafold_res: "<<splicemanager.template_from_res() -1<<",source scafold_res: "<< source_from_res() -1<<std::endl;
	utility::vector1< numeric::xyzVector< core::Real > > init_coords( coords( source_pose, pose_positions ) ), ref_coords( coords( pose/*this is the starting pose, the structure on which we want to graft the loop from the source protein*/, template_positions ));
	TR<<"template ref coords: "<<ref_coords[1][1]<<std::endl;
	TR<<"source ref coords: "<<init_coords[1][1]<<std::endl;

	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;

	superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );

	apply_superposition_transform( source_pose, rotation, to_init_center, to_fit_center );
	/// DEBUGGING
	pose.dump_pdb("pose_pdb.pdb");
	source_pose.dump_pdb("source_pose_pdb.pdb");
}

void SpliceOutTail::build_ideal_segment(core::pose::Pose & pose){
	using namespace core::chemical;
	using namespace core::conformation;
	ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose() );
	TR<<"build_ideal_segment:"<<splicemanager.pose_from_res()<<std::endl;
	TR<<"build_ideal_segment:"<<splicemanager.pose_to_res()<<std::endl;



	ResidueCOP new_res = ResidueFactory::create_residue( residue_set->name_map( name_from_aa( aa_from_oneletter_code( 'A' ) ) ) );
	core::Size new_segment_length(splicemanager.pose_to_res()-splicemanager.pose_from_res());

	if ( (boost::iequals(splicemanager.tail_seg(), "n")) ) {
		core::Size new_to_res = splicemanager.pose_to_res()-(splicemanager.pose_to_res()-splicemanager.pose_from_res());
		pose.delete_residue_range_slow( splicemanager.pose_from_res(), splicemanager.pose_to_res()-1);
		//pose.dump_pdb("after_delete_segment.pdb");
		for ( core::Size res=1; res<=new_segment_length; res++ ) {
			pose.conformation().safely_prepend_polymer_residue_before_seqpos(*new_res,new_to_res,true);
			//TR<<"Appending residue before residue: "<<res+new_to_res-1<<std::endl;
			//pose.dump_pdb("build_ideal_segment.pdb");
		}
	} else if ( (boost::iequals(splicemanager.tail_seg(), "c")) ) {
		pose.delete_residue_range_slow( splicemanager.pose_from_res()+1, splicemanager.pose_to_res());
		//pose.dump_pdb("after_delete_segment.pdb");
		for ( core::Size res=splicemanager.pose_from_res()+1; res<=splicemanager.pose_to_res(); res++ ) {
			pose.append_polymer_residue_after_seqpos(*new_res,res-1,true);
			if ( splicemanager.debug() ) {
				pose.dump_pdb(splicemanager.mover_name() +"build_ideal_segment.pdb");
			}
		}
	}

}
std::string SpliceOutTail::name_for_filter(){
	return "SpliceOutTail";
}


} //splice
} //protocols
