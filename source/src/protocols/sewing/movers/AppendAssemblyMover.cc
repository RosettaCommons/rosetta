// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/movers/AppendAssemblyMover.cc
/// @brief an AssemblyMover for adding to existing poses
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <protocols/sewing/movers/AppendAssemblyMover.hh>
#include <protocols/sewing/movers/AppendAssemblyMoverCreator.hh>
#include <protocols/sewing/scoring/AssemblyScorerCreator.hh>
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>
#include <protocols/sewing/requirements/AssemblyRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <protocols/sewing/hashing/EdgeMapGenerator.hh>
#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/hashing/BasisMapGenerator.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/moves/mover_schemas.hh>
// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/selection.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <numeric/random/random.hh>
static basic::Tracer TR( "protocols.sewing.movers.AppendAssemblyMover" );

namespace protocols {
namespace sewing {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
AppendAssemblyMover::AppendAssemblyMover():
	AssemblyMover()
{
	set_min_cycles( 10000 );
	set_max_cycles( 100000 );
	set_start_temperature( 0.6 );
	set_end_temperature( 0.6 );
	set_add_probability( 0.05 );
	set_delete_probability( 0.005 );
	set_modifiable_terminus('B');
	set_output_partner(true);
	set_start_node_vital_segments("all");
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
AppendAssemblyMover::AppendAssemblyMover( AppendAssemblyMover const & src ):
	AssemblyMover( src )
{
	if ( src.get_partner_pdb() != nullptr ) {
		partner_pdb_ = src.get_partner_pdb()->clone();
	}
	ligands_ = src.ligands_;
	required_resnums_ = src.get_required_resnums();
	required_selector_ = src.get_required_selector();
	TR << src.class_name() << std::endl;
	segments_from_dssp_ = src.get_segments_from_dssp();
	modifiable_terminus_ = src.get_modifiable_terminus();
	output_partner_ = src.get_output_partner();
	pose_segment_starts_string_ = src.get_pose_segment_starts_string();
	pose_segment_ends_string_ = src.get_pose_segment_ends_string();
	pose_segment_dssp_ = src.get_pose_segment_dssp();
	strict_dssp_changes_ = src.get_strict_dssp_changes();
	extend_mode_ = src.get_extend_mode();
	start_node_vital_segments_ = src.get_start_node_vital_segments();
}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
AppendAssemblyMover::apply( core::pose::Pose& pose){

	//Try to generate an assembly
	core::Size starttime = time(NULL);

	data_storage::SmartAssemblyOP assembly = set_up_assembly( pose );
	assembly->set_modifiable_terminus(this->get_modifiable_terminus());
	assembly->set_output_partner(this->get_output_partner());
	TR <<"Beginning assembly generation" << std::endl;
	TR << "Modifiable terminus set to " << this->get_modifiable_terminus() << std::endl;

	assembly->set_start_node_vital_segments(this->get_start_node_vital_segments());
	TR <<"Setting start node vital segments: " <<this->get_start_node_vital_segments() <<std::endl;

	if ( extend_mode_ ) {
		set_min_cycles(1);
		set_max_cycles(1);
	}
	if ( get_max_cycles() > 0 ) {
		AssemblyMover::generate_assembly(assembly,pose);
	}
	if ( extend_mode_ ) {
		assembly->trim_terminal_segments(this->get_modifiable_terminus(),2);
	}
	//If we failed, set to FAIL_RETRY
	if ( assembly == nullptr ) {
		TR << "Failed to generate an Assembly" << std::endl;
		set_last_move_status(protocols::moves::FAIL_RETRY);
		return;
	}

	core::Size endttime = time(NULL);
	TR << "Assembly successfully generated in " << endttime - starttime << " seconds" << std::endl;
	AssemblyMover::print_statistics( assembly );

	//Now convert the assembly to a pose
	pose = assembly->to_pose("fa_standard");

	add_motif_scorers_to_score_file( pose, assembly );

}


data_storage::SmartAssemblyOP
AppendAssemblyMover::set_up_assembly( core::pose::Pose & pose){

	data_storage::SmartAssemblyOP assembly;
	utility::vector1< data_storage::LigandDescription > expanded_ligands;
	std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
	if ( get_hashed() ) {
		data_storage::HashedSmartAssemblyOP hashed_assembly( new data_storage::HashedSmartAssembly( get_segment_vector() ) );
		hashed_assembly->set_basis_map_generator( get_basis_map_generator() );
		assembly = hashed_assembly;
		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector( pose, partner_pdb_, get_basis_map_generator(), ligands_, partner_ligands, expanded_ligands, required_resnums_, required_selector_, segments_from_dssp_ );
		assembly->pdb_segments( get_basis_map_generator()->pdb_segments() ); //This also adds them to the local segment vector
		assembly->set_partner(partner_pdb_);
		TR << "Set partner pdb" << std::endl;
		assembly->set_partner_ligands( partner_ligands );

		for ( data_storage::LigandDescription ligdes: expanded_ligands ) {
			assembly->load_initial_conformers( ligdes );
		}

	} else {
		assembly = data_storage::SmartAssemblyOP( new data_storage::SmartAssembly( this->get_segment_vector(), this->get_window_width() ) );
		std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs;
		pdbsegs.clear();



		//This version now needs to take pose_segment_starts_ and pose_segment_ends_--put them right after pdbsegs
		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector( pose, partner_pdb_, assembly->get_segment_vector(), pdbsegs, pose_segment_starts_string_, pose_segment_ends_string_, pose_segment_dssp_, ligands_, partner_ligands, expanded_ligands, required_resnums_, required_selector_, strict_dssp_changes_ );
		TR << "PDB segments added" << std::endl;
		assembly->pdb_segments( pdbsegs );
		assembly->set_partner(partner_pdb_);
		TR << "Set partner pdb" << std::endl;
		assembly->set_partner_ligands( partner_ligands );


		for ( data_storage::LigandDescription ligdes: expanded_ligands ) {
			assembly->load_initial_conformers( ligdes );
		}

		//Our starting segment will be one of the segments from our input pose (either the first or the last if it's not hashed

	}
	TR << "Set ligand conformers" << std::endl;
	data_storage::SmartSegmentOP starting_segment = find_starting_segment( assembly );
	assembly->set_start_node_vital_segments(this->get_start_node_vital_segments());
	assembly->set_starting_segment( starting_segment, assembly->get_start_node_vital_segments() );
	TR << "Starting segment set to: " << starting_segment->get_segment_id() << std::endl;
	return assembly;

}


data_storage::SmartSegmentOP
AppendAssemblyMover::find_starting_segment( data_storage::SmartAssemblyOP assembly ){
	data_storage::SmartSegmentOP starting_segment;
	if ( assembly->local_segments().size()<1 ) {
		utility_exit_with_message("local_segments is empty");
	}
	if ( get_hashed() ) {
		starting_segment = assembly->local_segments().at( assembly->get_segment_vector()->size() + numeric::random::rg().random_element( get_basis_map_generator()->alignment_settings().match_segments_  ) );
	} else {
		//We either want to get the N-terminal segment or the C-terminal segment from pdb_segments
		starting_segment = ( numeric::random::rg().uniform() < 0.5 ? assembly->local_segments().at( (assembly->pdb_segments().begin() )->second->get_segment_id() ) : assembly->local_segments().at( (assembly->pdb_segments().rbegin() )->second->get_segment_id() ) );
	}
	return starting_segment;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
AppendAssemblyMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

/// @brief Get the name of the Mover
std::string
AppendAssemblyMover::get_name() const {
	return AppendAssemblyMover::class_name();
}

core::pose::PoseOP
AppendAssemblyMover::get_partner_pdb() const{
	return partner_pdb_;
}

void
AppendAssemblyMover::set_strict_dssp_changes( bool input ){
	strict_dssp_changes_ = input;
}



void
AppendAssemblyMover::set_ligands( utility::vector1< data_storage::LigandDescription > input){
	ligands_ = input;
}

utility::vector1< data_storage::LigandDescription >
AppendAssemblyMover::get_ligands() const{
	return ligands_;
}

utility::vector1< data_storage::LigandDescription > &
AppendAssemblyMover::get_nonconst_ligands(){
	return ligands_;
}


bool
AppendAssemblyMover::get_segments_from_dssp() const{
	return segments_from_dssp_;
}


void
AppendAssemblyMover::set_segments_from_dssp( bool set ){
	segments_from_dssp_ = set;
}


void
AppendAssemblyMover::set_partner_pdb( core::pose::PoseOP partner ){
	partner_pdb_ = partner;
}

std::string
AppendAssemblyMover::get_required_resnums() const{
	return required_resnums_;
}

core::select::residue_selector::ResidueSelectorCOP
AppendAssemblyMover::get_required_selector() const{
	return required_selector_;
}

void
AppendAssemblyMover::set_required_resnums( std::string required ){
	required_resnums_ = required;
}

void
AppendAssemblyMover::set_required_selector( core::select::residue_selector::ResidueSelectorCOP select ){
	required_selector_ = select;
}

void
AppendAssemblyMover::set_modifiable_terminus(char modifiable_terminus){
	modifiable_terminus_ = modifiable_terminus;
}
char
AppendAssemblyMover::get_modifiable_terminus() const{
	return modifiable_terminus_;
}

void
AppendAssemblyMover::set_start_node_vital_segments( std::string start_node_vital_segments ) {
	start_node_vital_segments_ = start_node_vital_segments;
}

std::string
AppendAssemblyMover::get_start_node_vital_segments() const {
	return start_node_vital_segments_;
}

void
AppendAssemblyMover::set_output_partner(bool output_partner){
	output_partner_ = output_partner;
}
bool
AppendAssemblyMover::get_output_partner() const{
	return output_partner_;
}

std::string
AppendAssemblyMover::get_pose_segment_starts_string() const{
	return pose_segment_starts_string_;
}

std::string
AppendAssemblyMover::get_pose_segment_ends_string() const{
	return pose_segment_ends_string_;
}

void
AppendAssemblyMover::set_pose_segment_starts_string( std::string const & input ){
	pose_segment_starts_string_ = input;
}

void
AppendAssemblyMover::set_pose_segment_ends_string( std::string const & input ){
	pose_segment_ends_string_ = input;
}

std::string
AppendAssemblyMover::get_pose_segment_dssp() const{
	return pose_segment_dssp_;
}

bool
AppendAssemblyMover::get_strict_dssp_changes() const{
	return strict_dssp_changes_;
}

void
AppendAssemblyMover::set_pose_segment_dssp( std::string const & input ){
	pose_segment_dssp_ = input;
}

bool
AppendAssemblyMover::get_extend_mode() const{
	return extend_mode_;
}

void
AppendAssemblyMover::set_extend_mode( bool new_mode){
	extend_mode_ = new_mode;
}
////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
AppendAssemblyMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & filtermap,
	protocols::moves::Movers_map const & movermap,
	core::pose::Pose const & pose)
{
	AssemblyMover::parse_my_tag(tag, datamap, filtermap, movermap, pose);

	//get alignment settings
	if ( get_hashed() ) {
		hashing::AlignmentFileGeneratorMover::alignment_settings_from_tag( tag, datamap, filtermap, movermap, pose, get_basis_map_generator() );
	}
	pose_segment_starts_string_ = tag->getOption< std::string >( "pose_segment_starts", "" );
	pose_segment_ends_string_ = tag->getOption< std::string >( "pose_segment_ends", "" );
	pose_segment_dssp_ = tag->getOption< std::string >( "pose_segment_dssp", "" );
	strict_dssp_changes_ = tag->getOption< bool >( "strict_dssp_changes", false );
	//get partner pdb
	if ( tag->hasOption( "partner_pdb" ) ) {
		//pose_from_file returns a PoseOP
		partner_pdb_ = core::import_pose::pose_from_file( tag->getOption< std::string >( "partner_pdb") );
	}

	required_resnums_ = tag->getOption< std::string >( "required_resnums", "");
	//required_starting_residues_ = core::pose::get_resnum_list_ordered( str_required_res, pose );
	if ( tag->hasOption( "required_selector" ) ) {
		required_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "required_selector" ), datamap );
	}
	segments_from_dssp_ = tag->getOption< bool >( "set_segments_from_dssp", false );
	std::map< std::string, std::string > pdb_conformers;
	if ( tag->hasTag( "Ligands" ) ) {
		hashing::AlignmentFileGeneratorMover::parse_ligands_tag( tag->getTag( "Ligands" ), datamap, ligands_);
		hashing::AlignmentFileGeneratorMover::parse_ligand_conformers( ligands_ ); //This loads each of the ligand conformers as a ResidueCOP
	}
	if ( tag->hasOption( "modifiable_terminus" ) ) {
		modifiable_terminus_ = tag->getOption<std::string>("modifiable_terminus").at(0);
		TR << "Modifiable terminus set to " << modifiable_terminus_ << std::endl;
	}
	if ( tag->hasOption( "output_partner" ) ) {
		output_partner_ = tag->getOption<bool>("output_partner", true );
	}
	if ( tag->hasOption( "extend_mode" ) ) {
		extend_mode_ = tag->getOption<bool>("extend_mode", false );
	}

	if ( tag->hasOption( "start_node_vital_segments" ) ) {
		start_node_vital_segments_ = tag->getOption<std::string>("start_node_vital_segments", "all");
	}

}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AppendAssemblyMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AppendAssemblyMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AppendAssemblyMover::clone() const
{
	return protocols::moves::MoverOP( new AppendAssemblyMover( *this ) );
}

std::string
AppendAssemblyMover::class_name()
{
	return "AppendAssemblyMover";
}



////////////////////////////////////////////////////////////////////////////////
/// Creator ///
///////////////

/////////////// Creator ///////////////

protocols::moves::MoverOP
AppendAssemblyMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AppendAssemblyMover );
}

std::string
AppendAssemblyMoverCreator::keyname() const
{
	return AppendAssemblyMover::class_name();
}
std::string
AppendAssemblyMoverCreator::class_name()
{
	return AppendAssemblyMover::class_name();
}
std::string
AppendAssemblyMoverCreator::mover_name()
{
	return AppendAssemblyMover::class_name();
}


void
AppendAssemblyMover::attributes_for_append_assembly_mover( utility::tag::AttributeList & attributes ){
	using namespace utility::tag;

	AssemblyMover::define_generic_assembly_mover_attributes( attributes );
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "recursive_depth", xsct_non_negative_integer,  "How many nodes after the terminal node should we keep track of alignments for?", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "pose_segment_starts", xsct_int_cslist, "Residue numbers of the first residue in each segment in the input pose", ""  )
		+ XMLSchemaAttribute::attribute_w_default( "pose_segment_ends", xsct_int_cslist,  "Residue numbers of the last residue in each segment in the input pose. Length must match that of pose_segment_starts.", "" )
		+ XMLSchemaAttribute::attribute_w_default( "pose_segment_dssp", xs_string, "String indicating the secondary structure of user-specified segments, one character per segment (e.g. HLH for a helix-loop-helix motif). Length should match that of pose_segment_starts and pose_segment_ends if specified.", "" )
		+ XMLSchemaAttribute::attribute_w_default( "strict_dssp_changes", xsct_rosetta_bool, "Segments require at least a 2-residue change in DSSP to specify a new segment", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "set_segments_from_dssp", xsct_rosetta_bool, "Determine segment boundaries based on pose secondary structure", "false" )
		+ XMLSchemaAttribute( "match_segments", xsct_int_cslist,  "Which segments from the input pose should we be able to append onto? Defaults to exterior segments." )
		+ XMLSchemaAttribute( "partner_pdb", xs_string, "Name of PDB file containing binding partner for this assembly" )
		+ XMLSchemaAttribute( "required_resnums", xsct_refpose_enabled_residue_number_cslist, "Residue numbers of residues in the input structure that must be preserved" )
		+ XMLSchemaAttribute::attribute_w_default( "max_recursion", xsct_non_negative_integer, "How many alignments from the end nodes should be stored in memory?", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "modifiable_terminus", xs_string, "Which terminus of the starting node may be modified.","B")
		+ XMLSchemaAttribute::attribute_w_default( "output_partner", xsct_rosetta_bool, "Should the output pdb contain the partner?","true" )
		+ XMLSchemaAttribute::attribute_w_default( "extend_mode", xsct_rosetta_bool, "Should SEWING append only a single helix?","false" )
		+ XMLSchemaAttribute::attribute_w_default( "start_node_vital_segments", xs_string, "Which segments from starting node are vital? (terminal or all)", "all");
	core::select::residue_selector::attributes_for_parse_residue_selector( attributes, "required_selector",  "Residue selector specifying residues in the input structure that must be preserved" );
}




void
AppendAssemblyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	scoring::AssemblyScorerFactory::get_instance()->define_assembly_scorer_subtag( xsd );
	requirements::AssemblyRequirementFactory::get_instance()->define_assembly_requirement_subtag( xsd );


	//Define the Ligands subelement
	//Since this one has not just extra attributes but also extra subelements, I'll just get AssemblyMover's attributes and not its ct_gen.
	XMLSchemaSimpleSubelementList mover_subelements;
	mover_subelements
		.add_already_defined_subelement( "AssemblyScorers", & AssemblyMover::assembly_mover_subtag_ct_namer )
		.add_already_defined_subelement( "AssemblyRequirements", & AssemblyMover::assembly_mover_subtag_ct_namer );
	hashing::AlignmentFileGeneratorMover::append_ligands_subelement( mover_subelements, xsd, true );

	//Add attributes
	AttributeList attributes;
	attributes_for_append_assembly_mover( attributes );

	//Define complex type for the overall mover
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.element_name( AppendAssemblyMoverCreator::mover_name() )
		.add_attributes( attributes )
		.description( "Builds an assembly around the segment provided in the input PDB file" )
		.add_optional_name_attribute( "Name to identify this mover" )
		.set_subelements_single_appearance_optional( mover_subelements )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}

void
AppendAssemblyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{

	AppendAssemblyMover::provide_xml_schema( xsd );

}



////////////////////////////////////////////////////////////////////////////////
/// private methods ///
//////////////////////

std::ostream &
operator<<( std::ostream & os, AppendAssemblyMover const & mover )
{
	mover.show( os );
	return os;
}

} //protocols
} //sewing
} //movers

