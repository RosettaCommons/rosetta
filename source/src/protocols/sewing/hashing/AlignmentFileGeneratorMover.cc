// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/AlignmentFileGeneratorMover.cc
/// @brief Given a model file, edge file ,and one or more input structures, generates alignment files for use with AppendAssemblyMover.
/// @author guffysl (guffy@email.unc.edu)
/// @author minniel (minnie@email.unc.edu)

#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/hashing/AlignmentFileGeneratorMoverCreator.hh>
#include <protocols/sewing/hashing/AlignmentGenerator.hh>
#include <protocols/sewing/hashing/Hasher.hh>
#include <protocols/sewing/hashing/EdgeMapGenerator.hh>
#include <protocols/sewing/hashing/BasisMapGenerator.hh>
#include <protocols/sewing/hashing/LigandBindingResPlacer.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/LigandSegment.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/sewing/movers/AssemblyMover.hh>
//Temporary, functions will move later
//#include <protocols/sewing/util/io.hh >
#include <protocols/moves/mover_schemas.hh>
#include <core/util/metalloproteins_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/BondedResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/selection.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/AtomID.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/io/util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

static basic::Tracer TR( "protocols.sewing.AlignmentFileGeneratorMover" );


namespace protocols {
namespace sewing {
namespace hashing {

AlignmentFileGeneratorMover::AlignmentFileGeneratorMover():
	protocols::moves::Mover( "AlignmentFileGeneratorMover" )
{
	// initialize_from_options();
	basis_map_generator_ = BasisMapGeneratorOP( new BasisMapGenerator );
}


AlignmentFileGeneratorMover::AlignmentFileGeneratorMover( BasisMapGeneratorOP bmg ):
	protocols::moves::Mover( "AlignmentFileGeneratorMover" ),
	basis_map_generator_( bmg )
{
}


AlignmentFileGeneratorMover::AlignmentFileGeneratorMover(
	std::string model_file_name,
	std::string edge_file_name,
	utility::vector1< core::Size > match_segments,
	utility::vector1< core::Size > pose_segment_starts,
	utility::vector1< core::Size > pose_segment_ends,
	core::Size recursive_depth )
{

	basis_map_generator_ = BasisMapGeneratorOP( new BasisMapGenerator( model_file_name, edge_file_name, match_segments, pose_segment_starts, pose_segment_ends, recursive_depth) );

}

AlignmentFileGeneratorMover::AlignmentFileGeneratorMover(
	EdgeMapGeneratorOP edge_file_reader,
	std::string model_file_name,
	std::string alignment_file_name )
{
	basis_map_generator_ = BasisMapGeneratorOP(  new BasisMapGenerator( edge_file_reader, model_file_name, alignment_file_name ) );
}

AlignmentFileGeneratorMover::AlignmentFileGeneratorMover(
	EdgeMapGeneratorOP edge_file_reader,
	std::string model_file_name )
{
	basis_map_generator_= BasisMapGeneratorOP( new BasisMapGenerator( edge_file_reader, model_file_name ) );
}

AlignmentFileGeneratorMover::~AlignmentFileGeneratorMover(){}

AlignmentFileGeneratorMover::AlignmentFileGeneratorMover( AlignmentFileGeneratorMover const & src ):
	protocols::moves::Mover( src )
{
	basis_map_generator_ = src.basis_map_generator();
}

void
AlignmentFileGeneratorMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & filtermap,
	protocols::moves::Movers_map const & movermap,
	core::pose::Pose const & pose)
{
	//NOTE this currently allows users to specify options either in tag or on command line, but we may remove the command-line options at some point
	//model file
	if ( tag->hasOption( "model_file_name" ) ) {
		basis_map_generator_->set_model_file( tag->getOption<std::string>( "model_file_name", "" ) );
	} else {
		basis_map_generator_->model_file_from_options();
	}
	//edge file
	if ( tag->hasOption( "edge_file_name" ) ) {
		basis_map_generator_->set_edge_file( tag->getOption< std::string>( "edge_file_name", "" ) );
	} else {
		basis_map_generator_->edge_file_from_options();
	}
	if ( tag->hasTag( "Ligands" ) ) {
		parse_ligands_tag( tag->getTag( "Ligands" ), datamap, ligands_ );
		//We won't parse ideal contacts for this mover
		//It seems unnecessary to parse conformers, too
	}
	//Now get the information for the alignment settings
	alignment_settings_from_tag( tag, datamap, filtermap, movermap, pose, basis_map_generator_ );
	//Not currently compatible with providing these as command-line options
	required_resnums_ = tag->getOption< std::string >( "required_resnums", "");
	//required_starting_residues_ = core::pose::get_resnum_list_ordered( str_required_res, pose );
	if ( tag->hasOption( "required_selector" ) ) {
		required_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "required_selector" ), datamap );
	}
	set_segments_from_dssp_ = tag->getOption< bool >( "set_segments_from_dssp", false );
}

void
AlignmentFileGeneratorMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	AlignmentFileGeneratorMover::provide_xml_schema( xsd );
}

void
AlignmentFileGeneratorMover::parse_ligands_tag(
	utility::tag::TagCOP ligands_tag,
	basic::datacache::DataMap const & data,
	utility::vector1< data_storage::LigandDescription > & ligands
)
{
	//LigandContact( segid, resnum, resat, ligat )
	runtime_assert( ligands_tag->getName() == "Ligands" );
	core::Size counter = 1;
	for ( utility::tag::TagCOP ligand_tag: ligands_tag->getTags() ) {

		data_storage::LigandDescription new_ligand;

		new_ligand.partner_ligand = ligand_tag->getOption< bool >( "partner_ligand", false );
		if ( ligand_tag->hasOption( "ligand_resnum" ) ) {
			new_ligand.ligand_resnum_string = ligand_tag->getOption< std::string >( "ligand_resnum" );
			//We'll parse this later with core::pose::parse_resnum( string, pose, true)
		} else if ( ligand_tag->hasOption( "ligand_selector" ) ) {
			TR << "Selecting ligand(s) by residue selector" << std::endl;
			new_ligand.ligand_resnum_string = "selector_" + utility::to_string( counter );
			new_ligand.ligand_selector = core::select::residue_selector::get_residue_selector( ligand_tag->getOption< std::string >( "ligand_selector" ), data );
		} else {
			utility_exit_with_message( "You must specify either a residue number or a residue selector for your ligand!" );
		}
		new_ligand.auto_detect_contacts = ligand_tag->getOption< bool >( "auto_detect_contacts", true  );


		//Add the pdb_conformers
		if ( ligands_tag->hasOption( "pdb_conformers" ) ) {
			new_ligand.pdb_conformers_string = ligand_tag->getOption< std::string >( "pdb_conformers" );
			if ( !ligand_tag->hasOption( "alignment_atoms" ) ) {
				utility_exit_with_message( "If you provide PDB conformers for your ligand, you must also provide at least three atom names for alignment!" );
			}
			std::string alnatom_string = ligand_tag->getOption< std::string >( "alignment_atoms" );
			utility::vector1< std::string > alnatom_names = utility::string_split( alnatom_string, ',' );
			if ( alnatom_names.size() < 3 ) {
				utility_exit_with_message( "You must provide at least three atom names for alignment (separated by commas )!" );
			}
			new_ligand.alignment_atoms_str = alnatom_names;
		}

		//Loop through contact tags
		for ( utility::tag::TagCOP contact_tag: ligand_tag->getTags() ) {
			TR << "Found contact tags for ligand! If multiple ligands were specified using the same residue selector, contact tags are not recommended as they will likely not result in the correct behavior. Contacts will be added to the first ligand found." << std::endl;
			//Can have any number of Contact tags and any number of IdealCoordination tags
			if ( !( contact_tag->getName() ==  "Contact" ) ) {
				continue;
			}
			bool partner_contact = contact_tag->getOption< bool >( "partner_contact", false );
			//ligat = pose.residue( ligand_residues[ 1 ] ).atom_index( contact_tag->getOption< std::string >( "ligand_atom_name" ) );
			std::string ligat = contact_tag->getOption< std::string >( "ligand_atom_name", "" );
			std::string contact_resnum = contact_tag->getOption< std::string >( "contact_resnum" );
			//std::string contact_resnum = core::pose::get_resnum( contact_tag, pose, "contact_" );
			std::string resat = contact_tag->getOption< std::string >( "contact_atom_name", "" );


			data_storage::ContactDescription new_contact;
			new_contact.partner_contact = partner_contact;
			//new_contact.partner_ligand = partner_ligand;
			//new_contact.ligand_resnum_string=ligand_resnum;
			new_contact.contact_resnum_string = contact_resnum;
			new_contact.ligand_atom_name=ligat;
			new_contact.contact_atom_name=resat;
			new_ligand.ligand_contacts.push_back( new_contact );

		}
		ligands.push_back( new_ligand );
		++counter;
	}
}

std::string
AlignmentFileGeneratorMover::ligands_subtag_ct_namer( std::string tag_name ){
	return "ligands_subtag_" + tag_name + "_complex_type";
}

void
AlignmentFileGeneratorMover::append_ligands_subelement( utility::tag::XMLSchemaSimpleSubelementList & subs, utility::tag::XMLSchemaDefinition & xsd, bool include_coord ){
	using namespace utility::tag;
	AttributeList ligand_attributes;
	ligand_attributes
		+ XMLSchemaAttribute::attribute_w_default( "partner_ligand", xsct_rosetta_bool, "Is this ligand found in the partner PDB?", "false" )
		+ XMLSchemaAttribute( "pdb_conformers", xs_string, "Name of file containing a list of PDBs (or other Rosetta-compatible input files) containing alternate ligand conformations to sample" )
		+ XMLSchemaAttribute( "alignment_atoms", xs_string, "Comma-separated list of atom names to use when aligning ligand conformers to one another" )
		+ XMLSchemaAttribute::attribute_w_default( "auto_detect_contacts", xsct_rosetta_bool, "Should we automatically detect contacts that are joined to the ligand by inter-residue chemical bonds?", "true" )
		+ XMLSchemaAttribute( "ligand_resnum", xsct_refpose_enabled_residue_number, "Residue number of ligand in either PDB or Rosetta numbering");
	// + optional_name_attribute( "Optional name to identify this ligand" );
	core::select::residue_selector::attributes_for_parse_residue_selector( ligand_attributes, "ligand_selector", "Residue selector indicating ligand(s) covered in this tag" );
	AttributeList contact_attributes;
	//core::pose::attributes_for_get_resnum( contact_attributes, "contact_" );
	contact_attributes
		+ XMLSchemaAttribute::attribute_w_default( "partner_contact", xsct_rosetta_bool, "Does this tag specify a contact with the partner PDB?", "false" )
		+ XMLSchemaAttribute::required_attribute( "contact_resnum", xsct_refpose_enabled_residue_number, "Number of residue participating in this contact in PDB or Rosetta numbering" )
		+ XMLSchemaAttribute( "ligand_atom_name", xs_string, "Rosetta name for the ligand atom participating in the contact" )
		+ XMLSchemaAttribute( "contact_atom_name", xs_string, "Rosetta name for the protein atom participating in the contact" );
	XMLSchemaSimpleSubelementList ligand_subelements;
	ligand_subelements
		.add_simple_subelement( "Contact", contact_attributes, "Describes a contact between the ligand and the input pose" );
	if ( include_coord ) {
		LigandBindingResPlacer::define_ideal_contacts_subelement( ligand_subelements, xsd );
	}
	XMLSchemaComplexTypeGenerator ligand_ct_gen;
	ligand_ct_gen.complex_type_naming_func( & ligands_subtag_ct_namer )
		.element_name( "Ligand" )
		.add_attributes( ligand_attributes )
		.set_subelements_repeatable( ligand_subelements )
		.description( "Specifies the position of a ligand and the contacts that it forms with the input pose" )
		.write_complex_type_to_schema( xsd );


	XMLSchemaSimpleSubelementList ligands_subelements;
	ligands_subelements
		//.add_simple_subelement( "Ligand", ligand_attributes, "Specifies the position of a ligand and the protein residues in the pose with which it forms contacts" );
		.add_already_defined_subelement( "Ligand", & ligands_subtag_ct_namer );

	XMLSchemaComplexTypeGenerator ligands_ct_gen;
	ligands_ct_gen.complex_type_naming_func( & movers::AssemblyMover::assembly_mover_subtag_ct_namer )
		.element_name( "Ligands" )
		.set_subelements_repeatable( ligands_subelements )
		.description( "Subtags of this tag specify the ligands present in the input pose and their respective protein contacts." )
		.write_complex_type_to_schema( xsd );
	subs
		.add_already_defined_subelement( "Ligands", & movers::AssemblyMover::assembly_mover_subtag_ct_namer );
}
void
AlignmentFileGeneratorMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::required_attribute( "model_file_name", xs_string, "Name of input segment file for alignment file generation" )
		+ XMLSchemaAttribute::required_attribute( "edge_file_name", xs_string, "Name of input edge file for alignment file generation" )
		+ XMLSchemaAttribute::attribute_w_default( "recursive_depth", xsct_non_negative_integer, "How many nodes away from the starting segment for which to generate alignments", "1" )
		//comma-separated list of core::Size
		+ XMLSchemaAttribute::attribute_w_default( "match_segments", xsct_int_cslist, "Which segments in the input pose should be able to form chimaera?", "1" )
		//comma-separated lists of residue numbers
		+ XMLSchemaAttribute( "pose_segment_starts", xsct_int_cslist, "Comma-separated list of the first residue of each segment in the pose" )
		+ XMLSchemaAttribute( "required_resnums", xsct_refpose_enabled_residue_number_cslist, "Residue numbers of residues in the input structure that must be preserved" )
		+ XMLSchemaAttribute::attribute_w_default( "set_segments_from_dssp", xsct_rosetta_bool, "Determine segment boundaries based on pose secondary structure", "false" )
		+ XMLSchemaAttribute( "pose_segment_ends", xsct_int_cslist, "Comma-separated list of the last residue of each segment in the pose" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attributes, "required_selector",  "Residue selector specifying residues in the input structure that must be preserved" );
	//No need for score function or task operations with this one
	//No subelements
	//mover xsd_type_definition_w_attributes
	XMLSchemaSimpleSubelementList subelements;
	append_ligands_subelement( subelements, xsd, true );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen
		.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.element_name( AlignmentFileGeneratorMoverCreator::mover_name() )
		.set_subelements_single_appearance_optional( subelements )
		.add_attributes( attributes )
		.write_complex_type_to_schema( xsd );
}


protocols::moves::MoverOP
AlignmentFileGeneratorMover::clone() const{
	return protocols::moves::MoverOP( new AlignmentFileGeneratorMover( *this ) );
}


protocols::moves::MoverOP
AlignmentFileGeneratorMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AlignmentFileGeneratorMover );
}

std::string
AlignmentFileGeneratorMover::get_name() const {
	return AlignmentFileGeneratorMover::class_name();
}

std::string
AlignmentFileGeneratorMover::class_name() {
	return "AlignmentFileGeneratorMover";
}

void
AlignmentFileGeneratorMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, AlignmentFileGeneratorMover const &mover)
{
	mover.show(os);
	return os;
}


void
AlignmentFileGeneratorMover::apply( core::pose::Pose& pose){
	//A.Convert pose to segments
	//The new function should just go ahead and add it to the segment_vector in place

	//pose_model = create_model_from_pose( pose, alignment_settings_.segments_, new_model_id , alignment_settings_.match_segments_ );

	///*****************************EDITING************************************
	utility::vector1< data_storage::LigandDescription > expanded_ligands;
	std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands; //This will be filled in the function
	add_pose_segments_to_segment_vector( pose, nullptr, basis_map_generator_, ligands_, partner_ligands, expanded_ligands, required_resnums_, required_selector_, set_segments_from_dssp_ );




	//Append the new segments to the edge map (they should all be connected
	//Set up for initial recursion
	core::Size current_depth = 0;
	//Create a set of all segments to find alignments for in the model map
	std::set< data_storage::SmartSegmentCOP > current_segment_set;
	current_segment_set.clear();
	for ( std::pair< core::Size, data_storage::SmartSegmentOP > pdbseg: basis_map_generator_->pdb_segments() ) {
		basis_map_generator_->edge_file_reader()->append_to_edge_map( pdbseg.first );
		//Add the N- and C-terminal segments from the newly added model
		current_segment_set.insert( data_storage::SmartSegment::get_n_most_segment( pdbseg.second ,false ));//->get_n_most_segment( false ) );
		current_segment_set.insert( data_storage::SmartSegment::get_c_most_segment( pdbseg.second ,false ));//->get_c_most_segment( false ) );

	}
	// basis_map_generator_->edge_file_reader()->append_to_edge_map( basis_map_generator_->pdb_segments()



	basis_map_generator_->recurse_over_segments( current_depth, current_segment_set );
	std::string new_file_name = pose.pdb_info()->name() + "_" + pose.pdb_info()->modeltag() + "_" + utility::to_string( basis_map_generator_->version() ) + ".aln";
	basis_map_generator_->write_alignments_to_file( new_file_name );
}


//This version adds the segments to BMG's pdb_segments() (which returns a reference
core::Size
AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
	core::pose::Pose const & pose,
	core::pose::PoseOP partner_pose,
	BasisMapGeneratorOP bmg,
	utility::vector1< data_storage::LigandDescription > & ligands,
	std::map< core::Size, data_storage::LigandResidueCOP > & partner_ligands, //This will be filled in the function
	utility::vector1< data_storage::LigandDescription > & expanded_ligands,
	std::string required_resnums,
	core::select::residue_selector::ResidueSelectorCOP required_selector,
	bool set_segments_from_dssp
)
{
	core::Size added_segments = 0;
	core::Size current_ligand_id = 1;
	//Compute the pose's secondary structure now
	core::scoring::dssp::Dssp pose_dssp_object( pose );
	std::string secstruct = pose_dssp_object.get_dssp_reduced_IG_as_L_secstruct();
	core::chemical::AtomTypeSetCAP atom_types = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");

	//Step 0: Set up segments if getting from DSSP
	if ( set_segments_from_dssp ) {
		hashing::AlignmentSettings as = bmg->alignment_settings();
		utility::vector1< core::Size > pose_segment_starts;
		utility::vector1< core::Size > pose_segment_ends;
		char last_dssp = '-';
		if ( secstruct.size() == 1 ) { //loop won't execute
			last_dssp = secstruct[ 0 ];
			pose_segment_starts.push_back( 1 );
			pose_segment_ends.push_back( 1 );
		} else if ( pose.residue( 1 ).is_protein() ) {
			pose_segment_starts.push_back( 1 );
		}
		for ( core::Size i = 2; i <= secstruct.size(); ++i ) {
			if ( secstruct[ i - 1 ] != last_dssp && secstruct[ i - 2 ] == secstruct[ i - 1] ) {
				if ( secstruct[ i - 1 ] != '-' ) { //Ligands don't count as segments
					pose_segment_starts.push_back( i - 1 );
				}
				if ( i > 2 && ( pose_segment_starts.size() - pose_segment_ends.size() ) > 1  ) {
					pose_segment_ends.push_back( i - 2 );
				}
				last_dssp = secstruct[ i - 1 ];
			}
		}
		pose_segment_ends.push_back( secstruct.size() );
		//Set up the segments vector of pairs
		utility::vector1< std::pair< core::Size, core::Size > > segs;
		for ( core::Size j = 1; j <= pose_segment_starts.size(); ++j ) {
			if ( j > pose_segment_ends.size() ) {
				utility_exit_with_message( "Error initializing segments from DSSP!" );
			}
			segs.push_back( std::make_pair( pose_segment_starts.at( j ), pose_segment_ends.at( j ) ) );
		}
		as.segments_ = segs;
		bmg->alignment_settings( as );
	}
	//STEP 1: Parse all residue numbers/atom names from pose
	//1A: Initialize variables
	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), false );
	core::Size partner_resn = 1;
	if ( partner_pose ) {
		partner_resn = partner_pose->total_residue();
	}
	core::select::residue_selector::ResidueSubset partner_selection( partner_resn, false );
	std::set< core::Size > required_resnums_final;
	//1AA: parse required residue numbers
	//1AAi: parse required_resnums
	required_resnums_final = core::pose::get_resnum_list( required_resnums, pose );
	//1AAii: parse residue selector
	if ( required_selector != nullptr ) {
		selection = required_selector->apply( pose );
	}
	for ( core::Size resi = 1; resi <= pose.total_residue(); ++resi ) {
		if ( selection[resi] ) {
			required_resnums_final.insert( resi );
		}
	}

	std::set< core::Size > ligand_resnums_final; //To aid in skipping ligand residues during segment generation
	//1B: Parse ligand residue numbers
	for ( data_storage::LigandDescription & ligand_des: ligands ) {
		//utility::vector1< core::Size > resnums;
		utility::vector1< data_storage::LigandDescription > duplicate_ligand_des;
		if ( ligand_des.ligand_selector != nullptr ) {
			TR.Debug << "Ligand initialized from residue selector!" << std::endl;
			if ( ligand_des.partner_ligand ) {
				partner_selection = ligand_des.ligand_selector->apply( *partner_pose );
				for ( core::Size resi = 1, resimax = pose.total_residue() ; resi <= resimax; ++resi ) {
					if ( partner_selection[resi] && duplicate_ligand_des.size() == 0 ) {
						//resnums.push_back( resi );
						ligand_des.ligand_resnum = resi;
						duplicate_ligand_des.push_back( ligand_des );
					} else if ( partner_selection[ resi ] ) {
						//resnums.push_back( resi );
						data_storage::LigandDescription new_des( ligand_des );
						new_des.ligand_resnum = resi;
						duplicate_ligand_des.push_back( new_des );
					}
				}
			} else {
				selection = ligand_des.ligand_selector->apply( pose );
				for ( core::Size resi = 1, resimax = pose.total_residue() ; resi <= resimax; ++resi ) {
					if ( selection[resi] && duplicate_ligand_des.size() == 0 ) {
						//resnums.push_back( resi );
						ligand_des.ligand_resnum = resi;
						duplicate_ligand_des.push_back( ligand_des );
						ligand_resnums_final.insert( resi );
					} else if ( selection[ resi ] ) {
						//resnums.push_back( resi );
						data_storage::LigandDescription new_des( ligand_des );
						new_des.ligand_resnum = resi;
						duplicate_ligand_des.push_back( new_des );
						ligand_resnums_final.insert( resi );
					}
				}
			}
		} else if ( ligand_des.partner_ligand ) {
			core::Size resi = core::pose::parse_resnum ( ligand_des.ligand_resnum_string, *partner_pose, true );
			//resnums.push_back( resi );
			ligand_des.ligand_resnum = resi;
			duplicate_ligand_des.push_back( ligand_des );
		} else {
			core::Size resi = core::pose::parse_resnum( ligand_des.ligand_resnum_string, pose, true );
			//resnums.push_back( resi );
			ligand_des.ligand_resnum = resi;
			duplicate_ligand_des.push_back( ligand_des );
			ligand_resnums_final.insert( resi );
		}
		//From here on out, iterate over resnums and apply the corresponding values in the original data structures to all of them
		//1C: ligand_pdb_contacts
		for ( data_storage::LigandDescription & ligdes: duplicate_ligand_des ) {
			for ( data_storage::ContactDescription & contact_des: ligdes.ligand_contacts ) {
				if ( contact_des.partner_contact ) {
					contact_des.contact_resnum = core::pose::parse_resnum( contact_des.contact_resnum_string, *partner_pose, true );
					contact_des.contact_atom_num = partner_pose->residue( contact_des.contact_resnum ).atom_index( contact_des.contact_atom_name );
				} else {
					contact_des.contact_resnum = core::pose::parse_resnum( contact_des.contact_resnum_string, pose, true );
					contact_des.contact_atom_num = pose.residue( contact_des.contact_resnum ).atom_index( contact_des.contact_atom_name );
				}
				if ( ligdes.partner_ligand ) {
					contact_des.ligand_atom_num = partner_pose->residue( ligdes.ligand_resnum ).atom_index( contact_des.ligand_atom_name );
				} else {
					contact_des.ligand_atom_num = pose.residue( ligdes.ligand_resnum ).atom_index( contact_des.ligand_atom_name );
				}
			}

			//1D: alignment_atoms
			for ( std::string atm: ligdes.alignment_atoms_str ) {
				if ( ligdes.partner_ligand ) {
					ligdes.alignment_atoms_num.push_back( partner_pose->residue( ligdes.ligand_resnum ).atom_index( atm ) );
				} else {
					ligdes.alignment_atoms_num.push_back( pose.residue( ligdes.ligand_resnum ).atom_index( atm ) );
				}
			}

			//1E: ideal_contacts
			for ( std::pair< std::string, IdealContact > contact: ligdes.ideal_contacts_str ) {
				core::Size atom_no = 0;
				if ( ligdes.partner_ligand ) {
					atom_no = partner_pose->residue( ligdes.ligand_resnum ).atom_index( contact.first );
				} else {
					atom_no = pose.residue( ligdes.ligand_resnum ).atom_index( contact.first );
				}
				contact.second.ligand_atom = atom_no;
				ligdes.ideal_contacts_num[ atom_no ] = contact.second;
			}//end iterate over ideal contacts
			expanded_ligands.push_back( ligdes );
		}//end iterate over duplicate_ligand_des
		//expanded_ligands.push_back( duplicate_ligand_des );
	}
	//END STEP 1


	//segments is now part of AlignmentSettings
	//There is no model id, and the segment id will be determined by the segment's position in the vector
	//match_segments is also in AlignmentSettings


	SegmentVectorCOP segvec = bmg->segment_vector();
	core::Size starting_segvec_size = segvec->size();
	utility::vector1< core::conformation::Atom > res_atoms;
	core::Size residues_in_previous_segments = 0;
	for ( core::Size i=1; i<=bmg->alignment_settings().segments_.size(); ++i ) {
		// Begin iterate through segments
		core::Size current_seg_id = starting_segvec_size + i;
		data_storage::SmartSegmentOP seg ( new data_storage::SmartSegment( true ) );
		if ( ligand_resnums_final.size() > 0 ) {
			seg = data_storage::SmartSegmentOP(  new data_storage::LigandSegment( true )  );
		}
		//The segment id will be the original segvec size + i
		seg->set_segment_id( current_seg_id );
		//Get vector of residues
		utility::vector1< data_storage::SmartSewingResidueOP > seg_residues;
		for ( core::Size j=bmg->alignment_settings().segments_[i].first;
				j <= bmg->alignment_settings().segments_[i].second;
				++j ) {
			//Begin iterating through residues in segment i
			//if( j == ligand_index ){
			if ( std::find( ligand_resnums_final.begin(), ligand_resnums_final.end(), j ) != ligand_resnums_final.end() ) { //If j is one of the ligand resnums
				continue;
			}
			data_storage::SmartSewingResidueOP res( new data_storage::SmartSewingResidue );
			res->set_amino_acid_type( pose.residue( j ).type().base_name() );
			res->set_full_type_name( pose.residue( j ).type().name() );
			res->set_chis( pose.residue( j ).chi() );

			//Get the atom vector for this residue
			res_atoms = utility::vector1< core::conformation::Atom >( pose.residue( j ).atoms().begin(), pose.residue( j ).atoms().end() );
			utility::vector1< core::conformation::Atom > new_atoms;
			for ( core::Size ati = 1; ati <= res_atoms.size(); ++ati ) {
				TR << "Atom of type " <<  ( *atom_types.lock() )[ pose.residue( j ).atom( ati ).type() ].atom_type_name() << std::endl;
				if ( ( *atom_types.lock() )[ pose.residue( j ).atom( ati ).type() ].atom_type_name() != "VIRT" ) {
					new_atoms.push_back( res_atoms[ ati ] );
					TR <<  "Added atom" << std::endl;
				}
			}
			res->set_atom_vector( new_atoms );


			//If this residue is in required_resnums_final, set it to be vital
			if ( required_resnums_final.count( j ) != 0 ) {
				seg->add_vital_residue( j - residues_in_previous_segments );
				TR.Debug << "Adding required residue " << j << " as residue " << j - residues_in_previous_segments << "in segment " << current_seg_id << std::endl;
			}
			//Finally, get the secondary structure for this segment from the middle residue in that segment
			if ( j == (bmg->alignment_settings().segments_[i].first + bmg->alignment_settings().segments_[i].second)/2 ) {
				seg->set_dssp_code( secstruct[ j - 1] );
			}
			seg_residues.push_back( res );
		}//Done iterating through residues in segment i
		//We now have gone through all the residues in the segment
		seg->set_residue_vector( seg_residues );
		seg->set_length( seg_residues.size() );
		//Find if segment is supposed to be hashed
		//  if( std::find( basis_map_generator_->alignment_settings().match_segments_.begin(), basis_map_generator_->alignment_settings().match_segments_.end(), i ) == basis_map_generator_->alignment_settings().match_segments_.end() ) {
		//   seg->set_is_hashable( false );
		//  }
		//Next, we'll want to connect this segments to any connected segments
		//Connected segments will be any segments that we have added since beginning this function

		TR.Debug << "Vital residues for segment " << seg->get_segment_id() << std::endl;
		for ( core::Size vital: seg->get_vital_residues() ) {
			TR.Debug << vital << " ";
		}
		TR.Debug << std::endl;



		if ( i > 1 ) {
			data_storage::SmartSegment::link_to( bmg->pdb_segments().at( current_seg_id - 1 ), seg );
		}

		bmg->pdb_segments()[ current_seg_id ] = seg;

		++ added_segments;
		residues_in_previous_segments = bmg->alignment_settings().segments_[i].second;
	}//End iterate through segments





	//Now we can go back and add in the ligands, both partner and pose
	core::pose::Pose combined_pose;
	core::Size partner_res = 0;
	if ( partner_pose ) {
		partner_res = partner_pose->total_residue();
	}
	//Partner will be first, then pose
	if ( expanded_ligands.size() > 0 && partner_pose != nullptr ) {
		combined_pose = core::pose::Pose( *partner_pose );
		core::pose::append_pose_to_pose( combined_pose, pose, true );
	} else {
		combined_pose = pose; //No need to even bother making a copy
	}
	for ( data_storage::LigandDescription & current_ligand: expanded_ligands ) {

		data_storage::LigandResidueOP ligand( new data_storage::LigandResidue );
		ligand->set_ligand_id( current_ligand_id );
		current_ligand.ligand_id = current_ligand_id;
		++current_ligand_id;

		if ( current_ligand.partner_ligand ) {
			ligand->set_partner_ligand( true );
			if ( current_ligand.ligand_resnum > partner_res || current_ligand.ligand_resnum < 1 ) {
				TR.Debug << "Ligand index " << current_ligand.ligand_id << " is invalid!" << std::endl;
				utility_exit_with_message( "Invalid ligand resnum in partner!" );
			}
			ligand->set_amino_acid_type( partner_pose->residue( current_ligand.ligand_resnum ).type().base_name() );
			ligand->set_full_type_name( partner_pose->residue( current_ligand.ligand_resnum ).type().name() );
			//We want to include all of the ligand atoms
			res_atoms = utility::vector1< core::conformation::Atom >( partner_pose->residue( current_ligand.ligand_resnum ).atoms().begin(), partner_pose->residue( current_ligand.ligand_resnum ).atoms().end() );
			utility::vector1< core::conformation::Atom > lig_atoms;
			for ( core::Size ati = 1; ati <= res_atoms.size(); ++ati ) {
				TR << "Atom of type " <<  ( *atom_types.lock() )[ partner_pose->residue( current_ligand.ligand_resnum ).atom( ati ).type() ].atom_type_name() << std::endl;
				if ( ( *atom_types.lock() )[ partner_pose->residue( current_ligand.ligand_resnum ).atom( ati ).type() ].atom_type_name() != "VIRT" ) {
					lig_atoms.push_back( res_atoms[ ati ] );
				}
			}
			ligand->set_atom_vector( lig_atoms );
			ligand->set_type( data_storage::ligand );
			ligand->set_alignment_atoms( current_ligand.alignment_atoms_num );

			if ( current_ligand.ideal_contacts_num.size() != 0 ) {
				//The entry at this value is a map of ligand atom name to an IdealContact with all other info set
				for ( std::pair< core::Size, IdealContact > contact: current_ligand.ideal_contacts_num ) {
					//contact.second.ligand_atom = pose.residue( ligand_index ).atom_index( contact.first );
					ligand->add_ideal_contact( contact.second );
				}
			}
			//Initialize attached_ligand map
			//Has the ligand already been attached to THIS segment?
			std::map< core::Size, bool > attached_ligand;
			for ( std::pair< core::Size, data_storage::SmartSegmentOP > seg: bmg->pdb_segments() ) {
				attached_ligand[ seg.first ] = false;
			}

			//Begin auto detect contacts (if applicable)
			if ( current_ligand.auto_detect_contacts ) {
				core::conformation::Residue const & ligand_res = combined_pose.residue( current_ligand.ligand_resnum );

				//Begin bonded contacts
				//Iterate over all ligand atoms
				for ( core::Size ligand_atom_num = 1; ligand_atom_num <= ligand_res.natoms(); ++ligand_atom_num ) {
					if ( ( *atom_types.lock() )[ ligand_res.atom( ligand_atom_num ).type() ].atom_type_name() == "VIRT" ) {
						TR << "Skipping virtual ligand atom!" << std::endl;
						continue;
					}
					utility::vector1< core::id::AtomID > bonded_atoms = combined_pose.conformation().bonded_neighbor_all_res( core::id::AtomID( ligand_atom_num, current_ligand.ligand_resnum ) );

					for ( core::id::AtomID atom_id: bonded_atoms ) {
						TR << "Atom of type " <<  ( *atom_types.lock() )[ combined_pose.residue( atom_id.rsd() ).atom( atom_id.atomno() ).type() ].atom_type_name() << std::endl;
						if ( ( *atom_types.lock() )[ combined_pose.residue( atom_id.rsd() ).atom( atom_id.atomno() ).type() ].atom_type_name() != "VIRT" ) {

							//Make a new ContactDescription
							data_storage::ContactDescription new_contact;
							new_contact.ligand_atom_num = ligand_atom_num;
							new_contact.contact_atom_num = atom_id.atomno();
							if ( atom_id.rsd() <= partner_res ) {
								new_contact.partner_contact = true;
								new_contact.contact_resnum = atom_id.rsd();
							} else {
								new_contact.partner_contact = false;
								new_contact.contact_resnum = atom_id.rsd() - partner_res;
							}
							//We won't need atom names
							current_ligand.ligand_contacts.push_back( new_contact );
						}
					}//end for bonded atoms
				}//end for ligand atoms

				//Begin metal contacts
				std::map< core::Size, utility::vector1< core::id::AtomID > > metal_contacts;
				metal_contacts = core::util::find_metalbinding_atoms_for_complex( combined_pose, current_ligand.ligand_resnum, 1 );
				for ( std::pair< core::Size, utility::vector1< core::id::AtomID > > metal_atom_contacts: metal_contacts ) {
					for ( core::id::AtomID kk: metal_atom_contacts.second ) {
						//Make a new ContactDescription
						data_storage::ContactDescription new_contact;
						new_contact.ligand_atom_num = metal_atom_contacts.first;
						new_contact.contact_atom_num = kk.atomno();
						if ( kk.rsd() <= partner_res ) {
							new_contact.partner_contact = true;
							new_contact.contact_resnum = kk.rsd();
						} else {
							new_contact.partner_contact = false;
							new_contact.contact_resnum = kk.rsd() - partner_res;
						}
						//We won't need atom names
						current_ligand.ligand_contacts.push_back( new_contact );
					}//End for atom id
				}//end for metal contact
			}//end if auto_detect_contacts
			if ( current_ligand.ligand_contacts.size() == 0 ) {
				utility_exit_with_message( "No contacts found for ligand " + utility::to_string( current_ligand.ligand_id ) +  ". Did you forget to auto detect contacts?" );
			}
			//Begin adding ligands to assembly or to partner_ligands and attaching them to segments (if applicable)
			//Iterate over contacts
			for ( data_storage::ContactDescription & contact: current_ligand.ligand_contacts ) {
				//Internal contact
				//Resnums match and partner status matches
				if ( contact.contact_resnum == current_ligand.ligand_resnum && ((current_ligand.partner_ligand && contact.partner_contact ) || (!current_ligand.partner_ligand && !contact.partner_contact ) ) ) {
					TR.Debug << "Adding internal contact for ligand " << ligand->get_ligand_id() << std::endl;
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( 0, 0, contact.contact_atom_num, contact.ligand_atom_num ) ) );
					continue; //No need to try to attach it to a segment based on internal contacts
				}//End if internal contact
				//Non-internal partner contacts
				if ( contact.partner_contact ) {
					TR.Debug << "Adding partner contact for ligand " << ligand->get_ligand_id() << std::endl;
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( 0, contact.contact_resnum, contact.contact_atom_num, contact.ligand_atom_num ) ) );
					continue; //No need to try to attach it to a segment based on partner contacts
				}
				//Begin non-internal assembly contacts
				core::Size passed_residues = 0;
				for ( std::pair< core::Size, data_storage::SmartSegmentOP > seg: bmg->pdb_segments() ) {
					data_storage::LigandSegmentOP ligseg = std::dynamic_pointer_cast< data_storage::LigandSegment>( seg.second );
					TR.Debug << "Vital residues for segment " << ligseg << std::endl;
					for ( core::Size vital: ligseg->get_vital_residues() ) {
						TR.Debug << vital << " ";
					}
					TR.Debug << std::endl;
					if ( contact.contact_resnum > ligseg->get_length() + passed_residues ) {
						passed_residues += seg.second->get_length();
						continue; //Check the next segment
					}
					if ( contact.contact_resnum < passed_residues ) {
						utility_exit_with_message( "Error when adding ligand to pose segments! Check that your indices are all in range." );
					}
					//Otherwise this contact must be in this segment
					data_storage::SmartSewingResidueCOP contact_res = ligseg->get_residue( contact.contact_resnum - passed_residues );


					//Partner contacts will not have owners! They are free elves!
					if ( !attached_ligand.at( seg.first ) ) {
						ligseg->attach_ligand( ligand, false );
						attached_ligand[ seg.first ] = true;
					} //end if not attached
					runtime_assert( contact.contact_resnum > passed_residues );
					ligseg->add_ligand_contact( contact.contact_resnum - passed_residues );
					TR.Debug << "Vital residues for segment " << ligseg << std::endl;
					for ( core::Size vital: ligseg->get_vital_residues() ) {
						TR.Debug << vital << " ";
					}
					TR.Debug << std::endl;
					//Add the contact residue to the ligand
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( seg.first, (contact.contact_resnum - passed_residues ), contact.contact_atom_num, contact.ligand_atom_num ) ) );
					//Now add to partner ligands
					partner_ligands[ ligand->get_ligand_id() ] = ligand;
					break; //We don't need to keep looking through additional segments
				}//End iterate over pdb_segments
			}//end for current contacts
		} else { //end if partner ligand
			ligand->set_partner_ligand( false );
			if ( current_ligand.ligand_resnum > pose.total_residue() || current_ligand.ligand_resnum < 1 ) {
				TR.Debug << "Ligand index " << current_ligand.ligand_id << " is invalid!" << std::endl;
				utility_exit_with_message( "Invalid ligand resnum in pose!" );
				continue;
			}
			ligand->set_amino_acid_type( pose.residue( current_ligand.ligand_resnum ).type().base_name() );
			ligand->set_full_type_name( pose.residue( current_ligand.ligand_resnum ).type().name() );
			//We want to include all of the ligand atoms
			res_atoms = utility::vector1< core::conformation::Atom >( pose.residue( current_ligand.ligand_resnum ).atoms().begin(), pose.residue( current_ligand.ligand_resnum ).atoms().end() );
			utility::vector1< core::conformation::Atom > lig_atoms;
			for ( core::Size ati = 1; ati <= res_atoms.size(); ++ati ) {
				TR << "Atom of type " <<  ( *atom_types.lock() )[ pose.residue( current_ligand.ligand_resnum ).atom( ati ).type() ].atom_type_name() << std::endl;
				if ( ( *atom_types.lock() )[ pose.residue( current_ligand.ligand_resnum ).atom( ati ).type() ].atom_type_name() != "VIRT" ) {
					lig_atoms.push_back( res_atoms[ ati ] );
				}
			}
			ligand->set_atom_vector( lig_atoms );
			ligand->set_type( data_storage::ligand );
			ligand->set_alignment_atoms( current_ligand.alignment_atoms_num );

			if ( current_ligand.ideal_contacts_num.size() != 0 ) {
				//The entry at this value is a map of ligand atom name to an IdealContact with all other info set
				for ( std::pair< core::Size, IdealContact > contact: current_ligand.ideal_contacts_num ) {
					//contact.second.ligand_atom = pose.residue( ligand_index ).atom_index( contact.first );
					ligand->add_ideal_contact( contact.second );
				}
			}
			//Initialize attached_ligand map
			//Has the ligand already been attached to THIS segment?
			std::map< core::Size, bool > attached_ligand;
			for ( std::pair< core::Size, data_storage::SmartSegmentOP > seg: bmg->pdb_segments() ) {
				attached_ligand[ seg.first ] = false;
			}


			//Begin auto detect contacts ( if applicable)
			if ( current_ligand.auto_detect_contacts ) {
				core::conformation::Residue const & ligand_res = combined_pose.residue( current_ligand.ligand_resnum + partner_res );

				//Begin bonded contacts
				//Iterate over all ligand atoms
				for ( core::Size ligand_atom_num = 1; ligand_atom_num <= ligand_res.natoms(); ++ligand_atom_num ) {
					if ( ( *atom_types.lock() )[ ligand_res.atom( ligand_atom_num ).type() ].atom_type_name() == "VIRT" ) {
						TR << "Skipping virtual ligand atom!" << std::endl;
						continue;
					}
					utility::vector1< core::id::AtomID > bonded_atoms = combined_pose.conformation().bonded_neighbor_all_res( core::id::AtomID( ligand_atom_num, current_ligand.ligand_resnum + partner_res ) );

					for ( core::id::AtomID atom_id: bonded_atoms ) {
						TR << "Atom of type " <<  ( *atom_types.lock() )[ combined_pose.residue( atom_id.rsd() ).atom( atom_id.atomno() ).type() ].atom_type_name() << std::endl;
						if ( ( *atom_types.lock() )[ combined_pose.residue( atom_id.rsd() ).atom( atom_id.atomno() ).type() ].atom_type_name() != "VIRT" ) {

							//Make a new ContactDescription
							data_storage::ContactDescription new_contact;
							new_contact.ligand_atom_num = ligand_atom_num;
							new_contact.contact_atom_num = atom_id.atomno();
							if ( atom_id.rsd() <= partner_res ) {
								new_contact.partner_contact = true;
								new_contact.contact_resnum = atom_id.rsd();
							} else {
								new_contact.partner_contact = false;
								new_contact.contact_resnum = atom_id.rsd() - partner_res;
							}
							//We won't need atom names
							current_ligand.ligand_contacts.push_back( new_contact );
						}
					}//end for bonded atoms
				}//end for ligand atoms

				//Begin metal contacts
				std::map< core::Size, utility::vector1< core::id::AtomID > > metal_contacts;
				metal_contacts = core::util::find_metalbinding_atoms_for_complex( combined_pose, current_ligand.ligand_resnum + partner_res, 1 );
				for ( std::pair< core::Size, utility::vector1< core::id::AtomID > > metal_atom_contacts: metal_contacts ) {
					for ( core::id::AtomID kk: metal_atom_contacts.second ) {
						//Make a new ContactDescription
						data_storage::ContactDescription new_contact;
						new_contact.ligand_atom_num = metal_atom_contacts.first;
						new_contact.contact_atom_num = kk.atomno();
						if ( kk.rsd() <= partner_res ) {
							new_contact.partner_contact = true;
							new_contact.contact_resnum = kk.rsd();
						} else {
							new_contact.partner_contact = false;
							new_contact.contact_resnum = kk.rsd() - partner_res;
						}
						//We won't need atom names
						current_ligand.ligand_contacts.push_back( new_contact );
					}//End for atom id
				}//end for metal contact
			}//end if auto_detect_contacts
			if ( current_ligand.ligand_contacts.size() == 0 && !current_ligand.partner_ligand ) {
				utility_exit_with_message( "No contacts found for ligand " + utility::to_string( current_ligand.ligand_id ) +  ". Did you forget to auto detect contacts?" );
			}
			//Begin adding ligands to assembly or to partner_ligands and attaching them to segments (if applicable)
			//Iterate over contacts
			for ( data_storage::ContactDescription & contact: current_ligand.ligand_contacts ) {
				//Internal contact
				//Resnums match and partner status matches
				if ( contact.contact_resnum == current_ligand.ligand_resnum && ((current_ligand.partner_ligand && contact.partner_contact ) || (!current_ligand.partner_ligand && !contact.partner_contact ) ) ) {
					TR.Debug << "Adding internal contact for ligand " << ligand->get_ligand_id() << std::endl;
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( 0, 0, contact.contact_atom_num, contact.ligand_atom_num ) ) );
					continue; //No need to try to attach it to a segment based on internal contacts
				}//End if internal contact
				//Non-internal partner contacts
				if ( contact.partner_contact ) {
					TR.Debug << "Adding partner contact for ligand " << ligand->get_ligand_id() << std::endl;
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( 0, contact.contact_resnum, contact.contact_atom_num, contact.ligand_atom_num ) ) );
					continue; //No need to try to attach it to a segment based on partner contacts
				}
				//Begin non-internal assembly contacts
				core::Size passed_residues = 0;
				for ( std::pair< core::Size, data_storage::SmartSegmentOP > seg: bmg->pdb_segments() ) {
					data_storage::LigandSegmentOP ligseg = std::dynamic_pointer_cast< data_storage::LigandSegment>( seg.second );
					TR.Debug << "Vital residues for segment " << ligseg << std::endl;
					for ( core::Size vital: ligseg->get_vital_residues() ) {
						TR.Debug << vital << " ";
					}
					TR.Debug << std::endl;
					if ( contact.contact_resnum > ligseg->get_length() + passed_residues ) {
						passed_residues += seg.second->get_length();
						continue; //Check the next segment
					}
					if ( contact.contact_resnum < passed_residues ) {
						utility_exit_with_message( "Error when adding ligand to pose segments! Check that your indices are all in range." );
					}
					//Otherwise this contact must be in this segment
					data_storage::SmartSewingResidueCOP contact_res = ligseg->get_residue( contact.contact_resnum - passed_residues );
					if ( ligand->get_nonconst_owner_segment() == nullptr ) { //We're attaching the ligand at the first available contact.
						ligand->set_owner_segment( ligseg );
						//Add the contact residue to the ligand
						ligseg->attach_ligand( ligand, true );
						attached_ligand[ seg.first ] = true;
						//added_ligand = true;
					} else if ( !attached_ligand.at( seg.first ) ) { //end if no owner
						ligseg->attach_ligand( ligand, false );
						attached_ligand[ seg.first ] = true;
					} //end if not attached
					runtime_assert( contact.contact_resnum > passed_residues );
					ligseg->add_ligand_contact( contact.contact_resnum - passed_residues );
					TR.Debug << "Vital residues for segment " << ligseg << std::endl;
					for ( core::Size vital: ligseg->get_vital_residues() ) {
						TR.Debug << vital << " ";
					}
					TR.Debug << std::endl;
					//Add the contact residue to the ligand
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( seg.first, (contact.contact_resnum - passed_residues ), contact.contact_atom_num, contact.ligand_atom_num ) ) );
					break; //We don't need to keep looking through additional segments
				}//End iterate over pdb_segments
			}//end for current contacts
		}//End else (not partner ligand)
	}//End iterate over ligands





	return added_segments;
}//end add_pose_segments_to_segment_vector


//This method takes a reference to the Assembly's pdb_segments
core::Size
AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector(
	core::pose::Pose const & pose,
	core::pose::PoseOP partner_pose,
	SegmentVectorCOP segvec,
	std::map< core::Size, data_storage::SmartSegmentOP > & pdb_segments,
	std::string pose_segment_starts_string,
	std::string pose_segment_ends_string,
	std::string pose_segment_dssp,
	utility::vector1< data_storage::LigandDescription > & ligands,
	std::map< core::Size, data_storage::LigandResidueCOP > & partner_ligands, //This will be filled in the function
	utility::vector1< data_storage::LigandDescription > & expanded_ligands,
	std::string required_resnums,
	core::select::residue_selector::ResidueSelectorCOP required_selector,
	bool strict_dssp_changes
)
{
	core::Size added_segments = 0;
	core::Size current_ligand_id = 1;
	core::scoring::dssp::Dssp pose_dssp_object( pose );
	std::string secstruct = pose_dssp_object.get_dssp_reduced_IG_as_L_secstruct();
	core::chemical::AtomTypeSetCAP atom_types = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
	//STEP 0: Set up segments
	utility::vector1< core::Size > pose_segment_starts;
	utility::vector1< core::Size > pose_segment_ends;
	utility::vector1< std::pair< core::Size, core::Size > > segs;
	bool artificial_loop = false; //True if we artificially split a single segment into helix-loop-helix
	if ( pose_segment_starts_string != "" && pose_segment_ends_string != "" ) {
		TR << "Setting up User-specified pose segments" << std::endl;
		//Pose segments were manually specified
		pose_segment_starts = core::pose::get_resnum_list_ordered( pose_segment_starts_string, pose );
		pose_segment_ends = core::pose::get_resnum_list_ordered( pose_segment_ends_string, pose );
		for ( core::Size j = 1; j <= pose_segment_starts.size(); ++j ) {
			if ( j > pose_segment_ends.size() ) {
				utility_exit_with_message( "Error initializing segments from DSSP!" );
			}
			segs.push_back( std::make_pair( pose_segment_starts.at( j ), pose_segment_ends.at( j ) ) );
		}
	} else if ( pose_segment_starts_string != "" || pose_segment_ends_string != "" ) {
		utility_exit_with_message( "ERROR: Only one of pose_segment_starts/pose_segment_ends was specified!" );
	} else {
		TR << "Setting pose segments from DSSP" << std::endl;
		char last_dssp = '-';
		if ( strict_dssp_changes ) {
			if ( secstruct.size() == 1 ) { //loop won't execute
				last_dssp = secstruct[ 0 ];
				pose_segment_starts.push_back( 1 );
				pose_segment_ends.push_back( 1 );
			} else if ( pose.residue( 1 ).is_protein() ) {
				pose_segment_starts.push_back( 1 );
			}
			for ( core::Size i = 2; i <= secstruct.size(); ++i ) {
				if ( secstruct[ i - 1] != last_dssp && secstruct[ i - 2 ] == secstruct[ i - 1] ) {
					if ( secstruct[ i - 1 ] != '-' ) { //Ligands don't count as segments
						pose_segment_starts.push_back( i - 1 );
					}
					if ( i > 2  && ( pose_segment_starts.size() - pose_segment_ends.size() ) > 1 ) {
						pose_segment_ends.push_back( i - 2 );
					}
					last_dssp = secstruct[ i - 1 ];
				}
			}
			pose_segment_ends.push_back( secstruct.size() );
		} else {
			for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
				if ( secstruct[ i - 1 ] != last_dssp ) {
					pose_segment_starts.push_back( i );
					if ( i > 1 && pose_segment_starts.size() > 0 ) {
						pose_segment_ends.push_back( i - 1 );
					}
					last_dssp = secstruct[ i - 1];
				}
			}
			if ( secstruct[ pose.total_residue() - 1 ] == last_dssp ) {
				pose_segment_ends.push_back( pose.total_residue() );
			} else if ( pose_segment_starts.size() > pose_segment_ends.size() ) {
				pose_segment_ends.push_back( pose.total_residue() );
			}
		}
		//Set up the segments vector of pairs
		for ( core::Size j = 1; j <= pose_segment_starts.size(); ++j ) {
			if ( j > pose_segment_ends.size() ) {
				utility_exit_with_message( "Error initializing segments from DSSP!" );
			}
			segs.push_back( std::make_pair( pose_segment_starts.at( j ), pose_segment_ends.at( j ) ) );
		}

		//BEGIN REMOVE TERMINAL LOOPS
		if ( segs.size() == 1 ) {
			//This is really hacky and maybe should be removed, but it's good for now
			TR << "Single segment added. Forcing it to be helical." << std::endl;
			secstruct = std::string( secstruct.size(), 'H' );
		} else {
			TR << "Removing terminal segments that are loops:" << std::endl;
			char seg_dssp;
			core::Size res_to_check;
			//First check the first segment and delete if necessary
			res_to_check = std::min( (segs[ 1 ].second + segs[ 1 ].first) / 2, segs[ 1 ].second );
			seg_dssp = secstruct[ res_to_check - 1 ];
			if ( seg_dssp == 'L' ) {
				segs.erase( segs.begin());
			}
			//Then check the last segment and delete if necessary
			//res_to_check = std::max( ( segs[ segs.size() ].second + segs[ segs.size() ].first ) / 2, segs[ segs.size() ].first );
			res_to_check = segs[ segs.size() ].first;
			seg_dssp = secstruct[ res_to_check - 1 ];
			if ( seg_dssp == 'L' ) {
				segs.erase( segs.end() - 1);
			}
			//END REMOVE TERMINAL LOOPS
			if ( segs.size() == 1 ) {
				//This is really hacky and maybe should be removed, but it's good for now
				TR << "Single segment added. Forcing it to be helical." << std::endl;
				secstruct = std::string( secstruct.size(), 'H' );
			}
		}
	}
	/*
	if ( segs.size() == 1 ) {
	TR << "Single segment input detected. Splitting in half for chimerization" << std::endl;
	if ( pose.total_residue() < 3 ) {
	utility_exit_with_message( "Input pose too short!" );
	}
	//NOTE: We are going to the end of the stored segment, not necessarily the end of the pose
	//and starting at the beginning of the stored segment, not necessarily residue 1
	core::Size loop_res = ( segs[ 1 ].second -  segs[ 1].first + 1) / 2 + segs[ 1 ].first;
	//First segment begins at 1 and ends at loop_res - 1
	//Second segment begins and ends at loop_res
	//Final segment begins at loop_res + 1 and ends at total_residue
	artificial_loop = true;
	std::pair< core::Size, core::Size > seg1 = segs[ 1 ];
	segs.clear();
	segs.push_back( std::make_pair( seg1.first, loop_res - 1 ) );
	segs.push_back( std::make_pair( loop_res, loop_res ) );
	segs.push_back( std::make_pair( loop_res + 1, seg1.second ) );
	}
	*/
	TR << "Segments:" << std::endl;
	for ( core::Size i = 1; i <= segs.size(); ++i ) {
		TR << "Segment " << i << " " << segs[ i ].first << " " << segs[ i ].second << std::endl;
	}
	//STEP 1: Parse all residue numbers/atom names from pose
	//1A: Initialize variables
	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), false );
	core::Size partner_resn = 1;
	if ( partner_pose ) {
		partner_resn = partner_pose->total_residue();
	}
	core::select::residue_selector::ResidueSubset partner_selection( partner_resn, false );
	//utility::vector1< utility::vector1< data_storage::LigandDescription > > expanded_ligands;
	std::set< core::Size > required_resnums_final;
	std::set< core::Size > ligand_resnums_final; //To aid in skipping ligand residues during segment generation






	//1AA: parse required residue numbers
	//1AAi: parse required_resnums
	required_resnums_final = core::pose::get_resnum_list( required_resnums, pose );
	//1AAii: parse residue selector
	if ( required_selector != nullptr ) {
		selection = required_selector->apply( pose );
	}
	for ( core::Size resi = 1; resi <= pose.total_residue(); ++resi ) {
		if ( selection[resi] ) {
			required_resnums_final.insert( resi );
		}
	}




	//1B: Parse ligand residue numbers
	for ( data_storage::LigandDescription & ligand_des: ligands ) {
		//utility::vector1< core::Size > resnums;
		utility::vector1< data_storage::LigandDescription > duplicate_ligand_des;
		if ( ligand_des.ligand_selector != nullptr ) {
			TR.Debug << "Ligand initialized from residue selector!" << std::endl;
			if ( ligand_des.partner_ligand ) {
				partner_selection = ligand_des.ligand_selector->apply( *partner_pose );
				for ( core::Size resi = 1, resimax = pose.total_residue() ; resi <= resimax; ++resi ) {
					if ( partner_selection[resi] && duplicate_ligand_des.size() == 0 ) {
						//resnums.push_back( resi );
						ligand_des.ligand_resnum = resi;
						duplicate_ligand_des.push_back( ligand_des );
					} else if ( partner_selection[ resi ] ) {
						//resnums.push_back( resi );
						data_storage::LigandDescription new_des( ligand_des );
						new_des.ligand_resnum = resi;
						duplicate_ligand_des.push_back( new_des );
					}
				}
			} else {
				selection = ligand_des.ligand_selector->apply( pose );
				for ( core::Size resi = 1, resimax = pose.total_residue() ; resi <= resimax; ++resi ) {
					if ( selection[resi] && duplicate_ligand_des.size() == 0 ) {
						//resnums.push_back( resi );
						ligand_des.ligand_resnum = resi;
						duplicate_ligand_des.push_back( ligand_des );
						ligand_resnums_final.insert( resi );
					} else if ( selection[ resi ] ) {
						//resnums.push_back( resi );
						data_storage::LigandDescription new_des( ligand_des );
						new_des.ligand_resnum = resi;
						duplicate_ligand_des.push_back( new_des );
						ligand_resnums_final.insert( resi );
					}
				}
			}
		} else if ( ligand_des.partner_ligand ) {
			core::Size resi = core::pose::parse_resnum ( ligand_des.ligand_resnum_string, *partner_pose, true );
			//resnums.push_back( resi );
			ligand_des.ligand_resnum = resi;
			duplicate_ligand_des.push_back( ligand_des );
		} else {
			core::Size resi = core::pose::parse_resnum( ligand_des.ligand_resnum_string, pose, true );
			//resnums.push_back( resi );
			ligand_des.ligand_resnum = resi;
			duplicate_ligand_des.push_back( ligand_des );
			ligand_resnums_final.insert( resi );
		}
		//From here on out, iterate over resnums and apply the corresponding values in the original data structures to all of them
		//1C: ligand_pdb_contacts
		for ( data_storage::LigandDescription & ligdes: duplicate_ligand_des ) {
			for ( data_storage::ContactDescription & contact_des: ligdes.ligand_contacts ) {
				if ( contact_des.partner_contact ) {
					contact_des.contact_resnum = core::pose::parse_resnum( contact_des.contact_resnum_string, *partner_pose, true );
					contact_des.contact_atom_num = partner_pose->residue( contact_des.contact_resnum ).atom_index( contact_des.contact_atom_name );
				} else {
					contact_des.contact_resnum = core::pose::parse_resnum( contact_des.contact_resnum_string, pose, true );
					contact_des.contact_atom_num = pose.residue( contact_des.contact_resnum ).atom_index( contact_des.contact_atom_name );
				}
				if ( ligdes.partner_ligand ) {
					contact_des.ligand_atom_num = partner_pose->residue( ligdes.ligand_resnum ).atom_index( contact_des.ligand_atom_name );
				} else {
					contact_des.ligand_atom_num = pose.residue( ligdes.ligand_resnum ).atom_index( contact_des.ligand_atom_name );
				}
			}

			//1D: alignment_atoms
			for ( std::string atm: ligdes.alignment_atoms_str ) {
				if ( ligdes.partner_ligand ) {
					ligdes.alignment_atoms_num.push_back( partner_pose->residue( ligdes.ligand_resnum ).atom_index( atm ) );
				} else {
					ligdes.alignment_atoms_num.push_back( pose.residue( ligdes.ligand_resnum ).atom_index( atm ) );
				}
			}

			//1E: ideal_contacts
			for ( std::pair< std::string, IdealContact > contact: ligdes.ideal_contacts_str ) {
				core::Size atom_no = 0;
				if ( ligdes.partner_ligand ) {
					atom_no = partner_pose->residue( ligdes.ligand_resnum ).atom_index( contact.first );
				} else {
					atom_no = pose.residue( ligdes.ligand_resnum ).atom_index( contact.first );
				}
				contact.second.ligand_atom = atom_no;
				ligdes.ideal_contacts_num[ atom_no ] = contact.second;
			}//end iterate over ideal contacts
			expanded_ligands.push_back( ligdes );
		}//end iterate over duplicate_ligand_des
		//expanded_ligands.push_back( duplicate_ligand_des );
	}
	//END STEP 1




	//We will use expanded_ligands for all stages of ligand interaction from here on out

	TR.Debug << "pdb_segments has " << pdb_segments.size() << "segments" << std::endl;
	core::Size starting_segs = segvec->size() + pdb_segments.size() + 1;
	utility::vector1< core::conformation::Atom > res_atoms;

	core::Size residues_in_previous_segments = 0;
	core::Size last_added_seg_id;
	for ( core::Size i=1; i<=segs.size(); ++i ) {
		data_storage::SmartSegmentOP current_segment = data_storage::SmartSegmentOP(new data_storage::SmartSegment());
		if ( ligands.size() != 0 ) {
			current_segment = data_storage::SmartSegmentOP ( new data_storage::LigandSegment( true ) );
		}
		core::Size current_seg_id = starting_segs + i;
		last_added_seg_id = current_seg_id;
		current_segment->set_segment_id( current_seg_id );
		utility::vector1< data_storage::SmartSewingResidueOP > seg_residues;
		for ( core::Size j = segs[ i ].first; j <= segs[ i ].second; ++j ) {
			if ( std::find( ligand_resnums_final.begin(), ligand_resnums_final.end(), j ) != ligand_resnums_final.end() ) { //If j is one of the ligand resnums
				continue;
			}
			data_storage::SmartSewingResidueOP res( new data_storage::SmartSewingResidue );
			res->set_amino_acid_type( pose.residue( j ).type().base_name() );
			res->set_full_type_name( pose.residue( j ).type().name() );
			res->set_chis( pose.residue( j ).chi() );
			res_atoms = utility::vector1< core::conformation::Atom >( pose.residue( j ).atoms().begin(), pose.residue( j ).atoms().end() );
			utility::vector1< core::conformation::Atom > new_atoms;
			for ( core::Size ati = 1; ati <= res_atoms.size(); ++ati ) {
				TR.Debug << "Atom of type " <<  ( *atom_types.lock() )[ pose.residue( j ).atom( ati ).type() ].atom_type_name() << std::endl;
				if ( ( *atom_types.lock() )[ pose.residue( j ).atom( ati ).type() ].atom_type_name() != "VIRT" ) {
					new_atoms.push_back( res_atoms[ ati ] );
				}
			}
			res->set_atom_vector( new_atoms );
			if ( required_resnums_final.count( j ) != 0 ) {
				current_segment->add_vital_residue( j - residues_in_previous_segments );
				TR.Debug << "Adding required residue " << j << " as residue " << j - residues_in_previous_segments << "in segment " << current_seg_id << std::endl;
			}
			//Finally, get the secondary structure for this segment from the middle residue in that segment
			//REFACTOR!!!
			if ( j == core::Size( ( segs[i].first + segs[i].second )/2 ) ) {
				if ( artificial_loop && i == 2 ) {
					current_segment->set_dssp_code( 'L' );
				} else if ( pose_segment_starts_string != "" && pose_segment_dssp.size() >= i ) {
					current_segment->set_dssp_code( pose_segment_dssp[ i - 1 ] );
				} else {
					current_segment->set_dssp_code( secstruct[ j - 1 ] );
				}
			}
			seg_residues.push_back( res );
		}//Done iterating through residues in segment i
		current_segment->set_residue_vector( seg_residues );
		current_segment->set_length( seg_residues.size() );

		TR.Debug << "Vital residues for segment " << current_segment->get_segment_id() << std::endl;
		for ( core::Size vital: current_segment->get_vital_residues() ) {
			TR.Debug << vital << " ";
		}
		TR.Debug << std::endl;

		if ( i > 1 ) { //If this is not hte first segment, we need to link it to any previous segments
			data_storage::SmartSegment::link_to( pdb_segments.at( current_seg_id - 1 ), current_segment );
		}

		pdb_segments[ current_seg_id ] = current_segment;
		++ added_segments;
		residues_in_previous_segments = segs[i].second;
	}//End iterate through segments
	//Now we can go back and add in the ligands, both partner and pose
	core::pose::Pose combined_pose;
	core::Size partner_res = 0;
	if ( partner_pose ) {
		partner_res = partner_pose->total_residue();
	}
	//Partner will be first, then pose
	if ( expanded_ligands.size() > 0 && partner_pose != nullptr ) {
		combined_pose = core::pose::Pose( *partner_pose );
		core::pose::append_pose_to_pose( combined_pose, pose, true );
	} else {
		combined_pose = pose; //No need to even bother making a copy
	}
	for ( data_storage::LigandDescription & current_ligand: expanded_ligands ) {

		data_storage::LigandResidueOP ligand( new data_storage::LigandResidue );
		ligand->set_ligand_id( current_ligand_id );
		current_ligand.ligand_id = current_ligand_id;
		++current_ligand_id;

		if ( current_ligand.partner_ligand ) {
			ligand->set_partner_ligand( true );
			if ( current_ligand.ligand_resnum > partner_res || current_ligand.ligand_resnum < 1 ) {
				TR.Debug << "Ligand index " << current_ligand.ligand_id << " is invalid!" << std::endl;
				utility_exit_with_message( "Invalid ligand resnum in partner!" );
			}
			ligand->set_amino_acid_type( partner_pose->residue( current_ligand.ligand_resnum ).type().base_name() );
			ligand->set_full_type_name( partner_pose->residue( current_ligand.ligand_resnum ).type().name() );
			//We want to include all of the ligand atoms
			res_atoms = utility::vector1< core::conformation::Atom >( partner_pose->residue( current_ligand.ligand_resnum ).atoms().begin(), partner_pose->residue( current_ligand.ligand_resnum ).atoms().end() );
			utility::vector1< core::conformation::Atom > lig_atoms;
			for ( core::Size ati = 1; ati <= res_atoms.size(); ++ati ) {
				TR << "Atom of type " <<  ( *atom_types.lock() )[ partner_pose->residue( current_ligand.ligand_resnum ).atom( ati ).type() ].atom_type_name() << std::endl;
				if ( ( *atom_types.lock() )[ partner_pose->residue( current_ligand.ligand_resnum ).atom( ati ).type() ].atom_type_name() != "VIRT" ) {
					lig_atoms.push_back( res_atoms[ ati ] );
				}
			}
			ligand->set_atom_vector( lig_atoms );
			ligand->set_type( data_storage::ligand );
			ligand->set_alignment_atoms( current_ligand.alignment_atoms_num );

			if ( current_ligand.ideal_contacts_num.size() != 0 ) {
				//The entry at this value is a map of ligand atom name to an IdealContact with all other info set
				for ( std::pair< core::Size, IdealContact > contact: current_ligand.ideal_contacts_num ) {
					//contact.second.ligand_atom = pose.residue( ligand_index ).atom_index( contact.first );
					ligand->add_ideal_contact( contact.second );
				}
			}
			//Initialize attached_ligand map
			//Has the ligand already been attached to THIS segment?
			std::map< core::Size, bool > attached_ligand;
			for ( std::pair< core::Size, data_storage::SmartSegmentOP > seg: pdb_segments ) {
				attached_ligand[ seg.first ] = false;
			}

			//Begin auto detect contacts (if applicable)
			if ( current_ligand.auto_detect_contacts ) {
				core::conformation::Residue const & ligand_res = combined_pose.residue( current_ligand.ligand_resnum );

				//Begin bonded contacts
				//Iterate over all ligand atoms
				for ( core::Size ligand_atom_num = 1; ligand_atom_num <= ligand_res.natoms(); ++ligand_atom_num ) {
					if ( ( *atom_types.lock() )[ ligand_res.atom( ligand_atom_num ).type() ].atom_type_name() == "VIRT" ) {
						TR << "Skipping virtual ligand atom!" << std::endl;
						continue;
					}
					utility::vector1< core::id::AtomID > bonded_atoms = combined_pose.conformation().bonded_neighbor_all_res( core::id::AtomID( ligand_atom_num, current_ligand.ligand_resnum ) );

					for ( core::id::AtomID atom_id: bonded_atoms ) {
						TR << "Atom of type " <<  ( *atom_types.lock() )[ combined_pose.residue( atom_id.rsd() ).atom( atom_id.atomno() ).type() ].atom_type_name() << std::endl;
						if ( ( *atom_types.lock() )[ combined_pose.residue( atom_id.rsd() ).atom( atom_id.atomno() ).type() ].atom_type_name() != "VIRT" ) {

							//Make a new ContactDescription
							data_storage::ContactDescription new_contact;
							new_contact.ligand_atom_num = ligand_atom_num;
							new_contact.contact_atom_num = atom_id.atomno();
							if ( atom_id.rsd() <= partner_res ) {
								new_contact.partner_contact = true;
								new_contact.contact_resnum = atom_id.rsd();
							} else {
								new_contact.partner_contact = false;
								new_contact.contact_resnum = atom_id.rsd() - partner_res;
							}
							//We won't need atom names
							current_ligand.ligand_contacts.push_back( new_contact );
						}
					}//end for bonded atoms
				}//end for ligand atoms

				//Begin metal contacts
				std::map< core::Size, utility::vector1< core::id::AtomID > > metal_contacts;
				metal_contacts = core::util::find_metalbinding_atoms_for_complex( combined_pose, current_ligand.ligand_resnum, 1 );
				for ( std::pair< core::Size, utility::vector1< core::id::AtomID > > metal_atom_contacts: metal_contacts ) {
					for ( core::id::AtomID kk: metal_atom_contacts.second ) {
						//Make a new ContactDescription
						data_storage::ContactDescription new_contact;
						new_contact.ligand_atom_num = metal_atom_contacts.first;
						new_contact.contact_atom_num = kk.atomno();
						if ( kk.rsd() <= partner_res ) {
							new_contact.partner_contact = true;
							new_contact.contact_resnum = kk.rsd();
						} else {
							new_contact.partner_contact = false;
							new_contact.contact_resnum = kk.rsd() - partner_res;
						}
						//We won't need atom names
						current_ligand.ligand_contacts.push_back( new_contact );
					}//End for atom id
				}//end for metal contact
			}//end if auto_detect_contacts
			if ( current_ligand.ligand_contacts.size() == 0 ) {
				utility_exit_with_message( "No contacts found for ligand " + utility::to_string( current_ligand.ligand_id ) +  ". Did you forget to auto detect contacts?" );
			}
			//Begin adding ligands to assembly or to partner_ligands and attaching them to segments (if applicable)
			//Iterate over contacts
			for ( data_storage::ContactDescription & contact: current_ligand.ligand_contacts ) {
				//Internal contact
				//Resnums match and partner status matches
				if ( contact.contact_resnum == current_ligand.ligand_resnum && ((current_ligand.partner_ligand && contact.partner_contact ) || (!current_ligand.partner_ligand && !contact.partner_contact ) ) ) {
					TR.Debug << "Adding internal contact for ligand " << ligand->get_ligand_id() << std::endl;
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( 0, 0, contact.contact_atom_num, contact.ligand_atom_num ) ) );
					continue; //No need to try to attach it to a segment based on internal contacts
				}//End if internal contact
				//Non-internal partner contacts
				if ( contact.partner_contact ) {
					TR.Debug << "Adding partner contact for ligand " << ligand->get_ligand_id() << std::endl;
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( 0, contact.contact_resnum, contact.contact_atom_num, contact.ligand_atom_num ) ) );
					continue; //No need to try to attach it to a segment based on partner contacts
				}
				//Begin non-internal assembly contacts
				core::Size passed_residues = 0;
				for ( std::pair< core::Size, data_storage::SmartSegmentOP > seg: pdb_segments ) {
					data_storage::LigandSegmentOP ligseg = std::dynamic_pointer_cast< data_storage::LigandSegment>( seg.second );
					TR.Debug << "Vital residues for segment " << ligseg << std::endl;
					for ( core::Size vital: ligseg->get_vital_residues() ) {
						TR.Debug << vital << " ";
					}
					TR.Debug << std::endl;
					if ( contact.contact_resnum > ligseg->get_length() + passed_residues ) {
						passed_residues += seg.second->get_length();
						continue; //Check the next segment
					}
					if ( contact.contact_resnum < passed_residues ) {
						utility_exit_with_message( "Error when adding ligand to pose segments! Check that your indices are all in range." );
					}
					//Otherwise this contact must be in this segment
					data_storage::SmartSewingResidueCOP contact_res = ligseg->get_residue( contact.contact_resnum - passed_residues );


					//Partner contacts will not have owners! They are free elves!
					if ( !attached_ligand.at( seg.first ) ) {
						ligseg->attach_ligand( ligand, false );
						attached_ligand[ seg.first ] = true;
					} //end if not attached
					runtime_assert( contact.contact_resnum > passed_residues );
					ligseg->add_ligand_contact( contact.contact_resnum - passed_residues );
					TR.Debug << "Vital residues for segment " << ligseg << std::endl;
					for ( core::Size vital: ligseg->get_vital_residues() ) {
						TR.Debug << vital << " ";
					}
					TR.Debug << std::endl;
					//Add the contact residue to the ligand
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( seg.first, (contact.contact_resnum - passed_residues ), contact.contact_atom_num, contact.ligand_atom_num ) ) );
					//Now add to partner ligands
					partner_ligands[ ligand->get_ligand_id() ] = ligand;
					break; //We don't need to keep looking through additional segments
				}//End iterate over pdb_segments
			}//end for current contacts
		} else { //end if partner ligand
			ligand->set_partner_ligand( false );
			if ( current_ligand.ligand_resnum > pose.total_residue() || current_ligand.ligand_resnum < 1 ) {
				TR.Debug << "Ligand index " << current_ligand.ligand_id << " is invalid!" << std::endl;
				utility_exit_with_message( "Invalid ligand resnum in pose!" );
				continue;
			}
			ligand->set_amino_acid_type( pose.residue( current_ligand.ligand_resnum ).type().base_name() );
			ligand->set_full_type_name( pose.residue( current_ligand.ligand_resnum ).type().name() );
			//We want to include all of the ligand atoms
			res_atoms = utility::vector1< core::conformation::Atom >( pose.residue( current_ligand.ligand_resnum ).atoms().begin(), pose.residue( current_ligand.ligand_resnum ).atoms().end() );
			utility::vector1< core::conformation::Atom > lig_atoms;
			for ( core::Size ati = 1; ati <= res_atoms.size(); ++ati ) {
				TR << "Atom of type " <<  ( *atom_types.lock() )[ pose.residue( current_ligand.ligand_resnum ).atom( ati ).type() ].atom_type_name() << std::endl;
				if ( ( *atom_types.lock() )[ pose.residue( current_ligand.ligand_resnum ).atom( ati ).type() ].atom_type_name() != "VIRT" ) {
					lig_atoms.push_back( res_atoms[ ati ] );
				}
			}
			ligand->set_atom_vector( lig_atoms );
			ligand->set_type( data_storage::ligand );
			ligand->set_alignment_atoms( current_ligand.alignment_atoms_num );

			if ( current_ligand.ideal_contacts_num.size() != 0 ) {
				//The entry at this value is a map of ligand atom name to an IdealContact with all other info set
				for ( std::pair< core::Size, IdealContact > contact: current_ligand.ideal_contacts_num ) {
					//contact.second.ligand_atom = pose.residue( ligand_index ).atom_index( contact.first );
					ligand->add_ideal_contact( contact.second );
				}
			}
			//Initialize attached_ligand map
			//Has the ligand already been attached to THIS segment?
			std::map< core::Size, bool > attached_ligand;
			for ( std::pair< core::Size, data_storage::SmartSegmentOP > seg: pdb_segments ) {
				attached_ligand[ seg.first ] = false;
			}


			//Begin auto detect contacts ( if applicable)
			if ( current_ligand.auto_detect_contacts ) {
				core::conformation::Residue const & ligand_res = combined_pose.residue( current_ligand.ligand_resnum + partner_res );

				//Begin bonded contacts
				//Iterate over all ligand atoms
				for ( core::Size ligand_atom_num = 1; ligand_atom_num <= ligand_res.natoms(); ++ligand_atom_num ) {
					if ( ( *atom_types.lock() )[ ligand_res.atom( ligand_atom_num ).type() ].atom_type_name() == "VIRT" ) {
						TR << "Skipping virtual ligand atom!" << std::endl;
						continue;
					}
					utility::vector1< core::id::AtomID > bonded_atoms = combined_pose.conformation().bonded_neighbor_all_res( core::id::AtomID( ligand_atom_num, current_ligand.ligand_resnum + partner_res ) );

					for ( core::id::AtomID atom_id: bonded_atoms ) {
						TR << "Atom of type " <<  ( *atom_types.lock() )[ combined_pose.residue( atom_id.rsd() ).atom( atom_id.atomno() ).type() ].atom_type_name() << std::endl;
						if ( ( *atom_types.lock() )[ combined_pose.residue( atom_id.rsd() ).atom( atom_id.atomno() ).type() ].atom_type_name() != "VIRT" ) {

							//Make a new ContactDescription
							data_storage::ContactDescription new_contact;
							new_contact.ligand_atom_num = ligand_atom_num;
							new_contact.contact_atom_num = atom_id.atomno();
							if ( atom_id.rsd() <= partner_res ) {
								new_contact.partner_contact = true;
								new_contact.contact_resnum = atom_id.rsd();
							} else {
								new_contact.partner_contact = false;
								new_contact.contact_resnum = atom_id.rsd() - partner_res;
							}
							//We won't need atom names
							current_ligand.ligand_contacts.push_back( new_contact );
						}
					}//end for bonded atoms
				}//end for ligand atoms

				//Begin metal contacts
				std::map< core::Size, utility::vector1< core::id::AtomID > > metal_contacts;
				metal_contacts = core::util::find_metalbinding_atoms_for_complex( combined_pose, current_ligand.ligand_resnum + partner_res, 1 );
				for ( std::pair< core::Size, utility::vector1< core::id::AtomID > > metal_atom_contacts: metal_contacts ) {
					for ( core::id::AtomID kk: metal_atom_contacts.second ) {
						//Make a new ContactDescription
						data_storage::ContactDescription new_contact;
						new_contact.ligand_atom_num = metal_atom_contacts.first;
						new_contact.contact_atom_num = kk.atomno();
						if ( kk.rsd() <= partner_res ) {
							new_contact.partner_contact = true;
							new_contact.contact_resnum = kk.rsd();
						} else {
							new_contact.partner_contact = false;
							new_contact.contact_resnum = kk.rsd() - partner_res;
						}
						//We won't need atom names
						current_ligand.ligand_contacts.push_back( new_contact );
					}//End for atom id
				}//end for metal contact
			}//end if auto_detect_contacts
			if ( current_ligand.ligand_contacts.size() == 0 && !current_ligand.partner_ligand ) {
				utility_exit_with_message( "No contacts found for ligand " + utility::to_string( current_ligand.ligand_id ) +  ". Did you forget to auto detect contacts?" );
			}
			//Begin adding ligands to assembly or to partner_ligands and attaching them to segments (if applicable)
			//Iterate over contacts
			for ( data_storage::ContactDescription & contact: current_ligand.ligand_contacts ) {
				//Internal contact
				//Resnums match and partner status matches
				if ( contact.contact_resnum == current_ligand.ligand_resnum && ((current_ligand.partner_ligand && contact.partner_contact ) || (!current_ligand.partner_ligand && !contact.partner_contact ) ) ) {
					TR.Debug << "Adding internal contact for ligand " << ligand->get_ligand_id() << std::endl;
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( 0, 0, contact.contact_atom_num, contact.ligand_atom_num ) ) );
					continue; //No need to try to attach it to a segment based on internal contacts
				}//End if internal contact
				//Non-internal partner contacts
				if ( contact.partner_contact ) {
					TR.Debug << "Adding partner contact for ligand " << ligand->get_ligand_id() << std::endl;
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( 0, contact.contact_resnum, contact.contact_atom_num, contact.ligand_atom_num ) ) );
					continue; //No need to try to attach it to a segment based on partner contacts
				}
				//Begin non-internal assembly contacts
				core::Size passed_residues = 0;
				for ( std::pair< core::Size, data_storage::SmartSegmentOP > seg: pdb_segments ) {
					data_storage::LigandSegmentOP ligseg = std::dynamic_pointer_cast< data_storage::LigandSegment>( seg.second );
					TR.Debug << "Vital residues for segment " << ligseg << std::endl;
					for ( core::Size vital: ligseg->get_vital_residues() ) {
						TR.Debug << vital << " ";
					}
					TR.Debug << std::endl;
					if ( contact.contact_resnum > ligseg->get_length() + passed_residues ) {
						passed_residues += seg.second->get_length();
						continue; //Check the next segment
					}
					if ( contact.contact_resnum < passed_residues ) {
						utility_exit_with_message( "Error when adding ligand to pose segments! Check that your indices are all in range." );
					}
					//Otherwise this contact must be in this segment
					data_storage::SmartSewingResidueCOP contact_res = ligseg->get_residue( contact.contact_resnum - passed_residues );
					if ( ligand->get_nonconst_owner_segment() == nullptr ) { //We're attaching the ligand at the first available contact.
						ligand->set_owner_segment( ligseg );
						//Add the contact residue to the ligand
						ligseg->attach_ligand( ligand, true );
						attached_ligand[ seg.first ] = true;
						//added_ligand = true;
					} else if ( !attached_ligand.at( seg.first ) ) { //end if no owner
						ligseg->attach_ligand( ligand, false );
						attached_ligand[ seg.first ] = true;
					} //end if not attached
					runtime_assert( contact.contact_resnum > passed_residues );
					ligseg->add_ligand_contact( contact.contact_resnum - passed_residues );
					TR.Debug << "Vital residues for segment " << ligseg << std::endl;
					for ( core::Size vital: ligseg->get_vital_residues() ) {
						TR.Debug << vital << " ";
					}
					TR.Debug << std::endl;
					//Add the contact residue to the ligand
					ligand->add_contact( data_storage::LigandContactOP( new data_storage::LigandContact( seg.first, (contact.contact_resnum - passed_residues ), contact.contact_atom_num, contact.ligand_atom_num ) ) );
					break; //We don't need to keep looking through additional segments
				}//End iterate over pdb_segments
			}//end for current contacts
		}//End else (not partner ligand)
	}//End iterate over ligands

	if ( added_segments ==1 ) {
		TR << "Only one added segment. Forcing it to be helical." << std::endl;
		pdb_segments[last_added_seg_id]->set_dssp_code('H');
	}

	//Caveat--we don't want our terminal segments to be loops
	TR << "Removing terminal segments that are loops:" << std::endl;
	//First move backward from c terminus
	data_storage::SmartSegmentOP check_segment;
	while ( added_segments > 0 && pdb_segments[ last_added_seg_id ]->get_dssp_code() == 'L'  ) {

		TR.Debug << "Vital residues for terminal loop " << last_added_seg_id << std::endl;
		for ( core::Size vital: pdb_segments[ last_added_seg_id ]->get_vital_residues() ) {
			TR.Debug << vital << " ";
		}
		TR.Debug << std::endl;


		//If this segment contains the ligand and/or a vital ligand contact, break out and give a warning
		if ( ligand_resnums_final.size() != 0 && ( std::dynamic_pointer_cast< data_storage::LigandSegment >( pdb_segments.at( last_added_seg_id ) )->get_ligand_residues().size() != 0 || std::dynamic_pointer_cast< data_storage::LigandSegment >( pdb_segments.at( last_added_seg_id ) )->get_ligand_contact_indices().size() != 0 ) ) {
			TR << "WARNING: Vital terminal segment (containing ligand or ligand contacts) is a loop! Segments cannot be added to this terminus." << std::endl;
			break;
		}

		TR.Debug << "Found a terminal loop: ";
		std::map< core::Size, data_storage::SmartSegmentOP>::iterator it  = pdb_segments.find( last_added_seg_id );
		if ( it == pdb_segments.end() ) {
			utility_exit_with_message( "Can't find segment to delete!" );
		}
		pdb_segments.erase( it );
		TR.Debug << "Removed.";
		if ( pdb_segments.count( last_added_seg_id - 1 ) != 0 ) {
			pdb_segments[ last_added_seg_id - 1 ]->set_c_terminal_neighbor( nullptr );
		}
		if ( pdb_segments.count( last_added_seg_id + 1 ) != 0 ) {
			pdb_segments[ last_added_seg_id + 1 ]->set_n_terminal_neighbor( nullptr );
		}
		--added_segments;
		--last_added_seg_id;
	}
	//Then move forward from n terminus
	//Hopefully we took care of these above, but this is to double check
	core::Size first_segid = starting_segs + 1;
	while ( added_segments > 0 && pdb_segments[ first_segid ]->get_dssp_code() == 'L'  ) {

		TR.Debug << "Vital residues for terminal loop " << first_segid << std::endl;
		for ( core::Size vital: pdb_segments[ first_segid ]->get_vital_residues() ) {
			TR.Debug << vital << " ";
		}
		TR.Debug << std::endl;



		if ( ligand_resnums_final.size() != 0 && ( std::dynamic_pointer_cast< data_storage::LigandSegment >( pdb_segments.at( first_segid ) )->get_ligand_residues().size() != 0 || std::dynamic_pointer_cast< data_storage::LigandSegment >( pdb_segments.at( first_segid ) )->get_ligand_contact_indices().size() != 0 ) ) {
			TR << "WARNING: Vital terminal segment (containing ligand or ligand contacts) is a loop! Segments cannot be added to this terminus." << std::endl;
			break;
		}
		TR.Debug << "Found a terminal loop: SegID " << first_segid << std::endl;
		std::map< core::Size, data_storage::SmartSegmentOP>::iterator it  = pdb_segments.find( first_segid );
		if ( it == pdb_segments.end() ) {
			utility_exit_with_message( "Can't find segment to delete!" );
		}
		pdb_segments.erase( it );
		TR.Debug << "Removed." << std::endl;
		if ( pdb_segments.count( first_segid + 1 ) != 0 ) {
			check_segment = pdb_segments[ first_segid + 1 ];
			check_segment->set_n_terminal_neighbor( nullptr );
		}
		if ( pdb_segments.count( first_segid - 1 ) != 0 ) {
			pdb_segments[ first_segid - 1 ]->set_c_terminal_neighbor( nullptr );
		}
		--added_segments;
		++first_segid;
	}
	if ( added_segments ) {
		TR << "N terminal segment DSSP: " << pdb_segments.find( first_segid )->second->get_dssp_code() << std::endl;
		TR << "C terminal segment DSSP: " << pdb_segments.find( last_added_seg_id )->second->get_dssp_code() << std::endl;
	}
	return added_segments;
}//end add_pose_segments_to_segment_vector

//Getters
BasisMapGeneratorOP AlignmentFileGeneratorMover::basis_map_generator() const
{
	return basis_map_generator_;
}


std::string
AlignmentFileGeneratorMover::get_required_resnums() const{
	return required_resnums_;
}

core::select::residue_selector::ResidueSelectorCOP
AlignmentFileGeneratorMover::get_required_selector() const{
	return required_selector_;
}

//Setters
void
AlignmentFileGeneratorMover::set_required_resnums( std::string req){
	required_resnums_ = req;
}
void
AlignmentFileGeneratorMover::set_required_selector( core::select::residue_selector::ResidueSelectorCOP sel){
	required_selector_ = sel;
}

void
AlignmentFileGeneratorMover::basis_map_generator( BasisMapGeneratorOP bmg)
{
	basis_map_generator_ = bmg;
}

void
AlignmentFileGeneratorMover::alignment_settings_from_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose,
	BasisMapGeneratorOP bmg)
{
	//recursive_depth
	utility::vector1< core::Size> pose_segment_starts;
	utility::vector1< core::Size> pose_segment_ends;
	utility::vector1< core::Size> match_segments;
	core::Size recursive_depth = tag->getOption< core::Size>( "max_recursion", 1 );

	if ( tag->hasOption( "pose_segment_starts" ) ) {
		std::string str_seg_starts = tag->getOption< std::string >( "pose_segment_starts" );
		pose_segment_starts = core::pose::get_resnum_list_ordered( str_seg_starts, pose );
	} else {
		pose_segment_starts.clear();
		pose_segment_starts.push_back( 1 );
	}
	//Same for pose_segment_ends
	if ( tag->hasOption( "pose_segment_ends" ) ) {
		std::string str_seg_ends = tag->getOption< std::string >( "pose_segment_ends" );
		pose_segment_ends = core::pose::get_resnum_list_ordered( str_seg_ends, pose );
	} else {
		pose_segment_ends.clear();
		pose_segment_ends.push_back( pose.total_residue() );
	}
	//This one's a little different--not resnums, so we don't have a handy
	//built-in function available
	if ( tag->hasOption( "match_segments" ) ) {
		std::string str_ms = tag->getOption< std::string>( "match_segments" );
		utility::vector1< std::string > str_vec_ms( utility::string_split( str_ms, ',' ) );
		//Now we will loop through all of the strings and convert to core::Size
		for ( core::Size i = 1; i <= str_vec_ms.size(); ++i ) {
			match_segments.push_back( core::Size( std::stoi( str_vec_ms[ i ] ) ) );
		}
	} else {
		//If none provided, it will be the first and last segments
		match_segments.push_back( 1 );
		if ( pose_segment_starts.size() > 1 ) {
			match_segments.push_back( pose_segment_starts.size() );
		}
	}
	bmg->set_alignment_settings( match_segments, pose_segment_starts, pose_segment_ends, recursive_depth );

}


void
AlignmentFileGeneratorMover::parse_ligand_conformers(
	utility::vector1< data_storage::LigandDescription > & ligands
)
{
	for ( data_storage::LigandDescription & ligdes: ligands ) {
		if ( ligdes.pdb_conformers_string == "" ) {
			continue;
		}
		//Open the file
		//This function also removes leading/trailing whitespace, so this function should work
		utility::vector1< std::string > filenames = utility::io::get_lines_from_file_data( ligdes.pdb_conformers_string );
		//core::import_pose::read_all_poses takes a list of filenames and fills a reference vector1 of PoseOPs
		utility::vector1< core::pose::PoseOP > ligposes;
		core::import_pose::read_all_poses( filenames, ligposes );
		//Iterate through ligposes and get all the residues
		for ( core::pose::PoseOP ligpose: ligposes ) {
			for ( core::Size i = 1; i <= ligpose->total_residue(); ++i ) {
				ligdes.pdb_conformers.push_back( ligpose->residue( i ).clone() );
			}
		}
	}
}

/*
std::map< core::Size, std::string >
AlignmentFileGeneratorMover::map_ligand_id_to_resnum_string( core::pose::Pose const & pose, utility::vector1< std::string > const & resnum_strings,  std::map< std::string, core::select::residue_selector::ResidueSelectorCOP > const & ligand_selectors ){
core::select::residue_selector::ResidueSubset selection( pose.total_residue(), false );
std::map< core::Size, std::string > resnums;
std::map< core::Size, std::string > output;
for( std::string resnum_string: resnum_strings ){
if( ligand_selectors.count( resnum_string ) != 0 && ligand_selectors.find( resnum_string )->second != nullptr ){
selection = ligand_selectors.find( resnum_string )->second->apply( pose );
for ( core::Size resi = 1; resi <= pose.total_residue(); ++resi ) {
if ( selection[resi] ){
resnums[ resi ] = resnum_string;
}//end if
}//end for
}//end if selector
else{
core::Size resi = core::pose::parse_resnum( resnum_string, pose, true );
resnums[ resi ] = resnum_string;
}
}//end for
//Now fix output
core::Size ligID = 1;
for( std::pair< core::Size, std::string > ligpair: resnums ){
output[ ligID ] = ligpair.second;
++ligID;
}
return output;
}
*/


/////////////// Creator ///////////////

protocols::moves::MoverOP
AlignmentFileGeneratorMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AlignmentFileGeneratorMover );
}

std::string
AlignmentFileGeneratorMoverCreator::keyname() const {
	return AlignmentFileGeneratorMover::class_name();
}

std::string
AlignmentFileGeneratorMoverCreator::mover_name() {
	return "AlignmentFileGeneratorMover";
}

} //hashing
} //sewing
} //protocols

