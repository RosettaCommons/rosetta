// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ContactMap.cc
///
/// @brief Method code and full headers for ContactMapCreator, ContactMap, ContactPartner
///   and contact
///
/// @author Joerg Schaarschmidt

// unit headers
#include <protocols/contact_map/ContactMap.hh>
#include <protocols/contact_map/ContactMapCreator.hh>

// type headers
#include <core/types.hh>

// project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/ref_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/util.hh>
#include <protocols/jd2/util.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>

// ObjexxFCL headers

// C++ headers

#include <utility/vector0.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/conformation/Conformation.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace contact_map {

using namespace utility::tag;
using namespace core;

static basic::Tracer TR( "protocols.moves.ContactMap" );


///////////////////////////////  ContactMapCreator  ///////////////////////////////





///////////////////////////////  ContactMap  ///////////////////////////////

/// @brief Default constructor
ContactMap::ContactMap() :
	Mover(ContactMap::mover_name()),
	output_prefix_("contact_map_"),
	n_poses_(0),
	distance_cutoff_(10.0),
	models_per_file_(1),
	reset_count_(true),
	row_format_(false),
	distance_matrix_(false),
	region1_( utility::pointer::make_shared< core::select::residue_selector::TrueResidueSelector >() )
{
}

/// @brief Copy constructor
ContactMap::ContactMap(ContactMap const & ) = default;

/// @brief Destructor
ContactMap::~ContactMap() = default;


moves::MoverOP ContactMap::clone() const {
	return utility::pointer::make_shared< ContactMap >(*this);
}

moves::MoverOP ContactMap::fresh_instance() const {
	return utility::pointer::make_shared< ContactMap >();
}


/// @brief Processes options specified in xml-file and sets up the ContactMap
void ContactMap::parse_my_tag(TagCOP const tag, basic::datacache::DataMap & data
) {

	// 'distance_cutoff' option
	set_distance_cutoff(tag->getOption<core::Real>("distance_cutoff", distance_cutoff_));

	// 'output_prefix_' option
	set_output_prefix(tag->getOption<std::string>("prefix", output_prefix_));

	// 'reset_count' flag
	reset_count_ = tag->getOption<bool>("reset_count", true);

	// 'models_per_file' option
	models_per_file_ = tag->getOption<core::Size>("models_per_file", models_per_file_);

	// 'row_format_' flag
	row_format_ = tag->getOption<bool>("row_format", false);

	// 'distance_matrix_' flag
	distance_matrix_ = tag->getOption<bool>("distance_matrix", false);

	/////////////////////////////
	// 'region1', 'region2' and 'ligand' options

	reference_pose_ = protocols::rosetta_scripts::legacy_saved_pose_or_input( tag, data, mover_name(), /*use_native*/ false );

	// Parse 'region1' option if supplied
	// Defaults to complete pose if no option supplied.
	if ( tag->hasOption("region1") ) {
		region1_ = parse_region_string(tag->getOption<std::string>("region1"));
	}
	if ( tag->hasOption("region2") ) {
		// Parse 'region2' option if supplied and initialize ContactMap between 'region1' and 'region2'
		region2_ = parse_region_string(tag->getOption<std::string> ("region2"));
	} else if ( tag->hasOption("ligand") ) {
		// Parse 'ligand' option if supplied and initialize ContactMap between 'region1' and 'ligand'
		region2_ = parse_region_string(tag->getOption<std::string>("ligand"));
		region2_all_atom_ = true;
	}

}

core::select::residue_selector::ResidueSelectorCOP
ContactMap::parse_region_string(std::string const & region_def) const {

	// Check if region is defined by chain
	if ( region_def.size() == 1 && isalpha(region_def[0]) ) {
		return utility::pointer::make_shared< core::select::residue_selector::ChainSelector >(region_def[0]);
	}

	return utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( region_def );
}

/// @brief Initializes ContactMap within a single region
void ContactMap::fill_contacts(
	core::select::residue_selector::ResidueSelector const & region,
	Pose const & pose
) {

	ContactPartner p1, p2;

	utility::vector1< core::Size > residues = core::select::get_residues_from_subset( region.apply( pose ) );

	// Outer loop to assign ContactPartner1
	for ( core::Size ii(1); ii <= residues.size(); ++ii ) {
		core::Size res1( residues[ii] );
		std::string const & atom_name1 =  pose.residue_type(res1).name1() == 'G'  ?   "CA" : "CB";
		std::string const & resname1 = pose.residue(res1).name3();
		p1 = ContactPartner(res1, resname1, atom_name1);
		// Fill column and row names
		row_names_.push_back(p1.string_rep());
		column_names_.push_back(p1.string_rep());
		// Inner loop to assign ContactPartner2 and create contact
		for ( core::Size jj(ii+1); jj <= residues.size(); ++jj ) {
			core::Size res2( residues[jj] );
			std::string const & atom_name2 = pose.residue_type(res2).name1() == 'G'  ?   "CA" : "CB";
			std::string const & resname2 = pose.residue(res2).name3();
			p2 = ContactPartner(res2, resname2, atom_name2);
			// Create contact and add to contacts_
			Contact contact = Contact(p1,p2);
			contacts_.push_back(contact);
		}
	}
	// Add additional contact with same ContactPartner to fill identity contacts in Matrix
	contacts_.push_back(Contact(p1, p1));

	core::Size const length = residues.size();
	// Fill output_matrix with positions of corresponding contact in contacts_ vector
	for ( core::Size i = 1; i<=length; i++ ) {
		for ( core::Size j = 1; j<=length; j++ ) {
			if ( i==j ) {
				// Set matrix value to position of identity contact
				output_matrix_.push_back(contacts_.size());
			} else {
				// Set matrix value to corresponding position in contacts_ vector
				output_matrix_.push_back(i<j ?
					(2 * length - i) * (i - 1) / 2 - i + j :
					(2 * length - j) * (j - 1) / 2 - j + i);
			}
		}
	}
}

/// @brief Initializes ContactMap between two separate regions
void ContactMap::fill_contacts(
	core::select::residue_selector::ResidueSelector const & region1,
	core::select::residue_selector::ResidueSelector const & region2,
	Pose const & pose
) {
	core::Size matrix_position = 1;
	utility::vector1< core::Size > residues1 = core::select::get_residues_from_subset( region1.apply( pose ) );
	utility::vector1< core::Size > residues2 = core::select::get_residues_from_subset( region2.apply( pose ) );

	// Loop over residues in region1
	for ( core::Size res1: residues1 ) {
		// Assign ContactPartner1
		std::string atom_name1 = pose.residue_type(res1).name1() == 'G' ? "CA" : "CB";
		std::string resname1 = pose.residue(res1).name3();
		ContactPartner p1(res1, resname1, atom_name1);
		// Add ContactPartner1 string to 'row_names_'
		row_names_.push_back(p1.string_rep());
		// Loop over residues in region2
		for ( core::Size res2: residues2 ) {
			// Assign ContactPartner2
			std::string atom_name2 = pose.residue_type(res2).name1() == 'G' ? "CA" : "CB";
			std::string resname2 = pose.residue(res2).name3();
			ContactPartner p2(res2, resname2, atom_name2);
			// Create contact and add to 'contacts_'
			Contact contact = Contact(p1, p2);
			contacts_.push_back(contact);
			// Add ContactPartner2 string to 'column_names_'
			if ( res1 == residues1[1] ) {
				column_names_.push_back(p2.string_rep());
			}
			// Set matrix value to corresponding position in contacts_ vector
			output_matrix_.push_back(matrix_position++);
		}
	}
}

/// @brief Initializes ContactMap between a single region and a ligand
void ContactMap::fill_contacts_all_atom2(
	core::select::residue_selector::ResidueSelector const & region1,
	core::select::residue_selector::ResidueSelector const & region2,
	Pose const & pose
) {
	core::Size matrix_position = 1;

	utility::vector1< core::Size > residues1 = core::select::get_residues_from_subset( region1.apply( pose ) );
	utility::vector1< core::Size > residues2 = core::select::get_residues_from_subset( region2.apply( pose ) );
	if ( residues2.size() != 1 ) {
		utility_exit_with_message("In ContactMap mover, the ligand is more than one residue.");
	}

	core::Size ligand_seqpos = residues2[1];
	std::string ligand_resname = pose.residue(ligand_seqpos).name3();
	// Loop over residues in specified region
	for ( core::Size res1: residues1 ) {
		// Exclude ligand from residues
		if ( res1 == ligand_seqpos ) continue;
		std::string atom_name1 = pose.residue_type(res1).name1() == 'G' ? "CA" : "CB";
		std::string resname1 = pose.residue(res1).name3();
		ContactPartner p1(res1, resname1, atom_name1);
		// Add residue string to 'row_names_'
		row_names_.push_back(p1.string_rep());
		// Loop over atoms in ligand
		for ( core::Size jj = 1; jj <= pose.residue(ligand_seqpos).atoms().size(); ++jj ) {
			// Skip hydrogen atoms
			if ( pose.residue(ligand_seqpos).atom_is_hydrogen(jj) ) {
				continue;
			}
			// Initialize second partner, create contact and add to 'contacts_'
			std::string atom_name2 = pose.residue(ligand_seqpos).atom_name(jj);
			ContactPartner p2(ligand_seqpos, ligand_resname, atom_name2);
			Contact contact = Contact(p1, p2);
			contacts_.push_back(contact);
			// Add ligand atom string to 'column_names_'
			if ( res1 == residues1[1] ) {
				column_names_.push_back(p2.string_rep());
			}
			// Set matrix value to corresponding position in contacts_ vector
			output_matrix_.push_back(matrix_position++);
		}
	}
}


/// @brief Process supplied pose
void ContactMap::apply(Pose & pose) {

	if ( reference_pose_ != nullptr ) {
		runtime_assert( region1_ != nullptr );
		if ( region2_ == nullptr ) {
			fill_contacts( *region1_, *reference_pose_ );
		} else if ( region2_all_atom_ ) {
			fill_contacts_all_atom2( *region1_, *region2_, *reference_pose_ );
		} else {
			fill_contacts( *region1_, *region2_, *reference_pose_ );
		}

		reference_pose_ = nullptr; // Hacky "run-only-once" trick;
	}

	using namespace pose;
	n_poses_++;

	// Iterate over contacts
	for ( auto & contact : contacts_ ) {
		ContactPartner * p1(contact.partner1());
		ContactPartner * p2(contact.partner2());
		// Get coordinates of both contact partners and calculate distance
		numeric::xyzVector<core::Real> v1 = pose.residue(p1->seqpos()).atom(p1->atomname()).xyz();
		numeric::xyzVector<core::Real> v2 = pose.residue(p2->seqpos()).atom(p2->atomname()).xyz();
		core::Real distance = v1.distance(v2);
		// Add distance to contact if it's below the cutoff value
		if ( distance_matrix_ ) {
			contact.add_distance(distance);
			continue;
		}
		if ( distance <= distance_cutoff_ ) {
			contact.add_distance();
		}
	}

	// Output ContactMap to file if the number of processed poses since last output equals models_per_file_ variable
	if ( models_per_file_ > 0 && n_poses_ % models_per_file_ == 0 ) {
		// If the ContactMap is to be reset, create a unique output name and write the ContactMap to this file
		if ( reset_count_ ) {
			// Get the file name of the output structure.  This will be
			// something_XXXX.pdb if you're outputting to pdbs, something_XXXX if you're outputting to silent files
			std::string structure_output_name(protocols::jd2::current_output_name());
			// Append a prefix and a suffix to get a final filename for the contact map output
			std::string contact_map_file_name(output_prefix_ + structure_output_name + ".csv");
			// Call output function with the generated filename
			write_to_file(contact_map_file_name);
			// Reset ContactMap
			reset();
		} else {
			// If the ContactMap should just be updated, reconstruct the output name based on the current job name and
			// the specified prefix
			std::string current_job_tag(protocols::jd2::current_input_tag());
			// Append a prefix and a suffix to get a final filename for the contact map output
			std::string contact_map_file_name(output_prefix_ + current_job_tag + ".csv");
			// Call output function with the generated filename
			write_to_file(contact_map_file_name);
		}
	}
}

/// @brief Resets the movers n_poses_ variable and the counts of all contacts to 0
void ContactMap::reset(){
	n_poses_ = 0;
	for ( auto & contact : contacts_ ) {
		contact.reset_count();
	}
}

/// @brief Output function that writes the ContactMap to the specified file
void ContactMap::write_to_file(std::string filename) {
	if ( filename == "" ) {
		filename = output_prefix_ + ".csv";
	}
	TR.Info << "Writing ContactMap to '" << filename
		<< "', ContactMap includes " << n_poses_ << " structure(s)."
		<< std::endl;

	// Initialize output stream and write Header
	utility::io::ozstream output_stream;
	output_stream.open(filename, std::ios_base::out);
	write_to_stream( output_stream );

	// Finish up
	output_stream.close();
}

void
ContactMap::write_to_stream(std::ostream & output_stream) {
	output_stream << "# Number of Models:\t" << n_poses_ <<"\tDistance Cutoff:\t"<< distance_cutoff_<< std::endl
		<< std::endl;
	if ( row_format_ ) {
		for ( auto & contact : contacts_ ) {
			output_stream << contact.long_string_rep() << std::endl;
		}
	} else {

		// Print column header line
		for ( core::Size col = 1; col <= column_names_.size(); col++ ) {
			output_stream << "\t" << column_names_[col];
		}
		output_stream << std::endl;

		for ( core::Size row = 1; row <= row_names_.size(); row++ ) {
			// Print row header
			output_stream << row_names_[row];
			// Loop over columns and print corresponding fields
			for ( core::Size col = 1; col <= column_names_.size(); col++ ) {
				output_stream << "\t" << get_contact( row, col ).string_rep();
			}
			output_stream << std::endl;
		}
	}
}

Contact const &
ContactMap::get_contact( core::Size row, core::Size col) {
	core::Size pos_in_contacts = output_matrix_[(row - 1)
		* column_names_.size() + col];
	return contacts_.at(pos_in_contacts);
}

std::string ContactMap::get_name() const {
	return mover_name();
}

std::string ContactMap::mover_name() {
	return "ContactMap";
}

void ContactMap::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	//Format for region strings?
	//Option 1: a single alphabetical character (chain ID)
	//Option 3: a single residue number
	//Option 2: hyphen-separated residue numbers
	XMLSchemaRestriction region_string;
	region_string.name( "contact_map_region_string" );
	region_string.base_type( xs_string );
	std::string pattern = "[A-Za-z]|" + residue_number_string() + "|" + residue_number_string() + "[-]" + residue_number_string();
	region_string.add_restriction( xsr_pattern, pattern );
	xsd.add_top_level_element( region_string );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "distance_cutoff", xsct_real, "Maximum distance between two atoms that will be considered a contact" )
		+ XMLSchemaAttribute( "prefix", xs_string, "Prefix for output filenames" )
		+ XMLSchemaAttribute::attribute_w_default( "reset_count", xsct_rosetta_bool, "Should the count be reset to zero after outputting the contact map to a file?", "true" )
		+ XMLSchemaAttribute( "models_per_file", xsct_non_negative_integer, "Number of models to output per file" )
		+ XMLSchemaAttribute::attribute_w_default( "row_format", xsct_rosetta_bool, "Should the output be in row format instead of matrix format?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "distance_matrix", xsct_rosetta_bool, "Store distances of contacts in matrix", "false" )
		+ XMLSchemaAttribute( "region1", "contact_map_region_string", "Region definition for region1 of the contact map in format start-end or chainID" )
		+ XMLSchemaAttribute( "region2", "contact_map_region_string", "Region definition for region2 of the contact map" )
		+ XMLSchemaAttribute( "ligand", "contact_map_region_string", "Sequence position or chainID of a ligand" );

	core::pose::attributes_for_saved_reference_pose_w_description( attlist, "Reference pose to use when building the initial contact map description (defaults to input pose)." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Generates a contact map between specified regions", attlist );

}

std::string ContactMapCreator::keyname() const {
	return ContactMap::mover_name();
}

protocols::moves::MoverOP
ContactMapCreator::create_mover() const {
	return utility::pointer::make_shared< ContactMap >();
}

void ContactMapCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ContactMap::provide_xml_schema( xsd );
}


///////////////////////////////  ContactPartner  ///////////////////////////////

std::string ContactPartner::string_rep() const {
	std::ostringstream oss;
	oss << resname_ << seqpos_;
	//std::string test;
	// Append atomname if it's not the default
	if ( atomname_ != "CA" && atomname_ != "CB" ) {
		oss<< "-"<<atomname_;
	}
	return oss.str();
}


///////////////////////////////  Contact  ///////////////////////////////

/// @brief Adds distance to the contact
void Contact::add_distance(){
	++count_;
}

/// @brief Adds distance to the contact
void Contact::add_distance(core::Real dist){
	distance_ = dist;
}

/// @brief Resets count to 0
void Contact::reset_count(){
	count_ = 0;
}

/// @brief Returns string representation of the Contact
std::string Contact::string_rep() const {
	std::ostringstream oss;
	// oss << partner1_.string_rep() <<"|"<<partner2_.string_rep();
	if ( distance_ == 0.0 ) { oss << count_;}
	else {
		oss.precision(3);
		oss.setf( std::ios::fixed , std::ios::floatfield);
		oss << distance_;
	}
	return oss.str();
}

/// @brief Returns string representation of the Contact as percentage value
std::string Contact::string_rep(core::Size n_poses) const {
	std::ostringstream oss;
	oss.precision(3);
	oss.setf( std::ios::fixed , std::ios::floatfield);
	core::Real percentage = core::Real(count_) / core::Real(n_poses);
	if ( distance_ == 0.0 ) oss << percentage;
	else oss << distance_;
	return oss.str();
}


/// @brief Returns string representation of the Contact including partner names
std::string Contact::long_string_rep(core::Size n_poses) const {
	std::string stringrep;
	stringrep = partner1_.string_rep() + "\t"+ partner2_.string_rep() + "\t";
	stringrep = stringrep + ( (n_poses==0) ? string_rep() :string_rep(n_poses));
	return stringrep;
}

} // moves
} // protocols
