// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/chemical/ResidueType.hh>
#include <protocols/jd2/JobDistributor.hh>

// utility headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>

// ObjexxFCL headers

// C++ headers

#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/conformation/Conformation.hh>

namespace protocols {
namespace contact_map {

using namespace utility::tag;
using namespace core;

using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.moves.ContactMap" );


///////////////////////////////  ContactMapCreator  ///////////////////////////////

std::string ContactMapCreator::keyname() const {
	return ContactMapCreator::mover_name();
}

protocols::moves::MoverOP ContactMapCreator::create_mover() const {
	return protocols::moves::MoverOP( new ContactMap );
}

std::string ContactMapCreator::mover_name() {
	return "ContactMap";
}


///////////////////////////////  ContactMap  ///////////////////////////////

/// @brief Default constructor
ContactMap::ContactMap() :
	Mover(ContactMapCreator::mover_name()),
	output_prefix_("contact_map_"),
	n_poses_(0),
	distance_cutoff_(10.0),
	models_per_file_(1),
	reset_count_(true),
	row_format_(false),
	distance_matrix_(false)
{
}

/// @brief Copy constructor
ContactMap::ContactMap(ContactMap const & contact_map) :
	Mover(contact_map),
	contacts_(contact_map.contacts_),
	output_matrix_(contact_map.output_matrix_),
	column_names_(contact_map.column_names_),
	row_names_(contact_map.row_names_),
	output_prefix_(contact_map.output_prefix_),
	n_poses_(contact_map.n_poses_),
	distance_cutoff_(contact_map.distance_cutoff_),
	models_per_file_(contact_map.models_per_file_),
	reset_count_(contact_map.reset_count_),
	row_format_(contact_map.row_format_),
	distance_matrix_(contact_map.distance_matrix_)
{
}

/// @brief Destructor
ContactMap::~ContactMap() {
}


moves::MoverOP ContactMap::clone() const {
	return moves::MoverOP( new ContactMap(*this) );
}

moves::MoverOP ContactMap::fresh_instance() const {
	return moves::MoverOP( new ContactMap );
}

std::string ContactMap::get_name() const {
	return ContactMapCreator::mover_name();
}

/// @brief Processes options specified in xml-file and sets up the ContactMap
void ContactMap::parse_my_tag(TagCOP const tag, basic::datacache::DataMap &,
	protocols::filters::Filters_map const &, moves::Movers_map const &,
	Pose const & pose) {

	// 'distance_cutoff' option
	set_distance_cutoff(tag->getOption<core::Real>("distance_cutoff", distance_cutoff_));

	// 'output_prefix_' option
	set_output_prefix(tag->getOption<std::string>("prefix", output_prefix_));

	// 'reset_count' flag
	reset_count_ = tag->getOption<bool>("reset_count", 1);

	// 'models_per_file' option
	models_per_file_ = tag->getOption<core::Size>("models_per_file", models_per_file_);

	// 'row_format_' flag
	row_format_ = tag->getOption<bool>("row_format", 0);

	// 'distance_matrix_' flag
	distance_matrix_ = tag->getOption<bool>("distance_matrix", 0);

	// 'region1', 'region2' and 'ligand' options
	// Initialize 'region1' with complete pose in case no option is specified
	core::Size region1_begin = 1;
	core::Size region1_end = pose.n_residue();
	// Parse 'region1' option if supplied
	if ( tag->hasOption("region1") ) {
		parse_region_string(tag->getOption<std::string> ("region1"),
			region1_begin, region1_end, pose);
	}
	// Parse 'region2' option if supplied and initialize ContactMap between 'region1' and 'region2'
	if ( tag->hasOption("region2") ) {
		core::Size region2_begin, region2_end;
		parse_region_string(tag->getOption<std::string> ("region2"),
			region2_begin, region2_end, pose);
		fill_contacts(region1_begin, region1_end, region2_begin, region2_end, pose);
	} else if ( tag->hasOption("ligand") ) {
		// Parse 'ligand' option if supplied and initialize ContactMap between 'region1' and 'ligand'
		core::Size ligand_begin, ligand_end;
		parse_region_string(tag->getOption<std::string> ("ligand"),
			ligand_begin, ligand_end, pose);
		fill_contacts(region1_begin, region1_end, ligand_begin, pose);
	} else {
		// Initialize ContactMap with 'region1' if neither 'region2' nor 'ligand' were specified
		fill_contacts(region1_begin, region1_end, pose);
	}
}

/// @brief Parses region definition string end sets the boundaries accordingly
void ContactMap::parse_region_string(std::string region_def,
	core::Size & region_begin,
	core::Size & region_end,
	Pose const & pose) {

	// Check if region is defined by chain
	if ( region_def.size() == 1 && isalpha(region_def.at(0))  ) {
		// Retrieve chain ID
		core::Size const chain_id =
			core::pose::get_chain_id_from_chain(region_def, pose);
		// Set region according to chain boundaries
		region_begin = pose.conformation().chain_begin(chain_id);
		region_end = pose.conformation().chain_end(chain_id);
		//TR << "Region defined from " << *region_begin << " to " << *region_end <<"." << std::endl;
		return;
	}

	// Initialize stream and try to read in region_begin
	std::istringstream region_defstr(region_def);
	region_defstr >> region_begin;
	// Set region_end to region_begin if only one position is supplied
	if ( region_defstr.eof() ) {
		region_end = region_begin;
	} else {
		// Skip separator and read region end
		region_defstr.seekg(1, std::ios_base::cur);
		region_defstr >> region_end;
	}
	// Exit if something went wrong
	if ( region_defstr.fail() || ! region_defstr.eof() ) {
		utility_exit_with_message("Unable to parse region '" + region_def +"'!");
	}

	// Make sure region_begin is not greater than region end
	if ( region_begin > region_end ) {
		core::Size temp = region_begin;
		region_begin = region_end;
		region_end = temp;
	}

	// Check if region definition is within bounds of the pose
	if ( region_begin < 1  || region_end > pose.n_residue() ) {
		utility_exit_with_message("Specified region '" + region_def +"' is out of bounds!");
	}
	//TR << "Region defined from " << *region_begin << " to " << *region_end <<"." << std::endl;

}

/// @brief Initializes ContactMap within a single region
void ContactMap::fill_contacts(
	core::Size begin,
	core::Size end,
	Pose const & pose)
{
	ContactPartner p1, p2;
	// Outer loop to assign ContactPartner1
	for ( core::Size i = begin; i<=end; i++ ) {
		std::string atom_name1 =  pose.residue_type(i).name1() == 'G'  ?   "CA" : "CB";
		std::string resname1 = pose.residue(i).name3();
		p1 = ContactPartner(i, resname1, atom_name1);
		// Fill column and row names
		row_names_.push_back(p1.string_rep());
		column_names_.push_back(p1.string_rep());
		// Outer loop to assign ContactPartner1 and create contact
		for ( core::Size j = i+1; j<=end; j++ ) {
			std::string atom_name2 = pose.residue_type(j).name1() == 'G'  ?   "CA" : "CB";
			std::string resname2 = pose.residue(j).name3();
			p2 = ContactPartner(j, resname2, atom_name2);
			// Create contact and add to contacts_
			Contact contact = Contact(p1,p2);
			contacts_.push_back(contact);
		}
	}
	// Add additional contact with same ContactPartner to fill identity contacts in Matrix
	contacts_.push_back(Contact(p1, p1));

	// Set offset and length of region for calculation of array position
	core::Size offset = begin- 1;
	core::Size length = end - offset;

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

/// @brief Initializes ContactMap between a single region and a ligand
void ContactMap::fill_contacts(
	core::Size begin,
	core::Size end,
	core::Size ligand_seqpos,
	Pose const & pose)
{
	core::Size matrix_position = 1;
	std::string ligand_resname = pose.residue(ligand_seqpos).name3();
	// Make sure ligand_seqpos doesn't equal begin
	if ( ligand_seqpos == begin ) ++begin;
	// Loop over residues in specified region
	for ( core::Size i = begin; i <= end; i++ ) {
		// Exclude ligand from residues
		if ( i==ligand_seqpos ) continue;
		std::string atom_name1 = pose.residue_type(i).name1() == 'G' ? "CA" : "CB";
		std::string resname1 = pose.residue(i).name3();
		ContactPartner p1(i, resname1, atom_name1);
		// Add residue string to 'row_names_'
		row_names_.push_back(p1.string_rep());
		// Loop over atoms in ligand
		for ( core::Size j = 1; j <= pose.residue(ligand_seqpos).atoms().size(); j++ ) {
			// Skip hydrogen atoms
			if ( pose.residue(ligand_seqpos).atom_is_hydrogen(j) ) {
				continue;
			}
			// Initialize second partner, create contact and add to 'contacts_'
			std::string atom_name2 = pose.residue(ligand_seqpos).atom_name(j);
			ContactPartner p2(ligand_seqpos, ligand_resname, atom_name2);
			Contact contact = Contact(p1, p2);
			contacts_.push_back(contact);
			// Add ligand atom string to 'column_names_'
			if ( i == begin ) {
				column_names_.push_back(p2.string_rep());
			}
			// Set matrix value to corresponding position in contacts_ vector
			output_matrix_.push_back(matrix_position++);
		}
	}
}

/// @brief Initializes ContactMap between two separate regions
void ContactMap::fill_contacts(
	core::Size start1,
	core::Size end1,
	core::Size start2,
	core::Size end2,
	Pose const & pose)
{
	core::Size matrix_position = 1;
	// Loop over residues in region1
	for ( core::Size i = start1; i <= end1; i++ ) {
		// Assign ContactPartner1
		std::string atom_name1 = pose.residue_type(i).name1() == 'G' ? "CA"
			: "CB";
		std::string resname1 = pose.residue(i).name3();
		ContactPartner p1(i, resname1, atom_name1);
		// Add ContactPartner1 string to 'row_names_'
		row_names_.push_back(p1.string_rep());
		// Loop over residues in region2
		for ( core::Size j = start2; j <= end2; j++ ) {
			// Assign ContactPartner2
			std::string atom_name2 = pose.residue_type(j).name1() == 'G' ? "CA"
				: "CB";
			std::string resname2 = pose.residue(j).name3();
			ContactPartner p2(j, resname2, atom_name2);
			// Create contact and add to 'contacts_'
			Contact contact = Contact(p1, p2);
			contacts_.push_back(contact);
			// Add ContactPartner2 string to 'column_names_'
			if ( i == start1 ) {
				column_names_.push_back(p2.string_rep());
			}
			// Set matrix value to corresponding position in contacts_ vector
			output_matrix_.push_back(matrix_position++);
		}
	}
}

/// @brief Process supplied pose
void ContactMap::apply(Pose & pose) {
	using namespace pose;
	n_poses_++;

	// Iterate over contacts
	for ( utility::vector1<Contact>::iterator it = contacts_.begin(), end = contacts_.end(); it != end; ++it ) {
		ContactPartner * p1(it->partner1());
		ContactPartner * p2(it->partner2());
		// Get coordinates of both contact partners and calculate distance
		numeric::xyzVector<core::Real> v1 = pose.residue(p1->seqpos()).atom(p1->atomname()).xyz();
		numeric::xyzVector<core::Real> v2 = pose.residue(p2->seqpos()).atom(p2->atomname()).xyz();
		core::Real distance = v1.distance(v2);
		// Add distance to contact if it's below the cutoff value
		if ( distance_matrix_ ) {
			it->add_distance(distance);
			continue;
		}
		if ( distance <= distance_cutoff_ ) {
			it->add_distance();
		}
	}

	// Output ContactMap to file if the number of processed poses since last output equals models_per_file_ variable
	if ( models_per_file_ > 0 && n_poses_ % models_per_file_ == 0 ) {
		// If the ContactMap is to be reset, create a unique output name and write the ContactMap to this file
		if ( reset_count_ ) {
			// Get the file name of the output structure.  This will be
			// something_XXXX.pdb if you're outputting to pdbs, something_XXXX if you're outputting to silent files
			std::string structure_output_name(protocols::jd2::JobDistributor::get_instance()->current_output_name());
			// Append a prefix and a suffix to get a final filename for the contact map output
			std::string contact_map_file_name(output_prefix_ + structure_output_name + ".csv");
			// Call output function with the generated filename
			write_to_file(contact_map_file_name);
			// Reset ContactMap
			reset();
		} else {
			// If the ContactMap should just be updated, reconstruct the output name based on the current job name and
			// the specified prefix
			std::string current_job_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());
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
	for ( utility::vector1<Contact>::iterator it = contacts_.begin(), end = contacts_.end(); it != end; ++it ) {
		it->reset_count();
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
		for ( utility::vector1<Contact>::iterator it = contacts_.begin(), end = contacts_.end(); it != end; ++it ) {
			output_stream << it->long_string_rep() << std::endl;
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
