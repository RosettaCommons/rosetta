// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////////////////////////////////////////////////////
/// @file   ContactMap.hh
///
/// @brief  Mover to create a ContactMap
///
/// @details :
///  This mover takes a pose and creates a matrix with true(1) or false(0) values for
///  atom pairs depending on whether both atoms are within a specified cutoff distance.
///  For protein residues (specified via region1/2 tag) CB or CA (Gly) atoms and for ligands
///  (specified via ligand tag) all heavy atoms are used for distance calculations.
///  Depending on the options specified values may be added up to include multiple poses in
///  one ContactMap output file
///
/// @author Joerg Schaarschmidt
//////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_protocols_contact_map_ContactMap_hh
#define INCLUDED_protocols_contact_map_ContactMap_hh

// Unit headers
#include <protocols/contact_map/ContactMap.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace contact_map {

class ContactMap : public moves::Mover {

public:
	/// @brief Default constructor
	ContactMap();
	/// @brief Copy constructor
	ContactMap(ContactMap const &);

	/// @brief Destructor
	~ContactMap();

	virtual moves::MoverOP clone() const;
	virtual moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & );

	virtual void apply( Pose & pose );

	virtual std::string get_name() const;
	virtual void test_move( Pose & pose ){
		apply(pose);
	}

	/// @brief Setter functions for private class variables
	void set_output_prefix( std::string prefix){ output_prefix_ = prefix; }
	void set_distance_cutoff( core::Real cutoff ){ distance_cutoff_ = cutoff; }

	/// @brief :
	///  Parses region string end sets the boundaries accordingly
	/// @details:
	///  Possible formats for the region definition are "<int seqpos>" (e.g '12'),
	///  "<int begin>-<int end>" (e.g. '1-10') or "<char chainID>"(e.g. 'A')
	void parse_region_string(std::string,
		core::Size & ,
		core::Size &,
		Pose const &);

	/// @brief Resets the movers n_poses_ variable  and the counts of all contacts to 0
	void reset();

	/// @brief Output function that writes the ContactMap to the specified file
	void write_to_file(std::string filename="");
	void write_to_stream(std::ostream & output_stream);

	/// @brief Initializes ContactMap within a single region
	void fill_contacts(
		core::Size region_begin,
		core::Size region_end,
		Pose const & pose);

	/// @brief Initializes ContactMap between a single region and a ligand
	void fill_contacts(
		core::Size region_begin,
		core::Size region_end,
		core::Size seqpos_ligand,
		Pose const & pose);

	/// @brief Initializes ContactMap between two separate regions
	void fill_contacts(
		core::Size region1_begin,
		core::Size region1_end,
		core::Size region2_begin,
		core::Size region2_end,
		Pose const & pose);

	/// @brief Get the contact for a particular row and column
	/// The identity contact (the contact of a residue with itself)
	/// is shared between all identity contacts. (The last identity contact in region.)
	Contact const & get_contact( core::Size row, core::Size col);

	// utility::vector1<Contact> const & contacts() { return contacts_; }
	// utility::vector1<core::Size> const & output_matrix() { return output_matrix_; }
	// utility::vector1<std::string> const & column_names() { return column_names_; }
	// utility::vector1<std::string> const & row_names() { return row_names_; }

private:
	/// @brief Container for all AtomPairs included in the ContactMap
	utility::vector1<Contact> contacts_;
	/// @brief Vector that maps the matrix position to the index of the contact in contacts_
	/// @detail Vector of size columns*rows with field(row, column) stored at (row-1)*[#columns] + column
	utility::vector1<core::Size> output_matrix_;
	utility::vector1<std::string> column_names_;
	utility::vector1<std::string> row_names_;
	std::string output_prefix_;
	core::Size n_poses_;
	core::Real distance_cutoff_;
	/// @brief int value after how many models an output file will be generated
	core::Size models_per_file_;
	/// @brief bool value indicating whether the reset function will be called after file output
	bool reset_count_;
	bool row_format_;
	bool distance_matrix_;


}; //class ContactMap


/// @brief Simple class that holds the information on an atom involved in a contact
class ContactPartner{

public:
	/// @brief Default constructor
	ContactPartner(core::Size seqpos = 0, std::string resname = "", std::string aname= "") :
		seqpos_(seqpos),
		resname_(resname),
		atomname_(aname){}

	/// @brief Returns string representation of the ContactPartner
	std::string string_rep() const;

	/// @brief Accessor functions for private class variables
	core::Size seqpos() const {return seqpos_;}
	std::string resname() const {return resname_;}
	std::string atomname() const {return atomname_;}

private:
	core::Size seqpos_;
	std::string resname_;
	std::string atomname_;
}; // class ContactPartner


/// @brief Simple class representing a contact between two atoms
class Contact{

public:
	/// @brief Default constructor
	Contact(ContactPartner p1 = ContactPartner(), ContactPartner p2 = ContactPartner()) :
		partner1_(p1),
		partner2_(p2),
		count_(0),
		distance_(0.0){}

	/// @brief Adds distance to the contact
	void add_distance();

	/// @brief Adds distance to the contact
	void add_distance( core::Real distance );

	/// @brief Resets count to 0
	void reset_count();

	/// @brief Returns string representation of the Contact
	std::string string_rep() const;

	/// @brief Returns string representation of the Contact as percentage value
	std::string string_rep(core::Size n_poses) const;

	/// @brief Returns string representation of the Contact including partner names
	std::string long_string_rep(core::Size n_poses=0) const;

	/// @brief Accessor functions for private class variables
	ContactPartner * partner1() {return &partner1_;}
	ContactPartner * partner2() {return &partner2_;}

private:
	ContactPartner partner1_;
	ContactPartner partner2_;
	core::Size count_;
	core::Real distance_;


}; // class Contact


} // moves
} // protocols
#endif //INCLUDED_protocols_moves_ContactMap_HH
