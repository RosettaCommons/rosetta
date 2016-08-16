// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/PseudocontactShiftInput.hh
///
/// @brief Read input .npc input file
///
/// @details The following classes are responsable to read / parse the PCS input file (.npc format)
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


#ifndef INCLUDED_protocols_scoring_methods_pcs_PseudocontactShiftInput_hh
#define INCLUDED_protocols_scoring_methods_pcs_PseudocontactShiftInput_hh

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers
#include <string>
#include <map>

#include <utility/SingletonBase.hh>
#include <utility/vector1_bool.hh>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs {

///////////////////////////////////////////////////////////////////////////
/// @brief PCS_line_data class: hold a line of the input file information (.npc format)
/// One PCS_line_data per line in the input file
class PCS_line_data {

public:
	PCS_line_data();

	~PCS_line_data();

	PCS_line_data(PCS_line_data const & other);

	PCS_line_data &
	operator=( PCS_line_data const & other );

	PCS_line_data(
		core::Size const residue_num,
		std::string const atom_name,
		core::Real const PCS_experimental,
		core::Real const PCS_tolerance
	);

	core::Size
	residue_num() const;

	std::string
	atom_name() const;

	core::Real
	PCS_experimental() const;

	core::Real
	PCS_tolerance() const;

	friend std::ostream &
	operator<<(std::ostream& out, const PCS_line_data &PCS_l_d);

private:
	core::Size const residue_num_;
	std::string const atom_name_;
	core::Real const PCS_experimental_;
	core::Real const PCS_tolerance_;
};


//////////////////////////////////////////////////////////
/// @brief PCS_file_data contain all the information of a .npc file
/// one per lanthanide.
class PCS_file_data {
private:
	std::string const filename_;
	utility::vector1<PCS_line_data> PCS_data_line_all_;
	core::Real const weight_;

	void
	read_PCS_file();

public:
private:
	PCS_file_data();
public:
	~PCS_file_data();

	PCS_file_data(PCS_file_data const & other);

	PCS_file_data &
	operator=( PCS_file_data const & other );

	PCS_file_data(std::string const & filename, core::Real const my_weight );

	std::string
	get_filename() const;

	core::Real
	get_weight() const;

	utility::vector1<PCS_line_data> &
	get_PCS_data_line_all_reference();

	friend std::ostream &
	operator<<(std::ostream& out, const PCS_file_data &PCS_f_d);
};


//////////////////////////////////////////////////////////////
/// @brief PCS_data_input contain all the input information for the PCS.
/// This includes all the information from the .npc files
class PCS_data_input {
private:
	std::map< std::string, PCS_file_data > PCS_filename_and_data_;

public:
	PCS_data_input();

	~PCS_data_input();

	PCS_data_input(PCS_data_input const & other);

	PCS_data_input &
	operator=( PCS_data_input const & other );

	PCS_data_input(utility::vector1<std::string> const & filenames,  utility::vector1<core::Real> const & weight);

	std::map< std::string, PCS_file_data > &
	get_PCS_data_input_reference();

	friend std::ostream &
	operator<<(std::ostream& out, const PCS_data_input &PCS_d_i);

};

class PCS_data_input_manager : public utility::SingletonBase< PCS_data_input_manager >
{
public:
	friend class utility::SingletonBase< PCS_data_input_manager >;

private:
	PCS_data_input_manager();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static PCS_data_input_manager * create_singleton_instance();

private:

	std::map<std::string, PCS_data_input> file_2_data_map_;

public:

	PCS_data_input
	get_input_data(utility::vector1<std::string> const & filenames, utility::vector1<core::Real> const & vec_weight);
};


}//namespace pcs
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
