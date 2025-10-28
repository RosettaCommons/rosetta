// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/AntibodyNumberingParser.hh
/// @brief Read Antibody Numbering schemes and transform info from the database
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_AntibodyNumberingParser_hh
#define INCLUDED_protocols_antibody_AntibodyNumberingParser_hh

#include <protocols/antibody/AntibodyNumberingParser.fwd.hh>
#include <protocols/antibody/AntibodyEnumManager.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <utility/VirtualBase.hh> // AUTO IWYU For VirtualBase
#include <map> // AUTO IWYU For map
#include <core/types.hh> // AUTO IWYU For Size
#include <utility/vector1.hh> // AUTO IWYU For vector1

namespace protocols {
namespace antibody {

typedef utility::vector1< utility::vector1< core::Size > > Vector2DSize;

struct AntibodyNumbering {
	AntibodyNumberingSchemeEnum numbering_scheme;
	CDRDefinitionEnum cdr_definition;

	utility::vector1< utility::vector1< PDBLandmarkOP > > cdr_numbering;
	std::map< CDRDefinitionEnum , utility::vector1< utility::vector1< PDBLandmarkOP > > > cdr_definition_transform;

	std::map< AntibodyNumberingSchemeEnum, utility::vector1< PDBLandmarkOP > > numbering_scheme_transform;

};


/// @brief Class responsible for reading database Numbering Scheme definitions and their transforms from the database.
class AntibodyNumberingParser : public utility::VirtualBase {

public:

	/// @brief Default constructor.  We pass the enum manager which is constructed in AbInfo so we only should have one instance of it.
	/// Its not a singleton, so instead of global data, we are careful where and when we pass it around.
	AntibodyNumberingParser(AntibodyEnumManagerCOP enum_manager);

	~AntibodyNumberingParser() override;

	/// @brief Read numbering file and return AntibodyNumbering structure
	AntibodyNumbering
	get_antibody_numbering(AntibodyNumberingSchemeEnum const numbering_scheme, CDRDefinitionEnum const cdr_definition);

private:

	void
	read_cdr_definition_file(std::string const & file_path, AntibodyNumbering & numbering);

	void
	read_numbering_scheme_file(std::string const & file_path, AntibodyNumbering & numbering);


	/// @brief Reads lines defining start/end of each CDR and the relative transforms to the numbering schemes defined by the TRANSFORM line.
	void
	read_cdr_definition_numbering_line(utility::vector1< std::string> const & lineSP, AntibodyNumbering & numbering) const;

	/// @brief Reads line corresponding to TRANSFORM, which lists columns for which the transform to another numbering scheme is defined.
	void
	read_cdr_definition_transform_line(utility::vector1< std::string> const & lineSP, AntibodyNumbering & numbering);


	void
	read_scheme_numbering_line(utility::vector1< std::string > const & lineSP, AntibodyNumbering & numbering) const;

	void
	read_scheme_defines_line(utility::vector1< std::string > const & lineSP);


	/// @brief Check to make sure the path to the numbering scheme file is good.
	void
	check_path(std::string const & numbering_file_path) const;

	AntibodyNumberingSchemeEnum
	get_numbering_scheme_used_for_cdr_definition(CDRDefinitionEnum) const;

	/// @brief Gets equivalent landmark from that defined in landmark_to_match
	/// @details.  Ex.  landmark_to_match defines a CDR start point in Chothia_scheme.  Our numbering must be in Kabat.  What is the PDBLandmark for the same residue?
	PDBLandmarkOP
	get_equivalent_landmark(
		AntibodyNumbering & numbering,
		const AntibodyNumberingSchemeEnum scheme,
		PDBLandmark & landmark_to_match) const;

	void
	debug_print(AntibodyNumbering & numbering);

private:
	AntibodyEnumManagerCOP enum_manager_;
	std::string numbering_database_directory_; //Directory in database where numbering is stored
	std::string scheme_file_;
	std::string cdr_definition_file_;

	utility::vector1< CDRDefinitionEnum  > cdr_definitions_defined_; //A cdr definition needs a scheme to define it in.
	utility::vector1< AntibodyNumberingSchemeEnum >cdr_definitions_defined_using_;

	utility::vector1< AntibodyNumberingSchemeEnum > schemes_defined_;

};


/// @brief Class that was once a struct; Used for matching pdb information between numbering schemes and cdr definitions.
class PDBLandmark : public utility::VirtualBase {

public:
	PDBLandmark(std::string const & chain, core::Size resnum,  char insertion_code);

	/// @brief Alternative constructor to hold numbering scheme type as well.
	PDBLandmark(std::string const & chain, core::Size resnum, char insertion_code, AntibodyNumberingSchemeEnum scheme);

	core::Size
	resnum() const { return resnum_;}

	std::string const &
	chain() const { return chain_; }

	char
	insertion_code() const  {return insertion_code_;}

	AntibodyNumberingSchemeEnum
	numbering_scheme() const { return numbering_scheme_;}

	// Equivalent operators do not check scheme, since scheme is optional and we want to compare data with different schemes..
	bool
	operator==(const PDBLandmark & compare) const;

	bool
	operator!=(const PDBLandmark & compare) const;

	std::string
	get_string() const;

private:
	AntibodyNumberingSchemeEnum numbering_scheme_;

	core::Size resnum_; //PDB residue number
	std::string chain_;
	char insertion_code_;

};

} //antibody
} //protocols


#endif //INCLUDED_protocols_antibody_AntibodyNumberingParser.hh


