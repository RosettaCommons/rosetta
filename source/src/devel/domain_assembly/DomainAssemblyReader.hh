// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/olungu/
/// @brief  header of classes for domain assembly options


#ifndef INCLUDED_devel_domain_assembly_DomainAssemblyReader_hh
#define INCLUDED_devel_domain_assembly_DomainAssemblyReader_hh

// Unit Headers
#include <devel/domain_assembly/domain_assembly_setup.hh>

// Package Headers

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// STL Headers
#include <iosfwd>
#include <map>

#include <utility/vector1.hh>

namespace devel {
namespace domain_assembly {

class DomainAssemblyCommand;

typedef utility::pointer::shared_ptr< DomainAssemblyCommand > DomainAssemblyCommandOP;
typedef utility::pointer::shared_ptr< DomainAssemblyCommand const > DomainAssemblyCommandCOP;

//pares_domain_file will need to take vector of domain objects
//reaader function needs to pass one member of vector function to reader at a time
 ///////////////////////////////////////////////////
class DomainAssemblyCommand : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~DomainAssemblyCommand();
	//@brief domain_action is the activity code for the option
	virtual
	void domain_action(
		utility::vector1< std::string > const & tokens,
	  Size & which_token,
		DomainInfo & domain
	) const = 0;

};
//////////////////////////////////////////////////////////
///PDB class will input path to pdb(s) that will be read
///into domain class
class PDB : public DomainAssemblyCommand
{
public:

	virtual
	void domain_action(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	DomainInfo & domain
	) const;

	static std::string name() {return "PDB";}
};
	///////////////////////////////////////////////
//NTermLinker class will input one letter amino acid
///linker sequence after specified N terminus
class NTermLinker : public DomainAssemblyCommand
{
public:

	virtual
	void domain_action(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		DomainInfo & domain
		) const;

	static std::string name() {return "NTermLinker";}
};
	///////////////////////////////////////////////
//CTermLinker class will input one letter amino acid
///linker sequence after specified C terminus
class CTermLinker : public DomainAssemblyCommand
{
public:

	virtual
	void domain_action(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	DomainInfo & domain
) const;

	static std::string name() {return "CTermLinker";}
};
	///////////////////////////////////////////////
//NTrim class will remove specified number of amino
//acids from specified N terminus
class NTrim : public DomainAssemblyCommand
{
public:

	virtual
	void domain_action(
	utility::vector1< std::string > const & tokens,
	Size & which_token,
	DomainInfo & domain
	) const;

	static std::string name() {return "NTrim";}
};
	///////////////////////////////////////////////
//CTrim class will remove specified number of amino
//acids from specified C terminus
class CTrim : public DomainAssemblyCommand
{
public:

	virtual
	void domain_action(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		DomainInfo & domain
	) const;

	static std::string name() {return "CTrim";}
};

void
parse_da_option_file( utility::vector1< DomainInfo > & domains, std::string filename );

///////////utility functions for domain_assembly reader//////////
//add comment_line to resfile reader.hh
/// @brief utility for domain_assembly reader, commands MUST be entered into this hard-coded map
std::map< std::string, DomainAssemblyCommandOP >
create_command_map();

}
}


#endif
