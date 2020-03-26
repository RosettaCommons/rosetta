// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Steven Lewis
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_task_ResfileReader_hh
#define INCLUDED_core_pack_task_ResfileReader_hh

// Unit Headers
#include <core/pack/task/ResfileReader.fwd.hh>

// Package Headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pose/Pose.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/chemical/AA.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/VirtualBase.hh>
#include <utility/excn/Exceptions.hh>

// STL Headers
#include <iosfwd>
#include <list>
#include <map>
#include <string>


namespace core {
namespace pack {
namespace task {


class ResfileReaderException: public utility::excn::Exception
{
public:

	ResfileReaderException(char const *file, int line, std::string message ) :
		utility::excn::Exception(file, line, message )
	{}

	std::string get_message(){
		return msg();
	}

};

class ResfileContents : public utility::VirtualBase
{
public:

	// Moves fname, so by-value
	ResfileContents( pose::Pose const & pose, std::string const & fname, std::istream & istream );
	~ResfileContents() override;

	std::list< ResfileCommandCOP > const &
	default_commands() const;

	bool specialized_commands_exist_for_residue( Size resid ) const;

	std::list< ResfileCommandCOP > const &
	commands_for_residue( Size resid ) const;

	std::string const & fname_initialized_from() const;

private: // helper types and functions for parsing

	enum residue_identifier_type{
		SINGLE_RESID=1,
		RANGE_RESID,
		CHAIN_RESID
	};

	/// @brief Given a vector of strings corresponding to the words in a whitespace-separated line, look for
	/// deprecated commands and issue a suitable warning.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void check_for_deprecated_commands( utility::vector1< std::string > const & tokens ) const;

	void
	parse_header_line(
		utility::vector1< std::string > const & tokens,
		std::map< std::string, ResfileCommandOP > const & command_map,
		Size const lineno,
		bool & have_read_start_token);

	void
	parse_body_line(
		pose::Pose const & pose,
		utility::vector1< std::string > const & tokens,
		std::map< std::string, ResfileCommandOP > const & command_map,
		Size const lineno,
		utility::vector1< std::list< ResfileCommandOP > > & residue_range_commands,
		utility::vector1< std::list< ResfileCommandOP > > & residue_chain_commands);

	void
	parse_resid(
		Size & which_token,
		utility::vector1< std::string > const & tokens,
		Size const lineno,
		int & PDBnum,
		int & PDBnum_end, // only defined for RANGE id_type
		char & icode,
		char & icode_end, // only defined for RANGE id_type
		char & chain,
		residue_identifier_type & id_type) const;

	void
	parse_PDBnum_icode(
		std::string const & token,
		Size const lineno,
		int & PDBnum,
		char & icode) const;

	Size
	locate_resid(
		core::pose::Pose const & pose,
		char const chain,
		Size const PDBnum,
		char const icode,
		Size const lineno) const;

	ResfileCommandOP
	locate_command(
		Size const which_token,
		utility::vector1< std::string > const & tokens,
		std::map< std::string, ResfileCommandOP > const & command_map,
		Size const lineno) const;

private:
	std::string fname_initialized_from_;
	std::list< ResfileCommandCOP > default_commands_;
	utility::vector1< std::list< ResfileCommandCOP > > commands_;
};


/// @brief abstract/interface class for Resfile reader command objects
class ResfileCommand : public utility::VirtualBase
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from VirtualBase
	~ResfileCommand() override;
	virtual ResfileCommandOP clone() const = 0;

	// @brief Read the contents of the Resfile and store the state
	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) = 0;

	/// @brief Modify the packer task with the command that was read in
	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const = 0;

	//JAB - although this calls the name() function, the name function implemented in each command is static,
	// and the design of the ResfileReader requires it to be.
	// Since it can't be a virtual function in the base class,
	// we can't call it when we don't know the complete type of the class.
	//Hence this pure virtual.
	virtual std::string
	get_name() = 0;
};

///////////////////////////////////////////////////////////////////////
//In this section, list mode-style options (NATAA, etc)

/// @brief NATRO disables packing and designing at a position, the residue
///will be totally unchanged
class NATRO : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< NATRO >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "NATRO";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief NATAA allows repacking but no sequence changes (all rotamers are of the original residue)
class NATAA : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< NATAA >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "NATAA";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief ALLAA is deprecated; allows repacking and designing to any canonical residue (default state of PackerTask)
class ALLAA : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< ALLAA >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "ALLAA";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief ALLAAxc allows repacking and designing to any canonical noncysteine residue
class ALLAAxc : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< ALLAAxc >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "ALLAAXC";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief allows repacking and designing to any canonical residue (default state of PackerTask)
class ALLAAwc : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< ALLAAwc >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "ALLAAWC";}

	std::string
	get_name() override {
		return name();
	}


};

/// @brief PIKAA allows residues specifed in a following string.
/// @details In actuality, it is PROHIBITING any residue that is NOT in the
/// following string.  The string should be formatted as an all-caps string of
/// one-letter codes.  Noncanonical amino acids can be included using X[<full base name>].
/// For example, to allow tyrosine, threonine, tryptophan, and 2-aminoisobutyric acid,
/// you would use "PIKAA YTWX[AIB]".
/// @author Original author unknown.
/// @author Noncanonical pruning support added by Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class PIKAA : public ResfileCommand
{
public:
	/// @brief Default constructor.
	PIKAA() = default;

	/// @brief Default copy constructor.
	PIKAA( PIKAA const & /*src*/ ) = default;

	ResfileCommandOP clone() const override { return utility::pointer::make_shared< PIKAA >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	/// @brief Add a base name to the list of base names to keep.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	add_base_name_to_keep(
		std::string const & basename
	);

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "PIKAA";}

	std::string
	get_name() override {
		return name();
	}

private:
	bool initialized_ = false;
	utility::vector1< std::string > basenames_to_keep_;
};

/// @brief PIKNA allows nucleic acid residues specifed in a following string
class PIKNA : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< PIKNA >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "PIKNA";}

	std::string
	get_name() override {
		return name();
	}

private:
	utility::vector1< chemical::AA > keep_nas_;
};

/// @brief NOTAA disallows residues specified in a following string, and allows packing
class NOTAA : public ResfileCommand
{
public:
	/// @brief Default constructor.
	NOTAA() = default;

	ResfileCommandOP clone() const override { return utility::pointer::make_shared< NOTAA >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	/// @brief Add a base name to the list of base names to exclude.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void add_base_name_to_exclude( std::string const & basename );

	static std::string name() {return "NOTAA";}

	std::string
	get_name() override {
		return name();
	}

private:
	utility::vector1< std::string > basenames_to_exclude_;
	bool initialized_ = false;
};

/// @brief Allows designing on ANY residue type Property. (Only currently works with Cannonical AAs)
/// JAB
class PROPERTY : public ResfileCommand
{
public:
	/// @brief Constructor.
	PROPERTY() = default;

	ResfileCommandOP clone() const override { return utility::pointer::make_shared< PROPERTY >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "PROPERTY";}

	std::string
	get_name() override {
		return name();
	}

	core::chemical::ResidueProperty property_ = core::chemical::NO_PROPERTY;
};


/// @brief POLAR allows polar residues and packing
class POLAR : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< POLAR >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "POLAR";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief APOLAR allows nonpolar residues and packing
class APOLAR : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< APOLAR >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "APOLAR";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief APOLA is a deprecated version of APOLAR allows nonpolar residues and packing
class APOLA : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< APOLA >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "APOLA";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief CHARGED allows charged residues and packing
/// JAB
class CHARGED : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< CHARGED >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "CHARGED";}

	std::string
	get_name() override {
		return name();
	}

};


///@brief AROMATIC allows designing only aromatic residues.
/// JAB
class AROMATIC : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< AROMATIC >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "AROMATIC";}

	std::string
	get_name() override {
		return name();
	}

};

////////end of mode-style options////////////////////////

////////in this section list other options///////////////

/// @brief EX handles the various extrachi options
class EX : public ResfileCommand
{
public:
	EX():
		aro_specified_( false ), // Does this make sense?
		which_chi_( 0 ) // Does this make sense?
	{}

	ResfileCommandOP clone() const override { return utility::pointer::make_shared< EX >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;



	static std::string name() {return "EX";}

	std::string
	get_name() override {
		return name();
	}

private:
	bool           aro_specified_;
	Size           which_chi_;
	ExtraRotSample chi_sample_level_;
};

/// @brief EX_CUTOFF allows setting of the extrachi_cutoff (for determining burial for extra rotamers)
class EX_CUTOFF : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< EX_CUTOFF >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "EX_CUTOFF";}

	std::string
	get_name() override {
		return name();
	}

private:
	Size ex_cutoff_;
};

/// @brief USE_INPUT_SC turns on inclusion of the current rotamer for the packer
class USE_INPUT_SC : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< USE_INPUT_SC >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "USE_INPUT_SC";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief AUTO suggests that a packer can/should reconsider the design setting at a/each residue
/// @details This is a protocol-level flag to be used in non-vanilla packers. For example, one may want an ALLAA tag to be effective only if the residue is in an automatically-determined region of interest, without knowing which residues qualify a-priori
/// @author ashworth
class AUTO : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< AUTO >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "AUTO";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief SCAN suggests to some packing routines that if there are multiple type choices for this residue, then each of them should be considered explicitly in one way or another
/// @details This is a protocol-level flag to be used in non-vanilla packers
/// @author ashworth
class SCAN : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< SCAN >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "SCAN";}

	std::string
	get_name() override {
		return name();
	}
};

/// @brief TARGET flags the position as "targeted", and can optionally specify a "targeted" type
/// @details This is a protocol-level flag to be used in non-vanilla packers--positions flagged as "targeted" may be treated in a special fashion
///
///  The optional specification of a target type is be useful for multistate considerations:
///  multistate protocols need 1) rotamers and energies for all possible states, and 2) a target state
///  The target type must be a member of PackerTask's allowed_types_
///
/// @author ashworth
class TARGET : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< TARGET >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "TARGET";}

	std::string
	get_name() override {
		return name();
	}

private:
	std::string argstring_;
};

/// @brief NO_ADDUCTS will disable adducts, assuming they exist
/// @details This command exists because if adducts exist, then they are enabled by default for all residues.
/// @author ashworth
class NO_ADDUCTS : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< NO_ADDUCTS >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "NO_ADDUCTS";}

	std::string
	get_name() override {
		return name();
	}

};

/// @brief FIX_HIS_TAUTOMER: when a histidine is present when the PackerTask is initialized, this flag will fix its tautomer (whether its hydrogen is on ND1 or NE2.  Does nothing if not histidine at initialization (meaning if it mutates to histidine later this flag will have no effect).
class FIX_HIS_TAUTOMER : public ResfileCommand
{
public:
	ResfileCommandOP clone() const override { return utility::pointer::make_shared< FIX_HIS_TAUTOMER >(); }

	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) override;

	void residue_action(
		PackerTask &,
		Size resid
	) const override;

	static std::string name() {return "FIX_HIS_TAUTOMER";}

	std::string
	get_name() override {
		return name();
	}

};

///////////end of other options//////////////////////////
///////////utility functions for resfile reader//////////

/// @brief utility function to increment next token to be parsed
std::string
get_token(
	const Size which_token,
	const utility::vector1<std::string> & tokens,
	const bool make_upper_case = true);

/// @brief Split the line into vector of strings.
///   Each separated string is a 'token'
utility::vector1< std::string >
tokenize_line( std::istream & inputstream );

/// @brief Split the line into vector of strings.
///   Each separated string is a 'token'
utility::vector1< std::string >
tokenize_line( std::string const & line);

/// @brief utility for resfile reader, commands MUST be entered into this hard-coded map
std::map< std::string, ResfileCommandOP >
create_command_map();

/// @brief utility function for resfile reader (checks for a leading # signaling a comment)
bool
comment_begin( utility::vector1< std::string > const & tokens, Size which_token );

/// @brief changes the state of the given PackerTask according to the commands in the resfile at read in from the -pack:resfile option system.
void
parse_resfile(
	pose::Pose const & pose,
	PackerTask & the_task);

/// @brief changes the state of the given PackerTask according to the commands in the resfile at filename
void
parse_resfile(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string const & filename );

/// @brief changes the state of the given PackerTask according to the commands in the resfile.
/// @details This version calls the overloaded parse_resfile_string and just passes it a ResidueSubset
/// that's set to "true" for every residue of the pose.
void
parse_resfile_string(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string const & resfile_fname,
	std::string const & resfile_string
);


/// @brief changes the state of the given PackerTask according to the commands in the resfile.
/// @details This version accepts a ResidueSubset (a utility::vector1<bool>) that serves as a
/// mask.  Residues set to "false" don't have their packer behaviour altered by this TaskOperation.
void
parse_resfile_string(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string const & resfile_fname,
	std::string const & resfile_string,
	core::select::residue_selector::ResidueSubset const &mask

);

void
onError( std::string message );


}//namespace task
}//namespace pack
}//namespace core

#endif
