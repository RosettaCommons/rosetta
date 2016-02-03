// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/StructFileRep.cc
/// @brief  The representation of a structure file
/// @author Andy Watkins
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- conversion to class from struct
/// @author Labonte <JWLabonte@jhu.edu>

// Unit headers
#include <core/io/pdb/Field.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/StructFileRep.hh>

// Package headers
#include <core/io/StructFileRepOptions.hh>
//#include <core/io/pdb/file_data_fixup.hh>
#include <core/io/StructFileReaderOptions.hh>
#include <core/io/NomenclatureManager.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedAtomID_Map.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/cryst/util.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <boost/lexical_cast.hpp>

// C++ headers
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <utility>


namespace core {
namespace io {

using core::Size;
using core::SSize;

using core::chemical::chr_chains;

using basic::T;
using basic::Error;
using basic::Warning;

using ObjexxFCL::strip_whitespace;
using ObjexxFCL::stripped_whitespace;
using ObjexxFCL::rstripped_whitespace;
using namespace ObjexxFCL::format;

using std::string;
using std::iostream;

// Tracer instance for this file
static THREAD_LOCAL basic::Tracer TR( "core.io.StructFileRep" );

ResidueInformation::ResidueInformation() :
	// resid( "" ),
	resName_( "" ),
	chainID_( ' ' ),
	resSeq_( 0 ),
	iCode_( ' ' ),
	terCount_( 0 ),
	atoms_(),
	xyz_(),
	temps_(),
	segmentID_( "    " )
{}

ResidueInformation::ResidueInformation(
	AtomInformation const & ai
) :
	// resid( "" ),
	resName_( ai.resName ),
	chainID_( ai.chainID ),
	resSeq_( ai.resSeq ),
	iCode_( ai.iCode ),
	terCount_( ai.terCount ),
	atoms_(),
	xyz_(),
	temps_(),
	segmentID_( ai.segmentID )
{}

bool
ResidueInformation::operator==(
	ResidueInformation const & that) const {
	return
		resName_  == that.resName_ &&
		chainID_  == that.chainID_ &&
		resSeq_   == that.resSeq_  &&
		iCode_    == that.iCode_   &&
		terCount_ == that.terCount_;
}

bool
ResidueInformation::operator!=(
	ResidueInformation const & that) const {
	return !(*this == that);
}

String const & ResidueInformation::resName() const { return resName_; }
char ResidueInformation::chainID()  const { return chainID_; }
int  ResidueInformation::resSeq()   const { return resSeq_; }
char ResidueInformation::iCode()    const { return iCode_; }
int  ResidueInformation::terCount() const { return terCount_; }
std::string const & ResidueInformation::segmentID() const { return segmentID_; }

void ResidueInformation::resName(  String const & setting )
{
	resName_ = setting;
	for ( Size ii = 1 ; ii <= atoms_.size(); ++ii ) atoms_[ ii ].resName = setting;
}
void ResidueInformation::chainID(  char setting ) { chainID_ = setting; }
void ResidueInformation::resSeq(   int setting ) { resSeq_ = setting; }
void ResidueInformation::iCode(    char setting ) { iCode_ = setting; }
void ResidueInformation::terCount( int setting ) { terCount_ = setting; }
void ResidueInformation::set_xyz( std::string const &atomname, Vector const &vect ) { xyz_[atomname] = vect; }
void ResidueInformation::set_temp( std::string const &atomname, core::Real const &val ) { temps_[atomname] = static_cast<double>(val); }
void ResidueInformation::segmentID( std::string const & setting ) { segmentID_ = setting; }

utility::vector1< AtomInformation > const & ResidueInformation::atoms() const { return atoms_; }
void ResidueInformation::append_atom( AtomInformation const & new_atom ) {
	atoms_.push_back( new_atom );
	if ( ! xyz_.count( new_atom.name ) ) {
		xyz_[ new_atom.name ] = Vector( new_atom.x, new_atom.y, new_atom.z );
		temps_[ new_atom.name ] = new_atom.temperature;
	}
}

void ResidueInformation::rename_atom( std::string const & orig_name, std::string const & new_name )
{
	if ( xyz_.find( orig_name ) == xyz_.end() ) {
		utility_exit_with_message( "Failed to find atom \"" + orig_name + "\" in ResidueInformation when trying to rename it to \"" + new_name );
	}
	for ( Size ii = 1; ii <= atoms_.size(); ++ii ) {
		if ( atoms_[ ii ].name == orig_name ) {
			atoms_[ ii ].name = new_name;
			// Do not insert this atom into the xyz_ and temps_ maps if the
			// the atom would displace an older atom
			if ( ! xyz_.count( new_name ) ) {
				xyz_[ new_name ] = xyz_[ orig_name ];
				temps_[ new_name ] = atoms_[ ii ].temperature;
			}
			xyz_.erase( orig_name );
			temps_.erase( orig_name );
			return;
		}
	}
}

/// @details Append all the atoms from the source residue, updating their
/// resName so that they all have the same resName, and so that the xyz_
/// and temps_ maps are updated.
void ResidueInformation::append_atoms( ResidueInformation const & source )
{
	for ( core::Size ii = 1; ii <= source.atoms_.size(); ++ii ) {
		AtomInformation const & iiat = source.atoms_[ ii ];
		if ( ! xyz_.count( iiat.name ) ) {
			AtomInformation new_at( iiat );
			new_at.resName = resName_;
			atoms_.push_back( new_at );
			xyz_[ iiat.name ] = Vector( iiat.x, iiat.y, iiat.z );
			temps_[ iiat.name ] = iiat.temperature;
		}
	}
}

std::map< std::string, Vector > const & ResidueInformation::xyz() const {
	return xyz_;
}

//< map of names to B-factors;  redundant but used a lot in reader
std::map< std::string, double > const & ResidueInformation::temps() const {
	return temps_;
}

String ResidueInformation::resid() const {
	String buf;
	buf.resize(7);
	// This is horribly hacky. Is this necessary?
	sprintf(&buf[0], "%4d%c%c", resSeq_, iCode_, chainID_ );
	buf.resize(6);
	return buf;
}

///////////////////////////////////////////////////////////////////////////////
StructFileRep::StructFileRep() :
	utility::pointer::ReferenceCount(),
	filename_(""),
	modeltag_(""),
	header_(HeaderInformationOP( new HeaderInformation() ) ),
	remarks_(RemarksOP( new Remarks )),
	heterogen_names_(),
	residue_type_base_names_(),
	ssbond_map_(),
	link_map_(),
	crystinfo_(),
	chains_(),
	foldtree_string_(""),
	pdb_comments_(),
	additional_string_output_("")
{}

StructFileRep::~StructFileRep()
{}

/// @brief Clone operator.
/// @details Uses the compiler-default copy constructor.
StructFileRepOP
StructFileRep::clone() const {
	return StructFileRepOP( new StructFileRep( *this ) );
}

/// @details Debug/Info function.
/// Output StructFileRep object to TR like stream in human readable format.
std::ostream&
operator <<( std::ostream &os, StructFileRep const & sfr )
{
	os << "<StructFileRep>{";
	for ( Size i=0; i<sfr.chains().size(); i++ ) {
		os << "Chain<" << i << ">";
		for ( Size j=0; j<sfr.chains()[i].size(); j++ ) {
			os << "[" << j << ":" << sfr.chains()[i][j] << "]" << "\n";
		}
	}
	os << "}";
	return os;
}

/// @brief The name of the input file, if any.
std::string const & StructFileRep::filename() const { return filename_; }
std::string & StructFileRep::filename() { return filename_; }

/// @brief ????
std::string const & StructFileRep::modeltag() const { return modeltag_; }
std::string & StructFileRep::modeltag() { return modeltag_; }

/// @brief PDB Title Section
/// @details "header" is a misnomer, as it actually stores HEADER, TITLE, EXPDTA, KEYWDS, and COMPND records.
HeaderInformationCOP  StructFileRep::header() const { return header_; }
HeaderInformationOP & StructFileRep::header() { return header_; }

// Data for OBSLTE, SPLT, CAVEAT, NUMMDL, MDLTYP, AUTHOR, REVDAT, SPRSDE, and/or JRNL records should be declared
// here if ever implemented.
/// @brief PDB remarks.
///
RemarksCOP  StructFileRep::remarks() const { return remarks_; }
RemarksOP & StructFileRep::remarks() { return remarks_; }


// PDB Primary Structure Section //////////////////////////////////////////
// Data for DBREF, DBREF1, DBREF2, and/or SEQADV records should be declared here if ever implemented.
std::map< char, utility::vector1< std::string > > const & StructFileRep::chain_sequences() const { return chain_sequences_; }  // for storing SEQRES records
std::map< char, utility::vector1< std::string > > & StructFileRep::chain_sequences() { return chain_sequences_; }  // for storing SEQRES records

// map for storing MODRES records:
// key is 6-character resID of the modified residue
std::map< std::string, ModifiedResidueInformation > const & StructFileRep::modres_map() const { return modres_map_; }
std::map< std::string, ModifiedResidueInformation > & StructFileRep::modres_map() { return modres_map_; }

/// @brief list for storing HETNAM records:
/// @details first value of the pair is hetID\n
/// second value of the pair is the chemical name field
utility::vector1< std::pair< std::string, std::string > > const & StructFileRep::heterogen_names() const { return heterogen_names_; }
utility::vector1< std::pair< std::string, std::string > > & StructFileRep::heterogen_names() { return heterogen_names_; }

// map for storing HETSYN records:
//  key is hetID
//  value is the chemical synonym list field
std::map< std::string, utility::vector1< std::string > > const & StructFileRep::heterogen_synonyms() const { return heterogen_synonyms_; }
std::map< std::string, utility::vector1< std::string > > & StructFileRep::heterogen_synonyms() { return heterogen_synonyms_; }

// map for storing FORMUL records:
//  key is hetID
//  value is the chemical formula, including a potential asterisk character
std::map< std::string, std::string > const & StructFileRep::heterogen_formulae() const { return heterogen_formulae_; }
std::map< std::string, std::string > & StructFileRep::heterogen_formulae() { return heterogen_formulae_; }

/// @brief map for storing carbohydrate ResidueType base (non-variant) names; parsed from HETNAM records:
/// @details key is 6-character resID
std::map< std::string, std::string > const & StructFileRep::residue_type_base_names() const { return residue_type_base_names_; }
std::map< std::string, std::string > & StructFileRep::residue_type_base_names() { return residue_type_base_names_; }


// PDB Secondary Structure Section ////////////////////////////////////////
// Data for HELIX and/or SHEET records should be declared here if ever implemented.


// PDB Connectivity Annotation Section ////////////////////////////////////
// Data for SSBOND records should be declared here if ever implemented.

/// @brief Map for storing SSBOND records.
/// @details key is 6-character resID of 1st residue in ssbond\n
/// (A vector is needed because to futureproof if we ever handle weird disorder
/// situations.)
std::map< std::string, utility::vector1< SSBondInformation > > const & StructFileRep::ssbond_map() const { return ssbond_map_; }
std::map< std::string, utility::vector1< SSBondInformation > >       & StructFileRep::ssbond_map() { return ssbond_map_; }

/// @brief Map for storing LINK records.
/// @details Key is 6-character resID of 1st residue in link.\n
/// (A vector is needed because a single saccharide residue can have multiple branches.)
std::map< std::string, utility::vector1< LinkInformation > > const & StructFileRep::link_map() const { return link_map_; }
std::map< std::string, utility::vector1< LinkInformation > >       & StructFileRep::link_map() { return link_map_; }

// Map for storing CISPEP records.
/// @details Key is 6-character resID of 1st residue in the peptide bond.
/// @note    (A vector is needed because a single saccharide residue can have multiple branches.)
std::map< std::string, CisPeptideInformation > const & StructFileRep::cispep_map() const { return cispep_map_; }
std::map< std::string, CisPeptideInformation > & StructFileRep::cispep_map() { return cispep_map_; }


// PDB Miscellaneous Features Section /////////////////////////////////////
// Data for SITE records should be declared here if ever implemented.


// PDB Crystallographic and Coordinate Transformation Section /////////////
/// @brief Crystallographic information.
CrystInfo const & StructFileRep::crystinfo() const { return crystinfo_; }  //fpd:  CRYST1 line
CrystInfo       & StructFileRep::crystinfo() { return crystinfo_; }        //fpd:  CRYST1 line

// Data for MTRIX, ORIGX, and/or SCALE records should be declared here if ever implemented.


// PDB Coordinate Section /////////////////////////////////////////////////
/// @brief The actual atomic coordinates are stored here.
utility::vector0< ChainAtoms > const & StructFileRep::chains() const { return chains_; }
utility::vector0< ChainAtoms >       & StructFileRep::chains() { return chains_; }

// PDB Connectivity Section ///////////////////////////////////////////////
// Data for CONECT records should be declared here for consistency.

// Non-PDB Stuff //////////////////////////////////////////////////////////
/// @brief The foldtree, represented as a string.

/// @details Each file outputter must figure out how to write this out in its output format.
std::string const & StructFileRep::foldtree_string() const { return foldtree_string_; }
std::string       & StructFileRep::foldtree_string() { return foldtree_string_; }

/// @brief The PDB comments, represented as a map of string->string.
/// @details Each file outputter must figure out how to write this out in its output format.
std::map< std::string, std::string > const & StructFileRep::pdb_comments() const { return pdb_comments_; }
std::map< std::string, std::string >       & StructFileRep::pdb_comments() { return pdb_comments_; }

/// @brief A catch-all place to store additional data for output.
/// @details Each file outputter must figure out how to write this out in its output format.
std::string const & StructFileRep::additional_string_output() const { return additional_string_output_; }
std::string       & StructFileRep::additional_string_output() { return additional_string_output_; }

/// @brief Append more string data to the additional_string_output_ string in the SFR.
/// @details Retains all string data already added.
void
StructFileRep::append_to_additional_string_output(
	std::string const &input_string
) {
	additional_string_output_ = additional_string_output_ + input_string;
}


} // namespace io
} // namespace core
