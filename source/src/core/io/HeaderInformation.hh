// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/HeaderInformation.hh
///
/// @brief  Information stored in the HEADER record in the PDB format
/// @author Matthew O'Meara

#ifndef INCLUDED_core_io_HeaderInformation_hh
#define INCLUDED_core_io_HeaderInformation_hh

// Unit headers
#include <core/io/HeaderInformation.fwd.hh>
#include <core/io/pdb/Record.hh>

// Platform headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <list>
#include <utility>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>

#ifdef PYROSETTA
		#include <cereal/types/polymorphic.hpp>
#endif

#endif // SERIALIZATION

namespace core {
namespace io {

/// @brief Information stored in the header records
///http://www.wwpdb.org/documentation/format32/sect2.html
///HEADER    PEPTIDASE                               13-JAN-98   1A2Z
class HeaderInformation : public utility::pointer::ReferenceCount {
public:
	HeaderInformation();
	HeaderInformation(HeaderInformation const & src);
	~HeaderInformation() override;

	void
	store_record( pdb::Record & R );

	void
	finalize_parse();

	bool
	parse_in_progress() const;

	void
	fill_records( std::vector< pdb::Record > & VR ) const;

public: // Fields for HEADER Records
	void
	store_classification(std::string const & classification);

	std::string
	classification() const;

	/// @brief Store the deposition date using the format like "84-DEC-18"
	void
	store_deposition_date(std::string const & deposition_date);

	/// @brief Store the deposition date using the format of two digit numbers
	/// for the year, month and day of the month.
	void
	store_deposition_date(
		Size yy,
		Size mm,
		Size dd);

	std::string
	deposition_date() const;

	/// @brief get the deposition date using the format of two digit numbers
	/// for the year, month and day of the month
	void
	deposition_date(
		Size & yy,
		Size & mm,
		Size & dd) const;

	std::string
	idCode() const;

	void
	store_idCode(std::string const & id_code);

	void
	fill_header_record(std::vector< pdb::Record > & VR) const;

public: //  Fields for the TITLE Record

	/// @brief Append title record onto current title string
	void
	store_title(std::string const & title);

	void
	clear_title();

	std::string const &
	title() const;

	void
	fill_title_records(std::vector< pdb::Record > & VR) const;

public: //  Fields for KEYWDS Records

	typedef std::list< std::string > Keywords;

	void
	store_keywords(std::string const & keywords);

	void
	finalize_keyword_records();

	bool
	keyword_in_progress() const;

	void
	clear_keywords();

	Keywords const &
	keywords() const;

	// Undefined, commenting out to fix PyRosetta build  std::string keywords_string() const;

	void
	fill_keyword_records(std::vector< pdb::Record > & VR) const;

public: //  Fields for COMPND Records

	/*

	http://www.wwpdb.org/documentation/format23/sect2.html#COMPND

	The COMPND record describes the macromolecular contents of an
	entry. Each macromolecule found in the entry is described by a set of
	token: value pairs, and is referred to as a COMPND record
	component. Since the concept of a molecule is difficult to specify
	exactly, PDB staff may exercise editorial judgment in consultation
	with depositors in assigning these names.

	* In the general case the PDB tends to reflect the
	biological/functional view of the molecule. For example, the
	hetero-tetramer hemoglobin molecule is treated as a discrete
	component in COMPND.

	* In the case of synthetic molecules, e. g., hybrids, the depositor
	will provide the description.

	* No specific rules apply to the ordering of the tokens, except that
	the occurrence of MOL_ID or FRAGMENT indicates that the subsequent
	tokens are related to that specific molecule or fragment of the
	molecule.

	* Asterisks in nucleic acid names (in MOLECULE) are for ease of reading.

	* When insertion codes are given as part of the residue name, they
	must be given within square brackets, i.e., H57[A]N. This might
	occur when listing residues in FRAGMENT or OTHER_DETAILS.

	* For multi-chain molecules, e.g., the hemoglobin tetramer, a
	comma-separated list of CHAIN identifiers is used.

	* When non-blank chain identifiers occur in the entry, they must be specified.

	## Verification/Validation/Value Authority Control ##

	* CHAIN must match the chain identifiers(s) of the molecule(s). EC
	numbers are also checked

	## Relationships to Other Record Types ##

	* In the case of mutations, the SEQADV records will present
	differences from the reference molecule. REMARK records may
	further describe the contents of the entry. Also see verification
	above.
	*/

	enum CompoundToken {
		// Numbers each component; also used in SOURCE to associate the information.
		MOL_ID=1,

		// Name of the macromolecule.
		MOLECULE,

		// If the MOLECULE functions as part of a larger
		// biological unit, the entire functional unit may be
		// described.
		BIOLOGICAL_UNIT,

		// Comma-separated list of chain identifier(s).
		CHAIN,

		// Specifies a domain or region of the molecule.
		FRAGMENT,

		// Comma-separated list of synonyms for the MOLECULE.
		SYNONYM,

		// The Enzyme Commission number associated with the
		// molecule. If there is more than one EC number, they
		// are presented as a comma-separated list.
		EC,

		// Indicates that the molecule was produced using
		// recombinant technology or by purely chemical synthesis.
		ENGINEERED,

		// Indicates if there is a mutation.
		MUTATION,

		// Additional comments.
		OTHER_DETAILS,

		CompoundToken_max = OTHER_DETAILS
	};

	typedef utility::vector1< std::pair< CompoundToken, std::string > > Compounds;

	std::string
	static compound_token_to_string(CompoundToken token);

	CompoundToken
	static string_to_compound_token(std::string const & token);

	void
	store_compound(std::string const & compound);

	void
	store_compound(CompoundToken token, std::string const & value);

	Compounds const &
	compounds() const;

	void
	finalize_compound_records();

	bool
	compound_in_progress() const;

	void
	clear_compounds();

	void
	fill_compound_records(std::vector< pdb::Record > & VR) const;

public: /// Fields for the EXPDTA Record

	enum ExperimentalTechnique {
		// Experimental Techniques for spec version 3.3
		X_RAY_DIFFRACTION = 1,
		FIBER_DIFFRACTION,
		NEUTRON_DIFFRACTION,
		ELECTRON_CRYSTALLOGRAPHY,
		ELECTRON_MICROSCOPY,
		SOLID_STATE_NMR,
		SOLUTION_NMR,
		SOLUTION_SCATTERING,

		THEORETICAL_MODEL,
		ExperimentalTechnique_max_current = THEORETICAL_MODEL,

		// Obsolete technique codes
		ELECTRON_DEFRACTION,
		CRYO_ELECTRON_MICROSCOPY,
		SOLUTION_SCATTERING_THEORETICAL_MODEL,
		FLORECENCE_TRANSFER,
		NMR, // Note the qualifying information is parsed not stored

		ExperimentalTechnique_max = NMR
	};
	typedef std::list< ExperimentalTechnique > ExperimentalTechniques;

	std::string
	static experimental_technique_to_string(ExperimentalTechnique technique);

	ExperimentalTechnique
	static string_to_experimental_technique(std::string const &technique);

	/// @parse the list of techniques string following the technqiue
	/// field in the EXPDTA record of the pdb format and store the techniques
	void
	store_experimental_techniques(std::string const & exp);

	void
	store_experimental_technique(ExperimentalTechnique technique);

	ExperimentalTechniques const &
	experimental_techniques() const;

	void
	finalize_experimental_technique_records();

	bool
	experimental_technique_in_progress() const;

	void
	clear_experimental_techniques();

	bool
	is_experimental_technique(ExperimentalTechnique technique) const;

	void
	fill_experimental_technique_records(std::vector< pdb::Record > & VR) const;

private: // Helper functions

	/// @brief create enough records of <record_type> to express the
	/// <contents> string and save them into the Records vector
	void
	fill_wrapped_records(
		std::string const & record_type,
		std::string const & field_name,
		std::string const & contents,
		Size & line_no,
		std::vector< pdb::Record > & VR) const;

	void
	set_line_continuation(
		pdb::Record & R,
		Size const line_no) const;


private: // Data for HEADER Record

	/// @brief Possibly abbreviated classification type
	std::string classification_;

	/// @brief Deposition date DD-MON-YY
	Size dep_year_;
	Size dep_month_;
	Size dep_day_;

	/// @brief 4-character PDB unique identifier
	std::string idCode_;


private: // Data for TITLE Record

	std::string title_;

private: // Date for KEYWDS Record

	Keywords keywords_;

	bool keyword_in_progress_;

private: // Data for COMPND Record

	Compounds compounds_;

	bool compound_in_progress_;

private: // data for EXPDTA Record

	ExperimentalTechniques experimental_techniques_;

	std::string experimental_technique_in_progress_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace io
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_io_pdb_HeaderInformation )
#endif // SERIALIZATION


#endif // INCLUDED_core_io_pdb_HeaderInformation_HH
