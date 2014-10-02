// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/RotamerLibrary.hh
/// @brief  Rotamer Library classes and utility functions
/// @author Phil Bradley


#ifndef INCLUDED_core_pack_dunbrack_RotamerLibrary_hh
#define INCLUDED_core_pack_dunbrack_RotamerLibrary_hh

// Unit headers
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>


// Package headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueRotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Residue.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#ifdef WIN32 //VC++ needs full class declaration
	#include <core/graph/Graph.hh>
	#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
	#include <core/pack/dunbrack/SingleResidueRotamerLibrary.hh>
#endif

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>

// Numeric headers
#include <numeric/random/random.fwd.hh>

// C++ headers
#include <map>

#include <utility/vector1.hh>

namespace core {
namespace pack {
namespace dunbrack {

void
rotamer_from_chi(
	conformation::Residue const & rsd,
	RotVector & rot
);

void
rotamer_from_chi(
	chemical::ResidueType const & rsd_type,
	ChiVector const & chi,
	RotVector & rot
);

/// @brief Do not bother calling this function if you're interested in staying
/// up to date with changes in the rotamer library.  Call rotamer_from_chi instead.
/// It will eventually call this function iff dun10 is not set to true.
void
rotamer_from_chi_02(
	ChiVector const & chi,
	//ChiVector chi, //pass in by value and locally modified
	chemical::AA const res, // pass an AA, not an AA &
	RotVector & rot
);

void
rotamer_from_chi_02(
	Real4 const & chi,
	chemical::AA const res,
	Size nchi,
	Size4 & rot
);

/// @brief Find the difference in angles between two chi values using hard-coded
/// symmetry information for the symmetric amino acids.  Disappears for 2010 library.
Real
subtract_chi_angles(
	Real chi1,
	Real chi2,
	chemical::AA const & aa,
	int chino
);

///////////////////////////////////////////////////////////////////////////////
class RotamerLibrary : public utility::SingletonBase< RotamerLibrary >
{
public:
	friend class utility::SingletonBase< RotamerLibrary >;

	typedef chemical::AA AA;
	typedef conformation::Residue Residue;
	typedef chemical::ResidueType ResidueType;

public:
	virtual ~RotamerLibrary();

	Real
	rotamer_energy(
		Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

  Real
  best_rotamer_energy(
    Residue const & rsd,
    bool curr_rotamer_only,
    RotamerLibraryScratchSpace & scratch
  ) const;

	/// to do:
	Real
	rotamer_energy_deriv(
		Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const;

	///
	void
	add_residue_library(
		AA const & aa,
		SingleResidueRotamerLibraryCOP rot_lib
	) const;

	/// @brief Overrides the library for AA
	void
	add_residue_library(
		ResidueType const & rsd_type,
		SingleResidueRotamerLibraryCOP rot_lib
	) const;

	//XRW_B_T1
	//	RotamerLibraryOP coarsify(coarse::TranslatorSet const &map_set) const;
	//XRW_E_T1

	bool
	rsd_library_already_loaded( chemical::ResidueType const & rsd_type ) const;

	SingleResidueRotamerLibraryCOP
	get_rsd_library( chemical::ResidueType const & rsd_type ) const;

	SingleResidueRotamerLibraryCOP
	get_library_by_aa( chemical::AA const & aa ) const;

	//output dunbrack library  -- coarse rotamers can be cached on disk
	void
	write_to_file( std::string filename ) const;

	/// @brief Public interface to read in dunbrack libraries.  dun10 option checks
	/// are handled inside.
	void create_fa_dunbrack_libraries();

	/// @brief Get the rotamer library specified by the NCAA_ROTLIB flag in the residue type paramater file
	/// in the residue type paramater file if it exists
	SingleResidueRotamerLibraryCOP
	get_NCAA_rotamer_library( chemical::ResidueType const & rsd_type ) const;

	/// @brief Get the rotamer library specified by the PEPTOID_ROTLIB flag in the residue type paramater file
	/// in the residue type paramater file if it exists
	SingleResidueRotamerLibraryCOP
	get_peptoid_rotamer_library( chemical::ResidueType const & rsd_type ) const;


private:

	bool decide_read_from_binary() const;

	bool decide_read_from_binary_02() const;
	bool decide_read_from_binary_10() const;

	bool binary_is_up_to_date( utility::io::izstream & binlib ) const;

	bool binary_is_up_to_date_02( utility::io::izstream & binlib ) const;
	bool binary_is_up_to_date_10( utility::io::izstream & binlib ) const;

	bool decide_write_binary() const;

	bool decide_write_binary_02() const;
	bool decide_write_binary_10() const;

	std::string get_library_name_02() const;

	std::string get_binary_name() const;

	std::string get_binary_name_02() const;
	std::string get_binary_name_10() const;

	Size current_binary_format_version_id_02() const;
	Size current_binary_format_version_id_10() const;

	void create_fa_dunbrack_libraries_from_ASCII();

	void create_fa_dunbrack_libraries_from_binary();

	void create_fa_dunbrack_libraries_02_from_ASCII();
	void create_fa_dunbrack_libraries_02_from_binary();

	void create_centroid_rotamer_libraries_from_ASCII();

	void create_fa_dunbrack_libraries_10_from_ASCII();
	void create_fa_dunbrack_libraries_10_from_binary();

	void write_binary_fa_dunbrack_libraries() const;

	void write_binary_fa_dunbrack_libraries_02() const;
	void write_binary_fa_dunbrack_libraries_10() const;

	std::string random_tempname( std::string const & prefix ) const;

	template < Size T >
	void
	initialize_and_read_srsrdl(
		SemiRotamericSingleResidueDunbrackLibrary< T > & srsrdl,
		bool const nrchi_is_symmetric,
		Real const nrchi_start_angle,
		utility::io::izstream & rotamer_definitions,
		utility::io::izstream & regular_library,
		utility::io::izstream & continuous_minimization_bbdep
		// phasing out bbind sampling -- utility::io::izstream & rnchi_bbind_probabilities
	) const;

	template < Size T >
	void
	initialize_srsrdl(
		SemiRotamericSingleResidueDunbrackLibrary< T > & srsrdl,
		bool const nrchi_is_symmetric,
		Real const nrchi_start_angle
	) const;

	/// @brief Instantiate the appropriate RSRDL< T > library given the n_chi input and initialize
	/// the library from the input stream.  Return the string that the library returns after reading
	/// the from the stream (which may be the empty string, or the three letter code for the next AA
	/// in the input stream).  If the first three letter code for this amino acid has already been
	/// retrieved from the input stream, indicate so with the "first_three_letter_code_alread_read"
	/// boolean.
	SingleResidueDunbrackLibraryOP
	create_rotameric_dunlib(
		chemical::AA aa,
		Size const n_chi,
		utility::io::izstream & library,
		bool dun02,
		std::string & next_aa_in_library,
		bool first_three_letter_code_already_read
	) const;

	SingleResidueDunbrackLibraryOP
	create_rotameric_dunlib(
		chemical::AA aa,
		Size const n_chi,
		bool dun02
	) const;

	/// @brief Instantiate the appropriate SRSRDL< T > library given the nchi input.  nchi indicates
	/// the number of rotameric chi -- there are nchi+1 chi described by the library, the last of which
	/// is non-rotameric.  (E.g. Tyrosine has an nchi of 1, its chi2 is non-rotameric and its chi3 is not
	/// described by the library.)  The input streams are either read or not read depending on the
	/// boolean variables.
	SingleResidueDunbrackLibraryOP
	create_semi_rotameric_dunlib(
		chemical::AA aa,
		Size const nchi,
		bool const use_bbind_rnchi_scoring,
		bool const use_bbind_rnchi_sampling,
		bool const nrchi_is_symmetric,
		Real const nrchi_start_angle,
		utility::io::izstream & rotamer_definitions,
		utility::io::izstream & regular_library,
		utility::io::izstream & continuous_minimization_bbdep
		// phasing out bbind sampling -- utility::io::izstream & rnchi_bbind_probabilities
	) const;

	SingleResidueDunbrackLibraryOP
	create_semi_rotameric_dunlib(
		chemical::AA aa,
		Size const nchi,
		bool const use_bbind_rnchi_scoring,
		bool const use_bbind_rnchi_sampling,
		bool const nrchi_is_symmetric,
		Real const nrchi_start_angle
	) const;

	SingleResidueDunbrackLibraryOP
	create_srdl(
		chemical::AA aa_in
	) const;

public:
	static
	void
	initialize_dun10_aa_parameters(
		utility::vector1< chemical::AA > & rotameric_amino_acids,
		utility::vector1< Size > & rotameric_n_chi,
		utility::vector1< chemical::AA > & sraa,
		utility::vector1< Size > & srnchi,
		utility::vector1< bool > & scind,
		utility::vector1< bool > & sampind,
		utility::vector1< bool > & sym,
		utility::vector1< Real > & astr
	);

	static
	void
	initialize_dun02_aa_parameters(
		utility::vector1< chemical::AA > & rotameric_amino_acids,
		utility::vector1< Size > & rotameric_n_chi
	);

	static
	void
	initialize_dun02_aa_parameters(
		utility::vector1< Size > & nchi_for_aa,
		Size & n_rotameric_aas
	);

private:
	void
	write_to_binary( utility::io::ozstream & out ) const;

	void
	read_from_binary( utility::io::izstream & in );



	RotamerLibrary();
	RotamerLibrary( RotamerLibrary const & ); // unimplemented
	RotamerLibrary const & operator = ( RotamerLibrary const & ); // unimplemented

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static RotamerLibrary * create_singleton_instance();

private:

	typedef std::map< AA, SingleResidueRotamerLibraryCOP > LibraryMap;
	typedef LibraryMap::const_iterator library_iterator;

	// These entries take precedence over the entries for an amino acid.
	// The string used for the key is ResidueType.name()
	typedef std::map< std::string, SingleResidueRotamerLibraryCOP > ResLibraryMap;
	typedef ResLibraryMap::const_iterator reslibrary_iterator;

	// APL NOTE: All the mutable data in this class will need a ReadWriteMutex to
	// be made threadsafe.

	mutable ResLibraryMap reslibraries_;
	mutable LibraryMap libraries_;
	mutable utility::vector1< SingleResidueRotamerLibraryCOP > aa_libraries_;
	mutable utility::vector1< SingleResidueRotamerLibraryCOP > libraries_ops_;

	/// @brief Map that associates three letter codes with SRRLCOPs for the
	/// noncanonical alpha-amino acid sidechain rotamer libraries
	mutable std::map< std::string, SingleResidueRotamerLibraryCOP > ncaa_rotlibs_;

	/// @brief Map that associates three letter codes with SRRLCOPs for the
	/// peptoid sidechain rotamer libraries
	mutable std::map< std::string, SingleResidueRotamerLibraryCOP > peptoid_rotlibs_;

};

void
read_dunbrack_library(
	RotamerLibrary & rotamer_library
);


} // dunbrack
} // pack
} // core


#endif // INCLUDED_core_pack_dunbrack_RotamerLibrary_HH
