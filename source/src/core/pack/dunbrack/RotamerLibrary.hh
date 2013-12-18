// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/RotamerLibrary.hh
/// @brief  Rotamer Library classes
/// @author Phil Bradley


#ifndef INCLUDED_core_pack_dunbrack_RotamerLibrary_hh
#define INCLUDED_core_pack_dunbrack_RotamerLibrary_hh

// Unit headers
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>


// Package headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>
// AUTO-REMOVED #include <core/scoring/types.hh>

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
#include <core/graph/Graph.hh> // WIN32 INCLUDE
#endif

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/exit.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>

// Numeric headers
#include <numeric/random/random.fwd.hh>

// C++ headers
#include <map>

#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

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
// pure virtual base class
class SingleResidueRotamerLibrary : public utility::pointer::ReferenceCount
{
public:
	virtual
	~SingleResidueRotamerLibrary();

	virtual
	Real
	rotamer_energy_deriv(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const = 0;


	virtual
	Real
	rotamer_energy(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch
	) const = 0;

	/// @brief Returns the energy of the lowest-energy rotamer accessible to the given residue
	/// (based on e.g. its current phi and psi values).
	/// If curr_rotamer_only is true, then consider only the idealized version of the
	/// residue's current rotamer (local optimum); otherwise, consider all rotamers (global optimum).
	virtual
	Real
	best_rotamer_energy(
		conformation::Residue const & rsd,
		bool curr_rotamer_only,
		RotamerLibraryScratchSpace & scratch
	) const = 0;

	/// @brief Pick a rotamer for the input residue according to the rotamer probability
	/// distribution and assign chi angles to the input rsd.  If perturb_from_rotamer_center
	/// is true, then push the rotamer off from the center; for chi angles with a normal
	/// distribution, the perturbation is taken from a Gaussian random number with a standard
	/// deviation matching the chi angle's standard deviation.  For chi angles that are not
	/// normally distributed, the behavior is open to the derived classe's interpretation.
	virtual
	void
	assign_random_rotamer_with_bias(
		conformation::Residue const & rsd,
		RotamerLibraryScratchSpace & scratch,
		numeric::random::RandomGenerator & RG,
		ChiVector & new_chi_angles,
		bool perturb_from_rotamer_center
	) const = 0;

	virtual
	void
	fill_rotamer_vector(
		pose::Pose const & pose,
		scoring::ScoreFunction const & scorefxn,
		pack::task::PackerTask const & task,
		graph::GraphCOP packer_neighbor_graph,
		chemical::ResidueTypeCOP concrete_residue,
		conformation::Residue const& existing_residue,
		utility::vector1< utility::vector1< Real > > const & extra_chi_steps,
		bool buried,
		RotamerVector & rotamers
	) const = 0;

	//XRW_B_T1
	/*
	virtual
	SingleResidueRotamerLibraryOP
	coarsify(coarse::Translator const &map) const = 0;
	*/
	//XRW_E_T1

	virtual
	void
	write_to_file( utility::io::ozstream &out ) const = 0;

};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class RotamerLibrary // singleton -- no need to derive from RefCount : public utility::pointer::ReferenceCount
{

public:
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

	SingleResidueRotamerLibraryCAP
	get_rsd_library( chemical::ResidueType const & rsd_type ) const;

	SingleResidueRotamerLibraryCAP
	get_library_by_aa( chemical::AA const & aa ) const;

	//output dunbrack library  -- coarse rotamers can be cached on disk
	void
	write_to_file( std::string filename ) const;

	/// @brief Public interface to read in dunbrack libraries.  dun10 option checks
	/// are handled inside.
	void create_fa_dunbrack_libraries();

	pack::dunbrack::SingleResidueRotamerLibraryCAP
	get_NCAARotamerLibrary( chemical::ResidueType const & rsd_type ) const;

	static
	RotamerLibrary & get_instance();


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

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

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
	static RotamerLibrary * rotamer_library_;

private:

	// For thread safety, store the OPs in an unchanging and never-used vector;
	// keep copies of the APs in a map for actual use.
	// (OPs cannot be shared across threads because copying them onto the stack
	// results in increment / decrement operations on the reference count.)
	// If this doesn't make you weep for a language with real garbage collection,
	// I don't know what would.
	typedef std::map< AA, SingleResidueRotamerLibraryCAP > LibraryMap;
	typedef LibraryMap::const_iterator library_iterator;

	// These entries take precedence over the entries for an amino acid.
	// The string used for the key is ResidueType.name()
	typedef std::map< std::string, SingleResidueRotamerLibraryCAP > ResLibraryMap;
	typedef ResLibraryMap::const_iterator reslibrary_iterator;

	mutable ResLibraryMap reslibraries_;
	mutable LibraryMap libraries_;
	mutable utility::vector1< SingleResidueRotamerLibraryCAP > aa_libraries_;
	mutable utility::vector1< SingleResidueRotamerLibraryCOP > libraries_ops_;
	mutable std::map< std::string, pack::dunbrack::SingleResidueRotamerLibraryCOP > ncaa_rotlibs_;

};

// Unsure -- do the following methods really belong in this file?

// XRW if coarse is removed than this functions signature will have to chanhge
/* SingleResidueDunbrackLibraryOP
read_single_dunbrack_library(
	 utility::io::izstream &iunit,
	 bool coarse
); */


void
read_dunbrack_library(
	RotamerLibrary & rotamer_library
);


} // dunbrack
} // pack
} // core


#endif // INCLUDED_core_pack_dunbrack_RotamerLibrary_HH
