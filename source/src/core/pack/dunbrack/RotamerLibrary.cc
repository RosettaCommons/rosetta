// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

// Unit headers
#include <core/pack/dunbrack/RotamerLibrary.hh>

// Package headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh>

#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Project headers
#include <core/pack/task/PackerTask.hh>
#include <core/chemical/AA.hh>

// Basic headers
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/random/uniform.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

// External headers
#include <boost/lexical_cast.hpp>

// C++ Headers
#include <string>
#include <iostream>
#if defined(WIN32) || defined(__CYGWIN__)
#include <io.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include <sys/stat.h>

using basic::T;

using basic::Error;
using basic::Warning;


/*
rearrange the dunbrack arrays:

mapping from aa to single-residue-rotamer-library

to preserve current interpolation scheme, need ability to (quickly)
lookup the probability/chi-means/chi-sdevs for a given rotamer in
several torsion-angle bins.

so how to store prob/mean/sd for alt rotamers at different pp bins

quick conversion from rotno-tuple to index?

rotamer_data: probability, chi-mean, chi-sd
*/

namespace core {
namespace pack {
namespace dunbrack {

static THREAD_LOCAL basic::Tracer TR( "core.pack.dunbrack.RotamerLibrary" );

using namespace core::pack::rotamers;

/// @brief helper that'll determine the rotamer bins for a residue by asking its associated
/// rotamer library for that information.
void
rotamer_from_chi(
	conformation::Residue const & rsd,
	RotVector & rot
)
{
	debug_assert( rsd.aa() <= chemical::num_canonical_aas );

	SingleResidueRotamerLibraryCOP rotlib = rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() );
	if ( rotlib ) {
		SingleResidueDunbrackLibraryCOP dunlib( utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const > ( rotlib ));
		debug_assert( dunlib );
		dunlib->get_rotamer_from_chi( rsd.chi(), rot );
	} else { rot.resize( 0 ); }
}

/// @brief helper that'll determine the rotamer bins for a residue_type by asking its associated
/// rotamer library for that information. Takes chi vector as input parameter.
void
rotamer_from_chi(
	chemical::ResidueType const & rsd_type,
	ChiVector const & chi,
	RotVector & rot
)
{
	SingleResidueRotamerLibraryCOP rotlib = rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd_type );
	if ( rotlib ) {
		SingleResidueDunbrackLibraryCOP dunlib( utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const > ( rotlib ));
		debug_assert( dunlib );
		dunlib->get_rotamer_from_chi( chi, rot );
	} else { rot.resize( 0 ); }
}

void rotamer_from_chi_02(ChiVector const & chi, // pass by reference for speed
	chemical::AA const res, RotVector & rot) {
	rot.resize(chi.size());

	Size4 rot4;
	Real4 chi4;
	Size nchi(std::min<Size>(chi.size(), DUNBRACK_MAX_SCTOR));
	for ( Size ii = 1; ii <= nchi; ++ii ) {
		chi4[ii] = chi[ii];
	}

	rotamer_from_chi_02(chi4, res, nchi, rot4);

	for ( Size ii = 1; ii <= nchi; ++ii ) {
		rot[ii] = rot4[ii];
	}
	for ( Size ii = nchi + 1; ii <= rot.size(); ++ii ) {
		rot[ii] = 0;
	}
}

/// @brief DEPRECATED
/// convert between the real-valued chi dihedrals and the rotamer well indices.
///
/// @details This code comes directly from Roland Dunbrack.  In certain rare edge cases,
/// floating point comparison breaks down (e.g. ! x >= 120.0 && ! x < 120.0 ) and
/// Dunbrack's code leaves the rotamer well unassigned -- a value of zero. I'm modifying
/// the code to guarantee that the rotamer well is assigned, but so that
/// Dunbrack's original form is still recognizable.
/// e.g.
/// if (x) a; if (y) b; if (z) c; now reads as:
/// if (x) a; else if ( y ) b; else /*if (z)*/ c;
void rotamer_from_chi_02(Real4 const & chi, chemical::AA const res, Size nchi,
	Size4 & rot) {
	using namespace chemical;

	// This code assumes that we're dealing with a canonical aa - fail if not the case.
	// amw possibly temporarily commenting out so that this works for get_dihedral_b3aa
	//debug_assert( res <= num_canonical_aas );

	// default to 0
	// right now this means 0 for rot numbers larger than dunbrack's nchi
	// eg for H-chi or extra chi's

	for ( Size i = 1; i <= nchi; ++i ) {
		rot[i] = 0;
		Real chi_i = basic::periodic_range(chi[i], 360.0); // this saves us a vector construction/destruction
		if ( i == 1 ) { // chi 1 //KMa phospho_ser
			if ( res != aa_pro && res != aa_gly && res != aa_ala ) {
				// chi1 of all residues except P,G,A,RNA
				if ( (chi_i >= 0.0) && (chi_i <= 120.0) ) {
					rot[i] = 1;
				} else if ( std::abs(chi_i) >= 120.0 ) {
					rot[i] = 2;
				} else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/{
					rot[i] = 3;
				}
			} else if ( res == aa_pro ) {
				// chi1 of pro
				if ( chi_i >= 0.0 ) {
					rot[i] = 1;
				} else /*if ( chi_i <= 0.0 )*/{
					rot[i] = 2;
				}
			}

		} else if ( i == 2 ) { // chi 2
			if ( res == aa_arg ) {
				if ( (chi_i >= 0.0) && (chi_i <= 120.0) ) {
					rot[i] = 1;
				} else if ( std::abs(chi_i) >= 120.0 ) {
					rot[i] = 2;
				} else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/{
					rot[i] = 3;
				}

			} else if ( res == aa_glu || res == aa_his || res == aa_ile
					|| res == aa_lys || res == aa_leu || res == aa_met
					|| res == aa_gln ) {
				if ( (chi_i >= 0.0) && (chi_i <= 120.0) ) {
					rot[i] = 1;
				} else if ( std::abs(chi_i) >= 120.0 ) {
					rot[i] = 2;
				} else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/{
					rot[i] = 3;
				}

			} else if ( res == aa_asp ) {
				if ( ((chi_i >= 30.0) && (chi_i <= 90.0))
						|| ((chi_i <= -90.0) && (chi_i >= -150.0)) ) {
					rot[i] = 1;
				} else if ( ((chi_i >= -30.0) && (chi_i <= 30.0))
						|| (std::abs(chi_i) >= 150.0) ) {
					rot[i] = 2;
				} else {
					/*if ( ( ( chi_i >= -90.0 ) && ( chi_i <= -30.0 ) ) ||
					( ( chi_i >= 90.0 ) && ( chi_i <= 150.0 ) ) ) */
					rot[i] = 3;
				}

			} else if ( (res == aa_phe) || (res == aa_tyr) ) {
				if ( ((chi_i >= 30.0) && (chi_i <= 150.0))
						|| ((chi_i <= -30.0) && (chi_i >= -150.0)) ) {
					rot[i] = 1;
				} else { /*if ( ( ( chi_i >= -30.0 ) && ( chi_i <= 30.0 ) ) ||
					( std::abs(chi_i) >= 150.0 ) ) */
					rot[i] = 2;
				}

			} else if ( res == aa_trp ) {
				if ( (chi_i >= -180.0) && (chi_i <= -60.0) ) {
					rot[i] = 1;
				} else if ( (chi_i >= -60.0) && (chi_i <= 60.0) ) {
					rot[i] = 2;
				} else /*if ( ( chi_i >= 60.0 ) && ( chi_i <= 180.0 ) ) */{
					rot[i] = 3;
				}

			} else if ( res == aa_asn ) { // chi2 of asn
				// ctsa - note this is a special case of the new dunbrack rotamer set
				if ( rot[1] == 1 ) {
					if ( chi_i >= -150.0 && chi_i < -90.0 ) {
						rot[i] = 1;
					} else if ( chi_i >= -90.0 && chi_i < -30.0 ) {
						rot[i] = 2;
					} else if ( chi_i >= -30.0 && chi_i < 30.0 ) {
						rot[i] = 3;
					} else if ( chi_i >= 30.0 && chi_i < 90.0 ) {
						rot[i] = 4;
					} else if ( chi_i >= 90.0 && chi_i < 150.0 ) {
						rot[i] = 5;
					} else /*if ( std::abs(chi_i) >= 150.0 )*/{
						rot[i] = 6;
					}
				} else if ( rot[1] == 2 ) {
					if ( chi_i >= -180.0 && chi_i < -90.0 ) {
						rot[i] = 1;
					} else if ( chi_i >= -90.0 && chi_i < -45.0 ) {
						rot[i] = 2;
					} else if ( chi_i >= -45.0 && chi_i < 0.0 ) {
						rot[i] = 3;
					} else if ( chi_i >= 0.0 && chi_i < 45.0 ) {
						rot[i] = 4;
					} else if ( chi_i >= 45.0 && chi_i < 90.0 ) {
						rot[i] = 5;
					} else /*if ( chi_i >=  90.0 && chi_i <= 180.0 )*/{
						rot[i] = 6;
					}
				} else {
					if ( chi_i >= -180.0 && chi_i < -105.0 ) {
						rot[i] = 1;
					} else if ( chi_i >= -105.0 && chi_i < -45.0 ) {
						rot[i] = 2;
					} else if ( chi_i >= -45.0 && chi_i < 15.0 ) {
						rot[i] = 3;
					} else if ( chi_i >= 15.0 && chi_i < 60.0 ) {
						rot[i] = 4;
					} else if ( chi_i >= 60.0 && chi_i < 120.0 ) {
						rot[i] = 5;
					} else /*if ( chi_i >= 120.0 && chi_i <= 180.0 )*/{
						rot[i] = 6;
					}
				}
			} else if ( res == aa_pro ) {
				// apl proline has only two rotamers, distinguishable by chi1.  Chi2 is always 1.
				rot[i] = 1;
			}

		} else if ( i == 3 ) { // chi 3
			if ( res == aa_glu ) {
				if ( ((chi_i >= 30.0) && (chi_i <= 90.0))
						|| ((chi_i <= -90.0) && (chi_i >= -150.0)) ) {
					rot[i] = 1;
				} else if ( ((chi_i >= -30.0) && (chi_i <= 30.0))
						|| (std::abs(chi_i) >= 150.0) ) {
					rot[i] = 2;
				} else { /*if ( ( ( chi_i >= -90.0 ) && ( chi_i <= -30.0 ) ) ||
					( ( chi_i >= 90.0 ) && ( chi_i <= 150.0 ) ) )*/
					rot[i] = 3;
				}

			} else if ( res == aa_arg || res == aa_lys || res == aa_met ) {
				if ( (chi_i >= 0.0) && (chi_i <= 120.0) ) {
					rot[i] = 1;
				} else if ( std::abs(chi_i) > 120.0 ) {
					rot[i] = 2;
				} else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/{
					rot[i] = 3;
				}

			} else if ( res == aa_gln ) { // chi3 of gln
				// ctsa - note this is a special case of the new dunbrack rotamer set
				if ( rot[2] == 2 ) {
					if ( chi_i >= 135.0 || chi_i < -135.0 ) {
						rot[i] = 1;
					} else if ( chi_i >= -135.0 && chi_i < -45.0 ) {
						rot[i] = 2;
					} else if ( chi_i >= -45.0 && chi_i < 45.0 ) {
						rot[i] = 3;
					} else /*if ( chi_i >=  45.0 && chi_i < 135.0 )*/{
						rot[i] = 4;
					}
				} else {
					if ( chi_i >= -180.0 && chi_i < -90.0 ) {
						rot[i] = 1;
					} else if ( chi_i >= -90.0 && chi_i < 0.0 ) {
						rot[i] = 2;
					} else if ( chi_i >= 0.0 && chi_i < 90.0 ) {
						rot[i] = 3;
					} else /*if ( chi_i >=  90.0 && chi_i <= 180.0 )*/{
						rot[i] = 4;
					}
				}
			} else if ( res == aa_pro ) {
				// apl
				rot[i] = 1;
			}
		} else if ( i == 4 ) {  // chi 4
			if ( res == aa_arg || res == aa_lys ) {
				if ( (chi_i >= 0.0) && (chi_i <= 120.0) ) {
					rot[i] = 1;
				} else if ( std::abs(chi_i) > 120.0 ) {
					rot[i] = 2;
				} else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/{
					rot[i] = 3;
				}
			}
		}
	} // i=1, i<= nchi
}

///////////////////////////////////////////////////////////////////////////////
// taken from rotamer_functions -- total temporary hack
//
// removed pser and rna

RotamerLibrary::RotamerLibrary():
	aa_libraries_( chemical::num_canonical_aas ) // Resize with default (empty OP) constructor
{
	this->create_fa_dunbrack_libraries();
}

RotamerLibrary::~RotamerLibrary() {
}

Real
RotamerLibrary::rotamer_energy(
	Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	SingleResidueRotamerLibraryCOP const library( rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() ) );

	if ( library ) {
		return library->rotamer_energy(rsd, scratch);
	}
	return 0.0;
}

Real
RotamerLibrary::best_rotamer_energy(
	Residue const & rsd,
	bool curr_rotamer_only,
	RotamerLibraryScratchSpace & scratch
) const
{
	SingleResidueRotamerLibraryCOP const library( rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() ) );

	if ( library ) {
		if ( curr_rotamer_only ) {
			return library->best_rotamer_energy( rsd, true,scratch );
		} else {
			return library->best_rotamer_energy( rsd, false,scratch );
		}
	} else {
		return 0.0;
	}
}

Real
RotamerLibrary::rotamer_energy_deriv(
	Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	SingleResidueRotamerLibraryCOP const library( rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( rsd.type() ) );

	if ( library ) {
		return library->rotamer_energy_deriv(rsd, scratch);
	}

	Real5 & dE_dbb( scratch.dE_dbb() );
	Real4 & dE_dchi( scratch.dE_dchi() );

	// ensure that these guys are dimensioned
	std::fill( dE_dbb.begin(), dE_dbb.end(), 0 );
	std::fill( dE_dchi.begin(), dE_dchi.end(), 0 );
	return 0.0;

}

void
RotamerLibrary::add_residue_library(
	AA const & aa,
	SingleResidueRotamerLibraryCOP rot_lib
)
{
	if ( aa > chemical::num_canonical_aas ) {
		TR.Error << "Cannot add fullatom Dunbrack rotamer library of type " << aa << ": not a canonical amino acid." << std::endl;
		utility_exit_with_message("Cannot add a non-canonical fullatom Dunbrack library.");
	}
	if ( aa_libraries_[ aa ] != 0 ) {
		TR.Error << "Cannot add fullatom Dunbrack rotamer library of type " << aa << ": library already loaded." << std::endl;
		utility_exit_with_message("Can't add rsd library twice");
	}
	aa_libraries_[ aa ] = rot_lib;
}

SingleResidueRotamerLibraryCOP
RotamerLibrary::get_library_by_aa( chemical::AA const & aa ) const
{
	if ( (Size) aa <= aa_libraries_.size() ) {
		return aa_libraries_[ aa ];
	}
	TR.Error << "Cannot get fullatom Dunbrack rotamer library of type " << aa << ": not a canonical amino acid." << std::endl;
	utility_exit_with_message("Cannot get non-canonical fullatom Dunbrack library.");
	return SingleResidueRotamerLibraryCOP();
}

///////////////////////////////////////////////////////////////////////////////
//
// this is going to be wicked slow!!!
// apl -- I love Phil's comment above, and though it no longer makes sense to
// keep around, I have trouble deleting it.  Input is pretty fast now.
// (RM -- the Grade II* listed comments above were originally for a calling function,
// but were moved here as a heritage conservation effort.)

/// @details Generic decision tree for both 02 and 2010 libraries.  Requires
/// that low-level subroutines dispatch to 02 and 2010 versions based on
/// the dun10 flag.
void
RotamerLibrary::create_fa_dunbrack_libraries() {
	//MaximCode
	if ( basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable] ) {
		TR << "shapovalov_lib_fixes_enable option is true." << std::endl;
	} else {
		//MaximCode
		//Not to break some integration tests:
		//TR << "shapovalov_lib_fixes_enable option is false" << std::endl;
	}

	if ( decide_read_from_binary() ) {
		create_fa_dunbrack_libraries_from_binary();
	} else {
		create_fa_dunbrack_libraries_from_ASCII();
		if ( decide_write_binary() ) {
			write_binary_fa_dunbrack_libraries();
		}
	}
}

bool
RotamerLibrary::decide_read_from_binary() const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[in::file::no_binary_dunlib] ) {
		return false;
	}

	std::string binary_filename = get_binary_name();
	if ( binary_filename.size() == 0 ) {
		return false;
	}

	//MaximCode
	if (
			option[ corrections::shapovalov_lib_fixes_enable ] &&
			option[ corrections::shapovalov_lib::shap_dun10_enable ]
			) {
		std::string _smoothingRequsted = option[ corrections::shapovalov_lib::shap_dun10_smooth_level ];
		std::string _smoothingAsKeyword = "undefined";

		if ( _smoothingRequsted.compare("1")==0 ) {
			_smoothingAsKeyword = "lowest_smooth";
		} else if ( _smoothingRequsted.compare("2")==0 ) {
			_smoothingAsKeyword = "lower_smooth";
		} else if ( _smoothingRequsted.compare("3")==0 ) {
			_smoothingAsKeyword = "low_smooth";
		} else if ( _smoothingRequsted.compare("4")==0 ) {
			_smoothingAsKeyword = "average_smooth";
		} else if ( _smoothingRequsted.compare("5")==0 ) {
			_smoothingAsKeyword = "higher_smooth";
		} else if ( _smoothingRequsted.compare("6")==0 ) {
			_smoothingAsKeyword = "highest_smooth";
		} else {
			_smoothingAsKeyword = "unknown";
		}
		TR << "shapovalov_lib::shap_dun10_smooth_level of " << _smoothingRequsted << "( aka " << _smoothingAsKeyword << " )" << " got activated." << std::endl;
	}
	if ( option[ corrections::shapovalov_lib_fixes_enable ] ) {
		TR << "Binary rotamer library selected: " << binary_filename << std::endl;
	}

	utility::io::izstream binlib(binary_filename.c_str(),
		std::ios::in | std::ios::binary);

	if ( !binlib ) {
		TR.Info << "cannot find binary dunbrack library under name "
			<< binary_filename << std::endl;
		return false;
	}

	if ( !binary_is_up_to_date(binlib) ) {
		TR.Info << "dunbrack library is out-dated " << binary_filename
			<< std::endl;
		return false;
	}

	/// Look for specific objections for the 02 or 08 libraries (if any exist!)
	if ( option[corrections::score::dun10] ) {
		if ( !decide_read_from_binary_10() ) {
			return false;
		}
	} else {
		if ( !decide_read_from_binary_02() ) {
			return false;
		}
	}

	return true;
}

/// @details Put specific objections to reading the 02 library from binary here.
/// This function alone will not answer whether the binary file should be read;
/// rather, it is called by the more generic "decide_read_from_binary" function.
bool RotamerLibrary::decide_read_from_binary_02() const {
	/// TO DO: check time stamp of ascii file vs time stamp of binary file
	return true;
}

/// @details Put specific objections to reading the 10 library from binary here.
/// This function alone will not answer whether the binary file should be read;
/// rather, it is called by the more generic "decide_read_from_binary" function.
bool RotamerLibrary::decide_read_from_binary_10() const {
	/// TO DO: check time stamp of ascii files vs time stamp of binary file
	return true;
}

/// @details Dispatches binary file is up-to-date to appropriate library given
/// flag status (dun10 or dun08 yes/no).
bool RotamerLibrary::binary_is_up_to_date(
	utility::io::izstream & binlib) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[corrections::score::dun10] ) {
		return binary_is_up_to_date_10(binlib);
	} else {
		return binary_is_up_to_date_02(binlib);
	}
}

bool RotamerLibrary::binary_is_up_to_date_02(
	utility::io::izstream & binlib) const {
	boost::int32_t version(0);
	binlib.read((char*) &version, sizeof(boost::int32_t));
	return static_cast<Size>(version) == current_binary_format_version_id_02();
}

/// @details Binary file is up-to-date if all the following are true:
/// 1. Version # for binary file is equal to the current version #
/// 2. All the hard coded constants for the rotameric and semirotameric libraries
/// are the same from the previous hard coded values.
///
/// DANGER DANGER DANGER - True for dun10 as well? -JAB
/// When we reopen the binary file in create_fa_dunbrack_libraries_08_from_binary
/// we will want to skip past all this data; the code in that section assumes that
/// it will read 8 arrays: 2 representing the rotameric amino acids, and 6
/// representing the semirotameric amino acids.  The arithmetic it performs depends
/// on this basic structure not changing, therefore, if you change this structure,
/// you must also change the seekg arithmetic in the create...08_from_binary function.

bool RotamerLibrary::binary_is_up_to_date_10(
	utility::io::izstream & binlib) const {
	/// 1.
	boost::int32_t version(0);
	binlib.read((char*) &version, sizeof(boost::int32_t));
	if ( static_cast<Size>(version) != current_binary_format_version_id_10() ) {
		return false;
	}

	boost::int32_t nrotameric(0);
	binlib.read((char*) &nrotameric, sizeof(boost::int32_t));
	boost::int32_t nsemirotameric(0);
	binlib.read((char*) &nsemirotameric, sizeof(boost::int32_t));

	utility::vector1< chemical::AA > rotameric_amino_acids;
	utility::vector1< Size > rotameric_n_chi;
	utility::vector1< Size > rotameric_n_bb;

	utility::vector1< chemical::AA > sraa;
	utility::vector1< Size > srnchi;
	utility::vector1< Size > srnbb;
	utility::vector1< bool > scind;
	utility::vector1< bool > sampind;
	utility::vector1< bool > sym;
	utility::vector1< Real > astr;

	initialize_dun10_aa_parameters(
		rotameric_amino_acids, rotameric_n_chi, rotameric_n_bb,
		sraa, srnchi, srnbb, scind, sampind, sym, astr );

	/// 2.
	if ( static_cast< Size > (nrotameric) != rotameric_amino_acids.size() ) return false;
	if ( static_cast< Size > (nsemirotameric) != sraa.size() ) return false;

	boost::int32_t * rotaa_bin   = new boost::int32_t[ nrotameric ];
	boost::int32_t * rot_nchi_bin= new boost::int32_t[ nrotameric ];

	boost::int32_t * sraa_bin = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * srnchi_bin  = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * scind_bin   = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * sampind_bin = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * sym_bin  = new boost::int32_t[ nsemirotameric ];
	Real * astr_bin   = new Real[ nsemirotameric ];

	binlib.read( (char*) rotaa_bin, nrotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) rot_nchi_bin, nrotameric * sizeof( boost::int32_t  ));

	binlib.read( (char*) sraa_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) srnchi_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) scind_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) sampind_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) sym_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) astr_bin, nsemirotameric * sizeof( Real  ));

	bool good( true );
	for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
		if ( rotaa_bin[ ii - 1 ] != static_cast< boost::int32_t > ( rotameric_amino_acids[ ii ] ) ) good = false;
		if ( rot_nchi_bin[ ii - 1 ] != static_cast< boost::int32_t > ( rotameric_n_chi[ ii ] ) ) good = false;
	}
	for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
		if ( sraa_bin[ii - 1] != static_cast<boost::int32_t>(sraa[ii]) ) {
			good = false;
		}
		if ( srnchi_bin[ii - 1] != static_cast<boost::int32_t>(srnchi[ii]) ) {
			good = false;
		}
		if ( scind_bin[ii - 1] != static_cast<boost::int32_t>(scind[ii]) ) {
			good = false;
		}
		if ( sampind_bin[ii - 1] != static_cast<boost::int32_t>(sampind[ii]) ) {
			good = false;
		}
		if ( sym_bin[ii - 1] != static_cast<boost::int32_t>(sym[ii]) ) {
			good = false;
		}
		if ( astr_bin[ii - 1] != astr[ii] ) {
			good = false;
		}
	}

	delete[] rotaa_bin;
	delete[] rot_nchi_bin;
	delete[] sraa_bin;
	delete[] srnchi_bin;
	delete[] scind_bin;
	delete[] sampind_bin;
	delete[] sym_bin;
	delete[] astr_bin;

	return good;
}

/// @details First checks general conditions under which a binary library file should
/// not be writen. It relies on get_binary_name() and binary_is_up_to_date() to
/// dispatch to the 02 and 08 data.  If the binary file doesn't exist or isn't up-to-date
/// it then checks if there are any specific 02 or 08 objections to writing the
/// binary file.  If no one objects, then this function gives a green light.
bool RotamerLibrary::decide_write_binary() const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[in::file::no_binary_dunlib] ) {
		return false;
	}
	if ( option[out::file::dont_rewrite_dunbrack_database] ) {
		return false;
	}

	/// place other objections here, e.g. if we're running from a (genuine) database or on BOINC

#ifdef BOINC
	return false;
#endif

#if (defined WIN32) && (!defined PYROSETTA) // binary file handling doesnt appear to work properly on windows 64 bit machines.
	return false;
#endif

	std::string binary_filename = get_binary_name();
	if ( binary_filename.size() > 0 ) {
		utility::io::izstream binlib(binary_filename.c_str(),
			std::ios::in | std::ios::binary);

		if ( binlib && binary_is_up_to_date(binlib) ) {
			/// Do not overwrite a binary if it already exists and is up-to-date.
			return false;
		}
	}
	/// if we get this far, then the binary file either doesn't exist or isn't up-to-date.

	/// check if there are any objections to writing the library
	if ( option[corrections::score::dun10] ) {
		if ( !decide_write_binary_10() ) {
			return false;
		}
	} else {
		if ( !decide_write_binary_02() ) {
			return false;
		}
	}

	// otherwise, write it.
	return true;

}

/// @details Add specific objections to the writing of the 02 library here.
bool RotamerLibrary::decide_write_binary_02() const {
	return true;
}

/// @details Add specific objections to the writing of the 10 library here.
bool RotamerLibrary::decide_write_binary_10() const {
	return true;
}

std::string RotamerLibrary::get_library_name_02() const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	return basic::database::full_name(option[corrections::score::dun02_file]);
}

std::string RotamerLibrary::get_library_name_10() const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string dirname;
	if ( !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
			|| !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable] ) {
		dirname = basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir];
	} else {
		dirname = basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_dir];
	}
	// The dun10 directory can be a non-database directory
	if ( utility::file::file_exists( dirname ) && utility::file::is_directory( dirname ) ) {
		TR.Debug << "Using local directory '" << dirname << "' for dun10 library." << std::endl;
		return dirname;
	}
	return basic::database::full_name(dirname);
}

// Can return an empty string if the appropriate binary file does not exist.
std::string RotamerLibrary::get_binary_name() const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[corrections::score::dun10] ) {
		return get_binary_name_10();
	} else {
		return get_binary_name_02();
	}
}

std::string RotamerLibrary::get_binary_name_02(bool for_writing /*=false*/ ) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// Keep the binaries for multiple dun02 libraries seperate.
	return basic::database::full_cache_name( option[corrections::score::dun02_file] + ".Dunbrack02.lib.bin",
		get_library_name_02(),
		for_writing );
}

std::string RotamerLibrary::get_binary_name_10(bool for_writing /*=false*/ ) const {
	std::string dirname;
	if ( !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
			|| !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable] ) {
		dirname = basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir];
	} else {
		dirname = basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_dir];
	}

	return basic::database::full_cache_name( dirname + "/Dunbrack10.lib.bin",
		get_library_name_10() + "/Dunbrack10.lib", // Needed, as one layer will be stripped off.
		for_writing );
}

/// @details The older binary formats did not have a version number, instead,
/// a somewhat stupid time-stamp was hard coded and compared against the time
/// stamp of the binary file that existed... e.g. if a new binary format was
/// committed to trunk on 1/1/08 12:00 pm, and the time stamp on a particular binary
/// was from 12/1/07, then the old binary file would be overwritten.  However,
/// if the time stamp on some other binary file was from 2/1/08, the code assumed
/// that the binary file was up to date.  This is a poor assumption, since the
/// code from 12/1/07 might have been run on 2/1/08 to generate a binary file.
/// The time stamp would have falsely indicated that the 2/1/08 generated file
/// was of the same binary format imposed on 1/1/08.
///
/// Version numbering avoids this issue.  If the 12/1/07 code wrote version # 3,
/// the 1/1/08 code could look for version # 4 in a file generated on 2/1/08
/// and use the mismatch as an indication of an out-of-date binary format.
///
/// The challenge is to come up with a versioning scheme that handles the fact that
/// until now, there haven't been version #'s in the binary files.  How can a
/// binary file generated by a pre-versioning piece of code be detected?  The trick
/// is to look at the data in the older binary file at the position which will in
/// the future be used to hold a version number.  The first piece of data in the
/// old binary files was a count of the number of libraries in the binary file:
/// 18.  Therefore, any starting version number beyond 18 would distinguish the
/// pre-versioned files from the versioned files.
///
/// When changing the binary format, document the date of the change in the comment
/// section below and increment the version number in the function.  Do not delete
/// old comments from the comment block below.
///
/// Version 19: start vesion number.  8/9/08. Andrew Leaver-Fay
/// Version 20: Limit zero-probability rotamers to the resolution of the library (1e-4).
///    6/16/2009.  Andrew Leaver-Fay
/// Version 21: Write out bicubic spline parameters for -log(p(rot|phi,psi)).  3/28/2012.  Andrew Leaver-Fay
/// Version 22: Bicubic spline bugfix  3/29/2012.  Andrew Leaver-Fay
/// Version 23: Derivatives now are output in different order for better generality
///    with respect to number of backbone dihedrals. Also, resolution limit
///    repaired ( prob <= 1e-6 not prob == 0 ) 1/12/15 Andy Watkins
/// Version 24: Changed to templating on number of bbs, can't prove that this isn't needed 1/15/15 Andy Watkins
/// Version 25: Fixed input bug that was definitely having real effects 1/22/15 Andy Watkins
Size
RotamerLibrary::current_binary_format_version_id_02() const
{
	return 25;
}

/// @details Version number for binary format.  See comments for 02 version.
///
/// When changing the binary format, document the date of the change in the comment
/// section below and increment the version number in the function.  Do not delete
/// old comments from the comment block below.
///
/// Version 1: start vesion number.  6/3/11. Andrew Leaver-Fay
/// Version 2: Rotameric residues write out bicubic spline parameters for -log(p(rot|phi,psi)).  3/28/2012.  Andrew Leaver-Fay
/// Version 3: Tricubic interpolation data for semi-rotameric residues.  3/29/2012.  Andrew Leaver-Fay
/// Version 4: Generate and write out bicubic spline data for the rotameric portion of the semi-rotameric residues .  3/29/2012.  Andrew Leaver-Fay
/// Version 5: Derivatives now are output in different order for better generality
///    with respect to number of backbone dihedrals. Also, resolution limit
///    repaired ( prob <= 1e-6 not prob == 0 ) 1/12/15 Andy Watkins
/// Version 6: Changed to templating on number of bbs, can't prove that this isn't needed 1/15/15 Andy Watkins
/// Version 7: Fixed input bug that was definitely having real effects 1/22/15 Andy Watkins
/// Version 8: Update to beta_nov16 corrected library 3/2017 fpd
Size
RotamerLibrary::current_binary_format_version_id_10() const
{
	return 8;
}

void RotamerLibrary::create_fa_dunbrack_libraries_from_ASCII() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[corrections::score::dun10] ) {
		return create_fa_dunbrack_libraries_10_from_ASCII();
	} else {
		return create_fa_dunbrack_libraries_02_from_ASCII();
	}
}

void RotamerLibrary::create_fa_dunbrack_libraries_from_binary() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[corrections::score::dun10] ) {
		return create_fa_dunbrack_libraries_10_from_binary();
	} else {
		return create_fa_dunbrack_libraries_02_from_binary();
	}
}

void
RotamerLibrary::create_fa_dunbrack_libraries_02_from_ASCII()
{
	using namespace chemical;

	////////////////////////////////////////////
	//// DATA FOR THE ROTAMERIC AMINO ACIDS
	////////////////////////////////////////////
	utility::vector1< Size > nchi_for_aa( num_canonical_aas, 0 );
	utility::vector1< Size > nbb_for_aa( num_canonical_aas, 0 );
	Size nrot_aas;

	utility::vector1< Size > memory_use_rotameric( num_canonical_aas, 0 );

	initialize_dun02_aa_parameters( nchi_for_aa, nbb_for_aa, nrot_aas );


	/// Now read in the 02 library.
	clock_t starttime = clock();
	utility::io::izstream libstream(get_library_name_02());

	chemical::AA aan = chemical::aa_unk;
	std::string nextaa;
	libstream >> nextaa;

	//std::cout << "first aa: " << nextaa << std::endl;

	SingleResidueDunbrackLibraryOP rot_lib(0);
	Size count_libraries_read(0);
	while ( nextaa != "" ) {

		aan = chemical::aa_from_name( nextaa );
		SingleResidueDunbrackLibraryOP newlib = create_rotameric_dunlib(
			aan, nchi_for_aa[ aan ], nbb_for_aa[ aan ], libstream, true /*dun02*/, nextaa, true /*first aa string already read*/ );
		++count_libraries_read;

		add_residue_library( aan, newlib );
		memory_use_rotameric[ aan ] = newlib->memory_usage_in_bytes();

		//std::cout << "next aa: " << nextaa << " #" << count_libraries_read + 1 << std::endl;

	}
	libstream.close();

	if ( count_libraries_read != nrot_aas ) {
		std::cerr << "ERROR: expected to read " << nrot_aas
			<< " libraries from Dun02, but read " << count_libraries_read
			<< std::endl;
		utility_exit();
	}

	clock_t stoptime = clock();
	TR << "Dunbrack library took "
		<< ((double) stoptime - starttime) / CLOCKS_PER_SEC
		<< " seconds to load from ASCII" << std::endl;

	//TR << "Memory usage: " << std::endl;
	//Size total( 0 );
	//for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {
	// total += memory_use_rotameric[ ii ];
	//if ( memory_use_rotameric[ii] != 0 ) {
	//TR << chemical::name_from_aa( AA( ii ) ) << " with " << memory_use_rotameric[ ii ] << " bytes" << std::endl;
	//}
	//}
	//TR << "Total memory on Dunbrack Libraries: " << total << " bytes." << std::endl;

}

void RotamerLibrary::create_fa_dunbrack_libraries_02_from_binary() {
	clock_t starttime = clock();

	std::string binary_filename = get_binary_name_02();
	TR << "Using Dunbrack library binary file '" << binary_filename << "'." << std::endl;
	utility::io::izstream binlib(binary_filename.c_str(),
		std::ios::in | std::ios::binary);

	if ( !binlib ) {
		utility_exit_with_message(
			"Could not open binary Dunbrack02 file from database -- how did this happen?");
	}

	/// READ PREABMLE
	/// Even if binary file should change its structure in the future, version number should always be
	/// the first piece of data in the file.
	boost::int32_t version(0);
	binlib.read((char*) &version, sizeof(boost::int32_t));
	/// END PREABMLE

	debug_assert( static_cast< Size > (version) == current_binary_format_version_id_02() );

	read_from_binary( binlib );

	clock_t stoptime = clock();
	TR << "Dunbrack library took "
		<< ((double) stoptime - starttime) / CLOCKS_PER_SEC
		<< " seconds to load from binary" << std::endl;

}

void RotamerLibrary::create_fa_dunbrack_libraries_10_from_ASCII() {
	using namespace chemical;
	using namespace utility;

	TR << "Reading Dunbrack Libraries" << std::endl;

	std::string const dun10_dir( get_library_name_10() + "/" );

	clock_t starttime = clock();

	/// Hard code some stuff.
	/// Rotameric library file names are of the form
	/// <3lc>.Jun9.bbdep.regular.lib
	/// Semi-Rotameric library files names are of the form
	/// 1. <3lc>.Jun9.bbdep.regular.lib
	/// 2. <3lc>.Jun9.bbdep.contMin.lib
	/// 3. <3lc>.Jun9.bbind.chi<#>.Definitions.lib
	/// 4. <3lc>.Jun9.bbind.chi<#>.Probabities.lib

	//std::string const db_subdir = "dun10/";
	std::string const regular_lib_suffix = ".bbdep.rotamers.lib";
	std::string const bbdep_contmin = ".bbdep.densities.lib";
	std::string const bbind_midfix  = ".bbind.chi";
	std::string const bbind_defs = ".Definitions.lib";
	// Note: No Backbone-Independent NRCHI representation in '10 library
	// std::string const bbind_probs   = ".Probabilities.lib";

	////////////////////////////////////////////
	//// DATA FOR THE ROTAMERIC AMINO ACIDS
	////////////////////////////////////////////
	utility::vector1< chemical::AA > rotameric_amino_acids;
	utility::vector1< Size > rotameric_n_chi;
	utility::vector1< Size > rotameric_n_bb;
	utility::vector1< Size > memory_use_rotameric;

	////////////////////////////////////////////
	//// DATA FOR THE SEMI-ROTAMERIC AMINO ACIDS
	////////////////////////////////////////////

	/// AA code for entry i
	utility::vector1<chemical::AA> sraa;
	/// Number of rotameric chi for entry i; the index of the non-rotameric chi is +1 of the value stored.
	utility::vector1< Size > srnchi;
	/// N bb
	utility::vector1< Size > srnbb;
	/// Decision: SCore the nonrotameric chi in a backbone INDependent manner?
	utility::vector1<bool> scind;
	/// Decision: SAMPle the nonrotameric chi in a backbone INDependent manner?
	utility::vector1<bool> sampind;
	/// Is there symmetry about the rotameric chi?  If so, periodicity is 180 degrees, else, 360.
	/// Furthermore, if it's symmetric, then the bbdep sampling is at 5 degree intervals and
	/// the bbind is at 0.5 degree intervals; otherwise, bbdep is at 10 degrees, and bbind is at 1.
	utility::vector1<bool> sym;
	/// At what starting angle does the chi data come from?
	utility::vector1<Real> astr;

	utility::vector1<Size> memory_use_semirotameric;

	/// Initialize the data.

	initialize_dun10_aa_parameters(
		rotameric_amino_acids, rotameric_n_chi, rotameric_n_bb,
		sraa, srnchi, srnbb, scind, sampind, sym, astr
	); // use the same parameters as dun08 used to.

	for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
		std::string next_aa_in_library;

		chemical::AA ii_aa(rotameric_amino_acids[ii]);
		std::string ii_lc_3lc(chemical::name_from_aa(ii_aa));
		std::transform(ii_lc_3lc.begin(), ii_lc_3lc.end(), ii_lc_3lc.begin(),
			tolower);
		std::string library_name = dun10_dir + ii_lc_3lc + regular_lib_suffix;
		utility::io::izstream lib(library_name.c_str());
		if ( !lib.good() ) {
			utility_exit_with_message(
				"Unable to open database file for dun10 rotamer library: "
				+ library_name);
		}
		TR << "Reading " << library_name << std::endl;

		SingleResidueDunbrackLibraryOP newlib = create_rotameric_dunlib(
			ii_aa, rotameric_n_chi[ ii ], rotameric_n_bb[ ii ], lib, false, next_aa_in_library, false );
		if ( next_aa_in_library != "" ) {
			std::cerr << "ERROR: Inappropriate read from dun10 input while reading " << ii_aa
				<< ": \"" << next_aa_in_library << "\"" << std::endl;
			utility_exit();
		}

		memory_use_rotameric.push_back(newlib->memory_usage_in_bytes());

		add_residue_library( ii_aa, newlib );
	}

	for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
		chemical::AA ii_aa(sraa[ii]);
		std::string ii_lc_3lc(chemical::name_from_aa(ii_aa));
		std::transform(ii_lc_3lc.begin(), ii_lc_3lc.end(), ii_lc_3lc.begin(),
			tolower);
		Size const nrchi = srnchi[ii] + 1;

		std::string reg_lib_name = dun10_dir + ii_lc_3lc + regular_lib_suffix;
		std::string conmin_name = dun10_dir + ii_lc_3lc + bbdep_contmin;
		std::string rotdef_name = dun10_dir + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_defs;

		utility::io::izstream lib(  reg_lib_name.c_str() );
		utility::io::izstream contmin( conmin_name.c_str()  );
		utility::io::izstream defs( rotdef_name.c_str()  );
		//utility::io::izstream probs(   probs_name.c_str()   );

		if ( !lib.good() ) {
			utility_exit_with_message(
				"Unable to open database file for Dun10 rotamer library: "
				+ reg_lib_name);
		}
		if ( !contmin.good() ) {
			utility_exit_with_message(
				"Unable to open database file for Dun10 rotamer library: "
				+ conmin_name);
		}
		if ( !defs.good() ) {
			utility_exit_with_message(
				"Unable to open database file for Dun10 rotamer library: "
				+ rotdef_name);
		}

		TR << "Reading " << reg_lib_name << std::endl;
		TR << "Reading " << conmin_name << std::endl;
		TR << "Reading " << rotdef_name << std::endl;
		//TR << "Reading " << probs_name << std::endl;

		SingleResidueDunbrackLibraryOP newlib = create_semi_rotameric_dunlib(
			ii_aa, srnchi[ ii ], srnbb[ ii ],
			scind[ ii ], sampind[ ii ],
			sym[ ii ], astr[ ii ],
			defs, lib, contmin );

		memory_use_semirotameric.push_back(newlib->memory_usage_in_bytes());

		add_residue_library( ii_aa, newlib );
	}
	TR << "Finished reading Dunbrack Libraries" << std::endl;

	clock_t stoptime = clock();
	TR << "Dunbrack 2010 library took "
		<< ((double) stoptime - starttime) / CLOCKS_PER_SEC
		<< " seconds to load from ASCII" << std::endl;

	TR << "Memory usage: " << std::endl;
	Size total(0);
	for ( Size ii = 1; ii <= memory_use_rotameric.size(); ++ii ) {
		total += memory_use_rotameric[ii];
		TR << chemical::name_from_aa(rotameric_amino_acids[ii]) << " with "
			<< memory_use_rotameric[ii] << " bytes" << std::endl;
	}

	for ( Size ii = 1; ii <= memory_use_semirotameric.size(); ++ii ) {
		total += memory_use_semirotameric[ii];
		TR << chemical::name_from_aa(sraa[ii]) << " with "
			<< memory_use_semirotameric[ii] << " bytes" << std::endl;
	}
	TR << "Total memory on Dunbrack Libraries: " << total << " bytes."
		<< std::endl;

}

void RotamerLibrary::create_fa_dunbrack_libraries_10_from_binary() {
	clock_t starttime = clock();

	std::string binary_filename = get_binary_name_10();
	TR << "Using Dunbrack library binary file '" << binary_filename << "'." << std::endl;
	utility::io::izstream binlib(binary_filename.c_str(),
		std::ios::in | std::ios::binary);

	if ( !binlib ) {
		utility_exit_with_message(
			"Could not open binary Dunbrack10 file from database -- how did this happen?");
	}

	/// READ PREABMLE
	/// Even if binary file should change its structure in the future, version number should always be
	/// the first piece of data in the file.
	boost::int32_t version(0);
	binlib.read((char*) &version, sizeof(boost::int32_t));

	/// Version 1: write out hard coded parameters for 08 library.  Verify when reading in
	/// the binary file that the parameters have not changed.
	/// DANGER DANGER DANGER
	/// If you change the binary format and adjust the preamble to the binary file,
	/// change the arithmetic in this function so that the correct number of bytes
	/// are jumped past.

	boost::int32_t n_rotameric_aas(0);
	binlib.read((char*) &n_rotameric_aas, sizeof(boost::int32_t));

	boost::int32_t n_semirotameric_aas(0);
	binlib.read((char*) &n_semirotameric_aas, sizeof(boost::int32_t));

	Size const n_rotameric_arrays_read_of_int32s = 2;
	Size const n_semirot_arrays_read_of_int32s = 5;
	Size const n_semirot_arrays_read_of_Reals = 1;
	Size const bytes_to_skip = (n_rotameric_arrays_read_of_int32s
		* n_rotameric_aas
		+ n_semirot_arrays_read_of_int32s * n_semirotameric_aas)
		* sizeof(boost::int32_t)
		+ n_semirot_arrays_read_of_Reals * n_semirotameric_aas
		* sizeof(Real);

	//binlib.seekg( bytes_to_skip, std::ios::cur ); //seekg not available for izstream...

	char * temp_buffer = new char[bytes_to_skip];
	binlib.read(temp_buffer, bytes_to_skip);
	delete[] temp_buffer;
	/// END PREABMLE

	read_from_binary(binlib);

	clock_t stoptime = clock();
	TR << "Dunbrack 2010 library took "
		<< ((double) stoptime - starttime) / CLOCKS_PER_SEC
		<< " seconds to load from binary" << std::endl;
}

void RotamerLibrary::write_binary_fa_dunbrack_libraries() const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[corrections::score::dun10] ) {
		return write_binary_fa_dunbrack_libraries_10();
	} else {
		return write_binary_fa_dunbrack_libraries_02();
	}
}

void RotamerLibrary::write_binary_fa_dunbrack_libraries_02() const {
#ifndef __native_client__
	std::string binary_filename = get_binary_name_02(/*for_writing=*/ true );
	if ( binary_filename.size() == 0 ) {
		TR << "Unable to open temporary file for writing the binary version of the Dunbrack02 library." << std::endl;
		return;
	}
	std::string tempfilename = random_tempname(binary_filename, "dun02_binary");

	TR << "Opening file " << tempfilename << " for output." << std::endl;

	utility::io::ozstream binlib(tempfilename.c_str(),
		std::ios::out | std::ios::binary);
	if ( binlib ) {

		/// WRITE PREABMLE
		/// Even if binary file should change its structure in the future, version number should always be
		/// the first piece of data in the file.
		boost::int32_t version(current_binary_format_version_id_02());
		binlib.write((char*) &version, sizeof(boost::int32_t));
		/// END PREABMLE

		write_to_binary(binlib);
		binlib.close();

		// Move the temporary file to its permanent location
		TR << "Moving temporary file " << tempfilename.c_str() << " to "
			<< binary_filename.c_str() << std::endl;
		rename(tempfilename.c_str(), binary_filename.c_str());
#ifndef WIN32
		chmod(binary_filename.c_str(), S_IROTH | S_IWUSR | S_IRUSR | S_IRGRP);
#endif
	} else {
		TR
			<< "Unable to open temporary file in rosetta database for writing the binary version of the Dunbrack02 library."
			<< std::endl;
	}
#endif
}

void RotamerLibrary::write_binary_fa_dunbrack_libraries_10() const {
#ifndef __native_client__
	std::string binary_filename = get_binary_name_10(/*for_writing=*/ true);
	if ( binary_filename.size() == 0 ) {
		TR << "Unable to open temporary file for writing the binary version of the Dunbrack10 library." << std::endl;
		return;
	}
	std::string tempfilename = random_tempname(binary_filename, "dun10_binary");

	TR << "Opening file " << tempfilename << " for output." << std::endl;

	utility::io::ozstream binlib(tempfilename.c_str(),
		std::ios::out | std::ios::binary);
	if ( binlib ) {

		/// WRITE PREABMLE
		/// Even if binary file should change its structure in the future, version number should always be
		/// the first piece of data in the file.
		boost::int32_t version(current_binary_format_version_id_10());
		binlib.write((char*) &version, sizeof(boost::int32_t));

		utility::vector1< chemical::AA > rotameric_amino_acids;
		utility::vector1< Size > rotameric_n_chi;
		utility::vector1< Size > rotameric_n_bb;

		utility::vector1< chemical::AA > sraa;
		utility::vector1< Size > srnchi;
		utility::vector1< Size > srnbb;
		utility::vector1< bool > scind;
		utility::vector1< bool > sampind;
		utility::vector1< bool > sym;
		utility::vector1< Real > astr;

		initialize_dun10_aa_parameters(
			rotameric_amino_acids, rotameric_n_chi, rotameric_n_bb,
			sraa, srnchi, srnbb, scind, sampind, sym, astr );

		boost::int32_t nrotameric(
			static_cast<boost::int32_t>(rotameric_amino_acids.size()));
		boost::int32_t nsemirotameric(static_cast<boost::int32_t>(sraa.size()));

		binlib.write((char*) &nrotameric, sizeof(boost::int32_t));
		binlib.write((char*) &nsemirotameric, sizeof(boost::int32_t));

		/// 2.
		boost::int32_t * rotaa_bin   = new boost::int32_t[ nrotameric ];
		boost::int32_t * rot_nchi_bin= new boost::int32_t[ nrotameric ];

		boost::int32_t * sraa_bin = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * srnchi_bin  = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * scind_bin   = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * sampind_bin = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * sym_bin  = new boost::int32_t[ nsemirotameric ];
		Real * astr_bin = new Real[ nsemirotameric ];

		for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
			rotaa_bin[ ii - 1 ] = static_cast< boost::int32_t > ( rotameric_amino_acids[ ii ] );
			rot_nchi_bin[ ii - 1 ] = static_cast< boost::int32_t > ( rotameric_n_chi[ ii ] );
		}
		for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
			sraa_bin[ ii - 1 ] = static_cast< boost::int32_t > ( sraa[ ii ] );
			srnchi_bin[ ii - 1 ]  = static_cast< boost::int32_t > ( srnchi[ ii ] );
			scind_bin[ ii - 1 ]   = static_cast< boost::int32_t > ( scind[ ii ] );
			sampind_bin[ ii - 1 ] = static_cast< boost::int32_t > ( sampind[ ii ] );
			sym_bin[ ii - 1 ]  = static_cast< boost::int32_t > ( sym[ ii ] );
			astr_bin[ ii - 1 ]  =  astr[ ii ];
		}

		binlib.write((char*) rotaa_bin, nrotameric * sizeof(boost::int32_t));
		binlib.write((char*) rot_nchi_bin, nrotameric * sizeof(boost::int32_t));

		binlib.write((char*) sraa_bin, nsemirotameric * sizeof(boost::int32_t));
		binlib.write((char*) srnchi_bin,
			nsemirotameric * sizeof(boost::int32_t));
		binlib.write((char*) scind_bin,
			nsemirotameric * sizeof(boost::int32_t));
		binlib.write((char*) sampind_bin,
			nsemirotameric * sizeof(boost::int32_t));
		binlib.write((char*) sym_bin, nsemirotameric * sizeof(boost::int32_t));
		binlib.write((char*) astr_bin, nsemirotameric * sizeof(Real));

		delete[] rotaa_bin;
		delete[] rot_nchi_bin;
		delete[] sraa_bin;
		delete[] srnchi_bin;
		delete[] scind_bin;
		delete[] sampind_bin;
		delete[] sym_bin;
		delete[] astr_bin;

		/// END PREABMLE

		write_to_binary(binlib);
		binlib.close();

		// Move the temporary file to its permanent location
		TR << "Moving temporary file to " << binary_filename.c_str()
			<< std::endl;
		rename(tempfilename.c_str(), binary_filename.c_str());
#ifndef WIN32
		chmod(binary_filename.c_str(), S_IROTH | S_IWUSR | S_IRUSR | S_IRGRP);
#endif
	} else {
		TR
			<< "Unable to open temporary file in rosetta database for writing the binary version of the Dunbrack '10 library."
			<< std::endl;
	}
#endif
}

std::string RotamerLibrary::random_tempname(std::string const & same_dir_as, std::string const & prefix) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string dirname = utility::file::FileName( same_dir_as ).path();
	std::string tempfilename = utility::file::create_temp_filename( dirname, prefix );

	TR << "Random tempname will be: " << tempfilename << std::endl;
	return tempfilename;
}

SingleResidueDunbrackLibraryOP
RotamerLibrary::create_rotameric_dunlib(
	chemical::AA aa,
	Size const n_chi,
	Size const n_bb,
	utility::io::izstream & library,
	bool dun02,
	std::string & next_aa_in_library,
	bool first_three_letter_code_already_read
) const
{
	SingleResidueDunbrackLibraryOP rotlib;
	// scope the case statements
	switch ( n_chi ) {
	case 1 : {
		switch ( n_bb ) {
		case 1 : {
			RotamericSingleResidueDunbrackLibrary< ONE, ONE > * r1 =
				new RotamericSingleResidueDunbrackLibrary< ONE, ONE >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 2 : {
			RotamericSingleResidueDunbrackLibrary< ONE, TWO > * r1 =
				new RotamericSingleResidueDunbrackLibrary< ONE, TWO >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 3 : {
			RotamericSingleResidueDunbrackLibrary< ONE, THREE > * r1 =
				new RotamericSingleResidueDunbrackLibrary< ONE, THREE >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 4 : {
			RotamericSingleResidueDunbrackLibrary< ONE, FOUR > * r1 =
				new RotamericSingleResidueDunbrackLibrary< ONE, FOUR >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 5 : {
			RotamericSingleResidueDunbrackLibrary< ONE, FIVE > * r1 =
				new RotamericSingleResidueDunbrackLibrary< ONE, FIVE >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(n_bb) );
		}
	} break;
	case 2 : {
		switch ( n_bb ) {
		case 1 : {
			RotamericSingleResidueDunbrackLibrary< TWO, ONE > * r2 =
				new RotamericSingleResidueDunbrackLibrary< TWO, ONE >( aa, dun02 );
			next_aa_in_library = r2->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r2);
			break;
		}
		case 2 : {
			RotamericSingleResidueDunbrackLibrary< TWO, TWO > * r2 =
				new RotamericSingleResidueDunbrackLibrary< TWO, TWO >( aa, dun02 );
			next_aa_in_library = r2->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r2);
			break;
		}
		case 3 : {
			RotamericSingleResidueDunbrackLibrary< TWO, THREE > * r2 =
				new RotamericSingleResidueDunbrackLibrary< TWO, THREE >( aa, dun02 );
			next_aa_in_library = r2->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r2);
			break;
		}
		case 4 : {
			RotamericSingleResidueDunbrackLibrary< TWO, FOUR > * r2 =
				new RotamericSingleResidueDunbrackLibrary< TWO, FOUR >( aa, dun02 );
			next_aa_in_library = r2->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r2);
			break;
		}
		case 5 : {
			RotamericSingleResidueDunbrackLibrary< TWO, FIVE > * r1 =
				new RotamericSingleResidueDunbrackLibrary< TWO, FIVE >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(n_bb) );
		}
	} break;
	case 3 : {
		switch ( n_bb ) {
		case 1 : {
			RotamericSingleResidueDunbrackLibrary< THREE, ONE > * r3 =
				new RotamericSingleResidueDunbrackLibrary< THREE, ONE >( aa, dun02 );
			next_aa_in_library = r3->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r3);
			break;
		}
		case 2 : {
			RotamericSingleResidueDunbrackLibrary< THREE, TWO > * r3 =
				new RotamericSingleResidueDunbrackLibrary< THREE, TWO >( aa, dun02 );
			next_aa_in_library = r3->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r3);
			break;
		}
		case 3 : {
			RotamericSingleResidueDunbrackLibrary< THREE, THREE > * r3 =
				new RotamericSingleResidueDunbrackLibrary< THREE, THREE >( aa, dun02 );
			next_aa_in_library = r3->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r3);
			break;
		}
		case 4 : {
			RotamericSingleResidueDunbrackLibrary< THREE, FOUR > * r3 =
				new RotamericSingleResidueDunbrackLibrary< THREE, FOUR >( aa, dun02 );
			next_aa_in_library = r3->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r3);
			break;
		}
		case 5 : {
			RotamericSingleResidueDunbrackLibrary< THREE, FIVE > * r1 =
				new RotamericSingleResidueDunbrackLibrary< THREE, FIVE >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(n_bb) );
		}
	} break;
	case 4 : {
		switch ( n_bb ) {
		case 1 : {
			RotamericSingleResidueDunbrackLibrary< FOUR, ONE > * r4 =
				new RotamericSingleResidueDunbrackLibrary< FOUR, ONE >( aa, dun02 );
			next_aa_in_library = r4->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r4);
			break;
		}
		case 2 : {
			RotamericSingleResidueDunbrackLibrary< FOUR, TWO > * r4 =
				new RotamericSingleResidueDunbrackLibrary< FOUR, TWO >( aa, dun02 );
			next_aa_in_library = r4->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r4);
			break;
		}
		case 3 : {
			RotamericSingleResidueDunbrackLibrary< FOUR, THREE > * r4 =
				new RotamericSingleResidueDunbrackLibrary< FOUR, THREE >( aa, dun02 );
			next_aa_in_library = r4->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r4);
			break;
		}
		case 4 : {
			RotamericSingleResidueDunbrackLibrary< FOUR, FOUR > * r4 =
				new RotamericSingleResidueDunbrackLibrary< FOUR, FOUR >( aa, dun02 );
			next_aa_in_library = r4->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r4);
			break;
		}
		case 5 : {
			RotamericSingleResidueDunbrackLibrary< FOUR, FIVE > * r1 =
				new RotamericSingleResidueDunbrackLibrary< FOUR, FIVE >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(n_bb) );
		}
	} break;
	default :
		utility_exit_with_message( "ERROR: too many chi angles desired for Dunbrack library: " +
			boost::lexical_cast<std::string>(n_chi) );
		break;
	}
	return rotlib;

}

SingleResidueDunbrackLibraryOP
RotamerLibrary::create_rotameric_dunlib(
	chemical::AA aa,
	Size const n_chi,
	bool dun02
) const
{
	return create_rotameric_dunlib( aa, n_chi, 2, dun02 );
}

SingleResidueDunbrackLibraryOP
RotamerLibrary::create_rotameric_dunlib(
	chemical::AA aa,
	Size const n_chi,
	Size const n_bb,
	bool dun02
) const
{
	SingleResidueDunbrackLibraryOP rotlib;
	switch ( n_chi ) {
	case 1 : {
		switch ( n_bb ) {
		case 1 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< ONE, ONE >( aa, dun02 ) );
			break;
		}
		case 2 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< ONE, TWO >( aa, dun02 ) );
			break;
		}
		case 3 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< ONE, THREE >( aa, dun02 ) );
			break;
		}
		case 4 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< ONE, FOUR >( aa, dun02 ) );
			break;
		}
		case 5 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< ONE, FIVE >( aa, dun02 ) );
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(n_bb) );
		}
	} break;
	case 2 : {
		switch ( n_bb ) {
		case 1 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< TWO, ONE >( aa, dun02 ) );
			break;
		}
		case 2 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< TWO, TWO >( aa, dun02 ) );
			break;
		}
		case 3 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< TWO, THREE >( aa, dun02 ) );
			break;
		}
		case 4 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< TWO, FOUR >( aa, dun02 ) );
			break;
		}
		case 5 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< TWO, FIVE >( aa, dun02 ) );
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(n_bb) );
		}
	} break;
	case 3 : {
		switch ( n_bb ) {
		case 1 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< THREE, ONE >( aa, dun02 ) );
			break;
		}
		case 2 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< THREE, TWO >( aa, dun02 ) );
			break;
		}
		case 3 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< THREE, THREE >( aa, dun02 ) );
			break;
		}
		case 4 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< THREE, FOUR >( aa, dun02 ) );
			break;
		}
		case 5 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< THREE, FIVE >( aa, dun02 ) );
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(n_bb) );
		}
	} break;
	case 4 : {
		switch ( n_bb ) {
		case 1 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< FOUR, ONE >( aa, dun02 ) );
			break;
		}
		case 2 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< FOUR, TWO >( aa, dun02 ) );
			break;
		}
		case 3 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< FOUR, THREE >( aa, dun02 ) );
			break;
		}
		case 4 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< FOUR, FOUR >( aa, dun02 ) );
			break;
		}
		case 5 : {
			rotlib = SingleResidueDunbrackLibraryOP( new RotamericSingleResidueDunbrackLibrary< FOUR, FIVE >( aa, dun02 ) );
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(n_bb) );
		}
	} break;
	default :
		utility_exit_with_message( "ERROR: too many chi angles desired for Dunbrack library: " +
			boost::lexical_cast<std::string>(n_chi) );
		break;
	}
	return rotlib;

}


SingleResidueDunbrackLibraryOP
RotamerLibrary::create_semi_rotameric_dunlib(
	chemical::AA aa,
	Size const nchi,
	Size const nbb,
	bool const use_bbind_rnchi_scoring,
	bool const use_bbind_rnchi_sampling,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle,
	utility::io::izstream & rotamer_definitions,
	utility::io::izstream & regular_library,
	utility::io::izstream & continuous_minimization_bbdep
	// phasing out bbind sampling -- utility::io::izstream & rnchi_bbind_probabilities
) const
{
	SingleResidueDunbrackLibraryOP rotlib;
	// scope the case statements
	switch ( nchi ) {
	case 1 : {
		switch ( nbb ) {
		case 1 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 2 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 3 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 4 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 5 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, FIVE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, FIVE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(nbb) );
		}
	} break;
	case 2 : {
		switch ( nbb ) {
		case 1 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 2 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 3 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 4 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 5 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, FIVE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, FIVE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(nbb) );
		}
	} break;

	default :
		utility_exit_with_message(
			"ERROR: too many chi angles desired for semi-rotameric Dunbrack library: "
			+ boost::lexical_cast<std::string>(nchi));
		break;
	}
	return rotlib;

}

SingleResidueDunbrackLibraryOP
RotamerLibrary::create_semi_rotameric_dunlib(
	chemical::AA aa,
	Size const nchi,
	bool const use_bbind_rnchi_scoring,
	bool const use_bbind_rnchi_sampling,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle
) const
{
	return create_semi_rotameric_dunlib( aa,
		nchi,
		2,
		use_bbind_rnchi_scoring,
		use_bbind_rnchi_sampling,
		nrchi_is_symmetric,
		nrchi_start_angle);
}


SingleResidueDunbrackLibraryOP
RotamerLibrary::create_semi_rotameric_dunlib(
	chemical::AA aa,
	Size const nchi,
	Size const nbb,
	bool const use_bbind_rnchi_scoring,
	bool const use_bbind_rnchi_sampling,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle
) const
{
	SingleResidueDunbrackLibraryOP rotlib;
	switch ( nchi ) {
	case 1 : {
		switch ( nbb ) {
		case 1 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 2 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 3 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 4 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 5 : {
			SemiRotamericSingleResidueDunbrackLibrary< ONE, FIVE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE, FIVE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(nbb) );
		}
	} break;
	case 2 : {
		switch ( nbb ) {
		case 1 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 2 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 3 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 4 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		case 5 : {
			SemiRotamericSingleResidueDunbrackLibrary< TWO, FIVE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO, FIVE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = SingleResidueDunbrackLibraryOP(r1);
			break;
		}
		default :
			utility_exit_with_message( "ERROR: too many bb angles desired for Dunbrack library: " +
				boost::lexical_cast<std::string>(nbb) );
		}
	} break;

	default :
		utility_exit_with_message(
			"ERROR: too many chi angles desired for semi-rotameric Dunbrack library: "
			+ boost::lexical_cast<std::string>(nchi));
		break;
	}
	return rotlib;

}

/// @details Build and initialize some of the data in an SRDL in preparation for
/// binary file input -- that is, the SRDL does not have "read_from_file(s)" called
/// before it is returned to the calling function.  ATM not threadsafe.
SingleResidueDunbrackLibraryOP RotamerLibrary::create_srdl(
	chemical::AA aa_in) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool initialized(false);

	static utility::vector1< chemical::AA > rotameric_amino_acids;
	static utility::vector1< Size > rotameric_n_chi;
	static utility::vector1< Size > rotameric_n_bb;

	static utility::vector1< chemical::AA > sraa;
	static utility::vector1< Size > srnchi;
	static utility::vector1< Size > srnbb;
	static utility::vector1< bool > scind;
	static utility::vector1< bool > sampind;
	static utility::vector1< bool > sym;
	static utility::vector1< Real > astr;

	if ( !initialized ) {
		/// put a mutex here...
		if ( option[ corrections::score::dun10 ] ) {
			initialize_dun10_aa_parameters(
				rotameric_amino_acids, rotameric_n_chi, rotameric_n_bb,
				sraa, srnchi, srnbb, scind, sampind, sym, astr );
		} else {
			initialize_dun02_aa_parameters(
				rotameric_amino_acids, rotameric_n_chi, rotameric_n_bb );
		}
		initialized = true;
	}

	Size find_rot_aa = 1;
	while ( find_rot_aa <= rotameric_amino_acids.size() ) {
		if ( rotameric_amino_acids[ find_rot_aa ] == aa_in ) {
			return create_rotameric_dunlib( aa_in, rotameric_n_chi[ find_rot_aa ], rotameric_n_bb[ find_rot_aa ], ! (option[ corrections::score::dun10 ] )  );
		}
		++find_rot_aa;
	}

	Size find_semirot_aa = 1;
	while ( find_semirot_aa <= sraa.size() ) {
		if ( sraa[find_semirot_aa] == aa_in ) {
			Size i = find_semirot_aa;
			return create_semi_rotameric_dunlib( aa_in, srnchi[i], srnbb[i], scind[i], sampind[i], sym[i], astr[i] );
		}
		++find_semirot_aa;
	}
	return 0;
}


/// @details Hard code this data in a single place; adjust the booleans for the backbone
/// independent
void
RotamerLibrary::initialize_dun10_aa_parameters(
	utility::vector1< chemical::AA > & rotameric_amino_acids,
	utility::vector1< Size > & rotameric_n_chi,
	utility::vector1< Size > & rotameric_n_bb,
	utility::vector1< chemical::AA > & sraa,
	utility::vector1< Size > & srnchi,
	utility::vector1< Size > & srnbb,
	utility::vector1< bool > & scind,
	utility::vector1< bool > & sampind,
	utility::vector1< bool > & sym,
	utility::vector1< Real > & astr
)
{
	using namespace chemical;

	/// Rotameric specifics
	rotameric_amino_acids.reserve( 10 );   rotameric_n_chi.reserve( 10 );  rotameric_n_bb.reserve( 10 );
	rotameric_amino_acids.push_back( aa_cys );   rotameric_n_chi.push_back( 1 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_ile );   rotameric_n_chi.push_back( 2 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_lys );   rotameric_n_chi.push_back( 4 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_leu );   rotameric_n_chi.push_back( 2 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_met );   rotameric_n_chi.push_back( 3 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_pro );   rotameric_n_chi.push_back( 3 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_arg );   rotameric_n_chi.push_back( 4 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_ser );   rotameric_n_chi.push_back( 1 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_thr );   rotameric_n_chi.push_back( 1 ); rotameric_n_bb.push_back( 2 );
	rotameric_amino_acids.push_back( aa_val );   rotameric_n_chi.push_back( 1 ); rotameric_n_bb.push_back( 2 );


	/// Semi rotameric specifics
	Size const n_semi_rotameric_aas = 8;
	sraa.resize( n_semi_rotameric_aas );
	srnchi.resize( n_semi_rotameric_aas );
	srnbb.resize( n_semi_rotameric_aas );
	scind.resize( n_semi_rotameric_aas );
	sampind.resize( n_semi_rotameric_aas );
	sym.resize( n_semi_rotameric_aas );
	astr.resize( n_semi_rotameric_aas );


	Size i = 0;
	sraa[++i] = aa_asp; srnchi[i] = 1; srnbb[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=1; astr[i] =  -90;
	sraa[++i] = aa_glu; srnchi[i] = 2; srnbb[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=1; astr[i] =  -90;
	sraa[++i] = aa_phe; srnchi[i] = 1; srnbb[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=1; astr[i] =  -30;
	sraa[++i] = aa_his; srnchi[i] = 1; srnbb[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=0; astr[i] = -180;
	sraa[++i] = aa_asn; srnchi[i] = 1; srnbb[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=0; astr[i] = -180;
	sraa[++i] = aa_gln; srnchi[i] = 2; srnbb[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=0; astr[i] = -180;
	sraa[++i] = aa_trp; srnchi[i] = 1; srnbb[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=0; astr[i] = -180;
	sraa[++i] = aa_tyr; srnchi[i] = 1; srnbb[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=1; astr[i] =  -30;

}

void
RotamerLibrary::initialize_dun02_aa_parameters(
	utility::vector1< chemical::AA > & rotameric_amino_acids,
	utility::vector1< Size > & rotameric_n_chi,
	utility::vector1< Size > & rotameric_n_bb
)
{
	using namespace chemical;

	utility::vector1< Size > nchi_for_aa;
	utility::vector1< Size > nbb_for_aa;
	Size n_rotameric_aas;
	initialize_dun02_aa_parameters( nchi_for_aa, nbb_for_aa, n_rotameric_aas );
	rotameric_amino_acids.resize( 0 ); rotameric_amino_acids.reserve( n_rotameric_aas );
	rotameric_n_chi.resize( 0 );    rotameric_n_chi.reserve( n_rotameric_aas );
	rotameric_n_bb.resize( 0 );     rotameric_n_bb.reserve( n_rotameric_aas );

	for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {
		if ( nchi_for_aa[ ii ] != 0 ) {
			rotameric_amino_acids.push_back( AA(ii) );
			rotameric_n_chi.push_back( nchi_for_aa[ ii ] );
			rotameric_n_bb.push_back( nbb_for_aa[ ii ] );
		}
	}

}

void
RotamerLibrary::initialize_dun02_aa_parameters(
	utility::vector1< Size > & nchi_for_aa,
	utility::vector1< Size > & nbb_for_aa,
	Size & n_rotameric_aas
)
{
	using namespace chemical;
	nchi_for_aa.resize( num_canonical_aas );
	nbb_for_aa.resize( num_canonical_aas );
	std::fill( nchi_for_aa.begin(), nchi_for_aa.end(), 0 ); // overkill! ala, gly = 0...
	std::fill( nbb_for_aa.begin(), nbb_for_aa.end(), 0 ); // overkill! ala, gly = 0...
	n_rotameric_aas = 0;

	nchi_for_aa[ aa_cys ] = 1; nbb_for_aa[ aa_cys ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_asp ] = 2; nbb_for_aa[ aa_asp ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_glu ] = 3; nbb_for_aa[ aa_glu ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_phe ] = 2; nbb_for_aa[ aa_phe ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_his ] = 2; nbb_for_aa[ aa_his ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_ile ] = 2; nbb_for_aa[ aa_ile ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_lys ] = 4; nbb_for_aa[ aa_lys ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_leu ] = 2; nbb_for_aa[ aa_leu ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_met ] = 3; nbb_for_aa[ aa_met ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_asn ] = 2; nbb_for_aa[ aa_asn ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_pro ] = 3; nbb_for_aa[ aa_pro ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_gln ] = 3; nbb_for_aa[ aa_gln ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_arg ] = 4; nbb_for_aa[ aa_arg ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_ser ] = 1; nbb_for_aa[ aa_ser ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_thr ] = 1; nbb_for_aa[ aa_thr ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_val ] = 1; nbb_for_aa[ aa_val ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_trp ] = 2; nbb_for_aa[ aa_trp ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_tyr ] = 2; nbb_for_aa[ aa_tyr ] = 2; ++n_rotameric_aas;

}

/// @brief Generic "write out a library to binary" which has a straightforward
/// nlibraries (aa, aa_data)* format.  Works for both 02 and 08 libraries.  The
/// kind of aa_data written is left to the discression of the SinResDunLib
/// derived classes.
void RotamerLibrary::write_to_binary(utility::io::ozstream & out) const {

	Size count_libraries( 0 );
	for ( core::Size ii(1); ii <= aa_libraries_.size(); ++ii ) {
		SingleResidueDunbrackLibraryCOP srdl =
			utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const > ( aa_libraries_[ ii ] );

		if ( srdl ) ++count_libraries; /// write out only the dunbrack libraries.
	}

	/// 1. How many libraries?
	boost::int32_t const nlibraries = count_libraries;
	out.write((char*) &nlibraries, sizeof(boost::int32_t));

	for ( core::Size ii(1); ii <= aa_libraries_.size(); ++ii ) {
		SingleResidueDunbrackLibraryCOP srdl =
			utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const > ( aa_libraries_[ii] );

		if ( !srdl ) {
			continue; /// write out only the dunbrack libraries.
		}

		/// 2. Amino acid type of next library
		boost::int32_t const which_aa( srdl->aa() );
		out.write( (char*) &which_aa, sizeof( boost::int32_t ) );

		/// 3. Data for this amino acid type.
		srdl->write_to_binary(out);
	}

}

/// @details Generic binary file reader.  Intended to work for both 02 and 08 Libraries.
/// If a "preamble" starts the binary format, the input stream to this function must have already
/// been walked past it.  Somewhat circularly, the preamble is anything in the binary file that
/// leads the "nlibraries (aa, aa_data)*" format perscribed in this function.
void RotamerLibrary::read_from_binary(utility::io::izstream & in) {

	using namespace boost;

	/// 1. How many libraries?
	boost::int32_t nlibraries(0);
	in.read((char*) &nlibraries, sizeof(boost::int32_t));
	for ( Size ii = 1; ii <= Size(nlibraries); ++ii ) {

		/// 2.  Which amino acid is next in the binary file?
		boost::int32_t which_aa32(chemical::aa_unk);
		in.read((char*) &which_aa32, sizeof(boost::int32_t));
		AA which_aa(static_cast<AA>(which_aa32));

		//TR << "amw reading binlib for " << which_aa << std::endl;

		/// 3. Read the data associated with that amino acid.
		SingleResidueDunbrackLibraryOP single_lib = create_srdl( which_aa );
		if ( !single_lib ) {
			utility_exit_with_message( "Error reading single residue rotamer library for " + name_from_aa( which_aa ) );
		}
		single_lib->read_from_binary(in);
		add_residue_library(which_aa, single_lib);
	}
}

/// @brief Reload the Dunbrack Rotamer libraries from ASCII, and make sure that they match the ones loaded from binary.
/// NOTE WELL: This is *not* a const function, as reloading from ASCII modifies internals.
bool
RotamerLibrary::validate_dunbrack_binary() {

#if defined MULTI_THREADED
	TR.Error << "RotamerLibrary::validate_dunbrack_binary() is HIGHLY un-threadsafe. Don't call it from a multithreaded context. " << std::endl;
	return true; // Just return success
#endif

	std::string const binary_name( get_binary_name() ); // Handles Dun10/Dun02 decision internally
	// If the existing rotamer library wasn't loaded from binary, we don't need to validate it.
	if ( ! decide_read_from_binary() ) {
		TR.Warning << "Existing Dunbrack binary file '" << binary_name << "' either doesn't exist, or is known to be invalid. Skipping reloading check." << std::endl;
		return true; // Because any reads against it will get the valid ASCII-parsed version
	}

	// We need to offload the existing rotamer libraries, and then reload them from ASCII
	// Entries in libraries_ and aa_libraries_ will be over-written on reload
	// Need to make a copy of aa_libraries, and then zero out the originals so they can be reloaded
	utility::vector1 < SingleResidueRotamerLibraryCOP > aa_libraries_copy( aa_libraries_ );
	for ( core::Size aa( 1 ); aa <= core::chemical::num_canonical_aas; ++aa ) {
		aa_libraries_[ aa ] = 0;
	}

	TR << "Comparing ASCII version of Dunbrack library to binary version at " << binary_name << std::endl;

	// Now, we reload the libraries from ASCII specifically
	create_fa_dunbrack_libraries_from_ASCII(); // Handles Dun10/Dun02 decision internally

	// And compare
	if ( aa_libraries_copy.size() != aa_libraries_.size() ) {
		TR.Error << "AA Library sizes don't match! *SEVERE* logic error!" << std::endl;
		return false;
	}
	bool valid( true );
	for ( core::Size ii( 1 ); ii <= aa_libraries_copy.size(); ++ii ) {
		if ( aa_libraries_copy[ii] == 0 && aa_libraries_[ii] == 0 ) {
			// Both null? We're fine.
			TR << "Rotamer library for " << core::chemical::name_from_aa( core::chemical::AA( ii ) ) << " NULL in both ASCII and binary." << std::endl;
		} else if ( aa_libraries_copy[ii] == 0 ) {
			TR.Error << "Rotamer library for " << core::chemical::name_from_aa( core::chemical::AA( ii ) ) << " NULL in binary, but present in ASCII." << std::endl;
			valid = false;
		} else if ( aa_libraries_[ii] == 0 ) {
			TR.Error << "Rotamer library for " << core::chemical::name_from_aa( core::chemical::AA( ii ) ) << " NULL in ASCII, but present in binary." << std::endl;
			valid = false;
		} else if ( *(aa_libraries_[ii]) == *(aa_libraries_copy[ii]) ) { // ASCII first (as reference), binary second
			// Match. Good!
			TR << "Rotamer library for " << core::chemical::name_from_aa( core::chemical::AA( ii ) ) << " matches in both ASCII and binary." << std::endl;
		} else {
			// No match. BAD!
			TR.Error << "Rotamer library for " << core::chemical::name_from_aa( core::chemical::AA( ii ) ) << " DOES NOT MATCH for ASCII and binary." << std::endl;
			valid = false;
		}
	}
	return valid;
}

} // dunbrack
} // namespace scoring
} // namespace core
