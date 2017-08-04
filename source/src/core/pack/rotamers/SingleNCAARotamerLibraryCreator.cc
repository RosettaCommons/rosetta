// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleNCAARotamerLibraryCreator.hh
/// @brief  Class for instantiating a particular SingleResidueRotamerLibrary
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <core/pack/rotamers/SingleNCAARotamerLibraryCreator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Program header
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//External headers
#include <boost/lexical_cast.hpp>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace rotamers {

static THREAD_LOCAL basic::Tracer TR("core.pack.rotamers.SingleNCAARotamerLibraryCreator");

core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleNCAARotamerLibraryCreator::create( core::chemical::ResidueType const & restype) const {
	using namespace core::chemical::rotamers;
	using namespace core::pack::dunbrack;

	RotamerLibrarySpecificationCOP libspec( restype.rotamer_library_specification() );
	// If the factory system is sound, these two checks should work.
	debug_assert( libspec );
	NCAARotamerLibrarySpecificationCOP ncaa_libspec = utility::pointer::dynamic_pointer_cast< NCAARotamerLibrarySpecification const >(libspec);
	debug_assert( ncaa_libspec );

	chemical::AA aan( chemical::aa_unk );
	// dun02 if rotameric, dun10 densities if semirotameric
	bool dun02( !ncaa_libspec->semirotameric_ncaa_rotlib() );

	// Some basic error checking against restype
	Size pose_n_rotlib_chi( restype.nchi() - restype.n_proton_chi() - ( dun02 ? 0 : 1 ) );
	Size n_rotlib_chi( ncaa_libspec->ncaa_rotlib_n_rotameric_bins() );
	if ( n_rotlib_chi != pose_n_rotlib_chi ) {
		TR.Error << "Number of chi mismatch. Expected " << n_rotlib_chi << " rotatable heavy atom chis, found " << pose_n_rotlib_chi << std::endl;
		utility_exit_with_message("Number of chi mismatch in NCAA rotlib loading.");
	}
	if ( restype.aa() != aan ) {
		TR.Warning << "AA designation " << restype.aa() << " for NCAA rotamer library will not be obseved." << std::endl;
	}

	// AMW: Long term, the rotamer library should specify how many torsions
	// it depends on, and you should be able to tell your protocol how many
	// (and which) dependencies you want. For example, you should be able to
	// switch between alpha AA rotamers and specialized beta rotamers for the
	// betas!
	// AMW: Presently, we are compatible with THREE backbone torsion beta rotlibs.
	Size n_rotlib_bb( restype.mainchain_atoms().size() - 1 ); // 2 for alpha 3 for beta

	// AMW: This code sends a hideous shiver down my spine and it should do the same for yours.
	// The last passage was a regrettable tragedy--where this function had to figure something out about
	// the residue type besides whether its rotamer library was already known.
	// This passage requires that we even check for variant types. What's next???
	if ( restype.has_variant_type( chemical::ACETYLATED_NTERMINUS_VARIANT ) ) {
		--n_rotlib_bb;
	}

	using namespace utility::options;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// trust the params file to correctly name the rotamer library even if it's densities format
	std::string file_name = ncaa_libspec->ncaa_rotlib_path();
	std::string dir_name, full_path;
	bool tried(false);
	if ( ! file_name.size() ) {
		utility_exit_with_message("Unspecified NCAA rotlib path with residue type " + restype.name() );
	}

	// try priority path if available
	utility::options::PathVectorOption & priorityvec = option[ in::file::override_rot_lib_path ];
	Size pveci(1);
	utility::io::izstream rotlib_in;
	while ( pveci <= priorityvec.size() ) {
		tried=true;
		dir_name = priorityvec[ pveci ].name();
		full_path = dir_name + file_name;
		rotlib_in.open( full_path );

		if ( rotlib_in.good() ) break;

		// Try flat hierarchy too
		full_path = dir_name + file_name.substr( file_name.find_last_of( "/\\" )+1 );
		rotlib_in.open( full_path );
		if ( rotlib_in.good() ) break;

		pveci++;
	}

	if ( !rotlib_in.good() || !tried ) {
		// if the priority paths don't yield a hit, try the relative file name.
		tried=true;
		dir_name = "";
		full_path = file_name;
		rotlib_in.open( full_path );
	}

	if ( !rotlib_in.good() || !tried ) {
		tried=true;
		// if that doesn't work, try the relative file name in the database.
		dir_name = basic::database::full_name( "" );
		full_path = dir_name + file_name;
		rotlib_in.open( full_path );
	}

	if ( !rotlib_in.good() || !tried ) {
		tried=true;
		// next try the ncaa_rotlibs directory in the database.
		dir_name = basic::database::full_name( "/rotamer/ncaa_rotlibs/" );
		full_path = dir_name + file_name;
		rotlib_in.open( full_path );
	}

	if ( !rotlib_in.good() || !tried ) {
		// if we cannot open in regular path, try alternate paths
		utility::options::PathVectorOption & pvec = option[ in::file::extra_rot_lib_path ];
		pveci = 1;
		while ( pveci <= pvec.size() ) {
			// tried=true; // last time around - don't set as it will never be used.
			dir_name = pvec[ pveci ].name();
			full_path = dir_name + file_name;
			rotlib_in.open( full_path );

			if ( rotlib_in.good() ) break;

			// Try flat hierarchy too
			full_path = dir_name + file_name.substr( file_name.find_last_of( "/\\" )+1 );
			rotlib_in.open( full_path );
			if ( rotlib_in.good() ) break;

			pveci++;
		}
	}

	runtime_assert_string_msg( rotlib_in.good(), "Error!  Could not open file " + file_name + " for read." );

	TR << "Reading in rot lib " << full_path << "..." << std::endl;

	// get an instance of RotamericSingleResidueDunbrackLibrary, but need a RotamerLibrary to do it
	// this means that when ever you read in the NCAA libraries you will also read in the Dunbrack libraries
	// this may need to be a pointer to the full type and not just a SRRLOP
	// AMW: this is never used
	//std::string empty_string("");

	// this comes almost directally from RotmerLibrary.cc::create_rotameric_dunlib()
	SingleResidueRotamerLibraryOP ncaa_rotlib;

	bool use_shapovalov = basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable] && basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_dun10_enable];
	bool use_bicubic = basic::options::option[ basic::options::OptionKeys::corrections::score::use_bicubic_interpolation ]();
	bool entropy_correction = basic::options::option[ basic::options::OptionKeys::corrections::score::dun_entropy_correction ]();
	core::Real prob_buried = basic::options::option[ basic::options::OptionKeys::packing::dunbrack_prob_buried ]();
	core::Real prob_nonburied = basic::options::option[ basic::options::OptionKeys::packing::dunbrack_prob_nonburied ]();
	if ( !dun02 ) {
		// semirotameric

		// amw: I suppose this is sort of gross, but the residue type now reports whether its nrchi is symmetric
		// and where its angle starts
		// come to think of it why do we have a symmetric nrchi for C95? we shouldn't. TODO: make the right library.
		bool nrchi_is_symmetric = ncaa_libspec->nrchi_symmetric();
		Real nrchi_start_angle = ncaa_libspec->nrchi_start_angle();


		std::string def_suffix = "_definitions.rotlib";
		std::string dens_suffix = "_densities.rotlib";

		std::string full_def = dir_name + restype.name3() + def_suffix;
		std::string full_dens = dir_name + restype.name3() + dens_suffix;
		utility::io::izstream defstream( full_def );
		utility::io::izstream densstream( full_dens );

		// pretty sure we only have 1 and 2 rotameric chi implemented
		// because in the canonicals that's how many semirotameric rotlib chi we have
		// so stick to that for now...
		switch ( n_rotlib_chi ) {
		case 1 : {
			switch ( n_rotlib_bb ) {
			case 1 : {
				SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 2 : {
				SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 3 : {
				SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 4 : {
				SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 5 : {
				SemiRotamericSingleResidueDunbrackLibrary< ONE, FIVE > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< ONE, FIVE >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			default :
				utility_exit_with_message( "ERROR: too many bb angles desired for NCAA library: " +
					boost::lexical_cast<std::string>(n_rotlib_bb) );
				break;
			}

		} break;
		case 2 : {
			switch ( n_rotlib_bb ) {
			case 1 : {
				SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 2 : {
				SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 3 : {
				SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 4 : {
				SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 5 : {
				SemiRotamericSingleResidueDunbrackLibrary< TWO, FIVE > * r1 =
					new SemiRotamericSingleResidueDunbrackLibrary< TWO, FIVE >( aan, false, false, use_shapovalov, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle, defstream, rotlib_in, densstream );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			default :
				utility_exit_with_message( "ERROR: too many bb angles desired for semirotameric NCAA library: " +
					boost::lexical_cast<std::string>(n_rotlib_bb) );
				break;
			}
		} break;
		default :
			utility_exit_with_message( "ERROR: too many chi angles desired for semirotameric NCAA library: " +
				boost::lexical_cast<std::string>(n_rotlib_chi) );
			break;
		}
	} else {
		// dun02 style
		switch ( n_rotlib_chi ) {
		case 1 : {
			switch ( n_rotlib_bb ) {
			case 1 : {
				RotamericSingleResidueDunbrackLibrary< ONE, ONE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< ONE, ONE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 2 : {
				RotamericSingleResidueDunbrackLibrary< ONE, TWO > * r1 =
					new RotamericSingleResidueDunbrackLibrary< ONE, TWO >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 3 : {
				RotamericSingleResidueDunbrackLibrary< ONE, THREE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< ONE, THREE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 4 : {
				RotamericSingleResidueDunbrackLibrary< ONE, FOUR > * r1 =
					new RotamericSingleResidueDunbrackLibrary< ONE, FOUR >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 5 : {
				RotamericSingleResidueDunbrackLibrary< ONE, FIVE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< ONE, FIVE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			default :
				utility_exit_with_message( "ERROR: too many bb angles desired for NCAA library: " +
					boost::lexical_cast<std::string>(n_rotlib_bb) );
				break;
			}

		} break;
		case 2 : {
			switch ( n_rotlib_bb ) {
			case 1 : {
				RotamericSingleResidueDunbrackLibrary< TWO, ONE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< TWO, ONE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 2 : {
				RotamericSingleResidueDunbrackLibrary< TWO, TWO > * r1 =
					new RotamericSingleResidueDunbrackLibrary< TWO, TWO >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 3 : {
				RotamericSingleResidueDunbrackLibrary< TWO, THREE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< TWO, THREE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 4 : {
				RotamericSingleResidueDunbrackLibrary< TWO, FOUR > * r1 =
					new RotamericSingleResidueDunbrackLibrary< TWO, FOUR >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 5 : {
				RotamericSingleResidueDunbrackLibrary< TWO, FIVE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< TWO, FIVE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			default :
				utility_exit_with_message( "ERROR: too many bb angles desired for NCAA library: " +
					boost::lexical_cast<std::string>(n_rotlib_bb) );
				break;
			}
		} break;
		case 3 : {
			switch ( n_rotlib_bb ) {
			case 1 : {
				RotamericSingleResidueDunbrackLibrary< THREE, ONE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< THREE, ONE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 2 : {
				RotamericSingleResidueDunbrackLibrary< THREE, TWO > * r1 =
					new RotamericSingleResidueDunbrackLibrary< THREE, TWO >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 3 : {
				RotamericSingleResidueDunbrackLibrary< THREE, THREE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< THREE, THREE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 4 : {
				RotamericSingleResidueDunbrackLibrary< THREE, FOUR > * r1 =
					new RotamericSingleResidueDunbrackLibrary< THREE, FOUR >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 5 : {
				RotamericSingleResidueDunbrackLibrary< THREE, FIVE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< THREE, FIVE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			default :
				utility_exit_with_message( "ERROR: too many bb angles desired for NCAA library: " +
					boost::lexical_cast<std::string>(n_rotlib_bb) );
				break;
			};
		} break;
		case 4 : {
			switch ( n_rotlib_bb ) {
			case 1 : {
				RotamericSingleResidueDunbrackLibrary< FOUR, ONE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< FOUR, ONE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 2 : {
				RotamericSingleResidueDunbrackLibrary< FOUR, TWO > * r1 =
					new RotamericSingleResidueDunbrackLibrary< FOUR, TWO >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 3 : {
				RotamericSingleResidueDunbrackLibrary< FOUR, THREE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< FOUR, THREE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 4 : {
				RotamericSingleResidueDunbrackLibrary< FOUR, FOUR > * r1 =
					new RotamericSingleResidueDunbrackLibrary< FOUR, FOUR >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			case 5 : {
				RotamericSingleResidueDunbrackLibrary< FOUR, FIVE > * r1 =
					new RotamericSingleResidueDunbrackLibrary< FOUR, FIVE >( aan, dun02, use_bicubic, entropy_correction, prob_buried, prob_nonburied );
				r1->set_n_chi_bins( ncaa_libspec->ncaa_rotlib_n_bin_per_rot() );
				r1->read_from_file( rotlib_in, false );
				ncaa_rotlib = SingleResidueRotamerLibraryOP(r1);
				break;
			}
			default :
				utility_exit_with_message( "ERROR: too many bb angles desired for NCAA library: " +
					boost::lexical_cast<std::string>(n_rotlib_bb) );
				break;
			}
		} break;
		default :
			utility_exit_with_message( "ERROR: too many chi angles desired for NCAA library: " +
				boost::lexical_cast<std::string>(n_rotlib_chi) );
			break;
		}
	}

	return ncaa_rotlib;
}

std::string
SingleNCAARotamerLibraryCreator::keyname() const {
	return core::chemical::rotamers::NCAARotamerLibrarySpecification::library_name();
}

} //namespace rotamers
} //namespace pack
} //namespace core

