// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/cyclic_peptide/count_cycpep_sequences.cc
/// @brief An application to count the number of unique sequences for a cyclic peptide given
/// its symmetry.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
/// @author Math from Todd Yeates (yeates@mbi.ucla.edu).

// devel headers
#include <devel/init.hh>

// core headers
#include <core/types.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

// Boost headers
#include <boost/math/common_factor_rt.hpp>

// C++ headers
#include <cmath> //For std::pow

static basic::Tracer TR( "apps.public.cyclic_peptide.count_cycpep_sequences" );
static const std::string errmsg( "Error in count_cycpep_sequence application:  " );

#define COUNT_CYCPEP_SEQUENCE_INTTYPE signed long long

OPT_KEY( Boolean, mirror_symmetry )
OPT_KEY( Boolean, do_numerical_count )
OPT_KEY( Integer, symmetry_number )
OPT_KEY( Integer, peptide_length )
OPT_KEY( Integer, options_per_position )
OPT_KEY( Integer, achiral_options )
OPT_KEY( Boolean, write_out_sequences )
OPT_KEY( Boolean, SN_semi_analytical )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( symmetry_number, "The fold symmetry.  Set to 3 for C3 symmetry, for example.  The default, 1, is for the asymmetric case.", 1 );
	NEW_OPT( mirror_symmetry, "If true, the symmetry type is SN (S2, S4, S6, etc.).  If false, we consider simple cyclic symmetry (C2, C3, C4, etc.).  False by default.", false );
	NEW_OPT( do_numerical_count, "If true, we try to do the enumeration manually, rather than by using the expressions derived from Burnside's lemma.  Note that this scales poorly for longer peptides.  False by default.", false );
	NEW_OPT( peptide_length, "The number of residues in the peptide.  Defaults arbitrarily to 8.", 8 );
	NEW_OPT( options_per_position, "The number of options allowed at each position.  For example, if we were considering ABEGO bins, there would be five options; if we were considering the 20 canonical amino acids, there would be 20 options.  Default 20.", 20 );
	NEW_OPT( achiral_options, "The number of options that are their own mirror image (e.g. GLY, AIB), if we're considering mirror symmetry.  Note that the total number of options per position minus the number of asymmetric options must be an even number, if symmetry is enabled.  Default 0.", 0 );
	NEW_OPT( write_out_sequences, "If true, sequences enumerated numerically will be written out.  False by default, since this can produce astronomical numbers of sequences.", false );
	NEW_OPT( SN_semi_analytical, "If true, we use a semi-analytical method for SN (improper rotational) symmetry.  If false, we use a fully analytical method.  False by default.", false );
}

/// @brief Check that the user has set options sensibly.
void
do_options_checks() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	runtime_assert_string_msg( option[symmetry_number]() > 0, errmsg + "The -symmetry_number flag must be used to specifiy a natural number (i.e. greater than zero)." );
	runtime_assert_string_msg( option[peptide_length]() > 0, errmsg + "The peptide length (-peptide_length option) must be greater than zero.");
	runtime_assert_string_msg( option[peptide_length]() % option[symmetry_number]() == 0, errmsg + "The peptide length (-peptide_length option) must be an integer multiple of the fold symmetry (-symmetry_number option).");
	runtime_assert_string_msg( option[options_per_position]() > 0, errmsg + "There must be at least one option per position (-options_per_position flag)." );
	if ( option[mirror_symmetry].value() ) {
		runtime_assert_string_msg( option[symmetry_number]() % 2 == 0, errmsg + "If the -mirror_symmetry option is used, the fold symmetry (-symmetry_number option) must be even.");
		runtime_assert_string_msg( option[achiral_options]() >= 0, errmsg + "The number of achiral options must be greater than or equal to zero (-achiral_options flag)." );
		runtime_assert_string_msg( option[options_per_position]() >= option[achiral_options](), errmsg + "The total number of options per position (-options_per_position flag) must be greater than or equal to the number of achiral options (-achiral_options flag)." );
		runtime_assert_string_msg( (option[options_per_position]() - option[achiral_options]()) % 2 == 0, errmsg + "The number of options per position minus the number of achiral options must be an even number." );
	} else {
		runtime_assert_string_msg( !option[SN_semi_analytical].user(), errmsg + "The -SN_semi_analytical option was specified without the -mirror_symmetry option." );
		runtime_assert_string_msg( !option[achiral_options].user(), errmsg + "The -achiral_options option was used without the -mirror_symmetry option.");
	}
	if ( option[write_out_sequences].value() ) {
		runtime_assert_string_msg( option[do_numerical_count].value(), errmsg + "The -write_out_sequences option requires the -do_numerical_count option." );
	}

	if ( option[do_numerical_count].value() && ( option[peptide_length]() > 15 ) ) {
		TR.Warning << "Warning!  The -do_numerical_count option was used with a large peptide.  This could result in very slow performance!" << std::endl;
	}
}

/// @brief Circularly permute a vector by one.
/// @details Assumes vector is more than zero units long.
void
do_cyclic_permutation(
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > & vect
) {
	COUNT_CYCPEP_SEQUENCE_INTTYPE const first_element( vect[1] );
	for ( core::Size i(1); i<vect.size(); ++i ) {
		vect[i] = vect[i+1];
	}
	vect[vect.size()] = first_element;
}

/// @brief Circularly permute a vector until the least vector is found.
utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE >
find_min_permutation(
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > const &input_vector
) {
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > output_vector( input_vector ), curvector( input_vector );
	for ( core::Size i(2); i<=output_vector.size(); ++i ) {
		do_cyclic_permutation( curvector );
		if ( curvector < output_vector ) output_vector = curvector;
	}
	return output_vector;
}

/// @brief Increment the asymmetric unit.
/// @details If mirror symmetric, uses negative index values (excludes zero) for symmetric
/// possibilities.  Uses positive exclusively otherwise.
/// @returns False when the vector can no longer be incremented (i.e. all possibilities
/// have been exhausted).
bool
increment_asymm_unit_recursively (
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > & cursequence,
	core::Size const curindex,
	bool const mirror_symm,
	core::Size const maxval,
	core::Size const minval_if_symmetric
) {
	//TR << "curindex=" << curindex << std::endl;

	++cursequence[curindex];
	if ( cursequence[curindex] == 0 ) ++cursequence[curindex]; //Don't allow zeros.

	if ( cursequence[curindex] > static_cast<COUNT_CYCPEP_SEQUENCE_INTTYPE>( maxval ) ) {
		if ( curindex == 1 ) {
			return false; //We've reached the end of what we can increment
		} else {
			cursequence[curindex] = (mirror_symm ? (-1) * static_cast<COUNT_CYCPEP_SEQUENCE_INTTYPE>( minval_if_symmetric ) : 1);
			return increment_asymm_unit_recursively( cursequence, curindex - 1, mirror_symm, maxval, minval_if_symmetric );
		}
	}
	return true;
}

COUNT_CYCPEP_SEQUENCE_INTTYPE
invert_value(
	COUNT_CYCPEP_SEQUENCE_INTTYPE const val,
	core::Size const max_symm_val
) {
	if ( val <= static_cast< COUNT_CYCPEP_SEQUENCE_INTTYPE >( max_symm_val ) ) {
		return -val;
	}
	return val;
}

/// @brief Copy the first M residues, making up the asymmetric unit, to the other
/// lobes of the peptide, mirroring appropriately if we're considering mirror symmetry.
void
copy_asymm_unit_to_other_subunits(
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > & cursequence,
	core::Size const asymm_unit_size,
	bool const mirror_symm,
	core::Size const minval_if_symmetric
) {
	//TR << "aymm_unit_size=" << asymm_unit_size << std::endl;

	debug_assert( asymm_unit_size <= cursequence.size() );
	if ( asymm_unit_size == cursequence.size() ) return; //Nothing to do in the asymmetric case.

	for ( core::Size i(1); i<=asymm_unit_size; ++i ) {
		core::Size counter(i+asymm_unit_size);
		bool even( true );
		do{
			if ( mirror_symm && even ) {
				cursequence[counter] = invert_value( cursequence[i], minval_if_symmetric );
				even = false;
			} else {
				cursequence[counter] = cursequence[i];
				even = true;
			}
			counter += asymm_unit_size;
		} while( counter <= cursequence.size() );
	}
}

/// @brief Increment a vector, preserving symmetry.
/// @details Increments the asymmetric unit, then copies result to all other
/// symmetry copies.
/// @returns False when the vector can no longer be incremented (i.e. all possibilities
/// have been exhausted).
bool
increment_sequence(
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > & cursequence,
	core::Size const asymm_unit_size,
	bool const mirror_symm,
	core::Size const maxval,
	core::Size const minval_if_symmetric
) {
	if ( increment_asymm_unit_recursively( cursequence, asymm_unit_size, mirror_symm, maxval, minval_if_symmetric ) ) {
		copy_asymm_unit_to_other_subunits( cursequence, asymm_unit_size, mirror_symm, minval_if_symmetric );
	} else {
		return false;
	}
	return true;
}

/// @brief Compute the number of sequences numerically.
core::Size
count_numerically(
	core::Size const fold_symmetry,
	bool const mirror_symm,
	core::Size const nresidue,
	core::Size const noptions,
	core::Size const noptions_achiral,
	bool const print_all_sequences
) {
	TR << "Iterating through all sequences to count.  (Note that this can be extremely slow for large peptides.)" << std::endl;

	// The number of pairs of options that have a symmetry mate and their symmetry mates:
	if ( mirror_symm ) runtime_assert( (noptions-noptions_achiral) % 2 == 0 ); //Should be true from options checks.
	core::Size const noptions_symmpairs( mirror_symm ? (noptions - noptions_achiral) / 2 : 1 );

	// The size of a single lobe of the peptide:
	runtime_assert( nresidue % fold_symmetry == 0 ); //Should be true from options checks.
	core::Size const asymm_unit_size( nresidue / fold_symmetry );

	// The number of symmetry pairs of options plus the number of achiral options:
	core::Size const noptions_symmpairs_plus_achiral( mirror_symm ? noptions_symmpairs + noptions_achiral : noptions  );

	// The sequences that we have already seen, sorted.  Each is entered as the cyclic
	// permutation that is less than all others (based on std::less).
	std::set< utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > > observed_sequences;

	// The sequence vector, which will go through all permutations.
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > current_sequence( nresidue, noptions_symmpairs );
	if ( mirror_symm ) {
		for ( core::Size i(1); i<=current_sequence.size(); ++i ) {
			if ( ((i-1) / asymm_unit_size ) % 2 == 0 ) current_sequence[i] *= -1; //Flip every second lobe.
		}
	}

	// The number of unique sequences observed:
	core::Size count(0);

	do {
		utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > const min_permutation( find_min_permutation( current_sequence ) );
		if ( observed_sequences.count(min_permutation) == 0 ) {
			observed_sequences.insert( min_permutation );
			++count;
		}
		//TR << "CONSIDERED " << min_permutation << std::endl;
	} while( increment_sequence( current_sequence, asymm_unit_size, mirror_symm, mirror_symm ? noptions_symmpairs_plus_achiral : noptions, noptions_symmpairs ) );

	//TR << "Numerically, we found " << count << " unique sequences." << std::endl;

	debug_assert( count == observed_sequences.size() );

	if ( print_all_sequences ) {
		TR << "Observed sequences:" << std::endl;
		for ( std::set< utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > >::const_iterator it( observed_sequences.cbegin() ); it!=observed_sequences.cend(); ++it ) {
			for ( core::Size i(1), imax(it->size()); i<=imax; ++i ) {
				TR << (*it)[i];
				if ( i<imax ) TR << " ";
			}
			TR << std::endl;
		}
	}

	return count;
}

/// @brief Attempts to set a given entry in the pattern to
/// a given value.  Returns false if unsuccessful (the pattern
/// becomes inconsistent), true otherwise.
bool
set_pattern(
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > & pattern,
	core::Size const lobesize,
	core::Size const index,
	core::Size const setting,
	core::Size const offset
) {
	if ( pattern[index] != 0 ) return false;
	pattern[index] = setting;

	//Propagate the change:
	bool updated_one(false);
	do {
		updated_one = false;
		for ( core::Size i(1), imax(lobesize*2); i<=imax; ++i ) {
			if ( pattern[i] == 0 ) {
				if ( i > lobesize ) {
					if ( pattern[ i - lobesize ] != 0 ) {
						pattern[i] = -1*pattern[i-lobesize];
						updated_one = true;
						continue;
					}
				} else {
					if ( pattern[ i + lobesize ] != 0 ) {
						pattern[i] = -1*pattern[i+lobesize];
						updated_one = true;
						continue;
					}
				}
				core::Size const ioffset_plus( i + offset <= imax ? i + offset : i + offset - imax );
				if ( pattern[ ioffset_plus ] != 0 ) {
					pattern[i] = pattern[ ioffset_plus ];
					updated_one = true;
					continue;
				}
				core::Size const ioffset_minus(
					static_cast< signed long >( i ) - static_cast< signed long >( offset ) > 0 ?
					i - offset :
					imax + i - offset
				);
				if ( pattern[ ioffset_minus ] != 0 ) {
					pattern[i] = pattern[ ioffset_minus ];
					updated_one = true;
					continue;
				}
			}
		}
	} while( updated_one == true );

	//Check validity:
	for ( core::Size i(1), imax(lobesize*2); i<=imax; ++i ) {
		if ( i <= lobesize && pattern[ i + lobesize] != -1*pattern[i] ) return false;
		if ( i > lobesize && pattern[ i - lobesize] != -1*pattern[i] ) return false;
		core::Size const ioffset_minus(
			static_cast< signed long >( i ) - static_cast< signed long >( offset ) > 0 ?
			i - offset :
			imax + i - offset
		);
		if ( pattern[i] != pattern[ioffset_minus] ) return false;
		core::Size const ioffset_plus( i + offset <= imax ? i + offset : i + offset - imax );
		if ( pattern[i] != pattern[ioffset_plus] ) return false;
	}
	return true;
}

/// @brief Figure out how many possible sequences there are which
/// are invariant on a given cyclic permutation, given SN symmetry.
core::Size
invariants_for_cyclic_permutation(
	core::Size const asymm_unit_size,
	core::Size const offset
) {
	debug_assert( offset <= asymm_unit_size && offset > 0 );
	runtime_assert( asymm_unit_size % 2 == 0 );
	core::Size const lobesize( asymm_unit_size / 2 );
	utility::vector1< COUNT_CYCPEP_SEQUENCE_INTTYPE > pattern( asymm_unit_size, 0 );
	core::Size cur_count(0);
	for ( core::Size i(1); i<=lobesize; ++i ) {
		if ( pattern[i] == 0 ) {
			if ( !set_pattern( pattern, lobesize, i, ++cur_count, offset ) ) return 0;
		}
	}
	//TR << std::endl; TR << "OFFSET " << offset << "\t" << pattern << std::endl;
	return cur_count;
}

/// @brief Determine the number of sequences analytically, using expressions derived from
/// Burnside's lemma by Todd Yeates.  This version is for cyclic symmetry.
core::Size
count_analytically_cyclic(
	core::Size const fold_symmetry,
	core::Size const nresidue,
	core::Size const noptions
) {
	runtime_assert( nresidue % fold_symmetry == 0 ); //Should be true.
	core::Size const asymm_unit_size( nresidue / fold_symmetry );

	TR << "U = 1/" << asymm_unit_size << " * ( ";

	core::Size accumulator(0);
	for ( core::Size i(1); i<=asymm_unit_size; ++i ) {
		core::Size const gcf( boost::math::gcd( i, asymm_unit_size ) );
		accumulator += std::pow(noptions, gcf);
		TR << noptions << "^" << gcf;
		if ( i < asymm_unit_size ) TR << " + ";
	}

	TR << " )" << std::endl;

	if ( asymm_unit_size != 0 /*Silly -- to keep code quality tests happy.*/ ) {
		runtime_assert( accumulator % asymm_unit_size == 0 );
	} else {
		utility_exit_with_message( errmsg + "Division by zero in count_analytically_cyclic()!" );
	}
	return accumulator / asymm_unit_size;
}

/// @brief Determine the number of sequences analytically, using expressions derived from
/// Burnside's lemma by Todd Yeates.  This version is for improper rotational symmetry.
/// @details This version does some manual counting of the number of sequences compatible with
/// a given cyclic permutation.
core::Size
count_semi_analytically_improper_rotational(
	core::Size const fold_symmetry,
	core::Size const nresidue,
	core::Size const noptions
) {
	runtime_assert( nresidue % fold_symmetry == 0 ); //Should be true.
	core::Size const asymm_unit_size( nresidue / fold_symmetry * 2 );

	TR << "1/" << asymm_unit_size << " * ( ";

	core::Size accumulator( 0 );
	bool first( true );
	for ( core::Size i(1); i<=asymm_unit_size; ++i ) {
		core::Size const cur_invariants( invariants_for_cyclic_permutation( asymm_unit_size, i ) );
		if ( cur_invariants != 0 ) {
			accumulator += std::pow(noptions, cur_invariants);
			if ( first ) {
				first = false;
			} else {
				TR << " + ";
			}
			TR << noptions << "^" << cur_invariants;
		}
	}
	TR << " )" << std::endl;
	if ( asymm_unit_size != 0 /*Silly -- to keep code quality tests happy.*/ ) {
		runtime_assert( accumulator % asymm_unit_size == 0 );
	} else {
		utility_exit_with_message( errmsg + "Division by zero in count_semi_analytically_improper_rotational()!" );
	}
	return accumulator / asymm_unit_size;
}

/// @brief Determine the number of sequences analytically, using expressions derived from
/// Burnside's lemma by Todd Yeates.  This version is for improper rotational symmetry.
/// @details This version uses a fully analytic expression for the number of sequences compatible
/// with a given cyclic permutation.
core::Size
count_analytically_improper_rotational(
	core::Size const fold_symmetry,
	core::Size const nresidue,
	core::Size const noptions
) {
	runtime_assert( nresidue % fold_symmetry == 0 ); //Should be true.
	core::Size const lobe_size( nresidue/fold_symmetry );
	core::Size const asymm_unit_size( lobe_size * 2 );

	TR << "1/" << asymm_unit_size << " * ( ";

	core::Size accumulator( 0 );
	bool first( true );
	for ( core::Size i(1); i<=lobe_size; ++i ) {
		core::Size const gcf( boost::math::gcd( i, lobe_size ) );
		debug_assert( lobe_size % gcf == 0); //Must be true
		if ( (lobe_size / gcf) % 2 == 0 ) continue; //Skip evens.

		if ( first ) {
			first = false;
		} else {
			TR << " + ";
		}
		accumulator += std::pow(noptions, gcf);
		TR << noptions << "^" << gcf;
	}
	TR << " )" << std::endl;
	runtime_assert( accumulator % asymm_unit_size == 0 );
	return accumulator / asymm_unit_size;
}

/// @brief Determine the number of sequences analytically, using expressions derived from
/// Burnside's lemma by Todd Yeates.
core::Size
count_analytically(
	core::Size const fold_symmetry,
	bool const mirror_symm,
	bool const improper_rotational_semi_analytical,
	core::Size const nresidue,
	core::Size const noptions,
	core::Size const noptions_achiral
) {
	if ( noptions_achiral > 0 ) {
		TR.Warning << "Analytical evaluation does not currently work if there are achiral options!" << std::endl;
		return 0;
	}

	if ( mirror_symm ) {
		if ( improper_rotational_semi_analytical ) {
			return count_semi_analytically_improper_rotational( fold_symmetry, nresidue, noptions );
		} else {
			return count_analytically_improper_rotational( fold_symmetry, nresidue, noptions );
		}
	}
	return count_analytically_cyclic( fold_symmetry, nresidue, noptions );
}

/// @brief Entry point for program execution.
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		TR << "Starting count_cycpep_sequences application." << std::endl;
		TR << "Application created 13 August 2019 by Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org)." << std::endl;
		TR << "The associated manuscript is in prepration, and may be cited as \"Mulligan VK, Kang CS, Sawaya MR, Rettie S, Li X, Antselovich I, Craven T, Watkins A, Labonte JW, DiMaio F, Yeates TO, and Baker D. (2019).  Computational design of internally-symmetric peptide macrocycles with switchable folds.  Manuscript in prepration.\"." << std::endl;

		// Check and store options:
		do_options_checks();
		core::Size const fold_symmetry( option[symmetry_number]() );
		core::Size const nresidue( option[peptide_length]() );
		bool const mirror_symm( option[mirror_symmetry].value() );
		bool const do_numerical_method( option[do_numerical_count].value() );
		core::Size const noptions( option[options_per_position]() );
		core::Size const noptions_achiral( option[achiral_options]() );
		bool const print_all_sequences( option[write_out_sequences]() );
		bool const improper_rotational_semi_analytical( option[SN_semi_analytical]() );

		// The result from evaluating the count analytically:
		core::Size const analytical_result( count_analytically( fold_symmetry, mirror_symm, improper_rotational_semi_analytical, nresidue, noptions, noptions_achiral ) );

		// Write output from the analytical method:
		TR << "Analytically, counted " << analytical_result << " unique sequences for a " << nresidue << " residue cyclic peptide with ";
		if ( fold_symmetry > 1 ) {
			TR << (mirror_symm ? "S" : "C" ) << fold_symmetry;
		} else {
			TR << "no";
		}
		TR << " symmetry, with " << noptions << " possibilities at each position";
		if ( mirror_symm ) {
			TR << ", " << noptions_achiral << " of which are achiral";
		}
		TR << "." << std::endl;

		// The result from counting numerically:
		core::Size const numerical_result( do_numerical_method ? count_numerically( fold_symmetry, mirror_symm, nresidue, noptions, noptions_achiral, print_all_sequences ) : 0 );

		// Write output from numerical method:
		if ( do_numerical_method ) {
			TR << "Numerically, counted " << numerical_result << " unique sequences for a " << nresidue << " residue cyclic peptide with ";
			if ( fold_symmetry > 1 ) {
				TR << (mirror_symm ? "S" : "C" ) << fold_symmetry;
			} else {
				TR << "no";
			}
			TR << " symmetry, with " << noptions << " possibilities at each position";
			if ( mirror_symm ) {
				TR << ", " << noptions_achiral << " of which are achiral";
			}
			TR << "." << std::endl;
			runtime_assert_string_msg( numerical_result == analytical_result, errmsg + "The numerical and analytical results do not match!" );
		}

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "Caught exception " << e.msg() << std::endl;
		return -1;
	}
	TR << "Completed count_cycpep_sequences application.  Exiting with no errors (status zero)." << std::endl;
	return 0;
}
