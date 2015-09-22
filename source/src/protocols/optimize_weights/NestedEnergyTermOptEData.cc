// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/optimize_weights/NestedEnergyTermOptEData.cc
/// @author  Ron Jacak

#ifdef USEMPI
#include <mpi.h>
#endif

// Unit headers
#include <protocols/optimize_weights/OptEData.hh>
#include <protocols/optimize_weights/NestedEnergyTermOptEData.hh>

// Project headers
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

#include <utility/string_util.hh>
#include <utility/vector1.functions.hh> // to get arg_min()


// C++ headers
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>

// option key includes
#include <basic/options/keys/optE.OptionKeys.gen.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/options/option.hh>


using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::format;

namespace protocols {
namespace optimize_weights {

typedef utility::vector1< std::string > Strings;

static THREAD_LOCAL basic::Tracer TR( "NestedEnergyTermOptEData" );

#define CAP_FA_REP 1


// ------------------- NestedEnergyTermPNatAAOptEPositionData -----------------------//

///
NestedEnergyTermPNatAAOptEPositionData::NestedEnergyTermPNatAAOptEPositionData() {}

///
NestedEnergyTermPNatAAOptEPositionData::~NestedEnergyTermPNatAAOptEPositionData() {}

///
/// @brief
/// Does actual work for OptE minimization
/// Special implementation of get_score that includes logic to handle unfolded state energy calculation.
/// See header file and OptEData.hh for more information.
///
Real
NestedEnergyTermPNatAAOptEPositionData::get_score(
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const num_total_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const & score_list,
	ScoreTypes const & fixed_score_list
) const
{
	return process_score( TR, false, component_weights, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, score_list, fixed_score_list );
}

///
/// @brief
/// Special implementation of print_score that includes logic to handle unfolded state energy calculation.
///
void
NestedEnergyTermPNatAAOptEPositionData::print_score(
	std::ostream & ostr,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const num_total_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const & score_list,
	ScoreTypes const & fixed_score_list
) const
{
	process_score( ostr, true, component_weights, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, score_list, fixed_score_list );
}

///
/// @brief
/// One method to do the score processing which takes a boolean dictating whether to print to an ostream or not. With this function, changes
/// to how scoring works only need to be made in one place as opposed to two (when get_score() and print_score() both had scoring logic in them).
Real
NestedEnergyTermPNatAAOptEPositionData::process_score(
	std::ostream & ostr,
	bool print,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const & score_list,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core;
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	chemical::AA const this_native_aa( native_aa() );

	static Real const inv_kT( option[ optE::inv_kT_nataa ] );

	// contains the best energy for each amino acid at this position; init to a bad energy of 1000
	utility::vector1< Real > best_energy_by_aa( chemical::num_canonical_aas, 1000.0 );

	// containers for derivatives
	// all values are initialized to zero
	Multivec ref_deriv_weight( chemical::num_canonical_aas, 0.0 );
	utility::vector1< Real > vector_of_zeros( score_list.size(), 0.0 );  // assigned to unweighted_E_dof for bad rotamers
	utility::vector1< utility::vector1< Real > > unweighted_E_dof( chemical::num_canonical_aas, vector_of_zeros );

	process_rotamers( vars, num_energy_dofs, fixed_terms, score_list, fixed_score_list, chemical::num_canonical_aas,
		vector_of_zeros, best_energy_by_aa, unweighted_E_dof, ref_deriv_weight );

	//TR << "process_score(): best_energy_by_aa before: [ ";
	//for ( Size i=1; i <= best_energy_by_aa.size(); ++i ) {
	// TR << F(6,2,best_energy_by_aa[i]) << ", ";
	//}
	//TR << std::endl;

	// now do some special processing of the unfolded state energy
	// the parameters score_list and fixed_score_list are vectors of ScoreType objects; these should be in the same order
	// as the weights in the vars array.
	// this part has two steps: first we have to take the weights for the canonical score12 score terms and apply them
	// to the unweighted unfolded energies stored in the unfolded_energy_emap_vector. then, we have to take that entire
	// sum and multiply it by the unfolded term weight currently being evaluated.
	// If you neglect to multiply the total unfolded state energy of an aa by the current unfolded term weight, then
	// the unfolded term weight varies wildly and minimization doesn't do anything.

	for ( Size aa = 1; aa <= chemical::num_canonical_aas; ++aa ) {

		Real unfolded_energy_for_one_aa = 0.0;
		Real weighted_unfolded_energy_for_one_aa = 0.0;

		// Part 1a:
		// Variable-weighted energy terms
		//
		// Assume free params are fa_rep, solubility, and unfolded. The unfolded_energy_emap_vector only contains energies
		// for fa_rep; solubility and unfolded always return zeros.  Thus, even though we iterate through all the free terms
		// here, only the fa_rep will contribute to the unfolded_energy_for_one_aa.
		//
		// The neat thing about this setup is that if a score12 energy term is not included as a free or fixed param, then
		// it won't be included in the unfolded state energy either.
		//
		for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
			//if ( print ) {
			//TR << "process_score(): adding unfolded energy for aa '" << chemical::name_from_aa( (chemical::AA) aa )
			// << "' for unweighted free '" << name_from_score_type( score_list[ ii ] ) << "' energy: "
			// << unfolded_energy_emap_vector_[ aa ][ score_list[ ii ]] << " * '"
			// << name_from_score_type( score_list[ ii ] )
			// << "' weight: " << vars[ ii ]
			// << " to local unfolded_energy_for_one_aa variable." << std::endl;
			//}
			// the unweighted, unfolded energy for this ScoreType times the weight
			unfolded_energy_for_one_aa += (unfolded_energy_emap_vector_[ aa ][ score_list[ ii ]] * vars[ ii ]);

			// see comments in process_rotamers() commented out below for what's going on in the next block
			if ( ref_deriv_weight[ aa ] != 0.0 ) {

				// Subtract the unweighted, unfolded energy from unweighted_E_dof array (which is broken up by aa and eterm)
				// This is so that the minimizer knows how to adjust this weight during minimization.

				// aa is the amino acid, ii is the free weight; only free terms have a spot in the unweighted_E_dof array
				// so we only have to iterate over the array in this one section
				unweighted_E_dof[ aa ][ ii ] -= unfolded_energy_emap_vector_[ aa ][ score_list[ ii ] ];

				// What if the unfolded weight is free/floating? What do we use for the unweighted_E_dof for the unfolded
				// state energy? Well, the unfolded energy method returns 0.0 for the energy at this point of optE.
				// Therefore, the unweighted_E_dof should also be 0.0.  The unfolded_energy_emap_vector should contain a
				// 0.0 at the location for unfolded, so no special code is needed here to prevent problems if the unfolded
				// weight is variable.
			}
		}

		// Part 1b:
		// Fixed-weight energy terms
		for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
			//if ( print ) {
			//TR << "process_score(): adding unfolded energy for aa '" << chemical::name_from_aa( (chemical::AA) aa )
			// << "' unweighted fixed '" << name_from_score_type( fixed_score_list[ ii ] ) << "' energy: "
			// << unfolded_energy_emap_vector_[ aa ][ fixed_score_list[ ii ]] << " * '"
			// <<  name_from_score_type( fixed_score_list[ ii ] )
			// << "' weight: " << fixed_terms[ fixed_score_list[ ii ] ] << std::endl;
			//}
			unfolded_energy_for_one_aa += (unfolded_energy_emap_vector_[ aa ][ fixed_score_list[ ii ] ] * fixed_terms[ fixed_score_list[ ii ] ]);
		}

		// Part 1c:
		// Reference energy term
		// The unfolded state energy map should not have anything for a reference energy. Even if it does though, don't add anything
		// else here.

		// Part 2: Now use the current 'unfolded' term weight
		// 'unfolded' could be a free or fixed term so iterate over both lists; here it's necessary because this energy gets added
		// to the best energy by aa array regardless.
		for ( Size tt = 1; tt <= num_energy_dofs; ++tt ) {
			if ( name_from_score_type( score_list[ tt ] ) == "unfolded" ) {
				//if ( print ) {
				// TR << chemical::name_from_aa( (chemical::AA) aa )
				//  << " unweighted unfolded energy: " << unfolded_energy_for_one_aa << ", free '"
				//  << name_from_score_type( score_list[ tt ] ) << "' term weight: " << vars[tt]
				//  << ", weighted unfolded energy: " << unfolded_energy_for_one_aa * vars[ tt ] << std::endl;
				//}
				weighted_unfolded_energy_for_one_aa = unfolded_energy_for_one_aa * vars[ tt ];
			}
		}
		for ( Size tt = 1; tt <= fixed_score_list.size(); ++tt ) {
			if ( name_from_score_type( fixed_score_list[ tt ] ) == "unfolded" ) {
				//if ( print ) {
				// TR << chemical::name_from_aa( (chemical::AA) aa )
				//  << " unfolded energy: " << unfolded_energy_for_one_aa << ", FIXED '"
				//  << name_from_score_type( fixed_score_list[ tt ] ) << "' term weight: " << fixed_terms[ fixed_score_list[ tt ] ]
				//  << ", weighted unfolded energy: " << unfolded_energy_for_one_aa * fixed_terms[ fixed_score_list[ tt ] ] << std::endl;
				//}
				weighted_unfolded_energy_for_one_aa = unfolded_energy_for_one_aa * fixed_terms[ fixed_score_list[ tt ] ];
			}
		}

		// SUBTRACT this unfolded energy from the value for this aa in the best_energy_by_aa array. We subtract because only positive weights are
		// being used in the driver class, and unfolded energy should be subtracted from folded energy to get difference in free energy.
		// See the doxygen documentation for function optimize_weights() in IterativeOptEDriver. cc
		best_energy_by_aa[ aa ] -= weighted_unfolded_energy_for_one_aa;

		// do a check to make sure the energy doesn't exceed our cutoff (APL set at 300) or otherwise we get INF/NAN errors during minimization
		Real const cutoff( 300.0 );
		if ( best_energy_by_aa[ aa ] > cutoff ) {
			best_energy_by_aa[ aa ] = cutoff;
		} else if ( best_energy_by_aa[ aa ] < -1.0*cutoff ) {
			best_energy_by_aa[ aa ] = -1.0*cutoff;
		}

	}

	//TR << "process_score(): best_energy_by_aa after : [ ";
	//for ( Size i=1; i <= best_energy_by_aa.size(); ++i ) {
	// TR << F(6,2,best_energy_by_aa[i]) << ", ";
	//}
	//TR << std::endl;


	// now do the partition function analysis
	// all this should remain the same assuming the arrays were correctly set in process_rotamers()

	Real numerator(0.0), partition(0.0);
	Multivec dpartition( vars.size(), 0.0 ), dnumerator( vars.size(), 0.0 );

	for ( Size aa(1); aa <= chemical::num_canonical_aas; ++aa ) {

		Real const exp_term( std::exp( -1.0 * inv_kT * best_energy_by_aa[ aa ] ) );
		partition += exp_term;
		if ( aa == size_t( this_native_aa ) ) {
			numerator = exp_term;
		}

		// for reference energy derivatives, but don't assume the protocol is using them
		// note for derivatives: dE/dw( e^-(E*w+...) ) = -E * e^-(E*w+...) but as this is an energy, the 'weight' here is 0 or 1
		// reference energies do not have an unweighted energy.  really, the energy is just 1 or 0 depending on what the process_rotamers()
		// function decided. so d/dw[ e(-E*w) ] = -1 * (0|1) * e(-E*w)
		if ( num_ref_dofs != 0 ) {
			Real const ref_deriv_term( -1.0 * inv_kT * ref_deriv_weight[ aa ] * exp_term );
			dpartition[ num_energy_dofs + aa ] = ref_deriv_term;
			if ( aa == size_t(this_native_aa) ) {
				dnumerator[ num_energy_dofs + aa ] = ref_deriv_term;
			}
		}

		// partitions for energy derivatives
		// note for derivatives: dE/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
		// dE/dweight( e^(-unweightedE*weight + ...) ) = -1 * unweightedE * e^(-unweightedE*weight + ...)
		// and best_energy_by_aa is the same as (unweightedE * weight)

		// there's a potential problem here I'm not sure how it happens. If the best energy by aa is something really small (e.g. -700)
		// then the exponential e(700) is an extremely large number. the runtime just assigns it INF.  Then when you multiply it by
		// zero, you still get INF for some reason.
		for ( Size e_term = 1; e_term <= num_energy_dofs; ++e_term ) {
			Real e_term_deriv( -1.0 * inv_kT * unweighted_E_dof[ aa ][ e_term ] * exp_term );
			dpartition[ e_term ] += e_term_deriv;
			if ( aa == size_t( this_native_aa ) ) {
				dnumerator[ e_term ] = e_term_deriv;
			}
		}
	}

	// accumulate to passed-in derivative sums
	for ( Size dof(1); dof <= vars.size(); ++dof ) {
		dE_dvars[ dof ] += component_weights[ type() ] * ( dpartition[ dof ] / partition - dnumerator[ dof ] / numerator );

		if ( score_list[ dof ] == omega ) { dE_dvars[ dof ] = 0.0; }
		if ( score_list[ dof ] == hbond_lr_bb ) { dE_dvars[ dof ] = 0.0; }

		/*if ( tag() == "1nls_" && (int)this_native_aa == 3 ) {
		if ( score_list[ dof ] == omega ) {
		TR << "PNATAA " << tag() << X(1) << this_native_aa << "," << I(2, (int)this_native_aa) << X(1)
		<< "-lnp: " << F(6,4,-1.0 * std::log( numerator / partition ))
		<< ", dE_dvars[ omega ]: " << F(6,4, dE_dvars[ dof ])
		<< ", best_energy_by_aa: [ ";
		for ( Size i=1; i <= best_energy_by_aa.size(); ++i ) { TR << F(5,2,best_energy_by_aa[i]) << ", "; } TR << "], ";
		TR << ", dpart[ omega ]: " << F(7,3,dpartition[dof]) << ", part: " << F(7,3,partition)
		<< ", dnum[ omega ]: " << F(7,3,dnumerator[dof]) << ", num: " << F(7,3,numerator) << std::endl;
		}
		if ( score_list[ dof ] == unfolded ) {
		TR << "PNATAA " << tag() << X(1) << this_native_aa << "," << I(2, (int)this_native_aa) << X(1)
		<< " dpart[ unfolded ]: " << F(7,3,dpartition[dof]) << " part: " << F(7,3,partition)
		<< " dnum[ unfolded ]: " << F(7,3,dnumerator[dof]) << " num: " << F(7,3,numerator)
		<< " dE_dvars[ unfolded ]: " << F(6,4, dE_dvars[ dof ])
		<< " dE_dvars: [ ";
		for ( Size ii=1; ii <= vars.size(); ++ii ) { TR << F(5,2,dE_dvars[ii]) << ", "; }
		TR << "]" << std::endl;
		}
		}*/

	}

	if ( print ) {
		ostr << "PNATAA " << tag() << X(1) << this_native_aa << "," << I(2, (int)this_native_aa) << X(1)
			<< " nbs: " << I(2,neighbor_count())
			<< " num: " << F(7,3,numerator) << " part: " << F(7,3,partition)
			<< " p: " << F(7,5,numerator / partition)
			<< " -lnp: " << F(6,4,-1.0 * std::log( numerator / partition ))
			<< " -compwt_lnp: " << F(6, 4, component_weights[ type() ] * (-1.0 * std::log( numerator / partition )) )
			<< " best_energy_by_aa: [ ";

		for ( Size i=1; i <= best_energy_by_aa.size(); ++i ) {
			ostr << F(5,2,best_energy_by_aa[i]) << ", ";
		}
		ostr << "]" << std::endl;
	}

	return ( -1.0 * component_weights[ type() ] * std::log( numerator / partition ) );
}

///
/// @brief
/// To be sure we create the right types when writing/reading from files, need to add a special OptEPositionData
/// type that gets returned here.
///
OptEPositionDataType
NestedEnergyTermPNatAAOptEPositionData::type() const {
	return prob_native_amino_acid_with_unfolded_energy;
}

///
/// @brief
/// Add a special for loop to print out the unfolded state energy EnergyMap values.
///
void
NestedEnergyTermPNatAAOptEPositionData::write_to_file( std::ofstream & outfile ) const {

	outfile << "position " << position() << " "
		<< "nataa " << native_aa() << " "
		<< "neighbor_count " << neighbor_count() << " "
		<< "unfolded_energy" << std::endl;

	// print one line with all the emap info per aa
	for ( int aa = 1; aa < chemical::num_canonical_aas; ++aa ) {
		for ( int type = 1; type < scoring::n_score_types; ++type ) {
			outfile << name_from_score_type( ScoreType( type ) ) << "=" << (unfolded_energy_emap_vector_[ aa ])[ ScoreType(type) ] << " ";
		}
		outfile << std::endl;
	}
	outfile << std::endl;

	outfile << "nrots " << data().size() << "\n";
	for ( PNatAAOptERotamerDataOPs::const_iterator rot( rotamer_data_begin() ); rot != rotamer_data_end(); ++rot ) {
		outfile << *rot << std::endl;
	}

}

///
/// @brief
///
void
NestedEnergyTermPNatAAOptEPositionData::read_from_file( std::ifstream & infile ) {

	using namespace utility;

	// read first line with position, native aa, neighbor_count, and num_rotamers data
	std::string line;
	getline( infile, line );
	Strings words( string_split( line, ' ' ) );
	assert( words[ 1 ] == "position" );
	set_position( from_string( words[ 2 ], Size( 0 ) ) );
	assert( words[ 3 ] == "nataa" );
	set_native_aa ( chemical::aa_from_name( words[ 4 ] ) );
	assert( words[ 5 ] == "neighbor_count" );
	set_neighbor_count( from_string( words[ 6 ], Size( 0 ) ) );

	// extra logic to handle reading in the unfolded state energies into an EnergyMap
	assert( words[ 7 ] == "unfolded_energy" );
	utility::vector1 < EnergyMap > emap_vector;
	emap_vector.resize( chemical::num_canonical_aas );

	for ( int aa = 1; aa < chemical::num_canonical_aas; ++aa ) {
		getline( infile, line );
		Strings sections( string_split( line, ' ' ) );

		EnergyMap emap;
		for ( Strings::iterator section = sections.begin(); section != sections.end(); ++section ) {
			Strings pair( string_split( *section, '=' ) );
			ScoreType st = ScoreTypeManager::score_type_from_name( pair[1] );
			Real score;
			std::istringstream ss( pair[2] );
			ss >> score;
			emap[ st ] = score;
		}
		emap_vector[ aa ] = emap;
	}

	getline( infile, line );
	Strings rotamer_line_words( string_split( line, ' ' ) );
	assert( rotamer_line_words[ 1 ] == "nrots" );
	Size num_rotamers = from_string( rotamer_line_words[ 2 ], Size( 0 ) );

	for ( Size ii = 1; ii <= num_rotamers; ++ii ) {
		getline( infile, line );

		// rotamers for existing position: parse, append new OptERotamerDataOP to OptEPositionDataOP
		Strings sections( string_split( line, ',' ) );
		// sections:
		// 0 - rotnum, 1 - aa three-letter code, 2 - energies for fixed terms, 3 - energies for free terms
		Size rotnum;
		std::istringstream ss( sections[1] );
		ss >> rotnum;
		chemical::AA aa( chemical::aa_from_name( sections[2] ) );
		utility::vector1< Real > fixed_energies, energies;
		Strings fixed_vals( string_split( sections[3], ' ' ) ), free_vals( string_split( sections[4], ' ' ) );
		for ( Strings::iterator fixed_val( fixed_vals.begin() ); fixed_val != fixed_vals.end(); ++fixed_val ) {
			Real val;
			std::istringstream ss( *fixed_val );
			ss >> val;
			fixed_energies.push_back( val );
		}
		for ( Strings::iterator free_val( free_vals.begin() ); free_val != free_vals.end(); ++free_val ) {
			Real val;
			std::istringstream ss( *free_val );
			ss >> val;
			energies.push_back( val );
		}
		assert( !energies.empty() );
		PNatAAOptERotamerDataOP new_rot_data( new PNatAAOptERotamerData( aa, rotnum, energies, fixed_energies ) );
		add_rotamer_line_data( new_rot_data );
	}
}

///
/// @brief
/// Leaving this unimplemented since I don't feel like figuring out how to output an EnergyMap in binary
/// and also because reading/writing binary files is not being used in the optE protocol currently.
///
void
NestedEnergyTermPNatAAOptEPositionData::write_to_binary_file( std::ofstream & /* outfile */ ) const {}

///
/// @brief
/// Leaving this unimplemented since I don't feel like figuring out how to output an EnergyMap in binary
/// and also because reading/writing binary files is not being used in the optE protocol currently.
///
void
NestedEnergyTermPNatAAOptEPositionData::read_from_binary_file( std::ifstream & /* infile */ ) {}

///
Size
NestedEnergyTermPNatAAOptEPositionData::memory_use() const {

	Size total = sizeof( NestedEnergyTermPNatAAOptEPositionData ) + sizeof( PNatAAOptERotamerData ) * data().size();
	if ( data().size() > 0 ) {
		total +=  sizeof( Real ) * ( data()[ 1 ]->data().size() + data()[ 1 ]->fixed_data().size() ) * data().size();
	}
	// the emap vector uses space, too!
	total += sizeof( EnergyMap ) * unfolded_energy_emap_vector_.size();

	return total;
}


#ifdef USEMPI
///
///
void
NestedEnergyTermPNatAAOptEPositionData::send_to_node( int const destination_node, int const tag ) const {

	/// 1. Which position is this?
	int ii_pos = position();
	MPI_Send( & ii_pos, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 2. What is the native amino acid at this position?
	int ii_aa  = native_aa();
	MPI_Send( & ii_aa, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 3. How many neighbors did this position have?
	int ii_neighbor_count = neighbor_count();
	MPI_Send( & ii_neighbor_count, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 4. How many aa's are in the unfolded state energy map vector?
	/// This will probably always be 20, but send it anyway
	int ii_unfolded_energy_emap_vector_size = unfolded_energy_emap_vector_.size();
	MPI_Send( & ii_unfolded_energy_emap_vector_size, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 4b. The energies in the emap
	Real * unfolded_energies = new Real[ chemical::num_canonical_aas * scoring::n_score_types ];
	for ( int aa = 1; aa <= chemical::num_canonical_aas; ++aa ) {
		for ( Size ee = 1; ee <= scoring::n_score_types; ++ee ) {
			unfolded_energies[ (aa - 1) * scoring::n_score_types + ee - 1 ] = unfolded_energy_emap_vector_[ aa ][ (ScoreType) ee ];
		}
	}
	MPI_Send( unfolded_energies, chemical::num_canonical_aas * scoring::n_score_types, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );
	delete [] unfolded_energies;   unfolded_energies = 0;

	/// 5. The number of rotamers for this position
	Size ii_num_rotamers = size();
	MPI_Send( & ii_num_rotamers, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	if ( ii_num_rotamers == 0 )
		return;

	Size free_count = data()[1]->data().size();
	Size fixed_count = data()[1]->fixed_data().size();

	/// 6 The size of the free and fixed data, since that context is not available.
	MPI_Send( & free_count, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );
	MPI_Send( & fixed_count, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	int * ii_aa_types = new int[ ii_num_rotamers ];
	int * ii_rot_nums = new int[ ii_num_rotamers ];
	Real * free_data = new Real[ ii_num_rotamers * free_count ];
	Real * fixed_data = new Real[ ii_num_rotamers * fixed_count ];
	for ( Size jj = 1; jj <= ii_num_rotamers; ++jj ) {
		ii_aa_types[ jj - 1 ] = data()[ jj ]->this_aa();
		ii_rot_nums[ jj - 1 ] = data()[ jj ]->rot_number();
		for ( Size kk = 1; kk <= free_count; ++kk ) {
			free_data[ ( jj - 1 ) * free_count + kk - 1 ] = data()[ jj ]->data()[ kk ];
		}
		for ( Size kk = 1; kk <= fixed_count; ++kk ) {
			fixed_data[ ( jj - 1 ) * fixed_count + kk - 1 ] = data()[ jj ]->fixed_data()[ kk ];
		}
	}

	/// 7. All the amino acids for all rotamers at this position
	MPI_Send( ii_aa_types, ii_num_rotamers, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 8. All the rotamer indices for all rotamers at this position
	MPI_Send( ii_rot_nums, ii_num_rotamers, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 9. All the free data for all rotamers
	MPI_Send( free_data, ii_num_rotamers * free_count, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 10. All the fixed data for all rotamers
	MPI_Send( fixed_data, ii_num_rotamers * fixed_count, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	delete [] ii_aa_types; ii_aa_types = 0;
	delete [] ii_rot_nums; ii_rot_nums = 0;
	delete [] free_data;   free_data   = 0;
	delete [] fixed_data;  fixed_data  = 0;

	OptEPositionData::send_to_node( destination_node, tag );
}

///
void
NestedEnergyTermPNatAAOptEPositionData::receive_from_node( int const source_node, int const tag ) {

	MPI_Status stat;

	/// 1. Which sequence position is this?
	int ii_pos;
	MPI_Recv( & ii_pos, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	set_position( ii_pos );

	/// 2. What is the native amino acid at this position?
	int ii_aa;
	MPI_Recv( & ii_aa, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	set_native_aa( (chemical::AA) ii_aa );

	/// 3. How many neighbors did this position have?
	int ii_neighbor_count;
	MPI_Recv( & ii_neighbor_count, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	set_neighbor_count( ii_neighbor_count );

	/// 4. How many emaps are in the unfolded state energy map vector?
	/// Most likely 20, for the 20 canonical aas, but check anyway.
	int ii_unfolded_energy_emap_vector_size;
	MPI_Recv( & ii_unfolded_energy_emap_vector_size, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	utility::vector1 < EnergyMap > emap_vector;
	emap_vector.resize( ii_unfolded_energy_emap_vector_size );

	/// 4b. Now get the energies in the emap
	Real * unfolded_energies = new Real[ chemical::num_canonical_aas * scoring::n_score_types ];
	MPI_Recv( unfolded_energies, chemical::num_canonical_aas * scoring::n_score_types, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	for ( Size aa = 1; aa <= chemical::num_canonical_aas; ++aa ) {
		for ( Size ee = 1; ee <= scoring::n_score_types; ++ee ) {
			// be careful because the array is 0-based while the score type enum is 1-based
			emap_vector[ aa ][ (ScoreType) ee ] = unfolded_energies[ (aa-1) * scoring::n_score_types + ee - 1 ];
		}
	}
	set_unfolded_energy_emap_vector( emap_vector );
	delete [] unfolded_energies;   unfolded_energies = 0;


	/// 6. The number of rotamers for this position
	Size ii_num_rotamers;
	MPI_Recv( & ii_num_rotamers, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	if ( ii_num_rotamers == 0 )
		return;

	Size free_count(0);
	Size fixed_count(0);

	/// 6b The size of the free and fixed data, since that context is not available.
	MPI_Recv( & free_count, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	MPI_Recv( & fixed_count, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );


	int * ii_aa_types = new int[ ii_num_rotamers ];
	int * ii_rot_nums = new int[ ii_num_rotamers ];
	Real * free_data = new Real[ ii_num_rotamers * free_count ];
	Real * fixed_data = new Real[ ii_num_rotamers * fixed_count ];

	/// 7. All the amino acids for all rotamers at this position
	MPI_Recv( ii_aa_types, ii_num_rotamers, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 8. All the rotamer indices for all rotamers at this position
	MPI_Recv( ii_rot_nums, ii_num_rotamers, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 9. All the free data for all rotamers
	MPI_Recv( free_data, ii_num_rotamers * free_count, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 10. All the fixed data for all rotamers
	MPI_Recv( fixed_data, ii_num_rotamers * fixed_count, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	utility::vector1< Real > free_data_vect( free_count );
	utility::vector1< Real > fixed_data_vect( fixed_count );

	for ( Size jj = 1; jj <= ii_num_rotamers; ++jj ) {
		for ( Size kk = 1; kk <= free_count; ++kk ) {
			free_data_vect[ kk ] = free_data[ ( jj - 1 ) * free_count + kk - 1 ];
		}
		for ( Size kk = 1; kk <= fixed_count; ++kk ) {
			fixed_data_vect[ kk ] = fixed_data[ ( jj - 1 ) * fixed_count + kk - 1 ];
		}
		PNatAAOptERotamerDataOP jj_rotamer_data( new PNatAAOptERotamerData(
			(chemical::AA ) ii_aa_types[ jj - 1 ],
			ii_rot_nums[ jj - 1 ],
			free_data_vect,
			fixed_data_vect ) );
		add_rotamer_line_data( jj_rotamer_data );

	}

	delete [] ii_aa_types; ii_aa_types = 0;
	delete [] ii_rot_nums; ii_rot_nums = 0;
	delete [] free_data;   free_data   = 0;
	delete [] fixed_data;  fixed_data  = 0;

	OptEPositionData::receive_from_node( source_node, tag );

}
#endif


// ------------------- NestedEnergyTermDDGMutationOptEData -----------------------//

///
NestedEnergyTermDDGMutationOptEData::NestedEnergyTermDDGMutationOptEData() {}

///
NestedEnergyTermDDGMutationOptEData::~NestedEnergyTermDDGMutationOptEData() {}

///
/// @details
/// This get_score() method needs to contain some extra logic for the unfolded state energy term.  Right now, the value
/// under unfolded energy is 0.0, because that's what the EnergyMethod is coded to return.  But, as in the NatAA class
/// above, we need use the unweighted, unfolded energies and the current weight set to come up with a unfolded energy.
///
Real
NestedEnergyTermDDGMutationOptEData::get_score(
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const num_total_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list
) const
{
	return process_score( TR, false, component_weights, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, free_score_list, fixed_score_list );
}

///
void
NestedEnergyTermDDGMutationOptEData::print_score(
	std::ostream & ostr,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const num_total_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list
) const
{
	process_score( ostr, true, component_weights, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, free_score_list, fixed_score_list );
}

///
/// @brief
/// One method to do the score processing which takes a boolean dictating whether to print to an ostream or not. With this function, changes
/// to how scoring works only need to be made in one place as opposed to two (when get_score() and print_score() both had scoring logic in them).
Real
NestedEnergyTermDDGMutationOptEData::process_score(
	std::ostream & ostr,
	bool print,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const & score_list,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility;

	// if there are no structures to go through, return immediately
	if ( muts_.size() == 0 || wts_.size() == 0 ) return 0.0;

	// these vectors are sized to the number of structures there are for a wt name and mutant name;
	// they'll be used to determine which structure has the best energy
	utility::vector1< Real > wt_energies( wts_.size(), 0.0 );
	utility::vector1< Real > mut_energies( muts_.size(), 0.0 );

	// go through and come up with a total score for each structure in the wts_ and muts_ list.
	//
	// wts_ is a vector1 of SingleStructureData (SSD) objects. this for loop iterates over each free weight and
	// takes the unweighted energy for the current free term in wts_[jj] and multiplies it by the weight in vars.
	// so the term 'unfolded' will have an energy of zero, which is what we need to fix.  the SSD objects have
	// a vector1 of Reals accessible by free_data() and fixed_data() member functions.  These store the total
	// energies by score type for the entire pose.

	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {

			// cap the fa_rep term at some value - this at least keeps it around for most of the mutants
#ifdef CAP_FA_REP
			if ( ( score_list[ ii ] == fa_rep ) && ( vars[ ii ] * wts_[ jj ]->free_data()[ ii ] > 10 ) ) { wt_energies[ jj ] += 10; }
			else
#endif
			wt_energies[ jj ] += vars[ ii ] * wts_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
#ifdef CAP_FA_REP
			if ( ( score_list[ ii ] == fa_rep ) && ( vars[ ii ] * muts_[ jj ]->free_data()[ ii ] > 10 ) ) { mut_energies[ jj ] += 10; }
			else
#endif
			mut_energies[ jj ] += vars[ ii ] * muts_[ jj ]->free_data()[ ii ];
		}
	}
	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
#ifdef CAP_FA_REP
			if ( ( fixed_score_list[ ii ] == fa_rep ) && ( fixed_terms[ fixed_score_list[ ii ] ] * wts_[ jj ]->fixed_data()[ ii ] > 10 ) ) { wt_energies[ jj ] += 10; }
			else
#endif
			wt_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * wts_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
#ifdef CAP_FA_REP
			if ( ( fixed_score_list[ ii ] == fa_rep ) && ( fixed_terms[ fixed_score_list[ ii ] ] * muts_[ jj ]->fixed_data()[ ii ] > 10 ) ) { mut_energies[ jj ] += 10; }
			else
#endif
			mut_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * muts_[ jj ]->fixed_data()[ ii ];
		}
	}

	// I presume these are the reference energies that are being added in?
	// num_energy_dofs is the number of free, non-reference energy parameters in the run, so yes these are the refE's
	if ( num_ref_dofs != 0 ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
			wt_energies[ jj ] += vars[ num_energy_dofs + wt_aa_ ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			mut_energies[ jj ] += vars[ num_energy_dofs + mut_aa_ ];
		}
	}

	//TR << "process_score(): before unfolded wts_: [ ";
	//for ( Size jj = 1; jj <= wts_.size(); ++jj ) { TR << F(6,1,wt_energies[ jj ]) << ", "; }
	//TR << "]" << std::endl;

	// now do some special processing of the unfolded state energy
	// this part has two steps: first we have to take the weights for the canonical score12 score terms and apply them
	// to the unweighted unfolded energies stored in the member variable emap. then, we have to take that entire
	// sum of products and multiply it by the unfolded term weight currently being evaluated.
	Real wt_unweighted_unfolded_energy(0.0), wt_weighted_unfolded_energy(0.0);
	Real mut_unweighted_unfolded_energy(0.0), mut_weighted_unfolded_energy(0.0);

	// Part 1a: Variable-weighted energy terms
	//TR << "process_score(): weighted free unfolded energies: [ ";
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {

		//TR << name_from_score_type( score_list[ ii ] ) << ": "
		// << ( vars[ii] * wt_unfolded_energies_emap_[ score_list[ ii ] ] ) << ", "
		// << ( vars[ii] * mut_unfolded_energies_emap_[ score_list[ ii ] ] ) << " / ";

		//TR << "process_score(): adding unfolded '" << name_from_score_type( score_list[ ii ] )
		// << "' energy: " << wt_unfolded_energies_emap_[ score_list[ ii ] ]
		// << " * free weight: " << vars[ ii ]
		// << " = " << ( vars[ii] * wt_unfolded_energies_emap_[ score_list[ ii ] ] ) << std::endl;

		// Assume free params are fa_rep, solubility, and unfolded. The unfolded_energy_emap_vector only contains energies
		// for fa_rep; solubility and unfolded always return zeros.  Thus, even though we iterate through all the free terms
		// here, only the fa_rep will contribute to the unfolded_energy_for_one_aa.
		//
		// The neat thing about this setup is that if a score12 energy term is not included as a free or fixed param, then
		// it won't be included in the unfolded state energy either.
		wt_unweighted_unfolded_energy += ( vars[ii] * wt_unfolded_energies_emap_[ score_list[ ii ] ] );
		mut_unweighted_unfolded_energy += ( vars[ii] * mut_unfolded_energies_emap_[ score_list[ ii ] ] );
	}
	//TR << std::endl;

	// Part 1b: Fixed-weight energy terms
	//TR << "process_score(): weighted fixed unfolded energies: [ ";
	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {

		//TR << name_from_score_type( fixed_score_list[ ii ] ) << ": "
		// << ( fixed_terms[ fixed_score_list[ ii ] ] * wt_unfolded_energies_emap_[ fixed_score_list[ ii ] ] ) << ", "
		// << ( fixed_terms[ fixed_score_list[ ii ] ] * mut_unfolded_energies_emap_[ fixed_score_list[ ii ] ] ) << " / ";

		//TR << "process_score(): adding unfolded '" << name_from_score_type( fixed_score_list[ ii ] )
		// << "' energy: " << wt_unfolded_energies_emap_[ fixed_score_list[ ii ] ]
		// << " * fixed weight: " << fixed_terms[ fixed_score_list[ ii ] ]
		// << " = " << ( fixed_terms[ fixed_score_list[ ii ] ] * wt_unfolded_energies_emap_[ fixed_score_list[ ii ] ] ) << std::endl;

		wt_unweighted_unfolded_energy += ( fixed_terms[ fixed_score_list[ ii ] ] * wt_unfolded_energies_emap_[ fixed_score_list[ ii ] ] );
		mut_unweighted_unfolded_energy += ( fixed_terms[ fixed_score_list[ ii ] ] * mut_unfolded_energies_emap_[ fixed_score_list[ ii ] ] );
	}
	//TR << std::endl;

	// Part 2: Now use the current 'unfolded' term weight
	// 'unfolded' could be a free or fixed term so iterate over both lists;
	Real unfolded_weight = 0.0;
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		if ( name_from_score_type( score_list[ ii ] ) == "unfolded" ) {
			wt_weighted_unfolded_energy = vars[ii] * wt_unweighted_unfolded_energy;
			mut_weighted_unfolded_energy = vars[ii] * mut_unweighted_unfolded_energy;
			unfolded_weight = vars[ii];
			//TR << "process_score(): weighting unweighted unfolded energy: '" << wt_unweighted_unfolded_energy
			// << " by free unfolded term weight: " << vars[ii]
			// << " = " << vars[ii] * wt_unweighted_unfolded_energy << std::endl;
		}
	}

	// !!!!!!!!!!!!!!!!!  fixed_score_list[ ii ] gives you a ScoreType - that's not the weight.  fixed_terms[ fixed_score_list[ ii ] ]
	// is the actual weight.  But if you use the ScoreType in a multiplication, it will work just fine because ScoreTypes are
	// an enum.  They're position ~120 in the enum, so it'll be a nice large weight.

	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		if ( name_from_score_type( fixed_score_list[ ii ] ) == "unfolded" ) {
			wt_weighted_unfolded_energy = fixed_terms[ fixed_score_list[ ii ] ] * wt_unweighted_unfolded_energy;
			mut_weighted_unfolded_energy = fixed_terms[ fixed_score_list[ ii ] ] * mut_unweighted_unfolded_energy;
			unfolded_weight = fixed_terms[ fixed_score_list[ ii ] ];
			//TR << "process_score(): weighting unweighted unfolded energy: '" << wt_unweighted_unfolded_energy
			// << " by fixed unfolded term weight: " << fixed_terms[ fixed_score_list[ ii ] ]
			// << " = " << fixed_terms[ fixed_score_list[ ii ] ] * wt_unweighted_unfolded_energy << std::endl;
		}
	}

	// Now SUBTRACT this unfolded energy from the sums we have so far. See comments for optimize_weights() in IterativeOptEDriver.cc.
	// This weighted unfolded energy will be the same for every structure in the SSD vectors wts_ and muts_. So just iterate over
	// both of those and subtract out this unfolded energy.
	for ( Size ii = 1; ii <= wts_.size(); ++ii ) { wt_energies[ ii ] -= wt_weighted_unfolded_energy; }
	for ( Size ii = 1; ii <= muts_.size(); ++ii ) { mut_energies[ ii ] -= mut_weighted_unfolded_energy; }

	//TR << "process_score(): unf weight: " << unfolded_weight
	//  << ", unweighted energy: " << wt_unweighted_unfolded_energy << ", " << mut_unweighted_unfolded_energy
	// << "; weighted unfolded energy: " << F(7,2,wt_weighted_unfolded_energy) << ", "  << F(7,2,mut_weighted_unfolded_energy) << std::endl;

	//TR << "process_score(): after unfolded wts_: [ ";
	//for ( Size jj = 1; jj <= wts_.size(); ++jj ) { TR << F(6,1,wt_energies[ jj ]) << ", "; }
	//TR << "]" << std::endl;


	// This is where we branch on how the score is calculated.  The simplest approach is to take the minimum energy
	// of all the wts and all the muts and subtract them to get the ddG.  The mean-based approach use the difference
	// of the average of all muts and average of all wts. Finally, the boltzmann approach calculates a boltzmann
	// probability for the wts and muts to get a score.
	// The other ways have been removed.

	Real predicted_ddG( 0.0 );
	Real ddG_diff( 0.0 );

	// Do things the old-fashioned way: best energy mut - best energy wt
	Size const best_wt = arg_min( wt_energies );
	Size const best_mut = arg_min( mut_energies );

	Real const best_wt_energy = wt_energies[ best_wt ];
	Real const best_mut_energy = mut_energies[ best_mut ];

	predicted_ddG = best_mut_energy - best_wt_energy;
	ddG_diff = predicted_ddG - experimental_ddG_;

	if ( print ) {
		ostr << "DDG " << A( 20, tag() ) << X(1) << "pred: " << F(6,2,predicted_ddG) << " exp: " << F(6,2,experimental_ddG_)
			<< " diff^2: " << F(7,2,ddG_diff*ddG_diff) << " cmptwt_diff^2: " << F( 7,2,component_weights[ ddG_mutation_correlation_with_unfolded_energy ] * ddG_diff * ddG_diff )
			<< std::endl;

		TR << "process_score(): unf weight: " << F(5,3,unfolded_weight)
			<< ", raw unfE: " << wt_unweighted_unfolded_energy << ", " << mut_unweighted_unfolded_energy
			<< "; unfE: " << wt_weighted_unfolded_energy << ", "  << mut_weighted_unfolded_energy
			<< "; totalE: " << best_wt_energy << ", " << best_mut_energy
			<< "; pred: " << ObjexxFCL::format::F(4,2,predicted_ddG)
			<< ", exp'tal: " << experimental_ddG_ << ", ddG_diff^2: " << ObjexxFCL::format::F(6,3,ddG_diff * ddG_diff)
			<< ", tag: " << tag() << std::endl;

	} else {

		for ( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {

			if ( ( score_list[ e_dof ] == fa_rep ) && ( muts_[ best_mut ]->free_data()[ e_dof ] - wts_[ best_wt ]->free_data()[ e_dof ] ) > 10 ) {
				// deal with the really bad repulsive energy cases here
				dE_dvars[ e_dof ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff * 10;

			} else if ( score_list[ e_dof ] == unfolded ) {
				dE_dvars[ e_dof ] += 2 * component_weights[ ddG_mutation_correlation_with_unfolded_energy ] * ddG_diff *
					( mut_unweighted_unfolded_energy - wt_unweighted_unfolded_energy );

			} else {
				dE_dvars[ e_dof ] += 2 * component_weights[ ddG_mutation_correlation_with_unfolded_energy ] * ddG_diff *
					( muts_[ best_mut ]->free_data()[ e_dof ] - wts_[ best_wt ]->free_data()[ e_dof ] );
			}
		}

		if ( num_ref_dofs != 0 ) {
			dE_dvars[ num_energy_dofs + mut_aa_ ] += 2 * component_weights[ ddG_mutation_correlation_with_unfolded_energy ] * ddG_diff;
			dE_dvars[ num_energy_dofs + wt_aa_  ] -= 2 * component_weights[ ddG_mutation_correlation_with_unfolded_energy ] * ddG_diff;
		}
	}

	return component_weights[ ddG_mutation_correlation_with_unfolded_energy ] * ddG_diff * ddG_diff;
}

///
OptEPositionDataType
NestedEnergyTermDDGMutationOptEData::type() const {
	return ddG_mutation_correlation_with_unfolded_energy;
}

///
/// Only used for user feedback. Nothing in the code uses the result from this to allocate memory.
///
Size
NestedEnergyTermDDGMutationOptEData::memory_use() const {

	Size total = sizeof( DDGMutationOptEData ) +
		sizeof( SingleStructureData ) * wts_.size() +
		sizeof( SingleStructureData ) * muts_.size();
	if ( wts_.size() > 0 ) {
		total += sizeof( Real ) * ( wts_[ 1 ]->free_data().size() + wts_[ 1 ]->fixed_data().size() ) * wts_.size();
	}
	if ( muts_.size() > 0 ) {
		total += sizeof( Real ) * ( muts_[ 1 ]->free_data().size() + muts_[ 1 ]->fixed_data().size() ) * muts_.size();
	}

	// the emap uses some memory, too!  do I take the number of score types and multiply by sizeof(Real)?  Or can I just
	// do sizeof(emap)?  I think taking the size of the class is enough.  Emaps don't grow/shrink dynamically.  They're
	// essentially constant size, an array of Reals sized to n_score_types.
	total += sizeof( EnergyMap ) * 2; // there's the wt one and the mutant one

	return total;
}


#ifdef USEMPI
///
///
void
NestedEnergyTermDDGMutationOptEData::send_to_node( int const destination_node, int const tag ) const {

	/// 1. Experimental DDG, wt_aa, mut_aa
	int wt_aa( wt_aa_ ), mut_aa( mut_aa_ );
	Real experimental_ddG = experimental_ddG_; // stupid const pointer
	MPI_Send( & experimental_ddG, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );
	MPI_Send( & wt_aa, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );
	MPI_Send( & mut_aa, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	// 2a. n natives
	//std::cout << "sending nwts to node " << destination_node << std::endl;
	Size nwts = wts_.size();
	MPI_Send( & nwts, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 2b. n decoys
	//std::cout << "sending nmuts to node " << destination_node << std::endl;
	Size nmuts = muts_.size();
	MPI_Send( & nmuts, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	if ( nwts == 0 || nmuts == 0 )
		return;

	/// 3. n free
	Size n_free = muts_[ 1 ]->free_data().size();
	//std::cout << "sending n_free to node " << destination_node << " " << n_free << std::endl;
	MPI_Send( & n_free, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 4. n fixed
	Size n_fixed = muts_[ 1 ]->fixed_data().size();
	//std::cout << "sending n_fixed to node " << destination_node  << " " << n_fixed << std::endl;
	MPI_Send( & n_fixed, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// Send natives, then send decoys
	Real * free_data = new Real[ n_free * nwts ];
	Real * fixed_data = new Real[ n_fixed * nwts ];
	for ( Size ii = 1; ii <= nwts; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ] = wts_[ ii ]->free_data()[ jj ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ] = wts_[ ii ]->fixed_data()[ jj ];
		}
	}

	//std::cout << "sending native free_data to node " << destination_node << " " << free_data <<  std::endl;
	/// 5. native free data
	MPI_Send( free_data, nwts * n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//std::cout << "sending native fixed_data to node " << destination_node << " " << fixed_data <<  std::endl;
	/// 6. fixed data
	MPI_Send( fixed_data, nwts * n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//std::cout << "Sent -- about to delete data" << std::endl;

	/// now send decoys
	Real * decoy_free_data = new Real[ n_free * nmuts ];
	Real * decoy_fixed_data = new Real[ n_fixed * nmuts ];
	for ( Size ii = 1; ii <= nmuts; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			decoy_free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ] = muts_[ ii ]->free_data()[ jj ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			decoy_fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ] = muts_[ ii ]->fixed_data()[ jj ];
		}
	}
	/// 7. decoy free data
	//std::cout << "sending decoy free_data to node " << destination_node << std::endl;
	MPI_Send( decoy_free_data, nmuts * n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 8. decoy fixed data
	//std::cout << "sending decoy fixed_data to node " << destination_node << std::endl;
	MPI_Send( decoy_fixed_data, nmuts * n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	delete [] free_data;
	delete [] fixed_data;

	delete [] decoy_free_data;
	delete [] decoy_fixed_data;

	/// 9a. The energies in the unfolded emap, wt first
	Real * wt_unfolded_energies = new Real[ scoring::n_score_types ];
	for ( Size ee = 1; ee <= scoring::n_score_types; ++ee ) {
		wt_unfolded_energies[ ee - 1 ] = wt_unfolded_energies_emap_[ (ScoreType) ee ];
	}
	MPI_Send( wt_unfolded_energies, scoring::n_score_types, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );
	delete [] wt_unfolded_energies;
	wt_unfolded_energies = 0;

	/// 9b. mutant emap
	Real * mut_unfolded_energies = new Real[ scoring::n_score_types ];
	for ( Size ee = 1; ee <= scoring::n_score_types; ++ee ) {
		mut_unfolded_energies[ ee - 1 ] = mut_unfolded_energies_emap_[ (ScoreType) ee ];
	}
	MPI_Send( mut_unfolded_energies, scoring::n_score_types, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );
	delete [] mut_unfolded_energies;
	mut_unfolded_energies = 0;


	OptEPositionData::send_to_node( destination_node, tag );

}

///
void
NestedEnergyTermDDGMutationOptEData::receive_from_node( int const source_node, int const tag )
{
	MPI_Status stat;
	//TR << "PNatStructureOptEData::Recieving data from node... " << source_node << std::endl;

	/// 1. Experimental DDG, wt_aa, mut_aa
	int wt_aa( 0 ), mut_aa(0);
	MPI_Recv( & experimental_ddG_, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	MPI_Recv( & wt_aa, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	MPI_Recv( & mut_aa, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	wt_aa_  = static_cast< AA > ( wt_aa );
	mut_aa_ = static_cast< AA > ( mut_aa );

	/// 2a. n wts
	Size nwts( 0 );
	MPI_Recv( & nwts, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 2b. n decoys
	Size nmuts( 0 );
	MPI_Recv( & nmuts, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	if ( nwts == 0 || nmuts == 0 ) return;
	wts_.reserve( nwts );
	muts_.reserve( nmuts );

	/// 3. n free
	Size n_free( 0 );
	MPI_Recv( & n_free, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 4. n fixed
	Size n_fixed( 0 );
	MPI_Recv( & n_fixed, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// Recieve native data first, then decoys
	Real * free_data = new Real[ n_free * nwts ];
	Real * fixed_data = new Real[ n_fixed * nwts ];

	/// 5. free data
	MPI_Recv( free_data, nwts * n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 6. fixed data
	MPI_Recv( fixed_data, nwts * n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	utility::vector1< Real > free_data_v( n_free );
	utility::vector1< Real > fixed_data_v( n_fixed );
	for ( Size ii = 1; ii <= nwts; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data_v[ jj ] = free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data_v[ jj ] = fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ];
		}
		wts_.push_back( SingleStructureDataOP( new SingleStructureData( free_data_v, fixed_data_v ) ) );
	}


	delete [] free_data;// free_data = 0;
	delete [] fixed_data;// fixed_data = 0;

	//// Now receive decoy data
	free_data = new Real[ n_free * nmuts ];
	fixed_data = new Real[ n_fixed * nmuts ];

	/// 5. free data
	MPI_Recv( free_data, nmuts * n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 6. fixed data
	MPI_Recv( fixed_data, nmuts * n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	for ( Size ii = 1; ii <= nmuts; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data_v[ jj ] = free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data_v[ jj ] = fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ];
		}
		muts_.push_back( SingleStructureDataOP( new SingleStructureData( free_data_v, fixed_data_v ) ) );
	}

	delete [] free_data;
	delete [] fixed_data;


	/// 7a. Unfolded energies
	EnergyMap emap;

	Real * wt_unfolded_energies = new Real[ scoring::n_score_types ];
	MPI_Recv( wt_unfolded_energies, scoring::n_score_types, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	for ( Size ee = 1; ee <= scoring::n_score_types; ++ee ) {
		// be careful because the array is 0-based while the score type enum is 1-based
		emap[ (ScoreType) ee ] = wt_unfolded_energies[ ee - 1 ];
	}
	set_wt_unfolded_energies_emap( emap );

	/// 7b. Unfolded energies
	emap.zero();
	Real * mut_unfolded_energies = new Real[ scoring::n_score_types ];
	MPI_Recv( mut_unfolded_energies, scoring::n_score_types, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	for ( Size ee = 1; ee <= scoring::n_score_types; ++ee ) {
		emap[ (ScoreType) ee ] = mut_unfolded_energies[ ee - 1 ];
	}
	set_mut_unfolded_energies_emap( emap );

	delete [] wt_unfolded_energies;
	delete [] mut_unfolded_energies;

	wt_unfolded_energies = 0;
	mut_unfolded_energies = 0;


	OptEPositionData::receive_from_node( source_node, tag );

}
#endif


}
}
