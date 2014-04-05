// -*- Mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/optimization/OptEData.cc
/// @author  ashworth
/// @author  Jim Havranek
/// @author  Andrew Leaver-Fay

#ifdef USEMPI
#include <mpi.h>
#endif


// Unit headers
#include <protocols/optimize_weights/OptEData.hh>
#include <protocols/optimize_weights/NestedEnergyTermOptEData.hh>
#include <protocols/optimize_weights/PNatLigPoseOptEData.hh>
#include <protocols/optimize_weights/DGBindOptEData.hh>
#include <protocols/optimize_weights/DDGBindOptEData.hh>
#ifdef USEMPI
#include <protocols/optimize_weights/IterativeOptEDriver.hh>
#endif

// Project headers
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
// AUTO-REMOVED #include <basic/options/util.hh>

#include <utility/LexicographicalIterator.hh>
#include <utility/string_util.hh>
#include <utility/vector1.functions.hh>
#include <utility/exit.hh>

#include <ObjexxFCL/format.hh>

#include <numeric/numeric.functions.hh>
#include <numeric/statistics/functions.hh>

#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <cmath>

#include <basic/Tracer.hh>

// option key includes
#include <basic/options/keys/optE.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/options/option.hh>



using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace scoring;
using namespace ObjexxFCL::format;

namespace protocols {
namespace optimize_weights {

/// @details Auto-generated virtual destructor
OptEData::~OptEData() {}

/// @details Auto-generated virtual destructor
PNatRotOptERotamerData::~PNatRotOptERotamerData() {}

/// @details Auto-generated virtual destructor
PNatAAOptERotamerData::~PNatAAOptERotamerData() {}

static basic::Tracer TR("protocols.optimize_weights.OptEData");

///@begin (ostream operator for OptERotamerDataOP)
///@author ashworth
std::ostream & operator << ( std::ostream & os, PNatAAOptERotamerDataOP rd )
{
	os << rd->rot_number() << "," << rd->this_aa() << ",";

	utility::vector1< Real > const & fixed_data( rd->fixed_data() );
	for ( utility::vector1< Real >::const_iterator fd( fixed_data.begin() );
	      fd != fixed_data.end(); ++fd ) {
		if ( fd != fixed_data.begin() ) os << " ";
		os << *fd;
	}
	os << ",";

	utility::vector1< Real > const & data( rd->data() );
	for ( utility::vector1< Real >::const_iterator d( data.begin() );
	      d != data.end(); ++d ) {
		if ( d != data.begin() ) os << " ";
		os << *d;
	}
	return os;
}

//
// ------------------- OptEPositionData -----------------------//
//

OptEPositionData::OptEPositionData()
{}

OptEPositionData::~OptEPositionData()
{}

void
OptEPositionData::update_range(
	SingleStructureDataCOP structure,
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list,
	EnergyMap & lower_bound,
	EnergyMap & upper_bound
) const
{
	for ( Size ii = 1; ii <= free_score_list.size(); ++ii ) {
		Real ii_min = lower_bound[ free_score_list[ ii ] ];
		Real ii_max = upper_bound[ free_score_list[ ii ] ];
		if ( ii_min > structure->free_data()[ ii ] ) {
			ii_min = structure->free_data()[ ii ];
		}
		if ( ii_max < structure->free_data()[ ii ] ) {
			ii_max = structure->free_data()[ ii ];
		}
		lower_bound[ free_score_list[ ii ] ] = ii_min;
		upper_bound[ free_score_list[ ii ] ] = ii_max;
	}

	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		Real ii_min = lower_bound[ fixed_score_list[ ii ] ];
		Real ii_max = upper_bound[ fixed_score_list[ ii ] ];
		if ( ii_min > structure->fixed_data()[ ii ] ) {
			ii_min = structure->fixed_data()[ ii ];
		}
		if ( ii_max < structure->fixed_data()[ ii ] ) {
			ii_max = structure->fixed_data()[ ii ];
		}
		lower_bound[ fixed_score_list[ ii ] ] = ii_min;
		upper_bound[ fixed_score_list[ ii ] ] = ii_max;
	}
}

void
OptEPositionData::tag( std::string const & tag_in )
{
	tag_ = tag_in;
}

std::string const &
OptEPositionData::tag() const
{
	return tag_;
}

#ifdef USEMPI
void
OptEPositionData::send_to_node(
	int const destination_node,
	int const tag
) const
{
	int len( tag_.size() );
	MPI_Send( &len, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );
	MPI_Send( const_cast< char * > (tag_.c_str()), len, MPI_CHAR, destination_node, tag, MPI_COMM_WORLD );
}

void
OptEPositionData::receive_from_node(
	int const source_node,
	int const tag
)
{
	MPI_Status stat;

	int len( 0 );
	MPI_Recv( &len, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, & stat );
	char * str = new char[ len + 1 ];
	str[ len ] = '\0'; // ? do I need null terminated strings?
	MPI_Recv( str, len, MPI_CHAR, source_node, tag, MPI_COMM_WORLD, & stat );
	std::string tag_from_node( str, len );
	delete [] str;
	tag_ = tag_from_node;
}
#endif


//
// ------------------- PNatAAOptEPositionData -----------------------//
//

PNatAAOptEPositionData::PNatAAOptEPositionData()
{}

PNatAAOptEPositionData::~PNatAAOptEPositionData()
{}

/// Does actual work for OptE minimization
/// @details Determine the metric for optimization of energy weights.  Original source
/// is Brian Kuhlman's FORTRAN weight training program.  Since the goal is to
/// maximize site-by-site sequence recovery, the metric is decomposable by
/// positions, and this does the work for one position.  First, a pass is made
/// through all the rotamers at the position, and the lowest energy rotamer for
/// each amino acid is identified.  Next, a rough partition function is constructed.
/// The score to be minimized is the negative log of the probabilty of selecting
/// a native aa rotamer (needn't be the native _conformation_).
/// @details Limitations:
/// 1.  Assumes that the choices available at each position are the 20 canonical
/// amino acids - no DNA base recovery, no pTyr.
/// 2.  Doesn't have Chris Saunder's compositional constraints.  This would
/// require implementation of a direction-set ( aka powell ) minimizer, as his
/// optimization metric doesn't yield gradients.
//
Real
PNatAAOptEPositionData::get_score(
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

	Size const aa_range( chemical::num_canonical_aas );
	chemical::AA const this_native_aa( native_aa() );

	static Real const inv_kT( option[ optE::inv_kT_nataa ] );

	// contains the best energy for each amino acid at this position
	utility::vector1< Real > best_energy_by_aa( aa_range, 1000.0 );

	// containers for derivatives
	Multivec ref_deriv_weight( aa_range, 0.0 );
	utility::vector1< Real > dummy_set( score_list.size(), 0.0 );
	utility::vector1< utility::vector1< Real > > unweighted_E_dof( aa_range, dummy_set );

	process_rotamers(
		vars, num_energy_dofs, fixed_terms, score_list, fixed_score_list, aa_range, dummy_set,
		best_energy_by_aa, unweighted_E_dof, ref_deriv_weight );

	//Real const very_best = utility::min( best_energy_by_aa );
	//for ( Size ii = 1; ii <= best_energy_by_aa.size(); ++ii ) {
	//	best_energy_by_aa[ ii ] -= very_best;
	//}


	// now do the partition function analysis
	// reference energy deriv

	Real numerator(0.0), partition(0.0);
	Multivec dpartition( vars.size(), 0.0 ), dnumerator( vars.size(), 0.0 );

	for( Size aa(1); aa <= aa_range; ++aa ) {

		Real const exp_term( std::exp( -1.0 * inv_kT * best_energy_by_aa[ aa ] ) );
		partition += exp_term;
		if ( aa == size_t(this_native_aa) )
			numerator = exp_term;

		// for reference energy derivatives, but don't assume the protocol is using them
		// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...) but as this is an energy, the 'weight' here is 0 or 1
		if ( num_ref_dofs != 0 ) {
			Real const ref_deriv_term( -1.0 * inv_kT * ref_deriv_weight[ aa ] * exp_term );
			dpartition[ num_energy_dofs + aa ] = ref_deriv_term;
			if ( aa == size_t(this_native_aa) )
				dnumerator[ num_energy_dofs + aa ] = ref_deriv_term;
		}

		// partitions for energy derivatives
		// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
		for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
			Real e_dof_deriv( -1.0 * inv_kT * unweighted_E_dof[ aa ][ e_dof ] * exp_term );
			dpartition[ e_dof ] += e_dof_deriv;
			if ( aa == size_t(this_native_aa) )
				dnumerator[ e_dof ] = e_dof_deriv;
		}
	}

	// accumulate to passed-in derivative sums
	for ( Size dof(1); dof <= vars.size(); ++dof ) {
		dE_dvars[ dof ] += component_weights[ type() ] * ( dpartition[ dof ] / partition - dnumerator[ dof ] / numerator );
	}

	return ( -1.0 * component_weights[ type() ] * std::log( numerator / partition ) );
}

void
PNatAAOptEPositionData::print_score(
	std::ostream & ostr,
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
	//TR <<  "PNatAAOptEPositionData::print_score" << std::endl;
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Size const aa_range( chemical::num_canonical_aas );
	chemical::AA const this_native_aa( native_aa() );


	static Real const inv_kT( option[ optE::inv_kT_nataa ] );

	// contains the best energy for each amino acid at this position
	utility::vector1< Real > best_energy_by_aa( aa_range, 1000.0 );

	// containers for derivatives
	Multivec ref_deriv_weight( aa_range, 0.0 );
	utility::vector1< Real > dummy_set( score_list.size(), 0.0 );
	utility::vector1< utility::vector1< Real > > unweighted_E_dof( aa_range, dummy_set );

	process_rotamers(
		vars, num_energy_dofs, fixed_terms, score_list, fixed_score_list, aa_range, dummy_set,
		best_energy_by_aa, unweighted_E_dof, ref_deriv_weight );

	// now do the partition function analysis
	// reference energy deriv

	Real numerator(0.0), partition(0.0);
	Multivec dpartition( vars.size(), 0.0 ), dnumerator( vars.size(), 0.0 );

	for( Size aa(1); aa <= aa_range; ++aa ) {

		Real const exp_term( std::exp( -1.0 * inv_kT * best_energy_by_aa[ aa ] ) );
		partition += exp_term;
		if ( aa == size_t(this_native_aa) ) numerator = exp_term;

		// for reference energy derivatives, but don't assume the protocol is using them
		// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...) but as this is an energy, the 'weight' here is 0 or 1
		if ( num_ref_dofs != 0 ) {
			Real const ref_deriv_term( -1.0 * inv_kT * ref_deriv_weight[ aa ] * exp_term );
			dpartition[ num_energy_dofs + aa ] = ref_deriv_term;
			if ( aa == size_t(this_native_aa) ) dnumerator[ num_energy_dofs + aa ] = ref_deriv_term;
		}

		// partitions for energy derivatives
		for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
			// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
			Real e_dof_deriv( -1.0 * inv_kT * unweighted_E_dof[ aa ][ e_dof ] * exp_term );
			dpartition[ e_dof ] += e_dof_deriv;
			if ( aa == size_t(this_native_aa) ) dnumerator[ e_dof ] = e_dof_deriv;
		}
	}

	// accumulate to passed-in derivative sums
	for ( Size dof(1); dof <= vars.size(); ++dof ) {
		dE_dvars[ dof ] += dpartition[ dof ] / partition - dnumerator[ dof ] / numerator;
	}

	ostr << "PNATAA " << tag() << X(1) << this_native_aa << "," << I(2, (int)this_native_aa) << X(1)
			<< " nbs: " << I(2,neighbor_count())
			<< " num: " << F(7,3,numerator) << " part: " << F(7,3,partition)
			<< " p: " << F(7,5,numerator / partition)
			<< " -lnp: " << F(6,4,-1.0 * std::log( numerator / partition ))
			<< " -compwt_lnp: " << F(6, 4, component_weights[ type() ] * (-1.0 * std::log( numerator / partition )) )
			<< std::endl;

	//ostr << "PNATAA " << this_native_aa << " " << (int) this_native_aa << " ";
	//ostr << " nneighb: " << neighbor_count_ << " p: " << numerator / partition << " -lnp: ";
	//ostr << -1.0 * std::log( numerator / partition );
	//ostr << " -compwt_lnp: " << component_weights[ type() ] * -1.0 * std::log( numerator / partition );
	//ostr << " " << tag() << "\n";

	//TR << "PNatAAOptEPositionData::print_score done" << std::endl;

	return;// ( -1.0 * std::log( numerator / partition ) );
}

void
PNatAAOptEPositionData::range(
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list,
	EnergyMap & lower_bound,
	EnergyMap & upper_bound
) const
{
	for ( Size ii = 1; ii <= free_score_list.size(); ++ii ) {
		Real ii_min = lower_bound[ free_score_list[ ii ] ];
		Real ii_max = upper_bound[ free_score_list[ ii ] ];
		for( PNatAAOptERotamerDataOPs::const_iterator iter = rotamer_data_begin(),
				e_itr = rotamer_data_end() ; iter != e_itr; ++iter ) {
			if ( ii_min > (*iter)->data()[ ii ] ) {
				ii_min = (*iter)->data()[ ii ];
			}
			if ( ii_max < (*iter)->data()[ ii ] ) {
				ii_max = (*iter)->data()[ ii ];
			}
		}
		lower_bound[ free_score_list[ ii ] ] = ii_min;
		upper_bound[ free_score_list[ ii ] ] = ii_max;
	}

	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		Real ii_min = lower_bound[ fixed_score_list[ ii ] ];
		Real ii_max = upper_bound[ fixed_score_list[ ii ] ];
		for( PNatAAOptERotamerDataOPs::const_iterator iter = rotamer_data_begin(),
				e_itr = rotamer_data_end() ; iter != e_itr; ++iter ) {
			if ( ii_min > (*iter)->fixed_data()[ ii ] ) {
				ii_min = (*iter)->fixed_data()[ ii ];
			}
			if ( ii_max < (*iter)->fixed_data()[ ii ] ) {
				ii_max = (*iter)->fixed_data()[ ii ];
			}
		}
		lower_bound[ fixed_score_list[ ii ] ] = ii_min;
		upper_bound[ fixed_score_list[ ii ] ] = ii_max;
	}
}

void
PNatAAOptEPositionData::process_rotamers(
	optimization::Multivec const & vars,
	Size const num_energy_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const &,
	ScoreTypes const & fixed_score_list,
	Size const ,
	utility::vector1< Real > const & dummy_set,
	utility::vector1< Real > & best_energy_by_aa,
	utility::vector1< utility::vector1< Real > > & unweighted_E_dof,
	optimization::Multivec & ref_deriv_weight
) const
{
	for( PNatAAOptERotamerDataOPs::const_iterator itr = rotamer_data_begin(), e_itr = rotamer_data_end() ; itr != e_itr; ++itr ) {
		int const this_aa( (*itr)->this_aa() );

		Real weighted_energy( 0.0 );
		// Variable-weighted energy terms
		for( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
			weighted_energy += vars[ ii ] * ((**itr)[ ii ]);
		}
		// Fixed-weight energy terms
		for( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
			weighted_energy += fixed_terms[ fixed_score_list[ ii ] ] * (*itr)->fixed_data()[ ii ];
		}
		// Reference energy term
		// Don't count this unless we're definitely using reference energies. If the vars vector has more values in it than the
		// number of energy dofs, that means reference energies are in use.  Prob. a better way to do this, but this should work.
		if ( vars.size() > num_energy_dofs ) {
			weighted_energy += vars[ num_energy_dofs + this_aa ];
		}

		Real const cutoff( 300.0 );
		if( weighted_energy > cutoff ) weighted_energy = cutoff;
		if( weighted_energy < -1.0*cutoff ) weighted_energy = -1.0*cutoff;

		if( weighted_energy < best_energy_by_aa[ this_aa ] ) {

			best_energy_by_aa[ this_aa ] = weighted_energy;

			//if (true) {
			if( std::abs( weighted_energy ) < cutoff ) {
				unweighted_E_dof[ this_aa ] = (*itr)->data();
				ref_deriv_weight[ this_aa ] = 1.0;
			} else {
				// derivatives will be zero above and below cutoffs
				unweighted_E_dof[ this_aa ] = dummy_set;
				ref_deriv_weight[ this_aa ] = 0.0;
			}
		}
	} // done processing rotamers
}


OptEPositionDataType
PNatAAOptEPositionData::type() const
{
	return prob_native_amino_acid;
}


void
PNatAAOptEPositionData::write_to_file( std::ofstream & outfile ) const
{
	outfile
		<< "position " << position() << " "
		<< "nataa " << native_aa() << " "
		<< "neighbor_count " << neighbor_count() << " "
		<< "nrots " << data_.size() << "\n";
	for ( PNatAAOptERotamerDataOPs::const_iterator rot( rotamer_data_begin() );
			rot != rotamer_data_end(); ++rot ) {
		outfile << *rot << std::endl;
	}

}

void
PNatAAOptEPositionData::read_from_file( std::ifstream & infile )
{
	typedef utility::vector1< std::string > Strings;
	using namespace utility;

	// read first line with position, native aa, neighbor_count, and num_rotamers data
	std::string line;
	getline( infile, line );
	Strings words( string_split( line, ' ' ) );
	runtime_assert( words[ 1 ] == "position" );
	position_ = from_string( words[ 2 ], Size( 0 ) );
	runtime_assert( words[ 3 ] == "nataa" );
	native_aa_ = chemical::aa_from_name( words[ 4 ] );
	runtime_assert( words[ 5 ] == "neighbor_count" );
	neighbor_count_ = from_string( words[ 6 ], Size( 0 ) );
	runtime_assert( words[ 7] == "nrots" );
	Size num_rotamers = from_string( words[ 8 ], Size( 0 ) );

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
		Strings
			fixed_vals( string_split( sections[3], ' ' ) ),
			free_vals( string_split( sections[4], ' ' ) );
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
		runtime_assert( !energies.empty() );
		PNatAAOptERotamerDataOP new_rot_data = new PNatAAOptERotamerData( aa, rotnum, energies, fixed_energies );
		add_rotamer_line_data( new_rot_data );
	}
}

void
PNatAAOptEPositionData::write_to_binary_file( std::ofstream & outfile ) const
{
	typedef utility::vector1< Real > Energies;
	typedef Energies::const_iterator Energies_CItr;

	// compiler wanted these temporary variables instead of calls to class accessor functions,
	// which produced the error "non-lvalue in unary `&'"
	//Size const position( (position_ ),
	//			  neighbor_count( neighbor_count() );
	//chemical::AA native_aa( native_aa() );
	outfile.write( (char*) & position_, sizeof(Size) );
	outfile.write( (char*) & native_aa_, sizeof(chemical::AA) );
	outfile.write( (char*) & neighbor_count_, sizeof(Size) );

	Size const nrotamers( data().size() );
	outfile.write( (char*) &nrotamers, sizeof(Size) );
	for ( PNatAAOptERotamerDataOPs::const_iterator rot( data().begin() );
			rot != data().end(); ++rot ) {
		chemical::AA this_aa( (*rot)->this_aa() );
		Size rot_number( (*rot)->rot_number() );
		outfile.write( (char*) &this_aa, sizeof(chemical::AA) );
		outfile.write( (char*) &rot_number, sizeof(Size) );

		Size const nfree( (*rot)->data().size() );
		outfile.write( (char*) &nfree, sizeof(Size) );
		for ( Energies_CItr energy( (*rot)->data().begin() );
				energy != (*rot)->data().end(); ++energy ) {
			outfile.write( (char*) &(*energy), sizeof(Real) );
		}

		Size const nfixed( (*rot)->fixed_data().size() );
		outfile.write( (char*) &nfixed, sizeof(Size) );
		for ( Energies_CItr fixed_energy( (*rot)->fixed_data().begin() );
				fixed_energy != (*rot)->fixed_data().end(); ++fixed_energy ) {
			outfile.write( (char*) &(*fixed_energy), sizeof(Real) );
		}
	}


}


void
PNatAAOptEPositionData::read_from_binary_file( std::ifstream & infile )
{
	typedef utility::vector1< Real > Energies;
	//typedef Energies::const_iterator Energies_CItr;

	//Size position(0), neighbor_count(0);
	//chemical::AA nataa( chemical::aa_unk );
	infile.read( (char*) & position_, sizeof(Size) );
	infile.read( (char*) & native_aa_, sizeof(chemical::AA) );
	infile.read( (char*) & neighbor_count_, sizeof(Size) );
	//pos_data->set_position( position );
	//pos_data->set_neighbor_count( neighbor_count );
	//pos_data->set_native_aa( nataa );

	Size nrotamers(0);
	infile.read( (char*) &nrotamers, sizeof(Size) );
	for ( Size j(1); j <= nrotamers; ++j ) {

		chemical::AA rot_aa( chemical::aa_unk );
		Size rot_number(0);
		infile.read( (char*) &rot_aa, sizeof(chemical::AA) );
		infile.read( (char*) &rot_number, sizeof(Size) );

		Energies energies, fixed_energies;
		Size nfree(0);
		infile.read( (char*) &nfree, sizeof(Size) );
		for ( Size k(1); k <= nfree; ++k ) {
			Real energy(0.0);
			infile.read( (char*) &energy, sizeof(Real) );
			energies.push_back( energy );
		}

		Size nfixed(0);
		infile.read( (char*) &nfixed, sizeof(Size) );
		for ( Size k(1); k <= nfixed; ++k ) {
			Real fixed_energy(0.0);
			infile.read( (char*) &fixed_energy, sizeof(Real) );
			fixed_energies.push_back( fixed_energy );
		}
		runtime_assert( !energies.empty() );

		PNatAAOptERotamerDataOP rot_data =
			new PNatAAOptERotamerData( rot_aa, rot_number, energies, fixed_energies );
		add_rotamer_line_data( rot_data );
	}
}

Size
PNatAAOptEPositionData::memory_use() const
{
	Size total = sizeof( PNatAAOptEPositionData ) +
		sizeof( PNatAAOptERotamerData ) * data_.size();
	if ( data_.size() > 0 ) {
		total +=  sizeof( Real ) * ( data_[ 1 ]->data().size() + data_[ 1 ]->fixed_data().size() ) * data_.size();
	}
	return total;
}



#ifdef USEMPI

void
PNatAAOptEPositionData::send_to_node( int const destination_node, int const tag ) const
{

	//TR << "PNatAAOptEPositionData::Sending data to node " << destination_node << std::endl;

	/// 3. Which position is this?
	int ii_pos = position();
	MPI_Send( & ii_pos, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 4. What is the native amino acid at this position?
	int ii_aa  = native_aa();
	MPI_Send( & ii_aa, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 5. How many neighbors did this position have?
	int ii_neighbor_count = neighbor_count();
	MPI_Send( & ii_neighbor_count, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 6. The number of rotamers for this position
	Size ii_num_rotamers = size();
	MPI_Send( & ii_num_rotamers, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	if ( ii_num_rotamers == 0 ) return;

	Size free_count = data_[1]->data().size();
	Size fixed_count = data_[1]->fixed_data().size();

	/// 6b The size of the free and fixed data, since that context is not available.
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

	//TR << "...data sent to node " << destination_node << std::endl;
	OptEPositionData::send_to_node( destination_node, tag );
}


void
PNatAAOptEPositionData::receive_from_node( int const source_node, int const tag )
{
	MPI_Status stat;

	//TR << "Recieving data from node... " << source_node << std::endl;

	/// 3. Which sequence position is this?
	int ii_pos;
	MPI_Recv( & ii_pos, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	set_position( ii_pos );

	/// 4. What is the native amino acid at this position?
	int ii_aa;
	MPI_Recv( & ii_aa, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	set_native_aa( (chemical::AA) ii_aa );

	/// 5. How many neighbors did this position have?
	int ii_neighbor_count;
	MPI_Recv( & ii_neighbor_count, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	set_neighbor_count( ii_neighbor_count );

	/// 6. The number of rotamers for this position
	Size ii_num_rotamers;
	MPI_Recv( & ii_num_rotamers, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	if ( ii_num_rotamers == 0 ) return;

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
		PNatAAOptERotamerDataOP jj_rotamer_data = new PNatAAOptERotamerData(
			(chemical::AA ) ii_aa_types[ jj - 1 ],
			ii_rot_nums[ jj - 1 ],
			free_data_vect,
			fixed_data_vect );
		add_rotamer_line_data( jj_rotamer_data );

	}

	delete [] ii_aa_types; ii_aa_types = 0;
	delete [] ii_rot_nums; ii_rot_nums = 0;
	delete [] free_data;   free_data   = 0;
	delete [] fixed_data;  fixed_data  = 0;

	OptEPositionData::receive_from_node( source_node, tag );

	//TR << "...data received from node " << source_node << std::endl;
}

#endif


//
// ------------------- Position Specific Scoring Matrix -----------------------//
//

PSSMOptEPositionData::PSSMOptEPositionData() {}

PSSMOptEPositionData::~PSSMOptEPositionData() {}


Real
PSSMOptEPositionData::get_score(
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const & score_list,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static Real const inv_kT( option[ optE::inv_kT_nataa ] );

	Size const aa_range( chemical::num_canonical_aas );

	utility::vector1< Real > best_energy_by_aa( aa_range, 1000.0 );

	// cotainers for derivatives
	Multivec ref_deriv_weight( aa_range, 0.0 );
	utility::vector1< Real > dummy_set( score_list.size(), 0.0 );
	utility::vector1< utility::vector1< Real > > unweighted_E_dof( aa_range, dummy_set );

	process_rotamers(
		vars, num_energy_dofs, fixed_terms, score_list, fixed_score_list, aa_range, dummy_set,
		best_energy_by_aa, unweighted_E_dof, ref_deriv_weight );

	for ( Size ii = 1; ii <= aa_range; ++ii ) {
		if ( best_energy_by_aa[ ii ] == 1000.0 ) {
			//std::cerr << "Position " << position() << " with no energy represented for aa: " << chemical::AA( ii ) << std::endl;
			return 0.0;
		}
	}

	//Real const very_best = utility::min( best_energy_by_aa );
	//for ( Size ii = 1; ii <= best_energy_by_aa.size(); ++ii ) {
	//	best_energy_by_aa[ ii ] -= very_best;
	//}


	optimization::Multivec dpartition( vars.size(), 0.0 );
	optimization::Multivec dpssm_partition( vars.size(), 0.0 );
	Real partition( 0.0 ), pssm_partition( 0.0 );
	if ( option[ optE::sqrt_pssm ].user() ) {

		utility::vector1< Real > boltz_coeff( aa_range, 0.0 );
		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			Real const exp_term( std::exp( -1.0 *  inv_kT * best_energy_by_aa[ ii ] ) );
			partition += exp_term;
			boltz_coeff[ ii ] = exp_term;
		}
		Real const partition2 = partition * partition;
		Real dot_product = 0;
		utility::vector1< Real > energy_coeff( num_energy_dofs, 0.0 );
		//std::cout << "score: " << position() << " ";
		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			/// pssm_probabilities_ holds the squareroot of the probabilities
			dot_product += pssm_probabilities_[ ii ] * std::sqrt( boltz_coeff[ ii ] / partition );
			//std::cout << "( " << pssm_probabilities_[ ii ] << ", " << std::sqrt( boltz_coeff[ ii ] / partition ) << ") ";
			for ( Size jj = 1; jj <= num_energy_dofs; ++jj ) {
				energy_coeff[ jj ] += unweighted_E_dof[ ii ][ jj ] * boltz_coeff[ ii ];
			}

		}
		//std::cout << " total: " << -1.0 * dot_product << std::endl;
		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			Real const ii_negsqrt_prob = std::pow( boltz_coeff[ ii ] / partition, -0.5 );

			for ( Size jj = 1; jj <= num_energy_dofs; ++jj ) {
				dE_dvars[ jj ] += -0.5 * pssm_probabilities_[ ii ] *
					( ii_negsqrt_prob ) * component_weights[ type() ] *
					( boltz_coeff[ ii ] * energy_coeff[ jj ] -
					inv_kT * unweighted_E_dof[ ii ][ jj ] * boltz_coeff[ ii ] * partition ) /
					( partition2 );
			}
			for ( Size jj = 1; jj <= aa_range; ++jj ) {
				if ( ii == jj ) {
					dE_dvars[ num_energy_dofs + jj ] += -0.5 * pssm_probabilities_[ ii ] *
						( ii_negsqrt_prob ) * component_weights[ type() ] *
						( boltz_coeff[ ii ] * boltz_coeff[ ii ] -
						inv_kT * boltz_coeff[ ii ] * partition ) / ( partition2 );
				} else {
					dE_dvars[ num_energy_dofs + jj ] += -0.5 * pssm_probabilities_[ ii ] *
						( ii_negsqrt_prob ) * component_weights[ type() ] *
						( inv_kT * boltz_coeff[ ii ] * boltz_coeff[ jj ] ) / ( partition2 );
				}
			}
		}
		return -1.0 * component_weights[ type() ] * dot_product;
	} else {

		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			Real const exp_term( std::exp( -1.0 * inv_kT *  best_energy_by_aa[ ii ] ) );
			partition += exp_term;
			pssm_partition += pssm_probabilities_[ ii ] * exp_term;

			Real const ref_deriv_term( -1.0 * inv_kT * ref_deriv_weight[ ii ] * exp_term );
			dpartition[ num_energy_dofs + ii ] = ref_deriv_term;
			dpssm_partition[ num_energy_dofs + ii ] = pssm_probabilities_[ ii ] * ref_deriv_term;

			for ( Size jj = 1; jj <= num_energy_dofs; ++jj ) {
				Real e_dof_deriv( -1.0 * inv_kT * unweighted_E_dof[ ii ][ jj ] * exp_term );
				dpartition[ jj ] += e_dof_deriv;
				dpssm_partition[ jj ] += pssm_probabilities_[ ii ] * e_dof_deriv;
			}
		}
		for ( Size ii = 1; ii <= vars.size(); ++ii ) {
			dE_dvars[ ii ] += component_weights[ pssm_data ] * ( dpartition[ ii ] / partition - dpssm_partition[ ii ] / pssm_partition );
		}

		//std::cout << "Position data: ";
		//for ( Size ii = 1; ii <= aa_range; ++ii ) {
			//Real const exp_term( std::exp( -1.0 * best_energy_by_aa[ ii ] ) );
			//Real const ii_prob = exp_term / partition;
			//std::cout << "(" << ii_prob << ", " << pssm_probabilities_[ ii ] << ") ";
		//}
		//std::cout << std::endl << "score: " << -1 * std::log(  pssm_partition / partition ) << std::endl;

		return -1 * component_weights[ pssm_data ] * std::log(  pssm_partition / partition );
	}
}

void
PSSMOptEPositionData::print_score(
	std::ostream & ostr,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const & score_list,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//TR <<  "PSSMOptEPositionData::print_score" << std::endl;

	static Real const inv_kT( option[ optE::inv_kT_nataa ] );

	Size const aa_range( chemical::num_canonical_aas );

	utility::vector1< Real > best_energy_by_aa( aa_range, 1000.0 );

	// cotainers for derivatives
	Multivec ref_deriv_weight( aa_range, 0.0 );
	utility::vector1< Real > dummy_set( score_list.size(), 0.0 );
	utility::vector1< utility::vector1< Real > > unweighted_E_dof( aa_range, dummy_set );

	process_rotamers(
		vars, num_energy_dofs, fixed_terms, score_list, fixed_score_list, aa_range, dummy_set,
		best_energy_by_aa, unweighted_E_dof, ref_deriv_weight );

	for ( Size ii = 1; ii <= aa_range; ++ii ) {
		if ( best_energy_by_aa[ ii ] == 1000.0 ) {
			//std::cerr << "Position " << position() << " with no energy represented for aa: " << chemical::AA( ii ) << std::endl;
			return;// 0.0;
		}
	}

	optimization::Multivec dpartition( vars.size(), 0.0 );
	optimization::Multivec dpssm_partition( vars.size(), 0.0 );
	Real partition( 0.0 ), pssm_partition( 0.0 );
	if ( option[ optE::sqrt_pssm ].user() ) {

		utility::vector1< Real > boltz_coeff( aa_range, 0.0 );
		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			Real const exp_term( std::exp( -1.0 * inv_kT * best_energy_by_aa[ ii ] ) );
			partition += exp_term;
			boltz_coeff[ ii ] = exp_term;
		}
		Real const partition2 = partition * partition;
		Real dot_product = 0;
		utility::vector1< Real > energy_coeff( num_energy_dofs, 0.0 );
		//std::cout << "score: " << position() << " ";
		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			/// pssm_probabilities_ holds the squareroot of the probabilities
			dot_product += pssm_probabilities_[ ii ] * std::sqrt( boltz_coeff[ ii ] / partition );
			//std::cout << "( " << pssm_probabilities_[ ii ] << ", " << std::sqrt( boltz_coeff[ ii ] / partition ) << ") ";
			for ( Size jj = 1; jj <= num_energy_dofs; ++jj ) {
				energy_coeff[ jj ] += unweighted_E_dof[ ii ][ jj ] * boltz_coeff[ ii ];
			}

		}
		//std::cout << " total: " << -1.0 * dot_product << std::endl;
		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			Real const ii_negsqrt_prob = std::pow( boltz_coeff[ ii ] / partition, -0.5 );

			for ( Size jj = 1; jj <= num_energy_dofs; ++jj ) {
				dE_dvars[ jj ] += -0.5 * pssm_probabilities_[ ii ] *
					( ii_negsqrt_prob ) *
					( boltz_coeff[ ii ] * energy_coeff[ jj ] -
					inv_kT * unweighted_E_dof[ ii ][ jj ] * boltz_coeff[ ii ] * partition ) /
					( partition2 );
			}
			for ( Size jj = 1; jj <= aa_range; ++jj ) {
				if ( ii == jj ) {
					dE_dvars[ num_energy_dofs + jj ] += -0.5 * pssm_probabilities_[ ii ] *
						( ii_negsqrt_prob ) *
						( boltz_coeff[ ii ] * boltz_coeff[ ii ] -
						inv_kT * boltz_coeff[ ii ] * partition ) / ( partition2 );
				} else {
					dE_dvars[ num_energy_dofs + jj ] += -0.5 * pssm_probabilities_[ ii ] *
						( ii_negsqrt_prob ) *
						( inv_kT * boltz_coeff[ ii ] * boltz_coeff[ jj ] ) / ( partition2 );
				}
			}
		}

		ostr << "PSSM_SQRT " << native_aa() << " " << (int) native_aa();
		ostr << " nneighb: " << neighbor_count() << " dotpr: " << dot_product << " ";
		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			ostr << "([ " << ii << "] " << pssm_probabilities_[ ii ] << " , " << std::sqrt( boltz_coeff[ ii ] / partition ) << " ) ";
		}
		ostr << " " << tag() << "\n";

		//TR <<  "PSSMOptEPositionData::print_score done" << std::endl;

		return;// -1.0 * dot_product;
	} else {

		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			Real const exp_term( std::exp( -1.0 * inv_kT * best_energy_by_aa[ ii ] ) );
			partition += exp_term;
			pssm_partition += pssm_probabilities_[ ii ] * exp_term;

			Real const ref_deriv_term( -1.0 * inv_kT * ref_deriv_weight[ ii ] * exp_term );
			dpartition[ num_energy_dofs + ii ] = ref_deriv_term;
			dpssm_partition[ num_energy_dofs + ii ] = pssm_probabilities_[ ii ] * ref_deriv_term;

			for ( Size jj = 1; jj <= num_energy_dofs; ++jj ) {
				Real e_dof_deriv( -1.0 * inv_kT * unweighted_E_dof[ ii ][ jj ] * exp_term );
				dpartition[ jj ] += e_dof_deriv;
				dpssm_partition[ jj ] += pssm_probabilities_[ ii ] * e_dof_deriv;
			}
		}
		for ( Size ii = 1; ii <= vars.size(); ++ii ) {
			dE_dvars[ ii ] += dpartition[ ii ] / partition - dpssm_partition[ ii ] / pssm_partition;
		}

		//std::cout << "Position data: ";
		//for ( Size ii = 1; ii <= aa_range; ++ii ) {
			//Real const exp_term( std::exp( -1.0 * best_energy_by_aa[ ii ] ) );
			//Real const ii_prob = exp_term / partition;
			//std::cout << "(" << ii_prob << ", " << pssm_probabilities_[ ii ] << ") ";
		//}
		//std::cout << std::endl << "score: " << -1 * std::log(  pssm_partition / partition ) << std::endl;

		ostr << "PSSM_PART " << native_aa() << " " << (int) native_aa();
		ostr << " nneighb: " << neighbor_count() << " pssmpart: " << pssm_partition / partition << " ";
		ostr << " -lnpssmpart: " <<  -1.0 * std::log(  pssm_partition / partition ) << " ";
		ostr << " -compwt_lnpssm " << -1.0 * component_weights[ pssm_data ] * std::log(  pssm_partition / partition );
		for ( Size ii = 1; ii <= aa_range; ++ii ) {
			ostr << "([ " << ii << "] " << pssm_probabilities_[ ii ] << " , " << std::exp( -1.0 * best_energy_by_aa[ ii ] ) / partition << " ) ";
		}

		ostr << " " << tag() << "\n";

		//TR <<  "PSSMOptEPositionData::print_score done" << std::endl;

		return;// -1 * std::log(  pssm_partition / partition );
	}
}


void
PSSMOptEPositionData::set_pssm_probabilities(
	utility::vector1< Real > const & pssm_probs
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	/// runtime_assert something like that they sum to 1, or pretty close
	/// and that the size = num_canonincal_aas
	pssm_probabilities_ = pssm_probs;

	// Store the squareroot of the input probabilities if we're optimizing using the
	// dot product of unit vectors.
	if ( option[ optE::sqrt_pssm ].user() ) {
		for ( Size ii = 1; ii <= pssm_probs.size(); ++ii ) {
			pssm_probabilities_[ ii ] = std::sqrt( pssm_probs[ ii ] );
		}
	}
}

OptEPositionDataType
PSSMOptEPositionData::type() const
{
	return pssm_data;
}


void
PSSMOptEPositionData::write_to_file( std::ofstream & outfile ) const
{
	outfile << pssm_probabilities_.size() << " ";
	for ( Size ii = 1; ii <= pssm_probabilities_.size(); ++ii ) {
		outfile << pssm_probabilities_[ ii ] << " ";
	}
	outfile << "\n";
	parent::write_to_file( outfile );
}


void
PSSMOptEPositionData::read_from_file( std::ifstream & infile )
{
	Size nprobs;
	infile >> nprobs;
	pssm_probabilities_.resize( nprobs );
	for ( Size ii = 1; ii <= pssm_probabilities_.size(); ++ii ) {
		infile >> pssm_probabilities_[ ii ];
	}
	parent::read_from_file( infile );
}


void
PSSMOptEPositionData::write_to_binary_file( std::ofstream & /*outfile*/ ) const
{}


void
PSSMOptEPositionData::read_from_binary_file( std::ifstream & /*infile*/ )
{}

Size
PSSMOptEPositionData::memory_use() const
{
	Size total = pssm_probabilities_.size() * sizeof( Real );
	return total + parent::memory_use();
}


#ifdef USEMPI

void
PSSMOptEPositionData::send_to_node( int const destination_node, int const tag ) const
{
	//TR << "PSSMOptEPositionData::Sending data to node " << destination_node << std::endl;

	int n_probabilities = pssm_probabilities_.size();
	MPI_Send( & n_probabilities, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );
	if ( n_probabilities != 0 ) {
		Real * pssmprobs = new Real[ n_probabilities ];
		for ( int ii = 1; ii <= n_probabilities; ++ii ) {
			pssmprobs[ ii - 1 ] = pssm_probabilities_[ ii ];
		}
		MPI_Send( pssmprobs, n_probabilities, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );
		delete [] pssmprobs;
	}
	parent::send_to_node( destination_node, tag );
}


void
PSSMOptEPositionData::receive_from_node( int const source_node, int const tag )
{
	MPI_Status stat;
	int n_probabilities;
	MPI_Recv( & n_probabilities, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	if ( n_probabilities != 0 ) {
		Real * pssmprobs = new Real[ n_probabilities ];
		pssm_probabilities_.resize( n_probabilities );
		MPI_Recv( pssmprobs, n_probabilities, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
		for ( int ii = 1; ii <= n_probabilities; ++ii ) {
			pssm_probabilities_[ ii ] = pssmprobs[ ii - 1 ];
		}
		delete [] pssmprobs;
	}
	parent::receive_from_node( source_node, tag );

}

#endif //USEMPI


//
// ------------------- P Native Rotamer -----------------------//
//

PNatRotOptEPositionData::PNatRotOptEPositionData() : phi_(0.0), psi_(0.0) {}

PNatRotOptEPositionData::~PNatRotOptEPositionData() {}

Real
PNatRotOptEPositionData::get_score(
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const dummy1,
	int const dummy2,
	EnergyMap const & fixed_terms,
	ScoreTypes const & dummy3 /*free_score_list*/,
	ScoreTypes const & fixed_score_list
) const
{
	return process_score(
		TR, false, component_weights, vars,
		dE_dvars, num_energy_dofs, dummy1, dummy2,
		fixed_terms, dummy3, fixed_score_list );
}

/*
{
	using namespace core::optimization;
	using namespace utility;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static Real const inv_kT( option[ optE::inv_kT_natrot ] );


	utility::vector1< Real > best_energy_by_rotwell( n_wells_, 1000.0 );
	utility::vector1< Real > well_is_represented( n_wells_, false );

	// cotainers for derivatives
	utility::vector1< Real > dummy_set( score_list.size(), 0.0 );
	utility::vector1< utility::vector1< Real > > unweighted_E_dof( n_wells_, dummy_set );

	/// Zero the static data
	//for ( Size ii = 1; ii <= unweighted_E_dof.size(); ++ii ) {
	//	for ( Size jj = 1; jj <= unweighted_E_dof[ii].size(); ++jj ) {
	//		unweighted_E_dof[ ii ][ jj ] = 0.0;
	//	}
	//}

	Size const natrotwell = rotamer_index_2_well_id( native_rotamer_index_ );

	// HACK below
	Real best_energy( 0.0 );
	bool best_energy_is_native( false ), first_round( true );

	for( PNatRotOptERotamerDataOPs::const_iterator itr = rotamer_data_begin(),
			e_itr = rotamer_data_end() ; itr != e_itr; ++itr ) {
		utility::vector1< Size > const & this_rotindex( (*itr)->rotamer_index() );
		Size this_rotwell = rotamer_index_2_well_id( this_rotindex );

		if ( this_rotwell != natrotwell && count_rotamer_as_native( *itr ) ) {
			this_rotwell = natrotwell;
		}

		Real weighted_energy( 0.0 );
		// Variable-weighted energy terms
		for( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
			weighted_energy += vars[ ii ] * ((**itr)[ ii ]);
		}
		// Fixed-weight energy terms
		for( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
			weighted_energy += fixed_terms[ fixed_score_list[ ii ] ] * (*itr)->fixed_data()[ ii ];
		}

		Real const cutoff( 300.0 );
		if( weighted_energy > cutoff ) weighted_energy = cutoff;
		if( weighted_energy < -1.0*cutoff ) weighted_energy = -1.0*cutoff;

		if ( first_round || weighted_energy < best_energy ) {
			best_energy = weighted_energy;
			first_round = false;
			best_energy_is_native = this_rotwell == natrotwell;
		}

		if( weighted_energy < best_energy_by_rotwell[ this_rotwell ] ) {

			best_energy_by_rotwell[ this_rotwell ] = weighted_energy;
			well_is_represented[ this_rotwell ] = true;

			if( std::abs( weighted_energy ) < cutoff ) {
				unweighted_E_dof[ this_rotwell ] = (*itr)->free_data();
			} else {
				// derivatives will be zero above and below cutoffs
				unweighted_E_dof[ this_rotwell ] = dummy_set;
			}
		}
	} // done processing rotamers
	// now do the partition function analysis

	if ( ! well_is_represented[ rotamer_index_2_well_id( native_rotamer_index_ ) ] ) {
		/// This shouldn't happen, but it's not disasterous if it does.  ... er, it happens from time to time
		///std::cerr << "Skipping score in the absence of the native rotamer!" << std::endl;
		return 0.0;
	}

	Real numerator(0.0), partition(0.0);
	Multivec dpartition( vars.size(), 0.0 ), dnumerator( vars.size(), 0.0 );

	for( LexicographicalIterator iter( rotamer_well_counts_ ); ! iter.at_end(); ++iter ) {
		Size const iter_rotwell = rotamer_index_2_well_id( iter );
		if ( ! well_is_represented[ iter_rotwell ] ) continue;

		Real const exp_term( std::exp( -1.0 * inv_kT * best_energy_by_rotwell[ iter_rotwell ] ) );
		partition += exp_term;
		bool const native_rotamer( is_native_rotamer_well( iter ));
		if ( native_rotamer ) numerator = exp_term;

		// partitions for energy derivatives
		for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
			// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
			Real e_dof_deriv( -1.0 * inv_kT * unweighted_E_dof[ iter_rotwell ][ e_dof ] * exp_term );
			dpartition[ e_dof ] += e_dof_deriv;
			if ( native_rotamer ) dnumerator[ e_dof ] = e_dof_deriv;
		}
	}

	// accumulate to passed-in derivative sums
	for ( Size dof(1); dof <= vars.size(); ++dof ) {
		dE_dvars[ dof ] += component_weights[ prob_native_rotamer ] * ( dpartition[ dof ] / partition - dnumerator[ dof ] / numerator );
	}

	//return best_energy_is_native ? 0 : 1;
	return ( -1.0 * component_weights[ prob_native_rotamer ] * std::log( numerator / partition ) );

}*/

void
PNatRotOptEPositionData::print_score(
	std::ostream & ostr,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const dummy1,
	int const dummy2,
	EnergyMap const & fixed_terms,
	ScoreTypes const & dummy3,
	ScoreTypes const & fixed_score_list
) const
{
	process_score(
		ostr, true, component_weights, vars,
		dE_dvars, num_energy_dofs, dummy1, dummy2,
		fixed_terms, dummy3, fixed_score_list );
}



Real
PNatRotOptEPositionData::process_score(
	std::ostream & ostr,
	bool print,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const & score_list,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core::optimization;
	using namespace utility;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static Real const inv_kT( option[ optE::inv_kT_natrot ] );

	static bool const dun10active( option[ corrections::score::dun10 ] );

	//TR <<  "PNatRotOptEPositionData::print_score" << std::endl;

	utility::vector1< Real > best_energy_by_rotwell( n_wells_, 1000.0 );
	utility::vector1< PNatRotOptERotamerDataOP > best_rotamer_by_rotwell( n_wells_ );
	utility::vector1< Real > well_is_represented( n_wells_, false );

	// cotainers for derivatives
	utility::vector1< Real > dummy_set( score_list.size(), 0.0 );
	utility::vector1< utility::vector1< Real > > unweighted_E_dof( n_wells_, dummy_set );

	/// Zero the static data
	//for ( Size ii = 1; ii <= unweighted_E_dof.size(); ++ii ) {
	//	for ( Size jj = 1; jj <= unweighted_E_dof[ii].size(); ++jj ) {
	//		unweighted_E_dof[ ii ][ jj ] = 0.0;
	//	}
	//}

	// HACK below
	Real best_energyHACK( 0.0 );
	bool best_energy_is_native( false ), first_round( true );

	Size const natrotwell = rotamer_index_2_well_id( native_rotamer_index_ );
	Size nnearnative( 0 );
	Size ndecoys( 0 );

	for( PNatRotOptERotamerDataOPs::const_iterator itr = rotamer_data_begin(),
			e_itr = rotamer_data_end() ; itr != e_itr; ++itr ) {
		utility::vector1< Size > const & this_rotindex( (*itr)->rotamer_index() );
		Size this_rotwell = rotamer_index_2_well_id( this_rotindex );

		bool count_as_native( false );
		if ( this_rotwell != natrotwell && count_rotamer_as_native( *itr ) ) {
			this_rotwell = natrotwell;
			count_as_native = true;
			++nnearnative;
		} else {
			++ndecoys;
		}

		Real weighted_energy( 0.0 );
		if ( count_as_native && dun10active ) {
			using namespace core::chemical;
			/// Entropy bonus for the native rotamer; the semi-rotameric residues
			/// have fewer wells in the '02 library than they do in the '08 library.
			switch ( aa_ ) {
				case aa_asp : weighted_energy = std::log( 3.0f / 6.0f ); break;
				case aa_glu : weighted_energy = std::log( 3.0f / 6.0f ); break;
				case aa_phe : weighted_energy = std::log( 2.0f / 6.0f ); break;
				case aa_his : weighted_energy = std::log( 3.0f / 12.0f ); break;
				case aa_asn : weighted_energy = std::log( 6.0f / 12.0f ); break;
				case aa_gln : weighted_energy = std::log( 4.0f / 12.0f ); break;
				case aa_tyr : weighted_energy = std::log( 2.0f / 6.0f ); break;
				default:
				break;
			}
			///weighted_energy = 0; // TEMP TEMP TEMP deactivate entropy bonus.
		}
		// Variable-weighted energy terms
		for( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
			weighted_energy += vars[ ii ] * ((**itr)[ ii ]);
		}
		// Fixed-weight energy terms
		for( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
			weighted_energy += fixed_terms[ fixed_score_list[ ii ] ] * (*itr)->fixed_data()[ ii ];
		}

		Real const cutoff( 300.0 );
		if( weighted_energy > cutoff ) weighted_energy = cutoff;
		if( weighted_energy < -1.0*cutoff ) weighted_energy = -1.0*cutoff;

		if ( first_round || weighted_energy < best_energyHACK ) {
			best_energyHACK = weighted_energy;
			first_round = false;
			best_energy_is_native = this_rotwell == natrotwell;
		}

		if( weighted_energy < best_energy_by_rotwell[ this_rotwell ] ) {
			best_rotamer_by_rotwell[ this_rotwell ] = *itr;
			best_energy_by_rotwell[ this_rotwell ] = weighted_energy;
			well_is_represented[ this_rotwell ] = true;

			if( std::abs( weighted_energy ) < cutoff ) {
				unweighted_E_dof[ this_rotwell ] = (*itr)->free_data();
			} else {
				// derivatives will be zero above and below cutoffs
				unweighted_E_dof[ this_rotwell ] = dummy_set;
			}
		}
	} // done processing rotamers
	// now do the partition function analysis


	Real numerator(0.0), partition(0.0);
	Multivec dpartition( vars.size(), 0.0 ), dnumerator( vars.size(), 0.0 );
	static Real const SENTINEL = 1234;
	Real best_energy( SENTINEL ); Size best_rotind( 0 );
	utility::vector1< Size > best_rotwell( rotamer_well_counts_.size(), 0 );
	utility::vector1< Size > native_rotwell( rotamer_well_counts_.size(), 0 );


	for( LexicographicalIterator iter( rotamer_well_counts_ ); ! iter.at_end(); ++iter ) {
		Size const iter_rotwell = rotamer_index_2_well_id( iter );
		if ( ! well_is_represented[ iter_rotwell ] ) continue;

		Real const exp_term( std::exp( -1.0 * inv_kT * best_energy_by_rotwell[ iter_rotwell ] ) );
		partition += exp_term;
		bool const native_rotamer( is_native_rotamer_well( iter ));
		if ( native_rotamer ) numerator = exp_term;

		if ( best_energy == SENTINEL || best_energy > best_energy_by_rotwell[ iter_rotwell ] ) {
			best_energy = best_energy_by_rotwell[ iter_rotwell ];
			best_rotind = iter_rotwell;
			for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) {
				best_rotwell[ ii ] = iter[ ii ];
			}
		}

		if ( native_rotamer ) {
			for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) {
				native_rotwell[ ii ] = iter[ ii ];
			}
		}
		// partitions for energy derivatives
		for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
			// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
			Real e_dof_deriv( -1.0 * inv_kT * unweighted_E_dof[ iter_rotwell ][ e_dof ] * exp_term );
			dpartition[ e_dof ] += e_dof_deriv;
			if ( native_rotamer ) dnumerator[ e_dof ] = e_dof_deriv;
		}
	}


	if ( ! well_is_represented[ rotamer_index_2_well_id( native_rotamer_index_ ) ] ) {
		/// This shouldn't happen, but it's not disasterous if it does.  ... er, it happens from time to time
		/// std::cerr << "Skipping score in the absence of the native rotamer!" << std::endl;
		/// CHANGE OF HEART.  If the library has not built a near-native rotamer, assign a penalty!
		/// 100 seems like a decent penalty.
		if ( print ) {
			ostr << "PNATROT-NO-NEAR-NATIVE nnat:" << nnearnative << " ndec: " << ndecoys << " p: " << 0.0 << " -lnp: ";
			ostr << 100;
			ostr << " compwt_lnp: " << 100;
			ostr << " n<d? " << 0 << " " << tag() << " " << rotamer_well_counts_.size() << " n:";

			for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) ostr << " " << native_rotwell[ ii ] << " " << native_chi_[ ii ];
			ostr << " phi/psi:" << phi_ << " " << psi_ << "\n";

			for( PNatRotOptERotamerDataOPs::const_iterator itr = rotamer_data_begin(),
					e_itr = rotamer_data_end() ; itr != e_itr; ++itr ) {

				ostr << "NON-NATIVE-ROT " << tag() << " ";
				for ( Size ii = 1; ii <= (*itr)->rotamer_index().size(); ++ii ) {
					ostr << (*itr)->rotamer_index()[ ii ] << " "<< (*itr)->chi()[ ii ] << " ";
				}
				ostr << "\n";
			}
		}


		return 100;// 0.0;
	}


	// accumulate to passed-in derivative sums
	for ( Size dof(1); dof <= vars.size(); ++dof ) {
		dE_dvars[ dof ] += dpartition[ dof ] / partition - dnumerator[ dof ] / numerator;
	}

	if ( print ) {
		ostr << "PNATROT nnat:" << nnearnative << " ndec: " << ndecoys << " p: " << numerator / partition << " -lnp: ";
		ostr << -1.0 * std::log( numerator / partition );
		ostr << " compwt_lnp: " << -1.0 * component_weights[ prob_native_rotamer ] * std::log( numerator / partition );
		ostr << " n<d? " << best_energy_is_native << " " << tag() << " " << rotamer_well_counts_.size() << " n:";

		for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) ostr << " " << native_rotwell[ ii ] << " " << native_chi_[ ii ];
		ostr << " b:";
		for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) {
			ostr << " " << best_rotwell[ ii ] << " " << best_rotamer_by_rotwell[ best_rotind ]->chi()[ ii ] ;
		}
		ostr << " phi/psi:" << phi_ << " " << psi_ << "\n";
	}

	//TR <<  "PNatRotOptEPositionData::print_score done" << std::endl;

	return ( -1.0 * component_weights[ prob_native_rotamer ] * std::log( numerator / partition ) );

}

void
PNatRotOptEPositionData::range(
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list,
	EnergyMap & lower_bound,
	EnergyMap & upper_bound
) const
{
	for ( Size ii = 1; ii <= free_score_list.size(); ++ii ) {
		Real ii_min = lower_bound[ free_score_list[ ii ] ];
		Real ii_max = upper_bound[ free_score_list[ ii ] ];
		for( PNatRotOptERotamerDataOPs::const_iterator iter = rotamer_data_begin(),
				e_itr = rotamer_data_end() ; iter != e_itr; ++iter ) {
			if ( ii_min > (*iter)->free_data()[ ii ] ) {
				ii_min = (*iter)->free_data()[ ii ];
			}
			if ( ii_max < (*iter)->free_data()[ ii ] ) {
				ii_max = (*iter)->free_data()[ ii ];
			}
		}
		lower_bound[ free_score_list[ ii ] ] = ii_min;
		upper_bound[ free_score_list[ ii ] ] = ii_max;
	}

	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		Real ii_min = lower_bound[ fixed_score_list[ ii ] ];
		Real ii_max = upper_bound[ fixed_score_list[ ii ] ];
		for( PNatRotOptERotamerDataOPs::const_iterator iter = rotamer_data_begin(),
				e_itr = rotamer_data_end() ; iter != e_itr; ++iter ) {
			if ( ii_min > (*iter)->fixed_data()[ ii ] ) {
				ii_min = (*iter)->fixed_data()[ ii ];
			}
			if ( ii_max < (*iter)->fixed_data()[ ii ] ) {
				ii_max = (*iter)->fixed_data()[ ii ];
			}
		}
		lower_bound[ fixed_score_list[ ii ] ] = ii_min;
		upper_bound[ fixed_score_list[ ii ] ] = ii_max;
	}
}



Size
PNatRotOptEPositionData::size() const
{
	return data_.size();
}


OptEPositionDataType
PNatRotOptEPositionData::type() const
{
	return prob_native_rotamer;
}


void
PNatRotOptEPositionData::write_to_file( std::ofstream & /*outfile*/ ) const
{}


void
PNatRotOptEPositionData::read_from_file( std::ifstream & /*infile*/ )
{}


void
PNatRotOptEPositionData::write_to_binary_file( std::ofstream & /*outfile*/ ) const
{}


void
PNatRotOptEPositionData::read_from_binary_file( std::ifstream & /*infile*/ )
{}

Size
PNatRotOptEPositionData::memory_use() const
{
	Size total = sizeof( PNatRotOptEPositionData ) +
		sizeof( PNatRotOptERotamerData ) * data_.size() +
		sizeof( Size ) * native_rotamer_index_.size() +
		sizeof( Size ) * rotamer_well_counts_.size();
	if ( data_.size() > 0 ) {
		total += sizeof( Real ) * ( data_[ 1 ]->free_data().size() + data_[ 1 ]->fixed_data().size() ) * data_.size();
		total += sizeof( Size ) * data_[ 1 ]->rotamer_index().size() * data_.size();
	}
	return total;
}


#ifdef USEMPI
void
PNatRotOptEPositionData::send_to_node( int const destination_node, int const tag ) const
{
	//MPI_Status stat;

	//TR << "PNatRotOptEPositionData:: Sending data to node... " << destination_node << std::endl;

	/// 3. How many chi angles
	int n_chi( data_[ 1 ]->rotamer_index().size() );
	MPI_Send( & n_chi, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD);

	{
	/// 4. How many chi wells
	int * rotamer_well_counts = new int[ n_chi ];
	for ( int ii = 1; ii <= n_chi; ++ii ) {
		rotamer_well_counts[ ii - 1 ] = rotamer_well_counts_[ ii ];
	}
	MPI_Send( rotamer_well_counts, n_chi, MPI_INT, destination_node, tag, MPI_COMM_WORLD );
	delete [] rotamer_well_counts;
	}

	{
	/// 5. What is the native rotamer at this position?
	int * native_rotamer = new int[ n_chi ];
	for ( int ii = 1; ii <= n_chi; ++ii ) {
		native_rotamer[ ii -1 ] = native_rotamer_index_[ ii ];
	}
	MPI_Send( native_rotamer, n_chi, MPI_INT, destination_node, tag, MPI_COMM_WORLD );
	delete [] native_rotamer;
	}

	/// 6. The number of rotamers for this position
	Size ii_num_rotamers = data_.size();
	MPI_Send( & ii_num_rotamers, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	if ( ii_num_rotamers == 0 ) return;

	Size free_count(  data_[ 1 ]->free_data().size() );
	Size fixed_count( data_[ 1 ]->fixed_data().size() );

	/// 7 The size of the free and fixed data.
	MPI_Send( & free_count,  1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );
	MPI_Send( & fixed_count, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );


	Size * ii_rot_indices = new Size[ ii_num_rotamers * n_chi ];
	Real * chi = new Real[ ii_num_rotamers * n_chi ];
	Real * free_data = new Real[ ii_num_rotamers * free_count ];
	Real * fixed_data = new Real[ ii_num_rotamers * fixed_count ];

	for ( Size jj = 1; jj <= ii_num_rotamers; ++jj ) {
		utility::vector1< Size > const & rotamer_index_vect( data_[ jj ]->rotamer_index() );
		utility::vector1< Real > const & chi_vect( data_[ jj ]->chi() );
 		utility::vector1< Real > const & free_data_vect( data_[ jj ]->free_data() );
		utility::vector1< Real > const & fixed_data_vect( data_[ jj ]->fixed_data() );

		for ( int kk = 1; kk <= n_chi; ++kk ) {
			ii_rot_indices[ ( jj - 1 ) * n_chi + kk - 1 ] = rotamer_index_vect[ kk ];
			chi[ ( jj - 1 ) * n_chi + kk - 1 ] = chi_vect[ kk ];
		}
		for ( Size kk = 1; kk <= free_count; ++kk ) {
			free_data[ ( jj - 1 ) * free_count + kk - 1 ] = free_data_vect[ kk ];
		}
		for ( Size kk = 1; kk <= fixed_count; ++kk ) {
			fixed_data[ ( jj - 1 ) * fixed_count + kk - 1 ] = fixed_data_vect[ kk ];
		}

	}
	/// 8. All the amino acids for all rotamers at this position
	MPI_Send( ii_rot_indices, ii_num_rotamers * n_chi, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );
	MPI_Send( chi, ii_num_rotamers * n_chi, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 9. All the free data for all rotamers
	MPI_Send( free_data, ii_num_rotamers * free_count, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 10. All the fixed data for all rotamers
	MPI_Send( fixed_data, ii_num_rotamers * fixed_count, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );


	delete [] ii_rot_indices; ii_rot_indices = 0;
	delete [] chi; chi = 0;
	delete [] free_data;   free_data   = 0;
	delete [] fixed_data;  fixed_data  = 0;

	/// 11. Native chi data
	Size nnative_chi = native_chi_.size();
	MPI_Send( & nnative_chi, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );
	if ( nnative_chi != 0 ) {
		Real * native_chi = new Real[ nnative_chi ];
		for ( Size ii = 1; ii <= nnative_chi; ++ii ) native_chi[ ii - 1 ] = native_chi_[ ii ];
		MPI_Send( native_chi, nnative_chi, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );
		delete [] native_chi;
	}

	/// 12. Native chi periodicity data
	Size nnative_chi_periodicity = native_chi_periodicity_.size();
	MPI_Send( & nnative_chi_periodicity, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );
	if ( nnative_chi_periodicity != 0 ) {
		Real * native_chi_periodicity = new Real[ nnative_chi_periodicity ];
		for ( Size ii = 1; ii <= nnative_chi_periodicity; ++ii ) native_chi_periodicity[ ii - 1 ] = native_chi_periodicity_[ ii ];
		MPI_Send( native_chi_periodicity, nnative_chi_periodicity, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );
		delete [] native_chi_periodicity;
	}

	/// 13. phi/psi
	Real phi( phi_ ), psi( psi_ );
	MPI_Send( & phi, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );
	MPI_Send( & psi, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 14. aa_
	int aa( aa_ );
	MPI_Send( & aa, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );


	//TR << "... sent data to node " << destination_node << std::endl;
	OptEPositionData::send_to_node( destination_node, tag );
}


void
PNatRotOptEPositionData::receive_from_node( int const source_node, int const tag )
{
	MPI_Status stat;

	//TR << "PNatRotOptEPositionData::Recieving data from node... " << source_node << std::endl;

	/// 3. How many chi angles
	int n_chi;
	MPI_Recv( & n_chi, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );

	//TR << "PNatRotOptEPositionData:: received nchi = " << n_chi << std::endl;

	runtime_assert( n_chi > 0 );

	{ // scope
	/// 4. How many chi wells
	int * rotamer_well_counts = new int[ n_chi ];
	MPI_Recv( rotamer_well_counts, n_chi, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	utility::vector1< Size > rotamer_well_counts_v( n_chi );
	for ( int ii = 1; ii <= n_chi; ++ii ) {
		rotamer_well_counts_v[ ii ] = rotamer_well_counts[ ii - 1 ];
		//TR << "PNatRotOptEPositionData:: received rotamer_well_counts[  " << ii << "] "  << rotamer_well_counts[ ii - 1 ] << std::endl;
	}
	delete [] rotamer_well_counts;
	set_rotamer_well_counts( rotamer_well_counts_v );
	}

	{ // scope
	/// 5. What is the native rotamer at this position?
	int * native_rotamer = new int[ n_chi ];
	MPI_Recv( native_rotamer, n_chi, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	utility::vector1< Size > native_rotamer_v( n_chi );
	for ( int ii = 1; ii <= n_chi; ++ii ) {
		native_rotamer_v[ ii ] = native_rotamer[ ii -1 ];
		//TR << "PNatRotOptEPositionData:: received native_rotamer[  " << ii << "] "  << native_rotamer[ ii - 1 ] << std::endl;
	}
	set_native_rotamer_index( native_rotamer_v );
	delete [] native_rotamer;
	}

	/// 6. The number of rotamers for this position
	Size ii_num_rotamers;
	MPI_Recv( & ii_num_rotamers, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	//TR << "PNatRotOptEPositionData:: received ii_num_rotamers  " << ii_num_rotamers << std::endl;

	if ( ii_num_rotamers == 0 ) return;

	Size free_count(0);
	Size fixed_count(0);

	/// 7 The size of the free and fixed data.
	MPI_Recv( & free_count,  1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	MPI_Recv( & fixed_count, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );


	Size * ii_rot_indices = new Size[ ii_num_rotamers * n_chi ];
	Real * chi = new Real[ ii_num_rotamers * n_chi ];
	Real * free_data = new Real[ ii_num_rotamers * free_count ];
	Real * fixed_data = new Real[ ii_num_rotamers * fixed_count ];

	/// 8. All the amino acids for all rotamers at this position
	MPI_Recv( ii_rot_indices, ii_num_rotamers * n_chi, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	MPI_Recv( chi, ii_num_rotamers * n_chi, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 9. All the free data for all rotamers
	MPI_Recv( free_data, ii_num_rotamers * free_count, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 10. All the fixed data for all rotamers
	MPI_Recv( fixed_data, ii_num_rotamers * fixed_count, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	utility::vector1< Size > rotamer_index_vect( n_chi );
	utility::vector1< Real > chi_vect( n_chi );
 	utility::vector1< Real > free_data_vect( free_count );
	utility::vector1< Real > fixed_data_vect( fixed_count );

	for ( Size jj = 1; jj <= ii_num_rotamers; ++jj ) {
		for ( int kk = 1; kk <= n_chi; ++kk ) {
			rotamer_index_vect[ kk ] = ii_rot_indices[ ( jj - 1 ) * n_chi + kk - 1 ];
			chi_vect[ kk ] = chi[ ( jj - 1 ) * n_chi + kk - 1 ];
		}
		for ( Size kk = 1; kk <= free_count; ++kk ) {
			free_data_vect[ kk ] = free_data[ ( jj - 1 ) * free_count + kk - 1 ];
		}
		for ( Size kk = 1; kk <= fixed_count; ++kk ) {
			fixed_data_vect[ kk ] = fixed_data[ ( jj - 1 ) * fixed_count + kk - 1 ];
		}
		PNatRotOptERotamerDataOP jj_rotamer_data = new PNatRotOptERotamerData(
			rotamer_index_vect,
			chi_vect,
			free_data_vect,
			fixed_data_vect );
		add_rotamer_line_data( jj_rotamer_data );

	}

	delete [] ii_rot_indices; ii_rot_indices = 0;
	delete [] chi; chi = 0;
	delete [] free_data;   free_data   = 0;
	delete [] fixed_data;  fixed_data  = 0;

	/// 11. Native chi data
	Size nnative_chi;
	MPI_Recv( & nnative_chi, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	if ( nnative_chi != 0 ) {
		native_chi_.resize( nnative_chi );
		Real * native_chi = new Real[ nnative_chi ];
		MPI_Recv( native_chi, nnative_chi, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
		for ( Size ii = 1; ii <= nnative_chi; ++ii ) native_chi_[ ii ] = native_chi[ ii - 1 ];
		delete [] native_chi;
	}

	/// 12. Native chi periodicity data
	Size nnative_chi_periodicity = native_chi_periodicity_.size();
	MPI_Recv( & nnative_chi_periodicity, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	if ( nnative_chi_periodicity != 0 ) {
		Real * native_chi_periodicity = new Real[ nnative_chi_periodicity ];
		MPI_Recv( native_chi_periodicity, nnative_chi_periodicity, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
		native_chi_periodicity_.resize( nnative_chi_periodicity );
		for ( Size ii = 1; ii <= nnative_chi_periodicity; ++ii ) native_chi_periodicity_[ ii ] = native_chi_periodicity[ ii -1 ];
		delete [] native_chi_periodicity;
	}

	/// 13. phi/psi
	Real phi, psi;
	MPI_Recv( & phi, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	MPI_Recv( & psi, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	phi_ = phi;
	psi_ = psi;

	/// 14. aa_
	int aa;
	MPI_Recv( & aa, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	aa_ = core::chemical::AA( aa );

	//TR << "... received data from node " << source_node << std::endl;
	OptEPositionData::receive_from_node( source_node, tag );

}
#endif

void
PNatRotOptEPositionData::set_native_rotamer_index(
	utility::vector1< Size > const & native_rotamer_index
)
{
	runtime_assert( rotamer_well_counts_.size() != 0 ); // set well counts before native rotamer
	runtime_assert( native_rotamer_index.size() == rotamer_well_counts_.size() );
	native_rotamer_index_ = native_rotamer_index;
	for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) {
		runtime_assert( native_rotamer_index_[ ii ] > 0 );
		runtime_assert( native_rotamer_index_[ ii ] <= rotamer_well_counts_[ ii ] );
	}
}

void
PNatRotOptEPositionData::set_native_rotamer_chi(
	utility::vector1< Real > const & native_chi
)
{
	native_chi_ = native_chi;
}

void
PNatRotOptEPositionData::set_native_chi_periodicity(
	utility::vector1< Real > const & native_chi_periodicity
)
{
	native_chi_periodicity_ = native_chi_periodicity;
}

bool
PNatRotOptEPositionData::count_rotamer_as_native( PNatRotOptERotamerDataOP rotamer ) const
{
	static Real const TOLERANCE = 40.0;
	for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) {
		if ( std::abs( numeric::mod( rotamer->chi()[ ii ] - native_chi_[ ii ], native_chi_periodicity_[ ii ] )) > TOLERANCE ) {
			return false;
		}
	}
	//std::cout << "Counting as native:";
	//for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) {
	//	std::cout << " ( " << rotamer->chi()[ ii ] << " " << native_chi_[ ii ] << " " << native_chi_periodicity_[ ii ] << " )";
	//}
	//std::cout << std::endl;
	return true;
}

void
PNatRotOptEPositionData::set_rotamer_well_counts(
	utility::vector1< Size > const & rotamer_well_counts
)
{
	runtime_assert( rotamer_well_counts_.size() == 0 );
	runtime_assert( rotamer_well_counts.size() > 0 );
	rotamer_well_counts_ = rotamer_well_counts;
	n_wells_ = 1;
	for ( Size ii = 1; ii <= rotamer_well_counts_.size(); ++ii ) {
		n_wells_ *= rotamer_well_counts_[ ii ];
	}
}

void
PNatRotOptEPositionData::add_rotamer_line_data( PNatRotOptERotamerDataOP rot_in )
{
	data_.push_back( rot_in );
}

PNatRotOptERotamerDataOPs &
PNatRotOptEPositionData::data()
{
	return data_;
}

PNatRotOptERotamerDataOPs const &
PNatRotOptEPositionData::data() const
{
	return data_;
}

PNatRotOptERotamerDataOPs::const_iterator
PNatRotOptEPositionData::rotamer_data_begin() const
{
	return data_.begin();
}
PNatRotOptERotamerDataOPs::const_iterator
PNatRotOptEPositionData::rotamer_data_end() const
{
	return data_.end();
}

core::chemical::AA PNatRotOptEPositionData::aa() const { return aa_; }
core::chemical::AA & PNatRotOptEPositionData::aa() { return aa_; }

Real PNatRotOptEPositionData::phi() const { return phi_; }
Real PNatRotOptEPositionData::psi() const { return psi_; }
Real & PNatRotOptEPositionData::phi() { return phi_; }
Real & PNatRotOptEPositionData::psi() { return psi_; }


Size
PNatRotOptEPositionData::rotamer_index_2_well_id(
	utility::vector1< Size > const & rotamer_index
) const
{
	if ( rotamer_well_counts_.size() == 0 ) return 0;
	Size well_id = rotamer_index[ 1 ] - 1;
	for ( Size ii = 2; ii <= rotamer_well_counts_.size(); ++ii ) {
		well_id *= rotamer_well_counts_[ ii ];
		well_id += rotamer_index[ ii ] - 1;
	}
	return well_id + 1;
}

Size
PNatRotOptEPositionData::rotamer_index_2_well_id(
	utility::LexicographicalIterator const & lexiter
) const
{
	if ( rotamer_well_counts_.size() == 0 ) return 0;
	Size well_id = lexiter[ 1 ] - 1;
	for ( Size ii = 2; ii <= rotamer_well_counts_.size(); ++ii ) {
		well_id *= rotamer_well_counts_[ ii ];
		well_id += lexiter[ ii ] - 1;
	}
	return well_id + 1;
}

bool
PNatRotOptEPositionData::is_native_rotamer_well(
	utility::vector1< Size > const & rotamer_index
) const
{
	for ( Size ii = 1; ii <= native_rotamer_index_.size(); ++ii ) {
		if ( rotamer_index[ ii ] != native_rotamer_index_[ ii ] ) return false;
	}
	return true;
}

bool
PNatRotOptEPositionData::is_native_rotamer_well(
	utility::LexicographicalIterator const & lexiter
) const
{
	for ( Size ii = 1; ii <= native_rotamer_index_.size(); ++ii ) {
		if ( lexiter[ ii ] != native_rotamer_index_[ ii ] ) return false;
	}
	return true;
}


//
// ------------------- PNatStructureOptEData-----------------------//
//

Real const PNatStructureOptEData::high_entropy_rms_cutoff_( 4.0 );
Real PNatStructureOptEData::nativeness_rms_low_( 1.0 ); // Above this rms, nativeness starts to decline
Real PNatStructureOptEData::nativeness_rms_high_( 2.5 ); // Above this rms, nativeness is zero


PNatStructureOptEData::PNatStructureOptEData()
:
	total_residue_( 0 ),
	n_top_natives_to_score_( 1 ),
	n_high_entropy_decoys_( 0 ),
	normalize_decoy_stddev_( false ),
	initial_decoy_stddev_( 1.0 ),
	nativeness_sum_( 0.0 )
{}

PNatStructureOptEData::~PNatStructureOptEData()
{}


Real
PNatStructureOptEData::get_score(
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const dummy1,
	int const dummy2,
	EnergyMap const & fixed_terms,
	ScoreTypes const & dummy3 /*free_score_list*/,
	ScoreTypes const & fixed_score_list
) const
{
	return process_score(
		TR, false, component_weights, vars,
		dE_dvars, num_energy_dofs, dummy1, dummy2,
		fixed_terms, dummy3, fixed_score_list );
}

//{
//      using namespace core::optimization;
//      using namespace utility;
//      using namespace basic::options;
//      using namespace basic::options::OptionKeys;
//
//      static Real const inv_kT( option[ optE::inv_kT_natstruct ] );
//      Real score_factor = option[ optE::approximate_decoy_entropy ].user() ? 1 : total_residue_;
//
//      if ( decoys_.size() == 0 || natives_.size() == 0 ) return 0.0;
//
//      Real const ln_nhighentropy_decoys = n_high_entropy_decoys_ > 0 ? std::log( (Real) n_high_entropy_decoys_ ) : 0;
//
//      utility::vector1< Real > decoy_energies( decoys_.size(), 0.0 );
//      utility::vector1< Real > native_energies( natives_.size(), 0.0 );
//      for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
//              for ( Size jj = 1; jj <= natives_.size(); ++jj ) {
//                      native_energies[ jj ] += vars[ ii ]  * natives_[ jj ]->free_data()[ ii ];
//              }
//              for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
//                      decoy_energies[ jj ] += vars[ ii ] * decoys_[ jj ]->free_data()[ ii ];
//              }
//      }
//      for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
//              for ( Size jj = 1; jj <= natives_.size(); ++jj ) {
//                      native_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * natives_[ jj ]->fixed_data()[ ii ];
//              }
//              for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
//                      decoy_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * decoys_[ jj ]->fixed_data()[ ii ];
//              }
//      }
//
//      utility::vector1< Size > top_natives( n_top_natives_to_score_, 0 );
//      if ( n_top_natives_to_score_ == 1 ) {
//              top_natives[ 1 ] = arg_min( native_energies );
//      } else {
//              arg_least_several( native_energies, top_natives );
//      }
//  /*if ( top_natives.size() == 0  ) {
//              std::cerr << "No top natives?! "<< tag() << std::endl << "E:";
//    for ( Size ii = 1; ii <= native_energies.size(); ++ii ) {
//      std::cerr << " " << native_energies[ ii ];
//    }
//    std::cerr << std::endl;
//      }*/
//
//      if ( normalize_decoy_stddev_ ) {
//
//              Real sd = numeric::statistics::std_dev( decoy_energies.begin(), decoy_energies.end(), Real( 0.0 ));
//
//              /// Avoid dividing by zero -- but what is the appropriate scale factor if the SD is measured 0?
//              Real const scaling_factor = sd == 0.0 ? 1 : initial_decoy_stdev_ < sd : initial_decoy_stddev_ / sd : 1.0; /// APL TEMP -- do not prevent the compression of decoy stdev; only prevent the expansion.
//              for ( Size jj = 1; jj <= decoy_energies.size(); ++jj ) {
//                      decoy_energies[ jj ] *= scaling_factor;
//              }
//              for ( Size jj = 1; jj <= native_energies.size(); ++jj ) {
//                      native_energies[ jj ] *= scaling_factor;
//              }
//      }
//
//      /// Alpha is the expansion factor of conformation space as a function of the number of residues.
//      /// num states ~ alpha ^ N
//      /// therefore the entropy scales as N log alpha.
//      if ( option[ optE::approximate_decoy_entropy ].user() ) {
//              /// assume alpha constant over the course of an optE run
//              static Real const neg_lnalpha = -std::log( option[ optE::approximate_decoy_entropy ] );
//              Real const energy_offset = neg_lnalpha * total_residue_ + ln_nhighentropy_decoys;
//              Real const offset = ( energy_offset < 0 ) ? energy_offset : 0;
//              for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
//                      if ( decoys_[ jj ]->rms() > high_entropy_rms_cutoff_ ) {
//                              decoy_energies[ jj ] += offset;
//                      }
//              }
//      }
//
//      /// Now normalize for numerical stability by setting the best native or decoy to have an energy of 0
//      /// and set everything else to be a delta from there.
//      Real best_native_energy = native_energies[ top_natives[ 1 ] ];
//      Real const best_decoy_energy = min( decoy_energies );
//      Real const best_energy =  best_native_energy < best_decoy_energy ? best_native_energy : best_decoy_energy;
//
//      for ( Size ii = 1; ii <= top_natives.size(); ++ii ) {
//              native_energies[ top_natives[ ii ] ] -= best_energy;
//      }
//      for ( Size ii = 1; ii <= decoys_.size(); ++ii ) {
//              decoy_energies[ ii ] -= best_energy;
//      }
//
//
//      Real numerator(0.0), partition(0.0);
//      Multivec dpartition( vars.size(), 0.0 ), dnumerator( vars.size(), 0.0 );
//
//      for ( Size ii = 1; ii <= top_natives.size(); ++ii ) {
//              Real exp_term( std::exp( -1.0 * inv_kT *  native_energies[ top_natives[ ii ] ] ) );
//              // iwd Limit the improbability of each native to 1 in a million
//              // This prevents numerator ~ 0, which causes NANs and INFs in the derivatives
//              // It also limits the "force" that any one structure can exert on the minimization.
//              exp_term = std::max( 1e-6, exp_term );
//
//              partition += exp_term;
//              numerator += exp_term;
//              // partitions for energy derivatives
//              for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
//                      // note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
//                      Real e_dof_deriv( -1.0 * inv_kT * score_factor * natives_[ top_natives[ ii ] ]->free_data()[ e_dof ] * exp_term );
//                      dpartition[ e_dof ] += e_dof_deriv;
//                      dnumerator[ e_dof ] += e_dof_deriv;
//              }
//      }
//
//      for( Size ii(1); ii <= decoys_.size(); ++ii ) {
//
//              Real const exp_term( std::exp( -1.0 * inv_kT *  decoy_energies[ ii ] ) );
//              partition += exp_term;
//
//              // partitions for energy derivatives
//              for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
//                      // note for derivatives: d/dw( e^{-(E*w+...)/kT} ) = -E/kT * e^-(E*w+...)
//                      Real e_dof_deriv( -1.0 * inv_kT * score_factor * decoys_[ ii ]->free_data()[ e_dof ] * exp_term );
//                      dpartition[ e_dof ] += e_dof_deriv;
//              }
//      }
//
//      // accumulate to passed-in derivative sums -- excludes reference energies
//      //std::cout << "vars (dvars): ";
//      for ( Size dof(1); dof <= num_energy_dofs; ++dof ) {
//              dE_dvars[ dof ] += component_weights[ prob_native_structure ] * (dpartition[ dof ] / partition - dnumerator[ dof ] / numerator);
//              //std::cout << " " << vars[ dof ] << "(" << dpartition[ dof ] / partition - dnumerator[ dof ] / numerator << ")";
//      }
//
//      //std::cout << "struct score: " << best_native << " " << best_native_energy << " " ;
//      //std::cout << nat_exp_term << " " << partition << " : " <<  ( -1.0 * total_residue_ * std::log( numerator / partition ) ) << std::endl;
//
//
//      /*if ( partition == 0 ) {
//              std::cerr << tag() << " partition of zero: inv_kT " << inv_kT;
//    for ( Size ii = 1; ii <= decoys_.size(); ++ii ) {
//                      std::cerr << " d: " << decoy_energies[ ii ] << " " << std::exp( -1.0 * inv_kT *  decoy_energies[ ii ] ) << " ";
//                      for ( Size jj = 1; jj <= num_energy_dofs; ++jj ) {
//                              std::cerr << free_score_list[ jj ] << ": " << decoys_[ ii ]->free_data()[ jj ] << " ";
//                      }
//                      std::cerr << std::endl;
//              }
//              for ( Size ii = 1; ii <= top_natives.size(); ++ii ) {
//                      std::cerr << " top nat: " << native_energies[ top_natives[ ii ] ] << " " << std::exp( -1.0 * inv_kT *  native_energies[ top_natives[ ii ] ] ) << " ";
//                      for ( Size jj = 1; jj <= num_energy_dofs; ++jj ) {
//                              std::cerr << free_score_list[ jj ] << ": " << natives_[ top_natives[ ii ] ]->free_data()[ jj ] << " ";
//                      }
//                      std::cerr << std::endl;
//              }
//              std::cerr << std::endl;
//    for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
//                      std::cerr << " " << free_score_list[ ii ] << ": " << vars[ ii ] << std::endl;
//              }
//              return 0.0;
//      }*/
//
//
//      return ( -1.0 * component_weights[ prob_native_structure ] * score_factor * std::log( numerator / partition ) );
//}


void
PNatStructureOptEData::print_score(
	std::ostream & ostr,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const dummy1,
	int const dummy2,
	EnergyMap const & fixed_terms,
	ScoreTypes const & dummy3,
	ScoreTypes const & fixed_score_list
) const
{
	process_score(
		ostr, true, component_weights, vars,
		dE_dvars, num_energy_dofs, dummy1, dummy2,
		fixed_terms, dummy3, fixed_score_list );
}

Real
PNatStructureOptEData::process_score(
	std::ostream & ostr,
	bool print,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const &,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core::optimization;
	using namespace utility;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static Real const inv_kT( option[ optE::inv_kT_natstruct ] );
	Real score_factor = option[ optE::approximate_decoy_entropy ].user() ? 1 : total_residue_;

	if ( decoys_.size() == 0 || natives_.size() == 0 ) {
		if ( print ) {
			ostr << "PNATSTRUCT " << total_residue_ << " nnat: " << natives_.size() << " ndec: " << decoys_.size() << tag() << "\n";
		}
		return 0.0; // wtf?
	}

	//TR << "PNatStructureOptEData::print_score" << std::endl;

	/// Compute the energies for natives and decoys given the set of input weights;
	utility::vector1< Real > decoy_energies( decoys_.size(), 0.0 );
	utility::vector1< Real > native_energies( natives_.size(), 0.0 );
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		for ( Size jj = 1; jj <= natives_.size(); ++jj ) {
			native_energies[ jj ] += vars[ ii ]  * natives_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
			decoy_energies[ jj ] += vars[ ii ] * decoys_[ jj ]->free_data()[ ii ];
		}
	}
	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		for ( Size jj = 1; jj <= natives_.size(); ++jj ) {
			native_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * natives_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
			decoy_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * decoys_[ jj ]->fixed_data()[ ii ];
		}
	}


	// avoid allocating these arrays if you're not going to use them;
	utility::vector1< Real > noentropy_decoy_energies;
	utility::vector1< Real > noentropy_native_energies;
	if ( print ) {
		noentropy_decoy_energies = decoy_energies;
		noentropy_native_energies = native_energies;

		//// Normalize the "noentropy" energies
		Real const noentropy_best_native_energy = min( noentropy_native_energies );
		Real const noentropy_best_decoy_energy = min( noentropy_decoy_energies );
		bool const noentropy_native_beats_decoys( noentropy_best_native_energy < noentropy_best_decoy_energy );
		Real const noentropy_best_energy =  noentropy_native_beats_decoys ? noentropy_best_native_energy : noentropy_best_decoy_energy;

		for ( Size ii = 1; ii <= native_energies.size(); ++ii ) {
			noentropy_native_energies[ ii ] -= noentropy_best_energy;
		}
		for ( Size ii = 1; ii <= decoys_.size(); ++ii ) {
			noentropy_decoy_energies[ ii ] -= noentropy_best_energy;
		}

	}

	if ( normalize_decoy_stddev_ ) {

		Real sd = numeric::statistics::std_dev( decoy_energies.begin(), decoy_energies.end(), Real( 0.0 ));

		/// Avoid dividing by zero -- but what is the appropriate scale factor if the SD is measured 0?
		Real const scaling_factor = sd == 0.0 ? 1 : initial_decoy_stddev_ / sd;
		for ( Size jj = 1; jj <= decoy_energies.size(); ++jj ) {
			decoy_energies[ jj ] *= scaling_factor;
		}
		for ( Size jj = 1; jj <= native_energies.size(); ++jj ) {
			native_energies[ jj ] *= scaling_factor;
		}
	}

	/// Alpha is the expansion factor of conformation space as a function of the number of residues.
	/// num states ~= alpha ^ N
	/// therefore the entropy scales as N log( alpha ).

	if ( option[ optE::approximate_decoy_entropy ].user() ) {
		/// assume alpha is constant over the course of an optE run
		static Real const neg_lnalpha = -std::log( option[ optE::approximate_decoy_entropy ] );

		Real const energy_offset = neg_lnalpha * total_residue_ + std::log( (Real)n_high_entropy_decoys_ );
		Real const offset = ( energy_offset < 0 ) ? energy_offset : 0;
		for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
			if ( decoys_[ jj ]->rms() > high_entropy_rms_cutoff_ ) {
				decoy_energies[ jj ] += offset;
			}
		}
	}

	if ( option[ optE::ramp_nativeness ] ) {
		Real native_entropy_penalty = std::log( nativeness_sum_ );
		for ( Size jj = 1; jj <= native_energies.size(); ++jj ) {
			native_energies[ jj ] += native_entropy_penalty;
		}
	}

	utility::vector1< Size > top_natives( n_top_natives_to_score_, 0 );
	if ( n_top_natives_to_score_ == 1 ) {
		top_natives[ 1 ] = arg_min( native_energies );
	} else if ( option[ optE::approximate_decoy_entropy ].user() ) {
		/// ALL nativish structures are getting scored; it's only important
		/// that the best native structure is in position 1 of the top_natives array.
		top_natives[ 1 ] = arg_min( native_energies );
		for ( Size ii = 2; ii <= n_top_natives_to_score_; ++ii ) top_natives[ ii ] = ii;
		top_natives[ top_natives[ 1 ] ] = 1;
	} else {
		arg_least_several( native_energies, top_natives );
	}

	/// Normalize the scores so that the energy of the best structure
	/// is zero, and all other energies are positive.
	Real best_native_energy = native_energies[ top_natives[ 1 ] ];
	Real const best_decoy_energy = min( decoy_energies );
	bool const native_beats_decoys( best_native_energy < best_decoy_energy );
	Real const best_energy =  native_beats_decoys ? best_native_energy : best_decoy_energy;

	for ( Size ii = 1; ii <= native_energies.size(); ++ii ) native_energies[ ii ] -= best_energy;
	for ( Size ii = 1; ii <= decoys_.size(); ++ii ) decoy_energies[ ii ] -= best_energy;


	Real numerator(0.0), partition(0.0);
	Multivec dpartition( vars.size(), 0.0 ), dnumerator( vars.size(), 0.0 );

	for ( Size ii = 1; ii <= top_natives.size(); ++ii ) {
		Real exp_term( std::exp( -1.0 * inv_kT *  native_energies[ top_natives[ ii ] ] ) );
		// iwd Limit the improbability of each native to 1 in a million
		// This prevents numerator ~ 0, which causes NANs and INFs in the derivatives
		// It also limits the "force" that any one structure can exert on the minimization.
		exp_term = std::max( 1e-6, exp_term );

		Real iinativeness = nativeness( natives_[ top_natives[ ii ] ]->rms() );
		partition += exp_term;
		numerator += exp_term * iinativeness;
		// partitions for energy derivatives
		for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
			// note for derivatives: d/dw( e^-(E*w+...) ) = -E * e^-(E*w+...)
			Real e_dof_deriv( -1.0 * inv_kT * score_factor * natives_[ top_natives[ ii ] ]->free_data()[ e_dof ] * exp_term );
			dpartition[ e_dof ] += e_dof_deriv;
			dnumerator[ e_dof ] += iinativeness * e_dof_deriv ;
		}
	}

	for( Size ii(1); ii <= decoys_.size(); ++ii ) {

		Real const exp_term( std::exp( -1.0 * inv_kT * decoy_energies[ ii ] ) );
		partition += exp_term;

		// partitions for energy derivatives
		for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
			// note for derivatives: d/dw( e^{-(E*w+...)/kT} ) = -E/kT * e^-(E*w+...)
			Real e_dof_deriv( -1.0 * inv_kT * score_factor * decoys_[ ii ]->free_data()[ e_dof ] * exp_term );
			dpartition[ e_dof ] += e_dof_deriv;
		}
	}

	// accumulate to passed-in derivative sums -- excludes reference energies
	//std::cout << "vars (dvars): ";
	for ( Size dof(1); dof <= num_energy_dofs; ++dof ) {
		dE_dvars[ dof ] += dpartition[ dof ] / partition - dnumerator[ dof ] / numerator;
	}

	if ( print ) {
		ostr << "PNATSTRUCT " << total_residue_ << " nnat: " << natives_.size() << " ndec: " << decoys_.size();
		ostr << " n<d? " << native_beats_decoys << " p: " << numerator / partition;
		ostr << " -lnp: " << -1 * std::log( numerator / partition ) << " -wlnp " <<  -1.0 * score_factor * std::log( numerator / partition );
		ostr << " -compwt_lnp " << component_weights[ prob_native_structure ] *  -1.0 * score_factor * std::log( numerator / partition );
		ostr << " " << tag() << "\n";

		for ( Size ii = 1; ii <= natives_.size(); ++ii ) {
			ostr << "DECDISCRIM NATIVE " << tag() << " rms: " << natives_[ ii ]->rms() << " sc(w/o_etropy): " << noentropy_native_energies[ ii ] << " sc(w/entropy): "  << native_energies[ ii ] << " " << natives_[ ii ]->tag() << "\n";
		}
		for ( Size ii = 1; ii <= decoys_.size(); ++ii ) {
			ostr << "DECDISCRIM DECOY  " << tag() << " rms: " << decoys_[ ii ]->rms() << " sc(w/o_entropy): " << noentropy_decoy_energies[ ii ] << " sc(w/entropy): "  << decoy_energies[ ii ] << " " << decoys_[ ii ]->tag()  << "\n";
		}
	}
	return ( -1.0 * component_weights[ prob_native_structure ] * score_factor * std::log( numerator / partition ) );
}


void
PNatStructureOptEData::range(
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list,
	EnergyMap & lower_bound,
	EnergyMap & upper_bound
) const
{
	for ( Size jj = 1; jj <= natives_.size(); ++jj ) {
		update_range( natives_[ jj ], free_score_list, fixed_score_list, lower_bound, upper_bound );
	}
	for ( Size jj = 1; jj <= decoys_.size(); ++jj ) {
		update_range( decoys_[ jj ],  free_score_list, fixed_score_list, lower_bound, upper_bound );
	}
}

Size
PNatStructureOptEData::size() const
{
	return natives_.size() + decoys_.size();
}


OptEPositionDataType
PNatStructureOptEData::type() const
{
	return prob_native_structure;
}


void
PNatStructureOptEData::write_to_file( std::ofstream & /*outfile*/ ) const
{}


void
PNatStructureOptEData::read_from_file( std::ifstream & /*infile*/ )
{}


void
PNatStructureOptEData::write_to_binary_file( std::ofstream & /*outfile*/ ) const
{}


void
PNatStructureOptEData::read_from_binary_file( std::ifstream & /*infile*/ )
{}


Size
PNatStructureOptEData::memory_use() const
{
	Size total = sizeof( PNatStructureOptEData ) +
		sizeof( SingleStructureData ) * natives_.size() +
		sizeof( SingleStructureData ) * decoys_.size();
	if ( natives_.size() > 0 ) {
		total += sizeof( Real ) * ( natives_[ 1 ]->free_data().size() + natives_[ 1 ]->fixed_data().size() ) * natives_.size();
	}
	if ( decoys_.size() > 0 ) {
		total += sizeof( Real ) * ( decoys_[ 1 ]->free_data().size() + decoys_[ 1 ]->fixed_data().size() ) * decoys_.size();
	}
	return total;
}


#ifdef USEMPI

void
PNatStructureOptEData::send_to_node( int const destination_node, int const tag ) const
{
	//std::cout << "sending total_residue to node " << destination_node << std::endl;
	/// 1. total residue
	Size total_residue = total_residue_;
	MPI_Send( & total_residue, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 1b. n_top_natives_to_score_
	Size n_top = n_top_natives_to_score_;
	MPI_Send( & n_top, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	// 2a. n natives
	//std::cout << "sending nnatives to node " << destination_node << std::endl;
	Size nnatives = natives_.size();
	MPI_Send( & nnatives, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 2b. n decoys
	//std::cout << "sending ndecoys to node " << destination_node << std::endl;
	Size ndecoys = decoys_.size();
	MPI_Send( & ndecoys, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 2c. n_high_entropy_decoys_
	//std::cout << "sending ndecoys to node " << destination_node << std::endl;
	Size n_high_entropy_decoys = n_high_entropy_decoys_;
	MPI_Send( & n_high_entropy_decoys, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 2d. normalize_decoy_stddev_
	int normalize_decoy_stddev = static_cast< int > (normalize_decoy_stddev_);
	MPI_Send( & normalize_decoy_stddev, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );

	/// 2e. initial_decoy_stddev_
	Real initial_decoy_stddev = initial_decoy_stddev_;
	MPI_Send( & initial_decoy_stddev, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 2f. nativeness_sum_
	Real nativeness_sum( nativeness_sum_ );
	MPI_Send( & nativeness_sum, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );


	if ( nnatives == 0 || ndecoys == 0 ) {
		OptEPositionData::send_to_node( destination_node, tag );
		return;
	}

	/// 3. n free
	Size n_free = decoys_[ 1 ]->free_data().size();
	//std::cout << "sending n_free to node " << destination_node << " " << n_free << std::endl;
	MPI_Send( & n_free, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 4. n fixed
	Size n_fixed = decoys_[ 1 ]->fixed_data().size();
	//std::cout << "sending n_fixed to node " << destination_node  << " " << n_fixed << std::endl;
	MPI_Send( & n_fixed, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// Send natives, then send decoys
	Real * rms = new Real[ nnatives ];
	Real * free_data = new Real[ n_free * nnatives ];
	Real * fixed_data = new Real[ n_fixed * nnatives ];
	for ( Size ii = 1; ii <= nnatives; ++ ii ) {
		rms[ (ii-1) ] = natives_[ ii ]->rms();
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ] = natives_[ ii ]->free_data()[ jj ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ] = natives_[ ii ]->fixed_data()[ jj ];
		}
	}

	MPI_Send( rms, nnatives, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//std::cout << "sending native free_data to node " << destination_node << " " << free_data <<  std::endl;
	/// 5. native free data
	MPI_Send( free_data, nnatives * n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//std::cout << "sending native fixed_data to node " << destination_node << " " << fixed_data <<  std::endl;
	/// 6. fixed data
	MPI_Send( fixed_data, nnatives * n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//std::cout << "Sent -- about to delete data" << std::endl;

	/// now send decoys
	Real * decoy_rms = new Real[ ndecoys ];
	Real * decoy_free_data = new Real[ n_free * ndecoys ];
	Real * decoy_fixed_data = new Real[ n_fixed * ndecoys ];
	for ( Size ii = 1; ii <= ndecoys; ++ ii ) {
		decoy_rms[ (ii-1) ] = decoys_[ ii ]->rms();
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			decoy_free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ] = decoys_[ ii ]->free_data()[ jj ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			decoy_fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ] = decoys_[ ii ]->fixed_data()[ jj ];
		}
	}

	MPI_Send( decoy_rms, ndecoys, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 7. decoy free data
	//std::cout << "sending decoy free_data to node " << destination_node << std::endl;
	MPI_Send( decoy_free_data, ndecoys * n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 8. decoy fixed data
	//std::cout << "sending decoy fixed_data to node " << destination_node << std::endl;
	MPI_Send( decoy_fixed_data, ndecoys * n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	delete [] rms;
	delete [] free_data;
	delete [] fixed_data;

	delete [] decoy_rms;
	delete [] decoy_free_data;
	delete [] decoy_fixed_data;

	/// 9. structure tags
	/// 9a. natives
	for ( Size ii = 1; ii <= natives_.size(); ++ii ) {
		IterativeOptEDriver::send_string_to_node( destination_node, natives_[ ii ]->tag() );
	}
	/// 9b. decoys
	for ( Size ii = 1; ii <= decoys_.size(); ++ii ) {
		IterativeOptEDriver::send_string_to_node( destination_node, decoys_[ ii ]->tag() );
	}

	OptEPositionData::send_to_node( destination_node, tag );

}

void
PNatStructureOptEData::receive_from_node( int const source_node, int const tag )
{
	MPI_Status stat;
	//TR << "PNatStructureOptEData::Recieving data from node... " << source_node << std::endl;

	/// 1. total residue
	Size total_residue( 0 );
	MPI_Recv( & total_residue, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	total_residue_ = total_residue;

	/// 1b. n_top_natives_to_score_
	Size n_top;
	MPI_Recv( & n_top, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	n_top_natives_to_score_ = n_top;

	/// 2a. n natives
	Size nnatives( 0 );
	MPI_Recv( & nnatives, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 2b. n decoys
	Size ndecoys( 0 );
	MPI_Recv( & ndecoys, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 2c. n_high_entropy_decoys_
	//std::cout << "sending ndecoys to node " << destination_node << std::endl;
	Size n_high_entropy_decoys( 0 );
	MPI_Recv( & n_high_entropy_decoys, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	n_high_entropy_decoys_ = n_high_entropy_decoys;

	/// 2d. normalize_decoy_stddev_
	int normalize_decoy_stddev( 0 );
	MPI_Recv( & normalize_decoy_stddev, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
	normalize_decoy_stddev_ = static_cast< bool > (normalize_decoy_stddev);

	/// 2e. initial_decoy_stddev_
	Real initial_decoy_stddev( 1.0 );
	MPI_Recv( & initial_decoy_stddev, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	initial_decoy_stddev_ = initial_decoy_stddev;

	/// 2f. nativeness_sum_
	Real nativeness_sum( 0.0 );
	MPI_Recv( & nativeness_sum, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	nativeness_sum_ = nativeness_sum;

	if ( nnatives == 0 || ndecoys == 0 ) {
		OptEPositionData::receive_from_node( source_node, tag );
		return;
	}
	natives_.reserve( nnatives );
	decoys_.reserve( ndecoys );

	/// 3. n free
	Size n_free( 0 );
	MPI_Recv( & n_free, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 4. n fixed
	Size n_fixed( 0 );
	MPI_Recv( & n_fixed, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );

	/// Recieve native data first, then decoys
	Real * rms = new Real[ nnatives ];
	Real * free_data = new Real[ n_free * nnatives ];
	Real * fixed_data = new Real[ n_fixed * nnatives ];

	MPI_Recv( rms, nnatives, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 5. free data
	MPI_Recv( free_data, nnatives * n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 6. fixed data
	MPI_Recv( fixed_data, nnatives * n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	utility::vector1< Real > free_data_v( n_free );
	utility::vector1< Real > fixed_data_v( n_fixed );
	for ( Size ii = 1; ii <= nnatives; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data_v[ jj ] = free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data_v[ jj ] = fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ];
		}
		natives_.push_back( new SingleStructureData( free_data_v, fixed_data_v ) );
		natives_[ ii ]->rms( rms[ ii-1 ] );
		//std::cout << "recieved native with rms: " << rms[ (ii-1) ] << std::endl;
	}

	delete [] rms; rms = 0;
	delete [] free_data; free_data = 0;
	delete [] fixed_data; fixed_data = 0;

	//// Now receive decoy data
	rms = new Real[ ndecoys ];
	free_data = new Real[ n_free * ndecoys ];
	fixed_data = new Real[ n_fixed * ndecoys ];

	MPI_Recv( rms, ndecoys, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 5. free data
	MPI_Recv( free_data, ndecoys * n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 6. fixed data
	MPI_Recv( fixed_data, ndecoys * n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	for ( Size ii = 1; ii <= ndecoys; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data_v[ jj ] = free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data_v[ jj ] = fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ];
		}
		decoys_.push_back( new SingleStructureData( free_data_v, fixed_data_v ) );
		decoys_[ ii ]->rms( rms[ ii - 1 ] );
		//std::cout << "recieved decoy  with rms: " << rms[ (ii-1) ] << std::endl;
	}


	delete [] rms;
	delete [] free_data;
	delete [] fixed_data;

	/// 9. structure tags
	/// 9a. natives
	for ( Size ii = 1; ii <= natives_.size(); ++ii ) {
		natives_[ ii ]->tag( IterativeOptEDriver::receive_string_from_node( source_node ));
	}
	/// 9b. decoys
	for ( Size ii = 1; ii <= decoys_.size(); ++ii ) {
		decoys_[ ii ]->tag( IterativeOptEDriver::receive_string_from_node( source_node ));
	}


	OptEPositionData::receive_from_node( source_node, tag );

}
#endif


void
PNatStructureOptEData::add_native( SingleStructureDataOP native )
{
	natives_.push_back( native );
	if ( basic::options::option[ basic::options::OptionKeys::optE::ramp_nativeness ] ) {
		nativeness_sum_ += nativeness( native->rms() );
		n_top_natives_to_score_ = natives_.size();
	}
}

void
PNatStructureOptEData::add_decoy( SingleStructureDataOP decoy )
{
	decoys_.push_back( decoy );
	if ( decoy->rms() > high_entropy_rms_cutoff_ ) {
		++n_high_entropy_decoys_;
	}
}

void
PNatStructureOptEData::set_total_residue( Size total_residue ) {
	total_residue_ = total_residue;
	//total_residue_ = 1; // This line removes the total-residue reweighting.
}

void
PNatStructureOptEData::n_top_natives_to_score( Size n_top )
{
	if ( n_top < 1 ) {
		std::cerr << "ERROR: PNatStructureOptEData will not score fewer than one top native!" << std::endl;
		utility_exit();
	}
	n_top_natives_to_score_ = n_top;
}

Size
PNatStructureOptEData::n_top_natives_to_score() const
{
	return n_top_natives_to_score_;
}

void
PNatStructureOptEData::set_normalize_decoy_stddev( bool setting )
{
	normalize_decoy_stddev_ = setting;
}

void
PNatStructureOptEData::set_initial_decoy_stddev( Real setting )
{
	initial_decoy_stddev_ = setting;
}

Real
PNatStructureOptEData::nativeness( Real rms ) const
{
	static bool ramp_nativeness( basic::options::option[ basic::options::OptionKeys::optE::ramp_nativeness ] );
	return ! ramp_nativeness ? 1.0 :
		rms < nativeness_rms_low_ ? 1.0 : ( rms > nativeness_rms_high_ ? 0.0 :
		1 - ( rms - nativeness_rms_low_ ) / ( nativeness_rms_high_ - nativeness_rms_low_ ) );
}

void
PNatStructureOptEData::set_nativeness_low( Real nativeness_rms_low )
{
	nativeness_rms_low_ = nativeness_rms_low;
}

void
PNatStructureOptEData::set_nativeness_high( Real nativeness_rms_high )
{
	nativeness_rms_high_= nativeness_rms_high;

}

Real
PNatStructureOptEData::nativeness_low()
{
	return nativeness_rms_low_;
}

Real
PNatStructureOptEData::nativeness_high()
{
	return nativeness_rms_high_;
}


/////////////////////////////////////////////
//// ConstraintedOptimizationWeightFunc /////
/////////////////////////////////////////////

ConstraintedOptimizationWeightFunc::ConstraintedOptimizationWeightFunc() {}

ConstraintedOptimizationWeightFunc::ConstraintedOptimizationWeightFunc(
	ScoreTypes const & score_list
)
:
	free_terms_( score_list ),
	free_term_constraints_( score_list.size() )
{

}

ConstraintedOptimizationWeightFunc::~ConstraintedOptimizationWeightFunc()
{}


/// @details Constraint input file format:
/// ON <score_type> <min_val> <max_val> <spring_constant>
/// OFF <score_type>
/// Where score_type is a string; min_val, max_val, and spring_constant are reals
void
ConstraintedOptimizationWeightFunc::initialize_constraints_from_file( std::ifstream & infile )
{
	std::string on_off;
	std::string term_name;
	ScoreType term( ScoreType(1) );
	Real min_val;
	Real max_val;
	Real spring_constant;

	utility::vector1< bool > score_types_seen( core::scoring::n_score_types, false );
	while ( infile ) {
		infile >> on_off;
		if ( ! infile ) break;

		infile >> term_name;
		if ( core::scoring::ScoreTypeManager::is_score_type( term_name ) ) {
			term = core::scoring::ScoreTypeManager::score_type_from_name( term_name );
		} else {
			std::cerr << "ERROR reading '" << term_name << "' in weight-constraint input file as a ScoreType." << std::endl;
			utility_exit();
		}

		if ( score_types_seen[ term ] ) {
			std::cerr << "ERROR initializing weight-constraints from file: term '" << term_name << "' repeated in input file" << std::endl;
			utility_exit();
		}

		score_types_seen[ term ] = true;
		Size term_index( 0 );
		for ( Size ii = 1; ii <= free_terms_.size(); ++ii ) {
			if ( free_terms_[ ii ] == term ) {
				term_index = ii;
				break;
			}
		}

		if ( term_index == 0 ) {
			std::cerr << "Error initializing constraints for term: " << term_name << " which does not appear in the free_terms_ list";
			continue;
		}

		if ( on_off == "OFF" ) {
			free_term_constraints_[ term_index ].active_ = false;
		} else if ( on_off == "ON" ) {
			infile >> min_val >> max_val >> spring_constant;
			free_term_constraints_[ term_index ].active_ = true;
			free_term_constraints_[ term_index ].min_weight_ = min_val;
			free_term_constraints_[ term_index ].max_weight_ = max_val;
			free_term_constraints_[ term_index ].spring_constant_ = spring_constant;
		} else {
			std::cerr << "ERROR: In weight constraint file, expected 'ON' or 'OFF' but read '" << on_off << "' instead." << std::endl;
			utility_exit();
		}
	}
}

Real
ConstraintedOptimizationWeightFunc::get_score(
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const,
	int const,
	EnergyMap const &,
	ScoreTypes const &,
	ScoreTypes const &
) const
{
	using namespace core::optimization;
	Real total( 0 );
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		if ( ! free_term_constraints_[ ii ].active_ ) continue;

		WeightRangeConstraint wrc = free_term_constraints_[ ii ];

		if ( vars[ ii ] > wrc.max_weight_ ) {
			total += ( vars[ ii ] - wrc.max_weight_ )*( vars[ ii ] - wrc.max_weight_ ) * wrc.spring_constant_;
			dE_dvars[ ii ] += 2 * component_weights[ constrained_optimization_weight_func ] *
				( vars[ ii ] - wrc.max_weight_ ) * wrc.spring_constant_;
		} else if ( vars[ ii ] < wrc.min_weight_ ) {
			total += ( wrc.min_weight_ - vars[ ii ] )*( wrc.min_weight_ - vars[ ii ] ) * wrc.spring_constant_;
			dE_dvars[ ii ] += -2 * component_weights[ constrained_optimization_weight_func ] *
				( wrc.min_weight_ - vars[ ii ] ) * wrc.spring_constant_;
		}
	}
	return component_weights[ constrained_optimization_weight_func ] * total;
}

void
ConstraintedOptimizationWeightFunc::print_score(
	std::ostream & ostr,
	optimization::Multivec const & /*component_weights*/,
	optimization::Multivec const & vars,
	optimization::Multivec &,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const,
	int const,
	EnergyMap const &,
	ScoreTypes const & score_list,
	ScoreTypes const &
) const
{
	using namespace core::optimization;
	ostr << "CONSTR ";
	Real total( 0.0 );
	std::ostringstream sstr;
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {

		if ( ! free_term_constraints_[ ii ].active_ ) continue;

		WeightRangeConstraint wrc = free_term_constraints_[ ii ];

		if ( vars[ ii ] > wrc.max_weight_ ) {
			Real const pen( ( vars[ ii ] - wrc.max_weight_ )*( vars[ ii ] - wrc.max_weight_ ) * wrc.spring_constant_ );
			total += pen;
			sstr << score_list[ ii ] << ": " << vars[ ii ] << ", " << pen << "; ";
		} else if ( vars[ ii ] < wrc.min_weight_ ) {
			Real const pen( ( wrc.min_weight_ - vars[ ii ] )*( wrc.min_weight_ - vars[ ii ] ) * wrc.spring_constant_ );
			total += pen;
			sstr << score_list[ ii ] << " " << vars[ ii ] << ", " << pen << "; ";
		}
	}
	ostr << " total: " << total << " " << sstr.str() << "\n";
	return;
}


void
ConstraintedOptimizationWeightFunc::range(
	ScoreTypes const & ,
	ScoreTypes const & ,
	EnergyMap & ,
	EnergyMap &
) const
{
}

Size
ConstraintedOptimizationWeightFunc::size() const {
	return 1;
}



OptEPositionDataType
ConstraintedOptimizationWeightFunc::type() const
{
	return constrained_optimization_weight_func;
}


void
ConstraintedOptimizationWeightFunc::write_to_file( std::ofstream & /*outfile*/ ) const
{}


void
ConstraintedOptimizationWeightFunc::read_from_file( std::ifstream & /*infile*/ )
{}


void
ConstraintedOptimizationWeightFunc::write_to_binary_file( std::ofstream & /*outfile*/ ) const
{}


void
ConstraintedOptimizationWeightFunc::read_from_binary_file( std::ifstream & /*infile*/ )
{}


Size
ConstraintedOptimizationWeightFunc::memory_use() const
{
	Size total = sizeof( ConstraintedOptimizationWeightFunc );
	return total;
}


#ifdef USEMPI

void
ConstraintedOptimizationWeightFunc::send_to_node( int const /*destination_node*/, int const /*tag*/ ) const
{
	/// Stubbed out -- so far, ConstraintedOptimizationWeightFuncs are only used by node 0 and don't
	/// have to be passed around...

	utility_exit();
	//std::cout << "sending total_residue to node " << destination_node << std::endl;
	/// 1. max_weight_
	//double max_weight = max_weight_;
	//MPI_Send( & max_weight, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 2. min_weight_
	//double min_weight = min_weight_;
	//MPI_Send( & min_weight, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 3. spring_constant_
	//double spring_constant = spring_constant_;
	//MPI_Send( & spring_constant, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//OptEPositionData::send_to_node( destination_node, tag );

}

void
ConstraintedOptimizationWeightFunc::receive_from_node( int const /*source_node*/, int const /*tag*/ )
{
	utility_exit();
	//TR << "ConstraintedOptimizationWeightFunc::Recieving data from node... " << source_node << std::endl;

	//MPI_Status stat;

	/// 1. max_weight_
	//double max_weight;
	//MPI_Recv( & max_weight, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	//max_weight_ = max_weight;

	/// 2. min_weight_
	//double min_weight;
	//MPI_Recv( & min_weight, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	//min_weight_ = min_weight;

	/// 3. spring_constant_
	//double spring_constant;
	//MPI_Recv( & spring_constant, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	//spring_constant_ = spring_constant;

	//OptEPositionData::receive_from_node( source_node, tag );

}
#endif


//
// ------------------- DDGMutationOptEData-----------------------//
//

DDGMutationOptEData::DDGMutationOptEData()
:
	experimental_ddG_( 0 ),
	wt_aa_( core::chemical::aa_ala ),
	mut_aa_( core::chemical::aa_ala )
{}

DDGMutationOptEData::~DDGMutationOptEData()
{}

Real
DDGMutationOptEData::get_score(
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
	return process_score(
		TR, false, component_weights, vars,
		dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs,
		fixed_terms, free_score_list, fixed_score_list );
}


/*
Real
DDGMutationOptEData::get_score(
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const &,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core::optimization;
	using namespace utility;

	// if there are no structures to go through, return immediately
	if ( muts_.size() == 0 || wts_.size() == 0 ) return 0.0;

	// these vectors are sized to the number of structures there are for a wt name and mutant name;
	// they'll be used to determine which structure has the best energy
	utility::vector1< Real > wt_energies( wts_.size(), 0.0 );
	utility::vector1< Real > mut_energies( muts_.size(), 0.0 );

	// go through and come up with a total score for each structure in the wts_ and muts_ list
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
			wt_energies[ jj ] += vars[ ii ] * wts_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			mut_energies[ jj ] += vars[ ii ] * muts_[ jj ]->free_data()[ ii ];
		}
	}
	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
			wt_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * wts_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			mut_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * muts_[ jj ]->fixed_data()[ ii ];
		}
	}

	// I presume these are the reference energies that are being added in?
	// num_energy_dofs is the number of free, non-reference energy parameters in the run -ronj
	if ( num_ref_dofs != 0 ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
			wt_energies[ jj ] += vars[ num_energy_dofs + wt_aa_ ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			mut_energies[ jj ] += vars[ num_energy_dofs + mut_aa_ ];
		}
	}

	// Identify which wt and mut structure in the lists have the best energy
	// arg_min returns the index of the smallest value in a vector
	Size const best_wt = arg_min( wt_energies );
	Size const best_mut = arg_min( mut_energies );

	Real const best_wt_energy = wt_energies[ best_wt ];
	Real const best_mut_energy = mut_energies[ best_mut ];

	Real const predicted_ddG = best_mut_energy - best_wt_energy;
	Real const ddG_diff = predicted_ddG - experimental_ddG_;

	for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
		dE_dvars[ e_dof ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff * ( muts_[ best_mut ]->free_data()[ e_dof ] - wts_[ best_wt ]->free_data()[ e_dof ] );
	}

	if ( num_ref_dofs != 0 ) {
		dE_dvars[ num_energy_dofs + mut_aa_ ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff;
		dE_dvars[ num_energy_dofs + wt_aa_  ] -= 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff;
	}

	//TR << "DDGMutationOptEData::get_score(): predicted: " << ObjexxFCL::format::F(4,2,predicted_ddG)
	//	<< ", experimental: " << experimental_ddG_ << ", ddG_diff^2: " << ObjexxFCL::format::F(6,3,ddG_diff * ddG_diff)
	//	<< ", tag: " << tag() << std::endl;
	return component_weights[ ddG_mutation_correlation ] * ddG_diff * ddG_diff;
}*/

void
DDGMutationOptEData::print_score(
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
	process_score(
		ostr, true, component_weights, vars,
		dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs,
		fixed_terms, free_score_list, fixed_score_list );
}

Real
DDGMutationOptEData::process_score(
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

	utility::vector1< utility::vector1< Real > > wt_dG_energies( 20, utility::vector1 < Real > ( num_energy_dofs + fixed_score_list.size(), 0.0 ) );
	utility::vector1< utility::vector1< Real > > mut_dG_energies( 20, utility::vector1 < Real > ( num_energy_dofs + fixed_score_list.size(), 0.0 ) );

	#define CAP_FA_REP 1

	// go through and come up with a total score for each structure in the wts_ and muts_ list
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {

			// cap the fa_rep term at some value - this at least keeps it around for most of the mutants
		#ifdef CAP_FA_REP
			if ( ( score_list[ ii ] == fa_rep ) && ( vars[ ii ] * wts_[ jj ]->free_data()[ ii ] > 10 ) ) { wt_energies[ jj ] += 10; }
			else
		#endif
				wt_energies[ jj ] += vars[ ii ]  * wts_[ jj ]->free_data()[ ii ];

			wt_dG_energies[ jj ][ ii ] = vars[ ii ]  * wts_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			// cap the fa_rep term at some value - this at least keeps it around for most of the mutants
		#ifdef CAP_FA_REP
			if ( ( score_list[ ii ] == fa_rep ) && ( vars[ ii ] * muts_[ jj ]->free_data()[ ii ] > 10 ) ) { mut_energies[ jj ] += 10; }
			else
		#endif
				mut_energies[ jj ] += vars[ ii ] * muts_[ jj ]->free_data()[ ii ];

			mut_dG_energies[ jj ][ ii ] = vars[ ii ] * muts_[ jj ]->free_data()[ ii ];
		}
	}
	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
		#ifdef CAP_FA_REP
			if ( ( fixed_score_list[ ii ] == fa_rep ) && ( fixed_terms[ fixed_score_list[ ii ] ] * wts_[ jj ]->fixed_data()[ ii ] > 10 ) ) { wt_energies[ jj ] += 10; }
			else
		#endif
				wt_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * wts_[ jj ]->fixed_data()[ ii ];

			wt_dG_energies[ jj ][ num_energy_dofs + ii ] = fixed_terms[ fixed_score_list[ ii ] ] * wts_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
		#ifdef CAP_FA_REP
			if ( ( fixed_score_list[ ii ] == fa_rep ) && ( fixed_terms[ fixed_score_list[ ii ] ] * muts_[ jj ]->fixed_data()[ ii ] > 10 ) ) { mut_energies[ jj ] += 10; }
			else
		#endif
				mut_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * muts_[ jj ]->fixed_data()[ ii ];

			mut_dG_energies[ jj ][ num_energy_dofs + ii ] = fixed_terms[ fixed_score_list[ ii ] ] * muts_[ jj ]->fixed_data()[ ii ];
		}
	}

	// I presume these are the reference energies that are being added in?
	// num_energy_dofs is the number of free, non-reference energy parameters in the run -ronj
	if ( num_ref_dofs != 0 ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
			wt_energies[ jj ] += vars[ num_energy_dofs + wt_aa_ ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			mut_energies[ jj ] += vars[ num_energy_dofs + mut_aa_ ];
		}
	}

	//if ( print ) {
	//	for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
	//		ostr << "EXTRA " << tag() << " structure " << jj << " ddG by term: [ ";
	//		for ( Size ii = 1; ii <= num_energy_dofs; ++ii )
	//			ostr << name_from_score_type( score_list[ ii ] ) << ": " << ( mut_dG_energies[ jj ][ ii ] - wt_dG_energies[ jj ][ ii ] ) << ", ";
	//		for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii )
	//			ostr << name_from_score_type( fixed_score_list[ ii ] ) << ": " << ( mut_dG_energies[ jj ][ num_energy_dofs + ii ] - wt_dG_energies[ jj ][ num_energy_dofs + ii ] ) << ", ";
	//		ostr << "]" << std::endl;
	//	}
	//}


	Real predicted_ddG( 0.0 );
	Real ddG_diff( 0.0 );

	// This is where we branch on how the score is calculated.  The simplest approach is to take the minimum energy
	// of all the wts and all the muts and subtract them to get the ddG.  The mean-based approach use the difference
	// of the average of all muts and average of all wts. Finally, the boltzmann approach calculates a boltzmann
	// probability for the wts and muts to get a score.
	if ( option[ optE::optimize_ddGmutation_straight_mean ] ) {

		Real const wt_avg_energy = numeric::statistics::mean( wt_energies.begin(), wt_energies.end(), Real( 0.0 ) );
		Real const mut_avg_energy= numeric::statistics::mean( mut_energies.begin(), mut_energies.end(), Real( 0.0 ) );

		predicted_ddG = mut_avg_energy - wt_avg_energy;
		ddG_diff = predicted_ddG - experimental_ddG_;

		if ( print ) {
			ostr << "DDGMutationOptEData pred: " << predicted_ddG << "  exp: " << experimental_ddG_;
			ostr << " err: " << ddG_diff * ddG_diff << " cmptwt_err: ";
			ostr << component_weights[ ddG_mutation_correlation ] * ddG_diff * ddG_diff;
			ostr << " " << tag() << "\n";
		} else {
			Real const inv_nmuts = 1.0 / muts_.size();
			Real const inv_nwts  = 1.0 / wts_.size();
			for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
				for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
					dE_dvars[ e_dof ] += inv_nmuts * 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff * muts_[ jj ]->free_data()[ e_dof ];
				}
				for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
					dE_dvars[ e_dof ] += -1 * inv_nwts * 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff * wts_[ jj ]->free_data()[ e_dof ];
				}
			}

			if ( num_ref_dofs != 0 ) {
				dE_dvars[ num_energy_dofs + mut_aa_ ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff;
				dE_dvars[ num_energy_dofs + wt_aa_  ] -= 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff;
			}
		}

	} else if ( option[ optE::optimize_ddGmutation_boltzman_average ] ) {

		static Real const inv_kT( option[ optE::inv_kT_nataa ] );  // temp hack

		/// First normalize the energies of both the mutants and the wild type structures.
		Real best_wt = utility::min( wt_energies );
		Real best_mut = utility::min( mut_energies );
		Real best_of_best = std::min( best_wt, best_mut );
		for ( Size ii = 1; ii <= wt_energies.size(); ++ii ) wt_energies[ ii ] -= best_of_best;
		for ( Size ii = 1; ii <= mut_energies.size(); ++ii ) mut_energies[ ii ] -= best_of_best;

		/// Now compute the boltzman probabilities of each of the structures.

		Real wt_denom(0.0);
		Real mut_denom(0.0);
		utility::vector1< Real > wt_partition( wts_.size(), 0.0 );
		utility::vector1< Real > mut_partition( muts_.size(), 0.0 );
		for ( Size ii = 1; ii <= wt_partition.size(); ++ii ) wt_denom += wt_partition[ ii ] = std::exp( -1 * wt_energies[ ii ] * inv_kT );
		for ( Size ii = 1; ii <= mut_partition.size(); ++ii ) mut_denom += mut_partition[ ii ] = std::exp( -1 * mut_energies[ ii ] * inv_kT );

		Real wt_avg_energy( 0.0 ), mut_avg_energy( 0.0 );
		for ( Size ii = 1; ii <= wt_partition.size(); ++ii ) wt_avg_energy += wt_partition[ ii ] * wt_energies[ ii ];
		for ( Size ii = 1; ii <= mut_partition.size(); ++ii ) mut_avg_energy +=  mut_partition[ ii ] * mut_energies[ ii ];
		wt_avg_energy /= wt_denom;
		mut_avg_energy /= mut_denom;

		predicted_ddG = mut_avg_energy - wt_avg_energy;
		ddG_diff = predicted_ddG - experimental_ddG_;

		if ( print ) {
			ostr << "DDGMutationOptEData pred: " << predicted_ddG << "  exp: " << experimental_ddG_;
			ostr << " err: " << ddG_diff * ddG_diff << " cmptwt_err: ";
			ostr << component_weights[ ddG_mutation_correlation ] * ddG_diff * ddG_diff;
			ostr << " " << tag() << "\n";
		} else {

			for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
				/// quotient rule: d/dx f/g == (f'g - g'f)/(g*g).
				Real fmut(0.0), fprimemut(0.0), gmut( mut_denom ), gprimemut(0.0);
				for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
					fmut += mut_partition[ jj ] * mut_energies[ jj ];
					fprimemut += -1 * muts_[ jj ]->free_data()[ e_dof ] * inv_kT * mut_partition[ jj ] + muts_[ jj ]->free_data()[ e_dof ] * mut_partition[ jj ];
					gprimemut += -1 * muts_[ jj ]->free_data()[ e_dof ] * inv_kT * mut_partition[ jj ];
				}
				dE_dvars[ e_dof ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff * ( fprimemut * gmut - gprimemut * fmut ) / (gmut*gmut);

				Real fwt(0.0), fprimewt(0.0), gwt( wt_denom ), gprimewt(0.0);
				for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
					fwt += wt_partition[ jj ] * wt_energies[ jj ];
					fprimewt += -1 * wts_[ jj ]->free_data()[ e_dof ] * inv_kT * wt_partition[ jj ] + wts_[ jj ]->free_data()[ e_dof ] * wt_partition[ jj ];
					gprimewt += -1 * wts_[ jj ]->free_data()[ e_dof ] * inv_kT * wt_partition[ jj ];
				}
				dE_dvars[ e_dof ] -= 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff * ( fprimewt * gwt - gprimewt * fwt ) / (gwt*gwt);
			}

			if ( num_ref_dofs != 0 ) {
				dE_dvars[ num_energy_dofs + mut_aa_ ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff;
				dE_dvars[ num_energy_dofs + wt_aa_  ] -= 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff;
			}
		}


	} else {
		Size const best_wt = arg_min( wt_energies );
		Size const best_mut = arg_min( mut_energies );

		Real const best_wt_energy = wt_energies[ best_wt ];
		Real const best_mut_energy = mut_energies[ best_mut ];

		predicted_ddG = best_mut_energy - best_wt_energy;
		ddG_diff = predicted_ddG - experimental_ddG_;

		if ( print ) {

			ostr << "DDGMutationOptEData pred: " << predicted_ddG << "  exp: " << experimental_ddG_;
			ostr << " err: " << ddG_diff * ddG_diff << " cmptwt_err: ";
			ostr << component_weights[ ddG_mutation_correlation ] * ddG_diff * ddG_diff;
			ostr << " " << tag() << "\n";

		} else {

			for( Size e_dof(1); e_dof <= num_energy_dofs; ++e_dof ) {
				if ( ( score_list[ e_dof ] == fa_rep ) && ( muts_[ best_mut ]->free_data()[ e_dof ] - wts_[ best_wt ]->free_data()[ e_dof ] ) > 10 )
					dE_dvars[ e_dof ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff * 10;
				else
					dE_dvars[ e_dof ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff * ( muts_[ best_mut ]->free_data()[ e_dof ] - wts_[ best_wt ]->free_data()[ e_dof ] );
			}

			if ( num_ref_dofs != 0 ) {
				dE_dvars[ num_energy_dofs + mut_aa_ ] += 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff;
				dE_dvars[ num_energy_dofs + wt_aa_  ] -= 2 * component_weights[ ddG_mutation_correlation ] * ddG_diff;
			}

		}

	}

	return component_weights[ ddG_mutation_correlation ] * ddG_diff * ddG_diff;
}


/*
void
DDGMutationOptEData::print_score(
	std::ostream & ostr,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec &,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const &,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core::optimization;
	using namespace utility;

	if ( muts_.size() == 0 || wts_.size() == 0 ) return;

	utility::vector1< Real > wt_energies( wts_.size(), 0.0 );
	utility::vector1< Real > mut_energies( muts_.size(), 0.0 );
	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
			wt_energies[ jj ] += vars[ ii ]  * wts_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			mut_energies[ jj ] += vars[ ii ] * muts_[ jj ]->free_data()[ ii ];
		}
	}
	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
			wt_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * wts_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			mut_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * muts_[ jj ]->fixed_data()[ ii ];
		}
	}

	if ( num_ref_dofs != 0 ) {
		for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
			wt_energies[ jj ] += vars[ num_energy_dofs + wt_aa_ ];
		}
		for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
			mut_energies[ jj ] += vars[ num_energy_dofs + mut_aa_ ];
		}
	}

	Size const best_wt = arg_min( wt_energies );
	Size const best_mut = arg_min( mut_energies );

	Real const best_wt_energy = wt_energies[ best_wt ];
	Real const best_mut_energy = mut_energies[ best_mut ];

	Real const predicted_ddG = best_mut_energy - best_wt_energy;
	Real const ddG_diff = predicted_ddG - experimental_ddG_;

	TR << "DDGMutationOptEData::print_score(): predicted: " << predicted_ddG << ", experimental: " << experimental_ddG_
		<< ", ddG_diff^2: " << ddG_diff * ddG_diff << ", tag: " << tag() << std::endl;

	ostr << "DDGMutationOptEData pred: " << predicted_ddG << "  exp: " << experimental_ddG_;
	ostr << " err: " << ddG_diff * ddG_diff << " cmptwt_err: ";
	ostr << component_weights[ ddG_mutation_correlation ] * ddG_diff * ddG_diff;
	ostr << " " << tag() << "\n";

	//ostr << "DDGMutRaw: " << best_wt << " " << best_wt_energy << " " << best_mut << " " << best_mut_energy << " : ";
	//for ( Size ii = 1; ii <= wts_[ best_wt ]->free_data().size(); ++ii ) {
	//	ostr << " ( " << wts_[ best_wt ]->free_data()[ ii ] << " , " << muts_[ best_mut ]->free_data()[ ii ] << " ) ";
	//}
	//ostr << "\n";

	return;
}*/

void
DDGMutationOptEData::range(
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list,
	EnergyMap & lower_bound,
	EnergyMap & upper_bound
) const
{
	for ( Size jj = 1; jj <= wts_.size(); ++jj ) {
		update_range( wts_[ jj ],  free_score_list, fixed_score_list, lower_bound, upper_bound );
	}
	for ( Size jj = 1; jj <= muts_.size(); ++jj ) {
		update_range( muts_[ jj ], free_score_list, fixed_score_list, lower_bound, upper_bound );
	}
}

Size
DDGMutationOptEData::size() const
{
	return wts_.size() + muts_.size();
}


OptEPositionDataType
DDGMutationOptEData::type() const
{
	return ddG_mutation_correlation;
}


void
DDGMutationOptEData::write_to_file( std::ofstream & /*outfile*/ ) const
{}


void
DDGMutationOptEData::read_from_file( std::ifstream & /*infile*/ )
{}


void
DDGMutationOptEData::write_to_binary_file( std::ofstream & /*outfile*/ ) const
{}


void
DDGMutationOptEData::read_from_binary_file( std::ifstream & /*infile*/ )
{}


Size
DDGMutationOptEData::memory_use() const
{
	Size total = sizeof( DDGMutationOptEData ) +
		sizeof( SingleStructureData ) * wts_.size() +
		sizeof( SingleStructureData ) * muts_.size();
	if ( wts_.size() > 0 ) {
		total += sizeof( Real ) * ( wts_[ 1 ]->free_data().size() + wts_[ 1 ]->fixed_data().size() ) * wts_.size();
	}
	if ( muts_.size() > 0 ) {
		total += sizeof( Real ) * ( muts_[ 1 ]->free_data().size() + muts_[ 1 ]->fixed_data().size() ) * muts_.size();
	}
	return total;
}


#ifdef USEMPI

void
DDGMutationOptEData::send_to_node( int const destination_node, int const tag ) const
{
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

	if ( nwts == 0 || nmuts == 0 ) return;

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

	OptEPositionData::send_to_node( destination_node, tag );

}

void
DDGMutationOptEData::receive_from_node( int const source_node, int const tag )
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
		wts_.push_back( new SingleStructureData( free_data_v, fixed_data_v ) );
	}


	delete [] free_data; free_data = 0;
	delete [] fixed_data; fixed_data = 0;

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
		muts_.push_back( new SingleStructureData( free_data_v, fixed_data_v ) );
	}


	delete [] free_data;
	delete [] fixed_data;

	OptEPositionData::receive_from_node( source_node, tag );

}
#endif


void
DDGMutationOptEData::set_wt_aa( AA wt_aa )
{
	wt_aa_ = wt_aa;
}

void
DDGMutationOptEData::set_mut_aa( AA mut_aa )
{
	mut_aa_ = mut_aa;
}

void
DDGMutationOptEData::set_experimental_ddg( Real ddg )
{
	experimental_ddG_ = ddg;
}


void
DDGMutationOptEData::add_wt( SingleStructureDataOP wt )
{
	wts_.push_back( wt );
}

void
DDGMutationOptEData::add_mutant( SingleStructureDataOP mut )
{
	muts_.push_back( mut );
}


//
// ------------------- OptEData-----------------------//
//

core::Size
OptEData::num_rotamers() const
{
	core::Size num_rots(0);
	for ( OptEPositionDataOPs::const_iterator pos( data_.begin() );
	      pos != data_.end(); ++pos ) {
		num_rots += (*pos)->size();
	}
	return num_rots;
}


///@begin OptEData::read_from_file
///@brief slow
///@author ashworth
void
OptEData::read_from_file( std::string filename )
{
	runtime_assert( data_.empty() );
	data_.clear();
	ScoreTypes fixed_terms, free_terms;
	std::ifstream infile( filename.c_str() );
	std::string line;

	typedef utility::vector1< std::string > Strings;
	using utility::string_split;

	// first read the header to get the get list of fixed and free terms
	while( getline( infile, line ) ) {
		Strings words( string_split( line, ' ' ) );
		if ( line.substr(0,12) == "fixed terms:" ) {
			for ( Strings::iterator term( words.begin()+2 ); term != words.end(); ++term ) {
				fixed_terms.push_back( score_type_from_name( *term ) );
			}
		} else if ( line.substr(0,11) == "free terms:" ) {
			for ( Strings::iterator term( words.begin()+2 ); term != words.end(); ++term ) {
				free_terms.push_back( score_type_from_name( *term ) );
			}
			break; // (free terms expected to come after fixed terms)
		}
	}
	// assign fixed and free terms
	runtime_assert( !free_terms.empty() );
	runtime_assert( ( fixed_terms == fixed_energy_terms_ ) || fixed_energy_terms_.empty() );
	runtime_assert( ( free_terms == energy_terms_ ) || energy_terms_.empty() );
	fixed_energy_terms_ = fixed_terms;
	energy_terms_ = free_terms;

	while ( getline( infile, line ) ) {
		Strings words( string_split( line, ' ' ) );

		if ( words.size() == 0 ) { utility_exit_with_message("Bad line from optE data file"); }
		if ( words[ 1 ] != "optE_position_data_type" ) { utility_exit_with_message("Bad line from optE data file"); }

		OptEPositionDataType postype = (OptEPositionDataType) utility::from_string( words[2], int( 0 ) /* dummy */ );
		OptEPositionDataOP new_position_data = OptEPositionDataFactory::create_position_data( postype );
		new_position_data->read_from_file( infile );
		add_position_data( new_position_data );
	}

	infile.close();
	TR << "Read in " << num_rotamers() << " rotamers at "
	   << num_positions() << " positions from file " << filename << std::endl;
}

///@begin OptEData::write_to_file
///@brief human-readable
///@author ashworth
void
OptEData::write_to_file( std::string filename ) const
{
	std::ofstream outfile( filename.c_str() );
	// write explanatory header
	outfile << "# data line format:" << std::endl << "# rot_number, this_aa,"
	        << " [ energies for fixed terms ], [ energies for free terms ]"
	        << "\n";

	outfile << "fixed terms:";
	for ( ScoreTypes::const_iterator et( fixed_energy_terms_.begin() );
	      et != fixed_energy_terms_.end(); ++et ) {
		outfile << " " << name_from_score_type( *et );
	}
	outfile << "\n";

	outfile << "free terms:";
	for ( ScoreTypes::const_iterator et( energy_terms_.begin() );
	      et != energy_terms_.end(); ++et ) {
		outfile << " " << name_from_score_type( *et );
	}
	outfile << "\n";

	// write data
	for ( OptEPositionDataOPs::const_iterator pos( position_data_begin() );
	      pos != position_data_end(); ++pos ) {
		outfile << "optE_position_data_type " << (*pos)->type() << "\n";
		(*pos)->write_to_file( outfile );
	}
	outfile.close();
	TR << "Wrote " << num_rotamers() << " rotamers at "
	   << num_positions() << " positions to file " << filename << std::endl;
}


///@begin OptEData::write_to_binary_file
///@brief writes out the optE data to a binary file
///@author ashworth
void
OptEData::write_to_binary_file( std::string filename ) const
{
	std::ofstream outfile( filename.c_str(), std::ios::out | std::ios::binary );
	Size const npositions( data_.size() );
	outfile.write( (char*) &npositions, sizeof(Size) );
	for ( OptEPositionDataOPs::const_iterator pos( data_.begin() ); pos != data_.end(); ++pos ) {
		OptEPositionDataType pos_data_type = (*pos)->type();
		outfile.write( (char*) & pos_data_type, sizeof(OptEPositionDataType) );
		(*pos)->write_to_binary_file( outfile );
	}
	outfile.close();
	TR << "Wrote " << num_rotamers() << " rotamers at "
	   << num_positions() << " positions to binary file " << filename << std::endl;

}

///@begin OptEData::read_from_binary_file
///@brief binary I/O should be faster
///@author ashworth
void
OptEData::read_from_binary_file( std::string filename )
{

	std::ifstream infile( filename.c_str(), std::ios::in | std::ios::binary );
	Size npositions(0);
	// this while loop allows concatenation of multiple optE output files
	while ( infile.read( (char*) &npositions, sizeof(Size) ) ) {
		for ( Size i(1); i <= npositions; ++i ) {
			OptEPositionDataType pos_data_type;
			infile.read( (char*) &pos_data_type, sizeof( OptEPositionDataType ) );
			OptEPositionDataOP pos_data = OptEPositionDataFactory::create_position_data( pos_data_type );
			pos_data->read_from_binary_file( infile );
			add_position_data( pos_data );
		}
		npositions = 0;
	}
	infile.close();
	TR << "Read in " << num_rotamers() << " rotamers at "
	   << num_positions() << " positions from binary file " << filename << std::endl;

}


/// Initializers for static data
bool OptEPositionDataFactory::optE_type_name_map_initialized_( false );

utility::vector1< std::string >
OptEPositionDataFactory::optE_type_2_optE_type_name_;

std::map< std::string, OptEPositionDataType >
OptEPositionDataFactory::optE_type_name_map_;

OptEPositionDataOP
OptEPositionDataFactory::create_position_data( OptEPositionDataType const type )
{
	switch ( type ) {
		case prob_native_amino_acid :
			return new PNatAAOptEPositionData;
		case prob_native_amino_acid_with_unfolded_energy :
			return new NestedEnergyTermPNatAAOptEPositionData;
		case pssm_data :
			return new PSSMOptEPositionData;
		case prob_native_rotamer :
			return new PNatRotOptEPositionData;
		case prob_native_structure :
			return new PNatStructureOptEData;
		case prob_native_ligand_pose :
			return new PNatLigPoseOptEData;
		case dG_binding_correlation :
			return new DGBindOptEData;
		case ddG_mutation_correlation :
			return new DDGMutationOptEData;
		case ddG_mutation_correlation_with_unfolded_energy :
			return new NestedEnergyTermDDGMutationOptEData;
		case ddG_bind_correlation :
			return new DDGBindOptEData;
		case constrained_optimization_weight_func:
			return new ConstraintedOptimizationWeightFunc;
	}
	return 0;
}



std::string const &
OptEPositionDataFactory::optE_type_name( OptEPositionDataType const type )
{
	initialize_optE_type_name_map();
	runtime_assert( type > 0 && type <= n_optE_data_types );
	return optE_type_2_optE_type_name_[ type ];
}

bool
OptEPositionDataFactory::is_optE_type_name( std::string const & name )
{
	initialize_optE_type_name_map();
	return optE_type_name_map_.find( name ) != optE_type_name_map_.end();
}

OptEPositionDataType
OptEPositionDataFactory::optE_type_from_name( std::string const & name )
{
	initialize_optE_type_name_map();
	runtime_assert( is_optE_type_name( name ));
	return optE_type_name_map_[ name ];
}


void
OptEPositionDataFactory::initialize_optE_type_name_map()
{
	if ( optE_type_name_map_initialized_ ) return;

	optE_type_name_map_initialized_ = true;

	optE_type_2_optE_type_name_.resize( n_optE_data_types );
	optE_type_2_optE_type_name_[ prob_native_amino_acid ] = "prob_native_amino_acid";
	optE_type_2_optE_type_name_[ prob_native_amino_acid_with_unfolded_energy ] = "prob_native_amino_acid/unfolded_energy_scorefxn";
	optE_type_2_optE_type_name_[ pssm_data ] = "pssm_data";
	optE_type_2_optE_type_name_[ prob_native_rotamer ] = "prob_native_rotamer";
	optE_type_2_optE_type_name_[ prob_native_structure ] = "prob_native_structure";
	optE_type_2_optE_type_name_[ prob_native_ligand_pose ] = "prob_native_ligand_pose";
	optE_type_2_optE_type_name_[ dG_binding_correlation ] = "dG_binding_correlation";
	optE_type_2_optE_type_name_[ ddG_mutation_correlation ] = "ddG_mutation_correlation";
	optE_type_2_optE_type_name_[ ddG_mutation_correlation_with_unfolded_energy ] = "ddG_mutation_correlation/unfolded_energy_scorefxn";
	optE_type_2_optE_type_name_[ ddG_bind_correlation ] = "ddG_bind_correlation";
	optE_type_2_optE_type_name_[ constrained_optimization_weight_func ] = "constrained_optimization_weight_func";

	for ( Size ii = 1; ii <= n_optE_data_types; ++ii ) {
		optE_type_name_map_[ optE_type_2_optE_type_name_[ ii ] ] = static_cast< OptEPositionDataType > ( ii );
	}
}



}
}
