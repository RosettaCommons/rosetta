// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/NMerRefEnergy.hh
/// @brief  Reference energy method implementation
/// @author Chris King (dr.chris.king@gmail.com)

// Unit headers
#include <core/scoring/methods/NMerRefEnergy.hh>
#include <core/scoring/methods/NMerRefEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>

// C++ Headers
#include <string>
#include <vector>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>

static basic::Tracer TR( "core.scoring.methods.NMerRefEnergy" );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the NMerRefEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
NMerRefEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new NMerRefEnergy );
}

ScoreTypes
NMerRefEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( nmer_ref );
	return sts;
}

void
NMerRefEnergy::nmer_length( Size const nmer_length ){
	nmer_length_ = nmer_length;
	//nmer residue energy is attributed to position 1
	nmer_cterm_ = nmer_length_ - 1 ;
}

Size
NMerRefEnergy::n_tables() const
{
	return nmer_ref_energies_.size();
}

void
NMerRefEnergy::initialize_from_options()
{
	using namespace basic::options;
	NMerRefEnergy::nmer_length( option[ OptionKeys::score::nmer_ref_seq_length ]() );
}

NMerRefEnergy::NMerRefEnergy() :
	parent( methods::EnergyMethodCreatorOP( new NMerRefEnergyCreator ) )
{
	NMerRefEnergy::initialize_from_options();
	read_nmer_tables_from_options();
}

NMerRefEnergy::NMerRefEnergy( utility::vector1< std::map< std::string, core::Real > > const & nmer_ref_energies_in ):
	parent( methods::EnergyMethodCreatorOP( new NMerRefEnergyCreator ) )
{
	NMerRefEnergy::initialize_from_options();

	nmer_ref_energies_.clear();
	for ( Size i = 1; i <= nmer_ref_energies_in.size(); ++i ) {
		nmer_ref_energies_.push_back( nmer_ref_energies_in[ i ]  );
	}
}

NMerRefEnergy::NMerRefEnergy(
	core::Size const nmer_length,
	utility::vector1< std::string > const & fname_vec
) :
	parent( methods::EnergyMethodCreatorOP( new NMerRefEnergyCreator ) )
{
	NMerRefEnergy::nmer_length( nmer_length  );
	NMerRefEnergy::read_nmer_fname_vector( fname_vec );
}

NMerRefEnergy::~NMerRefEnergy() = default;

void NMerRefEnergy::read_nmer_tables_from_options() {

	using namespace basic::options;

	TR << "checking for NMerRefEnergy Ref list" << std::endl;

	//check for ref list file
	if ( option[ OptionKeys::score::nmer_ref_energies_list ].user() ) {
		std::string const ref_list_fname( option[ OptionKeys::score::nmer_ref_energies_list ] );
		NMerRefEnergy::read_nmer_table_list( ref_list_fname );
	}
	//use single ref file
	if ( option[ OptionKeys::score::nmer_ref_energies ].user() ) {
		std::string const ref_fname( option[ OptionKeys::score::nmer_ref_energies ] );
		NMerRefEnergy::read_nmer_table( ref_fname );
	}
}

//read energy table list
//entries from all lists just get added to the same map
void NMerRefEnergy::read_nmer_table_list( std::string const & ref_list_fname ) {
	TR << "reading NMerRefEnergy list from " << ref_list_fname << std::endl;
	utility::io::izstream in_stream;
	if ( utility::file::file_exists( ref_list_fname ) ) {
		in_stream.open( ref_list_fname );
	} else {
		in_stream.open( basic::database::full_name( ref_list_fname, false ) );
	}
	if ( !in_stream.good() ) {
		utility_exit_with_message( "[ERROR] opening NMerRefEnergy list file at " + ref_list_fname );
	}
	//now loop over all names in list
	std::string ref_fname;
	while ( getline( in_stream, ref_fname ) ) {
		utility::vector1< std::string > const tokens( utility::split( ref_fname ) );
		//skip comments
		if ( tokens[ 1 ][ 0 ] == '#' ) continue;
		NMerRefEnergy::read_nmer_table( ref_fname );
	}
}

//load tables from a vector of filenames
void
NMerRefEnergy::read_nmer_fname_vector( utility::vector1< std::string > const & fname_vec ) {
	//now loop over all names in vector
	for ( Size i = 1; i<= fname_vec.size(); ++i ) {
		std::string const & fname( fname_vec[ i ] );
		NMerRefEnergy::read_nmer_table( fname );
	}
}

// this now appends a new table (map) to the vector of tables (maps) every time new file is read
void NMerRefEnergy::read_nmer_table( std::string const & ref_fname ) {

	TR << "checking for NMerRefEnergy scores" << std::endl;

	utility::io::izstream in_stream;
	if ( utility::file::file_exists( ref_fname ) ) {
		in_stream.open( ref_fname );
	} else {
		in_stream.open( basic::database::full_name( ref_fname, false ) );
	}
	TR << "reading NMerRefEnergy scores from " << ref_fname << std::endl;

	if ( !in_stream.good() ) {
		utility_exit_with_message( "[ERROR] Error opening NMerRefEnergy file at " + ref_fname );
	}

	// single new empty map object to store nmer energies
	std::map< std::string, core::Real > nmer_ref_energy;

	std::string line;
	while ( getline( in_stream, line) ) {
		utility::vector1< std::string > const tokens ( utility::split( line ) );
		//skip comments
		if ( tokens[ 1 ][ 0 ] == '#' ) continue;
		if ( tokens.size() != 2 ) {
			utility_exit_with_message( "[ERROR] NMer ref energy database file "
				+ ref_fname + " does not have 2 entries at line " + line );
		}
		std::string const sequence( tokens[ 1 ] );
		if ( sequence.size() != nmer_length_ ) {
			utility_exit_with_message( "[ERROR] NMer ref energy database file "
				+ ref_fname + " has wrong length nmer at line " + line
				+ "\n\texpected: " + utility::to_string( nmer_length_ ) + " found: " + utility::to_string( sequence.size() ) );
		}
		//everything is cool! nothing is fucked!
		Real const energy( atof( tokens[ 2 ].c_str() ) );
		//Hmmmm... if we have duplicate entries in one table, what should we do? error, replace, or sum?
		//I think we should sum them; that is, all entries are correct, even if they say diff things about same sequence may be true statements that are combined in objective function
		if ( nmer_ref_energy.count( sequence ) ) {
			TR.Warning << "NMer ref energy database file "
				+ ref_fname + " has double entry for sequence " + sequence + " Summing with prev value..." << std::endl ;
			nmer_ref_energy[ sequence ] += energy;
		} else nmer_ref_energy[ sequence ] = energy;
	}
	// append map to vector of maps
	nmer_ref_energies_.push_back( nmer_ref_energy );
}


EnergyMethodOP
NMerRefEnergy::clone() const
{
	return EnergyMethodOP( new NMerRefEnergy( nmer_ref_energies_ ) );
}

void
NMerRefEnergy::get_residue_energy_by_table(
	pose::Pose const & pose,
	Size const & seqpos,
	Real & rsd_energy_sum,
	utility::vector1< Real > & rsd_table_energies
) const
{
	debug_assert( rsd_table_energies.size() == n_tables() );
	//for now, just assign all of the p1=seqpos frame's nmer_table energy to this residue
	//TODO: distribute frame's nmer_table energy evenly across the nmer
	//TODO: avoid wasting calc time by storing nmer_value, nmer_val_out_of_date in pose cacheable data
	Size p1_seqpos( seqpos );

	//get the nmer string
	//TODO: how deal w/ sequences shorter than nmer_length_?
	// this matters at both termini… maybe take max of all overlapping frames w/ missing res as 'X'?
	// go ahead and bail if we fall off the end of the chain
	if ( p1_seqpos + nmer_length_ - 1 <= pose.conformation().chain_end( pose.chain( p1_seqpos ) ) ) {
		//need p1 position in this chain's sequence, offset with index of first res
		std::string chain_sequence( pose.chain_sequence( pose.chain( p1_seqpos ) ) );

		rsd_energy_sum = 0.;
		for ( Size i = 1; i <= n_tables(); ++i ) {

			std::map < std::string, core::Real > nmer_ref_energy( nmer_ref_energies_[ i ] );
			if ( nmer_ref_energy.empty() ) return;
			//skip if not a NMer center
			if ( seqpos < 1 || seqpos > pose.size() - nmer_cterm_ ) return;
			//get the NMer centered on seqpos
			std::string sequence;
			for ( Size iseq = seqpos; iseq <= seqpos + nmer_cterm_; ++iseq ) {
				sequence += pose.residue( iseq ).name1();
			}
			//bail if seq not in table
			if ( !nmer_ref_energy.count( sequence ) ) continue;
			//must use find because this is a const function!
			//the [] operator will add the element to the map if key not found, which changes the state of this function
			Real const rsd_table_energy( nmer_ref_energy.find( sequence )->second );

			//store this tables rsd energy
			rsd_table_energies[ i ] = rsd_table_energy;
			//total score is just sum of all tables
			rsd_energy_sum += rsd_table_energy;
		}
	}
}

//retrieves ref energy of NMer centered on seqpos
void
NMerRefEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	if ( nmer_ref_energies_.empty() ) return;
	Size const seqpos( rsd.seqpos() );

	Real rsd_energy( 0. );
	utility::vector1< Real > rsd_table_energies( n_tables(), Real( 0. ) );
	get_residue_energy_by_table( pose, seqpos, rsd_energy, rsd_table_energies );
	emap[ nmer_ref ] += rsd_energy;

}


Real
NMerRefEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}

/// @brief NMerRefEnergy is context independent; indicates that no
/// context graphs are required
void
NMerRefEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
NMerRefEnergy::version() const
{
	return 1; // Initial versioning
}
} // methods
} // scoring
} // core
