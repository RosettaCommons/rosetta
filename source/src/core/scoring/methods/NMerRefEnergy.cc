// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
//// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/NMerRefEnergy.hh
/// @brief  Reference energy method implementation
/// @author Chris King (dr.chris.king@gmail.com)

// Unit headers
#include <core/scoring/methods/NMerRefEnergy.hh>
#include <core/scoring/methods/NMerRefEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
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
	return new NMerRefEnergy;
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

void
NMerRefEnergy::initialize_from_options()
{ 
  using namespace basic::options;
  NMerRefEnergy::nmer_length( option[ OptionKeys::score::nmer_ref_seq_length ]() );
}

NMerRefEnergy::NMerRefEnergy() :
	parent( new NMerRefEnergyCreator )
{
	NMerRefEnergy::initialize_from_options();
	read_nmer_tables_from_options();
}

NMerRefEnergy::NMerRefEnergy( std::map< std::string, core::Real > const & nmer_ref_energies_in ):
	parent( new NMerRefEnergyCreator )
{
	NMerRefEnergy::initialize_from_options();

	nmer_ref_energies_.clear();
	for ( std::map< std::string, Real >::const_iterator it = nmer_ref_energies_in.begin(); it != nmer_ref_energies_in.end(); ++it ) {
			nmer_ref_energies_.insert( *it );
	}
}

NMerRefEnergy::~NMerRefEnergy() {}

void NMerRefEnergy::read_nmer_tables_from_options() {

  using namespace basic::options;

  TR << "checking for NMerRefEnergy Ref list" << std::endl;

  //check for ref list file
  if ( option[ OptionKeys::score::nmer_ref_energies_list ].user() ) {
    std::string const ref_list_fname( option[ OptionKeys::score::nmer_ref_energies_list ] );
    NMerRefEnergy::read_nmer_table_list( ref_list_fname );
  }
  //use single ref file
  if( option[ OptionKeys::score::nmer_ref_energies ].user() ){
    std::string const ref_fname( option[ OptionKeys::score::nmer_ref_energies ] );
    NMerRefEnergy::read_nmer_table( ref_fname );
  }
}

//read energy table list
//entries from all lists just get added to the same map
void NMerRefEnergy::read_nmer_table_list( std::string ref_list_fname ) {
  TR << "reading NMerRefEnergy list from " << ref_list_fname << std::endl;
  if ( !utility::file::file_exists( ref_list_fname ) ) {
		ref_list_fname = basic::database::full_name( ref_list_fname, false );
	}
  utility::io::izstream in_stream( ref_list_fname );
  if (!in_stream.good()) {
    utility_exit_with_message( "[ERROR] opening NMerRefEnergy list file at " + ref_list_fname );
  }
  //now loop over all names in list
  std::string ref_fname;
  while( getline( in_stream, ref_fname ) ){
    utility::vector1< std::string > const tokens( utility::split( ref_fname ) );
    //skip comments
    if( tokens[ 1 ][ 0 ] == '#' ) continue;
    NMerRefEnergy::read_nmer_table( ref_fname );
  }
}

void NMerRefEnergy::read_nmer_table( std::string ref_fname ) {

	TR << "checking for NMerRefEnergy scores" << std::endl;

  if ( !utility::file::file_exists( ref_fname ) ) {
		ref_fname = basic::database::full_name( ref_fname, false );
	}
	TR << "reading NMerRefEnergy scores from " << ref_fname << std::endl;
	utility::io::izstream in_stream( ref_fname );
	if (!in_stream.good()) {
		utility_exit_with_message( "[ERROR] Error opening NMerRefEnergy file at " + ref_fname );
	}
	std::string line;
	while( getline( in_stream, line) ) {
		utility::vector1< std::string > const tokens ( utility::split( line ) );
		//skip comments
		if( tokens[ 1 ][ 0 ] == '#' ) continue;
		if( tokens.size() != 2 ) utility_exit_with_message( "[ERROR] NMer ref energy database file "
				+ ref_fname + " does not have 2 entries at line " + line );
		std::string const sequence( tokens[ 1 ] );
		if( sequence.size() != nmer_length_ ) utility_exit_with_message( "[ERROR] NMer ref energy database file "
				+ ref_fname + " has wrong length nmer at line " + line
				+ "\n\texpected: " + utility::to_string( nmer_length_ ) + " found: " + utility::to_string( sequence.size() ) );
		//everything is cool! nothing is fucked!
		Real const energy( atof( tokens[ 2 ].c_str() ) );
		//Hmmmm... if we have duplicate entries, say in multiple tables, what should we do? error, replace, or sum?
		//I think we should sum them; that is, all lists are correct, even if they say diff things about diff nmers	
//		if( nmer_ref_energies_.count( sequence ) ) utility_exit_with_message( "[ERROR] NMer ref energy database file "
//				+ ref_fname + " has double entry for sequence " + sequence );
		if( nmer_ref_energies_.count( sequence ) ){
			TR << "[WARNING]: NMer ref energy database file "
				+ ref_fname + " has double entry for sequence " + sequence + " Summing with prev value..." << std::endl ;
			nmer_ref_energies_[ sequence ] += energy;
		}
		else nmer_ref_energies_[ sequence ] = energy;
	}

}


EnergyMethodOP
NMerRefEnergy::clone() const
{
	return new NMerRefEnergy( nmer_ref_energies_ );
}


//retrieves ref energy of NMer centered on seqpos
void
NMerRefEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	using namespace chemical;

	if( nmer_ref_energies_.empty() ) return;
	Size const seqpos( rsd.seqpos() );
	//skip if not a NMer center 
	if( seqpos < 1 || seqpos > pose.total_residue() - nmer_cterm_ ) return;
	//get the NMer centered on seqpos
	std::string sequence;
	for( Size iseq = seqpos; iseq <= seqpos + nmer_cterm_; ++iseq ){
		sequence += pose.residue( iseq ).name1();
	}
	//bail if seq not in table
	if( !nmer_ref_energies_.count( sequence ) ) return;
	//must use find because this is a const function!
	//the [] operator will add the element to the map if key not found, which changes the state of this function
	Real const energy( nmer_ref_energies_.find( sequence )->second );
	emap[ nmer_ref ] += energy;

	return;

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

