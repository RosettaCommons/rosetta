// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/AtomVDW.hh
/// @brief
/// @author Phil Bradley


// Unit Headers
#include <core/scoring/AtomVDW.hh>

// Package headers


// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/io/izstream.hh>


#include <cmath>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace core {
namespace scoring {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.AtomVDW" );


/// @details ctor, reads data file. Need to configure to allow alternate tables/atom_sets
AtomVDW::AtomVDW( std::string const & atom_type_set_name_with_suffix ):vdw_suffix_("NULL") //atom_type_set_name_( atom_type_set_name )
{
	//init atom_type_set_name_ and vdw_suffix_
	get_atom_type_set_name(atom_type_set_name_with_suffix);

	std::string atom_type_set_name = atom_type_set_name_;

	// get the relevant atomset
	chemical::AtomTypeSet const & atom_set
		( *chemical::ChemicalManager::get_instance()->atom_type_set( atom_type_set_name ) );//chemical::CENTROID ) );

	// initialize the values
	Size const n_atom_types( atom_set.n_atomtypes() );
	atom_vdw_.clear();
	atom_vdw_.resize( n_atom_types );
	for ( Size j=1; j<= n_atom_types; ++j ) {
		atom_vdw_[j].resize( n_atom_types, 0.0 ); // default value is 0.0 --> no clashes
	}

	utility::io::izstream stream;
	if ( vdw_suffix_.compare("NULL")==0 ) {
		stream.open( atom_set.directory() + "/atom_vdw.txt" );
	} else {
		stream.open( atom_set.directory() + "/" + vdw_suffix_ + ".txt" );
		TR << "Openning alternative vdw file: " << atom_set.directory() + "/" + vdw_suffix_ + ".txt" << std::endl;
	}

	if ( !stream.good() ) {
		stream.clear();
		if ( atom_type_set_name != chemical::CENTROID ) {
			utility_exit_with_message( "cant find atom_vdw.txt in directory: "+ atom_set.directory() );
		}
		TR.Warning << "YOU NEED TO UPDATE YOUR DATABASE!" << std::endl;
		basic::database::open( stream, "scoring/AtomVDW/atom_vdw.txt" );
		if ( !stream.good() ) {
			utility_exit_with_message( "Unable to open scoring/AtomVDW/atom_vdw.txt!" );
		}
	}

	// read the entire file and figure out what atom_types are present and in what order
	utility::vector1< std::string > lines;
	utility::vector1< int > atom_type_index; // mapping from index in the file to index in atom_set

	std::string line, atom_type_name;
	while ( getline( stream, line ) ) {
		lines.push_back(line);
		std::istringstream l(line);
		l >> atom_type_name;
		atom_type_index.push_back( atom_set.atom_type_index( atom_type_name ) ); // will fail if name is not recognized
	}

	Size const natoms( lines.size() );
	for ( Size i=1; i<= natoms; ++i ) {
		std::istringstream l( lines[i] );
		l >> atom_type_name;

		Size const i_atom_type_index( atom_type_index[i] );

		for ( Size j=1; j<= natoms; ++j ) {
			Size const j_atom_type_index( atom_type_index[j] );
			l >> atom_vdw_[ i_atom_type_index ][ j_atom_type_index ];
		}
		if ( l.fail() ) utility_exit_with_message( "bad format in atom_vdw.txt ");
	}

	setup_approximate_vdw_radii( atom_type_index, atom_set );
}

// tool for parseing alternative vdw parameter file
void AtomVDW::get_atom_type_set_name( std::string ats_suff )
{
	std::vector<std::string> elems;
	std::stringstream ss(ats_suff);
	std::string item;
	while ( std::getline(ss, item, ':') ) {
		elems.push_back(item);
	}
	Size n = elems.size();
	if ( n==0 ) utility_exit_with_message( "Unable to parse the atom_type_set name" );
	atom_type_set_name_ = elems[0];
	if ( n==2 ) vdw_suffix_ = elems[1];
	if ( n>2 ) utility_exit_with_message( "Error format of atom_type_set name" );
}

/// @details  Calculates approximation to a single per-atom radius using the pairwise data from the file
/// For atoms not present in the file uses the lj_radius from atom_type_set

void
AtomVDW::setup_approximate_vdw_radii(
	utility::vector1< int > const & atom_type_index,
	chemical::AtomTypeSet const & atom_type_set
) {
	using namespace ObjexxFCL::format;

	Size const natoms_full( atom_vdw_.size() );
	Size const natoms( atom_type_index.size() ); // # atoms in the db file
	debug_assert( natoms_full == Size( atom_type_set.n_atomtypes() ) );

	// initialize
	utility::vector1< Real > R( natoms_full );
	for ( Size i=1; i<= natoms_full; ++i ) {
		R[i] = atom_type_set[i].lj_radius();
	}

	// just for the guys in the file
	for ( Size i=1; i<= natoms; ++i ) {
		Size const ii( atom_type_index[i] );
		R[ii] = std::sqrt( atom_vdw_[ii][ii] );
	}

	Size const niter( 10 );
	for ( Size r=1; r<= niter; ++r ) {

		utility::vector1< Real > new_R( R );
		// calculate deviations
		Real total_dev( 0.0 );

		for ( Size i=1; i<= natoms; ++i ) {
			Size const ii( atom_type_index[i] );

			Real dev(0.0), abs_dev(0.0);
			for ( Size j=1; j<= natoms; ++j ) {
				Size const jj( atom_type_index[j] );

				Real const actual( std::sqrt( atom_vdw_[ii][jj] ) );
				Real const expected( R[ii] + R[jj] );
				dev += actual - expected;
				total_dev += ( actual - expected ) * ( actual - expected );
				abs_dev += ( actual - expected ) * ( actual - expected );
			}
			dev /= natoms;
			abs_dev = std::sqrt( abs_dev / natoms );
			new_R[ii] = R[ii] + dev/2;
			//std::cout << "R " << I(4,i) << F(9,3,R[ii]) <<F(9,3,dev) << F(9,3,abs_dev ) << std::endl;
		}

		R = new_R;
		//std::cout << "round: " << r << ' ' << std::sqrt( total_dev / ( natoms * natoms ) ) << std::endl;
	}

	/// now set the long-name guy
	approximate_vdw_radii_ = R;
}

std::string
AtomVDW::atom_type_set_name() const
{
	return atom_type_set_name_;
}

} // namespace scoring
} // namespace core

