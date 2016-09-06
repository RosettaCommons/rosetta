// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Oliver Lange


// Unit Headers
#include <protocols/noesy_assign/MethylNames.hh>

// Package Headers
#include <protocols/noesy_assign/Exceptions.hh>

// Project Headers
#include <core/chemical/AA.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

static THREAD_LOCAL basic::Tracer tr( "protocols.noesy_assign.crosspeaks" );

// Singleton instance and mutex static data members
namespace utility {

using protocols::noesy_assign::MethylNameLibrary;

#if defined MULTI_THREADED
template <> std::mutex utility::SingletonBase< MethylNameLibrary >::singleton_mutex_{};
template <> std::atomic< MethylNameLibrary * > utility::SingletonBase< MethylNameLibrary >::instance_( 0 );
#else
template <> MethylNameLibrary * utility::SingletonBase< MethylNameLibrary >::instance_( 0 );
#endif

}

namespace protocols {
namespace noesy_assign {

using namespace core;

MethylNameLibrary *
MethylNameLibrary::create_singleton_instance()
{
	return new MethylNameLibrary;
}

MethylNameLibrary::MethylNameLibrary() {
	load_database_table();
}

MethylNames const& MethylNameLibrary::operator[]( chemical::AA aa ) const {
	MethylNameTable::const_iterator it = methyl_names_.find( aa );
	if ( it == methyl_names_.end() ) {
		utility_exit_with_message( name_from_aa( aa )+ " is not a known aminoacid-type in NMR proton Database" );
	}
	return it->second;
}

void MethylNameLibrary::load_database_table() {
	utility::io::izstream db_stream;
	basic::database::open(db_stream, "chemical/residue_type_sets/fa_standard/nmr_nameing_conventions.txt" );
	std::string line;
	std::string last_aa( "X" );
	std::string aa_name;
	MethylNameTable::iterator current;
	while ( getline( db_stream, line ) ) {
		if ( line[0] == '#' ) continue;
		std::istringstream line_stream( line );
		line_stream >> aa_name;
		if ( !line_stream.good() && aa_name.size() != 1 ) continue;
		if ( last_aa != aa_name ) {
			last_aa = aa_name;
			chemical::AA aa( chemical::aa_from_oneletter_code( aa_name[0] ) );
			methyl_names_[ aa ] = MethylNames( aa );
			current = methyl_names_.find( aa );
		}
		std::string nmr, rosetta;
		line_stream >> nmr >> rosetta;
		if ( nmr.size() && rosetta.size() ) {
			current->second.add_proton( nmr, rosetta);
		}

		std::string tag, methyl;
		line_stream >> tag >> methyl;
		tr.Trace << "read: " << tag << " " << methyl << std::endl;
		if ( !line_stream.bad() && tag == ">" ) {
			current->second.add_methyl( methyl, rosetta );
		}
		if ( line_stream.good() ) {
			line_stream >> methyl;
			if ( !line_stream.bad() ) {
				current->second.add_methyl( methyl, rosetta );
			}
		}
	}
}

MethylNames::MethylNames() : aa_( chemical::aa_unk ) {}

MethylNames::MethylNames( chemical::AA aa ) : aa_( aa ) {}

void MethylNames::add_proton( std::string const& nmr, std::string const& rosetta ) {
	tr.Trace << aa_ << " add proton " << nmr << " " << rosetta << std::endl;
	rosetta2nmr_[ rosetta ] = nmr;
	nmr2rosetta_[ nmr ] = AtomList( 1, rosetta );
}

void MethylNames::add_methyl( std::string const& methyl, std::string const& rosetta ) {
	tr.Trace << aa_  << " add methyl " << methyl << " " << rosetta << std::endl;

	NameTable::iterator rit = rosetta2methyl_.find( rosetta );
	if ( rit != rosetta2methyl_.end() ) {
		rit->second.push_back( methyl );
	} else {
		rosetta2methyl_[ rosetta ] = AtomList( 1, methyl );
	}

	NameTable::iterator it =  nmr2rosetta_.find( methyl );
	if ( it != nmr2rosetta_.end() ) {
		it->second.push_back( rosetta );
	} else {
		nmr2rosetta_[ methyl ] = AtomList( 1, rosetta );
	}
}

std::string MethylNames::aa_name() const {
	return name_from_aa( aa_ );
}

std::string const& MethylNames::rosetta2nmr( std::string const& proton ) const {
	std::map< std::string, std::string >::const_iterator it = rosetta2nmr_.find( proton );
	if ( it == rosetta2nmr_.end() ) {
		throw EXCN_UnknownAtomname("proton_name " + proton + " not recognized for aminoacid " + aa_name() );
	}
	return it->second;
}

MethylNames::AtomList const& MethylNames::nmr2rosetta( std::string const& proton ) const {
	NameTable::const_iterator it =  nmr2rosetta_.find( proton );
	if ( it == nmr2rosetta_.end() ) {
		throw EXCN_UnknownAtomname("proton_name " + proton + " not recognized for aminoacid " + aa_name() );
	}
	return it->second;
}

MethylNames::AtomList const& MethylNames::rosetta2methyl( std::string const& proton ) const {
	NameTable::const_iterator it =  rosetta2methyl_.find( proton );
	if ( it == rosetta2methyl_.end() ) {
		throw EXCN_UnknownAtomname("proton_name " + proton + " not recognized for aminoacid " + aa_name() );
	}
	return it->second;
}

core::Size MethylNames::proton_index( std::string const& proton ) const {
	Size ct( 1 );
	for ( NameTable::const_iterator it = begin(); it != end(); ++it, ++ct ) {
		if ( it->first == proton ) return ct;
	}
	throw EXCN_UnknownAtomname("proton_name " + proton + " not recognized for aminoacid " + aa_name() );
	return 0;
}

}
}
