// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/FormFactorManager.cc
/// @brief Loads and manages atomic form factors
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <core/scoring/saxs/FormFactorManager.hh>
#include <core/scoring/saxs/FormFactor.hh>

// utility headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <basic/database/open.hh>

// C++
#include <string>
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace saxs {

static THREAD_LOCAL basic::Tracer trFormFactorManager( "core.scoring.saxs.FormFactorManager" );

FormFactorManager* FormFactorManager::singleton_;

FormFactorManager* FormFactorManager::get_manager() {

	if ( singleton_==0 ) {
		singleton_ = new FormFactorManager();
	}

	return singleton_;
}

void FormFactorManager::load_ff(std::string file_name) {

	utility::io::izstream input(file_name.c_str());
	std::string line;

	while ( getline( input, line ) ) {
		if ( line.substr(0,1) == "#" ) continue;
		std::istringstream line_stream( line );
		trFormFactorManager.Debug << "line: " << line << std::endl;
		std::string n,f;
		bool g;
		line_stream >> n >> f >> g;
		FormFactorOP c( new FormFactor(n,f) );
		c->is_glob(g);
		trFormFactorManager.Debug << "Form factor for atom >"<<n<<"< loaded from a file: "<<f<<std::endl;
		register_ff( n,c);
	}
}

void FormFactorManager::load_ff_from_db(std::string file_name) {

	utility::io::izstream input(file_name.c_str());
	std::string line;

	while ( getline( input, line ) ) {
		if ( line.substr(0,1) == "#" ) continue;
		std::istringstream line_stream( line );
		trFormFactorManager.Debug << "line " << line << std::endl;
		std::string n,f;
		bool g;
		line_stream >> n >> f >> g;
		FormFactorOP c( new FormFactor( n, basic::database::full_name("scoring/score_functions/saxs/"+f) ) );
		c->is_glob(g);
		trFormFactorManager.Warning << "Form factor for atom >"<<n<<"< loaded from a minirosetta database file: "<<f<<std::endl;
		register_ff( n,c);
	}
}


void FormFactorManager::register_ff(std::string atom_name,FormFactorOP new_ff) {

	Size n = known_atoms_.size() + 1;
	new_ff->id_ = n;
	if ( ff_map_.find(atom_name) == ff_map_.end() ) {
		ff_map_.insert( std::pair<std::string,FormFactorOP> (atom_name,new_ff) );
		known_atoms_.push_back( atom_name );
		names_to_indexes_.insert( std::pair<std::string,Size> (atom_name, n) );
		ff_vector_.push_back(new_ff);
	}
}

/// @brief returns true if the manager has form factor function for a given atom
bool FormFactorManager::is_known_atom(std::string atom_name) {

	debug_assert( atom_name.size() > 0 );
	if ( ff_map_.find(atom_name) != ff_map_.end() ) {
		return true;
	}
	return false;
}

/// @brief returns form factor function for a given atom
FormFactorOP FormFactorManager::get_ff(std::string atom_name) {

	if ( ff_map_.find(atom_name) != ff_map_.end() ) {
		// trFormFactorManager.Debug << "Returning FF for this atom: "<<atom_name<<std::endl;
		return ff_map_.find(atom_name)->second;
	}

	trFormFactorManager.Debug << "The manager knows nothing about this atom: "<<atom_name<<std::endl;
	return 0;
}


} // saxs
} // scoring
} // core

