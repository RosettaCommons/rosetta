// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/noesy_assign/MethylNames.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_MethylNames_HH
#define INCLUDED_protocols_noesy_assign_MethylNames_HH


// Unit Headers
//#include <protocols/noesy_assign/MethylNames.fwd.hh>

// Package Headers
#include <protocols/noesy_assign/FragsToAtomDist.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>

// Utility headers
#include <utility/SingletonBase.hh>

// C++ headers
#include <string>

namespace protocols {
namespace noesy_assign {

class MethylNames {
public:
  typedef utility::vector1< std::string > AtomList;
  typedef std::map< std::string, AtomList > NameTable;
  typedef NameTable::const_iterator const_iterator;
  MethylNames( );
  MethylNames( core::chemical::AA aa );

	//which NMR proton name belongs to the ROSETTA proton ?
	//1HG1 ---> return HG11
  std::string const& rosetta2nmr( std::string const& proton ) const;

	//return which NMR methyl name belongs to the Rosetta proton ?
	//1HG1 --> returns  QG1, QQG
  AtomList const& rosetta2methyl( std::string const& proton ) const;

	//return list of rosetta protons that belong to a proton or methyl NMR name
  AtomList const& nmr2rosetta( std::string const& proton ) const;

	//iterate over protons/methyls
  const_iterator begin() const { return nmr2rosetta_.begin(); }
  const_iterator end() const { return nmr2rosetta_.end(); }

	//the amino-acid type
  core::chemical::AA aa() const { return aa_; }

	//return the amino-acid name (as from chemical::name_from_aa(x) )
  std::string aa_name() const;

	//return index of nmr-name (using the same sequence as begin()...end() implies
	core::Size proton_index( std::string const& ) const;

	//add a nmr / rosetta proton pair
  void add_proton( std::string const& nmr, std::string const& rosetta );
	//add a nmr-methyl / rosetta proton pair
  void add_methyl( std::string const& rosetta, std::string const& methyl );

private:
  core::chemical::AA aa_;
  std::map< std::string, std::string > rosetta2nmr_;
  NameTable rosetta2methyl_;
  NameTable nmr2rosetta_;
};

class MethylNameLibrary : public utility::SingletonBase< MethylNameLibrary > {
public:
	friend class utility::SingletonBase< MethylNameLibrary >;

private:
  //Singleton Class
  MethylNameLibrary();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static MethylNameLibrary * create_singleton_instance();

public:

  MethylNames const& operator[]( core::chemical::AA ) const;

private:
  void load_database_table();

  typedef std::map< core::chemical::AA, MethylNames > MethylNameTable;

  MethylNameTable methyl_names_;
};

}
}

#endif
