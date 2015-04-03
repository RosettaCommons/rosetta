// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/SS_Killhairpins_Info.hh
/// @brief  Scoring manager class header
/// @author Robert Vernon (rvernon@u.washington.edu)


#ifndef INCLUDED_core_scoring_SS_Killhairpins_Info_hh
#define INCLUDED_core_scoring_SS_Killhairpins_Info_hh

/// Unit headers

/// Package headers
#include <core/types.hh>
#include <basic/datacache/CacheableData.hh>

// utility headers

/// Utility headers

/// Numeric headers

// ObjexxFCL Headers

// C++ headers
#include <iosfwd>

#include <utility/vector1.hh>
#include <utility/io/izstream.fwd.hh>


namespace core {
namespace scoring {

//////////////////////////////////////////////////////////////////////////////////////////////////////

struct Hairpin {

	std::pair< std::pair< core::Size, core::Size>, std::pair< core::Size, core::Size> > range_pair_;

	Hairpin();

	Hairpin( core::Size s1_1, core::Size s1_2, core::Size s2_1, core::Size s2_2);

	~Hairpin();

 	core::Size s1_start() const;

	core::Size s1_end() const;

	core::Size s2_start() const;

	core::Size s2_end() const;

	/// @brief copy assignment
	Hairpin const &
	operator =( Hairpin const & s );

	friend
	std::ostream &
	operator<< ( std::ostream & out, Hairpin const & s );

};

struct Hairpins {

	utility::vector1< Hairpin > hairpin_list_;

	/// @brief default constructor
	Hairpins();

	/// @brief copy constructor
	Hairpins(
		Hairpins const & s
	);

	/// @brief default destructor
	~Hairpins();

	/// @brief copy assignment
	Hairpins const &
	operator =( Hairpins const & s );

	void append_hairpin( core::Size s1_1, core::Size s1_2, core::Size s2_1, core::Size s2_2);

	void clear();

	utility::vector1< Hairpin > list() const;

	core::Size size() const;

	friend
	std::ostream &
	operator<< ( std::ostream & out, Hairpins const & s );

};


//////////////////////////////////////////////////////////////////////////////////////////////////////
class SS_Killhairpins_Info : public basic::datacache::CacheableData {
public:

	SS_Killhairpins_Info();

	SS_Killhairpins_Info( SS_Killhairpins_Info const & src );


	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new SS_Killhairpins_Info( *this ) );
	}

	inline
	Hairpins const &
	hairpins() const
	{
		return hairpins_;
	}

	inline
	Hairpins &
	hairpins()
	{
		return hairpins_;
	}


	//Returns true if the two supplied residues numbers fall within a "killed" hairpin
	bool
	check_hairpin( core::Size const & strand1_res, core::Size const & strand2_res );

	void
	setup_from_psipred(utility::io::izstream & input_file);

	void
	setup_from_kill_hairpins_file(utility::io::izstream & input_file);

	void
	setup_killhairpins();

	inline
	bool kill_parallel() const {
		return kill_parallel_;
	}

	inline
	bool kill_antiparallel() const {
		return kill_antiparallel_;
	}


private:

	bool kill_parallel_;

	bool kill_antiparallel_;

	Hairpins hairpins_;

};

} // ns scoring
} // ns core

#endif
