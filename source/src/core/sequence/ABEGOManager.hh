// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/sequence/ABEGOManager.hh
/// @brief header file for class of ABEGO plus
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_core_sequence_ABEGOManager_hh
#define INCLUDED_core_sequence_ABEGOManager_hh

#include <core/sequence/ABEGOManager.fwd.hh>

// package headers
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <core/types.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace sequence {

/// @brief struct
struct Line {


	/// @brief default constructor
	inline
	Line() :
		slope( 0.0 ),
		intercept( 0.0 ),
		region( true ) // true means region above the line is indicated
	{}


	/// @brief value constructor
	inline
	Line(
		Real const r1,
		Real const r2,
		bool const b // true means region above the line is indicated
	) :
		slope( r1 ),
		intercept( r2 ),
		region( b )
	{}

	/// @brief copy constructor
	inline
	Line( Line const & ) = default;


	/// @brief default destructor
	inline
	~Line() = default;


	/// @brief copy assignment
	inline
	Line & operator =( Line const & rval ) {
		if ( this != &rval ) {
			slope = rval.slope;
			intercept = rval.intercept;
			region = rval.region;
		}
		return *this;
	}

	/// @brief slope of line
	Real slope;

	/// @brief intercept of line
	Real intercept;

	/// @brief region to be selected, true: above the line, false below the line
	bool region;

};

/// @brief abego elments
class ABEGO {
public:

	typedef core::Real Real;
	typedef std::string String;

public:

	/// @brief default constructor
	ABEGO() :
		name_( 'X' ),
		phi_min_( 0.0 ),
		phi_max_( 0.0 ),
		psi_min_( 0.0 ),
		psi_max_( 0.0 ),
		cis_omega_( false )
	{}

	/// @brief value constructor
	ABEGO( char const & name,
		Real phi_min,
		Real phi_max,
		Real psi_min,
		Real psi_max,
		bool cis_omega ) :
		name_( name ),
		phi_min_( phi_min ),
		phi_max_( phi_max ),
		psi_min_( psi_min ),
		psi_max_( psi_max ),
		cis_omega_( cis_omega )
	{
		runtime_assert( phi_min_ <= phi_max_ );
		runtime_assert( psi_min_ <= psi_max_ );
	}

	/// @brief destrurctor
	~ABEGO() = default;

public: // accessor

	inline char name()    { return name_; }
	inline Real phi_max() { return phi_max_; }
	inline Real phi_min() { return phi_min_; }
	inline Real psi_max() { return psi_max_; }
	inline Real psi_min() { return psi_min_; }
	inline bool cis_omega() { return cis_omega_; }

	void add_line( Real const slope, Real const intercept, bool const region );

public: //

	/// @brief check input torsion angle are in a given abego region
	bool check_rama2( Real const & phi, Real const & psi );

	/// @brief check input torsion angle are in a given abego region
	bool check_rama( Real const & phi, Real const & psi, Real const & omega );

private: // data

	char name_;
	Real phi_min_;
	Real phi_max_;
	Real psi_min_;
	Real psi_max_;
	bool cis_omega_;

	utility::vector1< Line > lines_;

};


/// @brief manager for abego
class ABEGOManager : public utility::pointer::ReferenceCount {
public: // typedef


	typedef core::Size Size;
	typedef core::Real Real;
	typedef std::string String;
	typedef core::pose::Pose Pose;

public:

	/// @brief default constructor
	ABEGOManager();

	/// @brief value constructor
	~ABEGOManager() override ; // auto-removing definition from header{}

	/// @brief copy constructor
	ABEGOManager( ABEGOManager const & rval );


public:


	/// @brief total number of abego definition
	Size total_number_abego() { return totnum_abego_; }


public:


	/// @brief initialize
	void initialize();

	/// @brief check input torsion angle are in a given abego region
	bool check_rama( char const & symbol, Real const & phi, Real const & psi, Real const & omega );

	/// @brief get abego index from torsion angles
	Size torsion2index( Real const phi, Real const psi, Real const omega, Size const level=1 );

	/// @brief get abego index from torsion angles at level 1
	Size torsion2index_level1( Real const phi, Real const psi, Real const omega );

	/// @brief get abego index from torsion angles at level 2
	Size torsion2index_level2( Real const phi, Real const psi, Real const omega );

	/// @brief get abego index from torsion angles at level 3
	Size torsion2index_level3( Real const phi, Real const psi, Real const omega );

	/// @brief get abego index from torsion angles at level 3
	Size torsion2index_level4( Real const phi, Real const psi, Real const omega );

	/// @brief all output level in current setup
	Size alllevel() { return 3; }

	/// @brief transform abego symbol to index
	Size symbol2index( char const & symbol );

	/// @brief transform abego symbol string to base5 index
	Size symbolString2base5index( std::string symbolString);

	/// @brief transform abego string to abego base5 index
	std::string base5index2symbolString( Size base5index,Size length);

	/// @brief transform abego index to symbol
	char index2symbol( Size const & idx );

	/// @brief get abego sequence from pose
	utility::vector1< String > get_symbols( Pose const & pose, Size const level=1 );

	/// @brief get abego sequence from pose
	utility::vector1< String > get_symbols( Pose const & pose, Size const begin, Size const end, Size const level );

	/// @brief get abego string
	String get_abego_string( utility::vector1< String > abego );


private: // data

	/// @brief total number of abego symbols
	Size totnum_abego_;

	/// @brief map relating the index to ABEGO class
	std::map< Size, ABEGO > name2abego_;


};

/// @brief utility for getting abego
utility::vector1< std::string >
get_abego( core::pose::Pose const & pose, core::Size const level=1 );


} // namespace util
} // namespace core

#endif /* INCLUDED_core_sequence_hh */
