// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoringManager.hh
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_SS_Info_hh
#define INCLUDED_core_scoring_SS_Info_hh

/// Unit headers
#include <core/scoring/SS_Info.fwd.hh>

/// Package headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/CacheableData.hh>

/// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>

/// Numeric headers
#include <numeric/xyzVector.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <string>
#include <iosfwd>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

//////////////////////////////////////////////////////////////////////////////////////////////////////
class BB_Pos {

public:


	void
	resize( int const nres );

	void
	take_coordinates_from_pose( pose::Pose const & pose );

	/// @details accessor for N's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	N( int const i ) const
	{
		return N_[i];
	}

	/// @details accessor for CA's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	CA( int const i ) const
	{
		return CA_[i];
	}

	/// @details accessor for CB's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	CB( int const i ) const
	{
		return CB_[i];
	}


	/// @details accessor for C's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	C( int const i ) const
	{
		return C_[i];
	}

	/// @details accessor for O's coordinate on residue i, requires take_coordinates_from_pose have been recently called.
	Vector const &
	O( int const i ) const
	{
		return O_[i];
	}


private:

	bool bbindices_up_to_date( pose::Pose const & pose ) const;
	void update_indices( pose::Pose const & pose );

private: // DATA
	utility::vector1< Vector > N_;
	utility::vector1< Vector > CA_;
	utility::vector1< Vector > CB_;
	utility::vector1< Vector > C_;
	utility::vector1< Vector > O_;

	/// Residue types must match those of the pose for the indices
	/// to match.
	utility::vector1< chemical::ResidueType const * > residue_types_;
	utility::vector1< Size > N_index_;
	utility::vector1< Size > CA_index_;
	utility::vector1< Size > CB_index_;
	utility::vector1< Size > C_index_;
	utility::vector1< Size > O_index_;

};



//////////////////////////////////////////////////////////////////////////////////////////////////////

struct Strands {

	typedef ObjexxFCL::FArray1D< int > FArray1D_int;
	typedef ObjexxFCL::FArray2D< int > FArray2D_int;

	/// @brief default constructor
	Strands();

	/// @brief total residue constructor
	Strands(
		int const & total_residue
	);

	/// @brief copy constructor
	Strands(
		Strands const & s
	);

	/// @brief default destructor
	~Strands();

	///
	void
	resize( Size const nres );

	void clear();

	/// @brief copy assignment
	Strands const &
	operator =( Strands const & s );

	friend
	std::ostream &
	operator<< ( std::ostream & out, Strands const & s );

	/// @brief number of strand dimers
	int total_SS_dimer;

	/// @brief residue number of strand dimer i
	FArray1D_int SS_resnum;

	/// @brief number of strands
	int total_strands;

	/// @brief strand number containing SS_dimer i
	FArray1D_int SS_strand;

	/// @brief dimer number starting with position i
	FArray1D_int SS_dimer;

	/// @brief residue number of first non-E res
	FArray2D_int SS_strand_end;

	/// @brief two neighbors, used for determining sheets
	mutable FArray2D_int dimer_neighbor; // hack to allow updating in sspair

	/// @brief strand-strand score array, resized for each calculation:
	// need to re-add this logic if it turns out to be useful
	//ObjexxFCL::FArray2D< Real > strand_strand_score;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////


struct Helices {

	typedef ObjexxFCL::FArray1D< int > FArray1D_int;
	typedef ObjexxFCL::FArray2D< int > FArray2D_int;

	/// @brief default constructor
	Helices();

	/// @brief total residue constructor
	Helices(
		int const & total_residue
	);

	///
	void
	resize( int const nres );

	/// @brief copy constructor
	Helices(
		Helices const & h
	);

	/// @brief default destructor
	~Helices();

	/// @brief copy assignment
	Helices const &
	operator = ( Helices const & h );

	friend
	std::ostream &
	operator<< ( std::ostream & out, Helices const & s );

	/// @brief number of helices
	int total_helices;

/// PRIVATE: !

	// variables
	int total_HH_dimer;
	FArray1D_int HH_resnum;
	FArray2D_int HH_helix_end;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////
class SS_Info : public basic::datacache::CacheableData {
public:

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new SS_Info( *this ) );
	}

	void
	resize( int const nres )
	{
		bb_pos_.resize( nres );
		strands_.resize( nres );
		helices_.resize( nres );
	}

	///
	inline
	BB_Pos const &
	bb_pos() const
	{
		return bb_pos_;
	}

	///
	inline
	BB_Pos &
	bb_pos()
	{
		return bb_pos_;
	}

	///
	inline
	Strands const &
	strands() const
	{
		return strands_;
	}

	///
	inline
	Strands &
	strands()
	{
		return strands_;
	}

	///
	inline
	Helices const &
	helices() const
	{
		return helices_;
	}

	///
	inline
	Helices &
	helices()
	{
		return helices_;
	}

private:

	BB_Pos bb_pos_;

	Strands strands_;

	Helices helices_;

};

} // ns scoring
} // ns core

#endif
