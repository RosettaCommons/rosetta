// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/fldgsn/BetaAlphaBetaMotif.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_BetaAlphaBetaMotif_hh
#define INCLUDED_protocols_fldsgn_topology_BetaAlphaBetaMotif_hh

// unit headers
#include <protocols/fldsgn/topology/BetaAlphaBetaMotif.fwd.hh>

// project headers
#include <core/types.hh>
#include <protocols/fldsgn/topology/Sheet.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <utility/assert.hh>

#include <utility/vector1_bool.hh>
#include <numeric/xyzVector.hh>

namespace protocols {
namespace fldsgn {
namespace topology {

class BetaAlphaBetaMotif : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef core::Real Real;
	typedef std::string String;
	typedef core::Vector Vector;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;


public:// construct/destruct


	/// @brief default constructor
	BetaAlphaBetaMotif();

	/// @brief value constructor
	BetaAlphaBetaMotif(
		Size const & strand1,
		Size const & strand2,
		Size const & helix,
		Size const & cross_over );

	/// @brief copy constructor
	BetaAlphaBetaMotif( BetaAlphaBetaMotif const & s );

	/// @brief destructor
	virtual ~BetaAlphaBetaMotif();


public:// operator


	/// @brief IO Operator
	friend std::ostream & operator<<(std::ostream & out, const BetaAlphaBetaMotif & s );


public:// accessor


	String name() const;

	inline Size helix() const { return helix_; }

	inline Size strand1() const { return strand1_; }

	inline Size strand2() const { return strand2_; }

	inline Size cross_over() const { return cross_over_; }

	inline bool is_lefthanded() const { return left_handed_; }

	inline Real hsheet_dist() const { return hs_dist_; }

	inline Real hs1_dist() const { return hs1_dist_; }

	inline Real hs2_dist() const { return hs2_dist_; }

	inline Real hs_angle() const { return hs_angle_; }

	inline Real hsheet_elev_angle() const { return hsheet_elev_angle_; }

	inline utility::vector1< Size > helix_cycle() const { return helix_cycle_; }

	String helix_cycle_as_string() const;

private:// mutator


	inline void left_handed( bool const v ) { left_handed_ = v; }

	void calc_helix_cycle( SS_Info2_COP const ssinfo );

public:


	Size calc_inout( SS_Info2_COP const ssinfo, Size const resi ) const;

	void calc_geometry( SS_Info2_COP const ssinfo,  SheetSetCOP const sheet_set );


private: /// data

	/// @brief id of strand
	Size strand1_;

	/// @brief id of strand
	Size strand2_;

	/// @brief id of helix
	Size helix_;

	/// @brief number of strands crossed over by bab-motif
	Size cross_over_;

	/// @brief is this left handed ?
	bool left_handed_;

	/// @brief vector to define sheet plane
	Vector sheet_plane_;

	/// @brief one positional vector in beta sheet plane
	Vector sheet_pos_;

	/// @brief distance between helix and sheet
	Real hs_dist_;

	/// @brief angle between helix projected on sheet and strands
	Real hs_angle_;

	/// @brief distance between helix and mid point of strand1_
	Real hs1_dist_;

	/// @brief distance between helix and mid point of strand2_
	Real hs2_dist_;

	/// @brief elevation angle of helix respect to sheet
	Real hsheet_elev_angle_;

	/// @brief cylce of helix against sheet
	utility::vector1< Size > helix_cycle_;

	bool geometry_is_initialized_;


};

////////////////////////////////////////////////////////////////////////////////////////////////////
class BetaAlphaBetaMotifSet : public utility::pointer::ReferenceCount {
public:


	typedef core::Real Real;
	typedef core::Size Size;
	typedef protocols::fldsgn::topology::SheetSetCOP SheetSetCOP;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;


public:// constructor/destructor


	/// @brief default constructor
	BetaAlphaBetaMotifSet();

	/// @brief value constructor
	BetaAlphaBetaMotifSet( BetaAlphaBetaMotifs const & bab_motifs );

	/// @brief value constructor
	BetaAlphaBetaMotifSet( SS_Info2_COP const ssinfo, SheetSetCOP const sheet_set );

	/// @brief copy constructor
	BetaAlphaBetaMotifSet( BetaAlphaBetaMotifSet const & s );

	/// @brief destructor
	virtual ~BetaAlphaBetaMotifSet();


public://


	/// @brief
	void push_back( BetaAlphaBetaMotifOP const bop );

	/// @brief
	void clear();

	/// @brief
	Size size() const { return bab_motifs_.size(); }


public:// accessor


	/// @brief
	BetaAlphaBetaMotifs const & bab_motifs() const;

	/// @brief
	BetaAlphaBetaMotifOP bab_motif( Size const & i ) const;

	/// @brief
	friend std::ostream & operator<<( std::ostream & out, const BetaAlphaBetaMotifSet & s );


public:


	/// @brief
	void set_babmotifs( SS_Info2_COP const ssinfo,  SheetSetCOP const sheet_set );

	/// @brief
	void calc_geometry( SS_Info2_COP const ssinfo,  SheetSetCOP const sheet_set );


private:

	BetaAlphaBetaMotifs bab_motifs_;

};

} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
