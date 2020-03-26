// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_DimerPairing_hh
#define INCLUDED_protocols_fldsgn_topology_DimerPairing_hh

#include <core/types.hh>
#include <utility/VirtualBase.hh>
#include <protocols/fldsgn/topology/DimerPairing.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace fldsgn {
namespace topology {


class DimerPairing  : public utility::VirtualBase {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from VirtualBase
	~DimerPairing() override;


	typedef core::Size Size;
	typedef core::Real Real;


public:

	DimerPairing (
		core::Size const res1,
		core::Size const res2,
		Real const dist,
		Real const phi,
		Real const theta,
		Real const sigma,
		Real const dp,
		core::Size const sign1,
		core::Size const sign2,
		Real const score );

	core::Size
	res1() const
	{
		return res1_;
	}

	core::Size
	res2() const
	{
		return res2_;
	}

	Real
	dist() const
	{
		return dist_;
	}

	Real
	phi() const
	{
		return phi_;
	}

	Real
	theta() const
	{
		return theta_;
	}

	Real
	dp() const
	{
		return dp_;
	}

	Real
	sigma() const
	{
		return sigma_;
	}

	core::Size
	sign1() const
	{
		return sign1_;
	}

	core::Size
	sign2() const
	{
		return sign2_;
	}

	Real
	score() const
	{
		return score_;
	}

	char
	orient() const
	{
		return orient_;
	}

	bool
	valid() const
	{
		return valid_;
	}

	void
	valid( bool const v )
	{
		valid_ = v;
	}


public:


	bool is_parallel( Real const phi, Real const theta );


public:


	friend
	std::ostream & operator<<( std::ostream & out, const DimerPairing &dp );


private: // data


	core::Size res1_;
	core::Size res2_;
	Real dist_;
	Real phi_;
	Real theta_;
	Real sigma_;
	Real dp_;
	core::Size sign1_;
	core::Size sign2_;
	Real score_;
	char orient_;
	bool valid_;

};

class DimerPairings : public utility::vector1< DimerPairingOP > {
public:

	typedef core::Size Size;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::DimerPairing DimerPairing;

public:

	void finalize( SS_Info2 const & ss_info );

	friend
	std::ostream & operator<<( std::ostream& out, const DimerPairings &dps );

};


} // ns topology
} // ns fldsgn
} // ns protocols

#endif
