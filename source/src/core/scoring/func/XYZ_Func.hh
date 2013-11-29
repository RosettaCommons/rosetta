// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#ifndef INCLUDED_core_scoring_constraints_XYZ_Func_hh
#define INCLUDED_core_scoring_constraints_XYZ_Func_hh

#include <core/scoring/func/XYZ_Func.fwd.hh>

#include <core/types.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/id/AtomID.fwd.hh>

//Utility Headers

// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
/// helper classes to reuse code
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

/// just a simple class that takes an AtomID and returns a Vector (position)
class XYZ_Func {
public:
	typedef id::AtomID AtomID;
	typedef conformation::Residue Residue;
	typedef conformation::Conformation Conformation;

public:

	virtual
	Vector const &
	operator()( AtomID const & id ) const = 0;


	virtual
	Residue const &
	residue( Size seqpos ) const = 0;

	virtual
	~XYZ_Func();
};


class ResidueXYZ : public XYZ_Func {
public:

	ResidueXYZ( Residue const & rsd_in );

	virtual
	Vector const &
	operator()( AtomID const & id ) const;

	virtual
	Residue const &
	residue( Size seqpos ) const;

	virtual
	~ResidueXYZ();

private:
	Residue const & rsd_;

};

class ResiduePairXYZ : public XYZ_Func {
public:

	ResiduePairXYZ( Residue const & rsd1_in, Residue const & rsd2_in );

	virtual
	Vector const &
	operator()( AtomID const & id ) const;

	virtual
	Residue const &
	residue( Size seqpos ) const;

	virtual
	~ResiduePairXYZ();

private:
	Residue const & rsd1_;
	Residue const & rsd2_;

};

class ConformationXYZ : public XYZ_Func {
public:

	ConformationXYZ( Conformation const & conformation_in );

	virtual
	Vector const &
	operator()( AtomID const & id ) const;

	virtual
	Residue const &
	residue( Size seqpos ) const;

	virtual
	~ConformationXYZ();

private:
	Conformation const & conformation_;

};

// class BBTorsionXYZ : public ResidueXYZ {
//
// public:
//
// 	virtual
// 	~BBTorsionXYZ;
//
// private:
// 	Residue const & rsd_;
// };


} // constraints
} // scoring
} // core

#endif
