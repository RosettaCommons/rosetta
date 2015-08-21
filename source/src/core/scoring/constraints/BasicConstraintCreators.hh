// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/ConstraintCreator.hh
/// @brief  Base class for ConstraintCreators for the Constraint load-time factory registration scheme
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_constraints_BasicConstraintCreators_hh
#define INCLUDED_core_scoring_constraints_BasicConstraintCreators_hh

// Unit Headers
#include <core/scoring/constraints/ConstraintCreator.hh>

// c++ headers

namespace core {
namespace scoring {
namespace constraints {

/// @brief Mover creator for the AtomPairConstraint constraint
class AtomPairConstraintCreator : public ConstraintCreator
{
public:
	AtomPairConstraintCreator();
	virtual ~AtomPairConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the AngleConstraint constraint
class AngleConstraintCreator : public ConstraintCreator
{
public:
	AngleConstraintCreator();
	virtual ~AngleConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the DihedralConstraint constraint
class DihedralConstraintCreator : public ConstraintCreator
{
public:
	DihedralConstraintCreator();
	virtual ~DihedralConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Constraint creator for DihedralPairConstraint
class DihedralPairConstraintCreator : public ConstraintCreator
{
public:
	DihedralPairConstraintCreator();
	virtual ~DihedralPairConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the BigBinConstraint constraint
class BigBinConstraintCreator : public ConstraintCreator
{
public:
	BigBinConstraintCreator();
	virtual ~BigBinConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the MultiConstraint constraint
class MultiConstraintCreator : public ConstraintCreator
{
public:
	MultiConstraintCreator();
	virtual ~MultiConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the AmbiguousConstraint constraint
class AmbiguousConstraintCreator : public ConstraintCreator
{
public:
	AmbiguousConstraintCreator();
	virtual ~AmbiguousConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the KofNConstraint constraint
class KofNConstraintCreator : public ConstraintCreator
{
public:
	KofNConstraintCreator();
	virtual ~KofNConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the CoordinateConstraint constraint
class CoordinateConstraintCreator : public ConstraintCreator
{
public:
	CoordinateConstraintCreator();
	virtual ~CoordinateConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the LocalCoordinateConstraint constraint
class LocalCoordinateConstraintCreator : public ConstraintCreator
{
public:
	LocalCoordinateConstraintCreator();
	virtual ~LocalCoordinateConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the AmbiguousNMRDistanceConstraint constraint
class AmbiguousNMRDistanceConstraintCreator : public ConstraintCreator
{
public:
	AmbiguousNMRDistanceConstraintCreator();
	virtual ~AmbiguousNMRDistanceConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};
/// @brief Mover creator for the AmbiguousNMRConstraint constraint
class AmbiguousNMRConstraintCreator : public ConstraintCreator
{
public:
	AmbiguousNMRConstraintCreator();
	virtual ~AmbiguousNMRConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the SiteConstraint constraint
class SiteConstraintCreator : public ConstraintCreator
{
public:
	SiteConstraintCreator();
	virtual ~SiteConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the SiteConstraintResidues constraint
class SiteConstraintResiduesCreator : public ConstraintCreator
{
public:
	SiteConstraintResiduesCreator();
	virtual ~SiteConstraintResiduesCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Mover creator for the FabConstraint constraint
class FabConstraintCreator : public ConstraintCreator
{
public:
	FabConstraintCreator();
	virtual ~FabConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

/// @brief Constraint creator for the NamedAngleConstraint
class NamedAngleConstraintCreator : public ConstraintCreator
{
public:
	NamedAngleConstraintCreator();
	virtual ~NamedAngleConstraintCreator();

	virtual ConstraintOP create_constraint() const;
	virtual std::string keyname() const;
};

} //namespace constraints
} //namespace scoring
} //namespace core

#endif
