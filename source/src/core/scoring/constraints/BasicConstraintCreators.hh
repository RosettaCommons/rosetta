// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

/// @brief Constraint creator for the AtomPairConstraint constraint
class AtomPairConstraintCreator : public ConstraintCreator
{
public:
	AtomPairConstraintCreator();
	~AtomPairConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the NamedAtomPairConstraint constraint
class NamedAtomPairConstraintCreator : public ConstraintCreator
{
public:
	NamedAtomPairConstraintCreator();
	~NamedAtomPairConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the BasePairConstraint constraint
class BasePairConstraintCreator : public ConstraintCreator
{
public:
	BasePairConstraintCreator();
	~BasePairConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the AngleConstraint constraint
class AngleConstraintCreator : public ConstraintCreator
{
public:
	AngleConstraintCreator();
	~AngleConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the DihedralConstraint constraint
class DihedralConstraintCreator : public ConstraintCreator
{
public:
	DihedralConstraintCreator();
	~DihedralConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for DihedralPairConstraint
class DihedralPairConstraintCreator : public ConstraintCreator
{
public:
	DihedralPairConstraintCreator();
	~DihedralPairConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the BigBinConstraint constraint
class BigBinConstraintCreator : public ConstraintCreator
{
public:
	BigBinConstraintCreator();
	~BigBinConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the MultiConstraint constraint
class MultiConstraintCreator : public ConstraintCreator
{
public:
	MultiConstraintCreator();
	~MultiConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the AmbiguousConstraint constraint
class AmbiguousConstraintCreator : public ConstraintCreator
{
public:
	AmbiguousConstraintCreator();
	~AmbiguousConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the KofNConstraint constraint
class KofNConstraintCreator : public ConstraintCreator
{
public:
	KofNConstraintCreator();
	~KofNConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the CoordinateConstraint constraint
class CoordinateConstraintCreator : public ConstraintCreator
{
public:
	CoordinateConstraintCreator();
	~CoordinateConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the LocalCoordinateConstraint constraint
class LocalCoordinateConstraintCreator : public ConstraintCreator
{
public:
	LocalCoordinateConstraintCreator();
	~LocalCoordinateConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the AmbiguousNMRDistanceConstraint constraint
class AmbiguousNMRDistanceConstraintCreator : public ConstraintCreator
{
public:
	AmbiguousNMRDistanceConstraintCreator();
	~AmbiguousNMRDistanceConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};
/// @brief Constraint creator for the AmbiguousNMRConstraint constraint
class AmbiguousNMRConstraintCreator : public ConstraintCreator
{
public:
	AmbiguousNMRConstraintCreator();
	~AmbiguousNMRConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the SiteConstraint constraint
class SiteConstraintCreator : public ConstraintCreator
{
public:
	SiteConstraintCreator();
	~SiteConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the SiteConstraintResidues constraint
class SiteConstraintResiduesCreator : public ConstraintCreator
{
public:
	SiteConstraintResiduesCreator();
	~SiteConstraintResiduesCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the FabConstraint constraint
class FabConstraintCreator : public ConstraintCreator
{
public:
	FabConstraintCreator();
	~FabConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

/// @brief Constraint creator for the NamedAngleConstraint
class NamedAngleConstraintCreator : public ConstraintCreator
{
public:
	NamedAngleConstraintCreator();
	~NamedAngleConstraintCreator() override;

	ConstraintOP create_constraint() const override;
	std::string keyname() const override;
};

} //namespace constraints
} //namespace scoring
} //namespace core

#endif
