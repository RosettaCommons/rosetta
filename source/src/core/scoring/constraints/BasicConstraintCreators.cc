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

// Unit Headers
#include <core/scoring/constraints/BasicConstraintCreators.hh>

// Package Headers
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BigBinConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/DihedralPairConstraint.hh>
//#include <core/scoring/constraints/DunbrackConstraint.hh>
#include <core/scoring/constraints/KofNConstraint.hh>
#include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/NamedAngleConstraint.hh>
#include <core/scoring/constraints/SiteConstraint.hh>
#include <core/scoring/constraints/SiteConstraintResidues.hh>
#include <core/scoring/constraints/FabConstraint.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <boost/bind.hpp>

namespace core {
namespace scoring {
namespace constraints {

AmbiguousConstraintCreator::AmbiguousConstraintCreator() {}
AmbiguousConstraintCreator::~AmbiguousConstraintCreator() {}

ConstraintOP AmbiguousConstraintCreator::create_constraint() const {
	return ConstraintOP( new AmbiguousConstraint );
}

std::string AmbiguousConstraintCreator::keyname() const
{
	return "AmbiguousConstraint";
}

AmbiguousNMRConstraintCreator::AmbiguousNMRConstraintCreator() {}
AmbiguousNMRConstraintCreator::~AmbiguousNMRConstraintCreator() {}

ConstraintOP AmbiguousNMRConstraintCreator::create_constraint() const {
	return ConstraintOP( new AmbiguousNMRConstraint );
}

std::string AmbiguousNMRConstraintCreator::keyname() const
{
	return "AmbiguousNMRConstraint";
}

AmbiguousNMRDistanceConstraintCreator::AmbiguousNMRDistanceConstraintCreator() {}
AmbiguousNMRDistanceConstraintCreator::~AmbiguousNMRDistanceConstraintCreator() {}

ConstraintOP AmbiguousNMRDistanceConstraintCreator::create_constraint() const {
	return ConstraintOP( new AmbiguousNMRDistanceConstraint );
}

std::string AmbiguousNMRDistanceConstraintCreator::keyname() const
{
	return "AmbiguousNMRDistance";
}

AngleConstraintCreator::AngleConstraintCreator() {}
AngleConstraintCreator::~AngleConstraintCreator() {}

ConstraintOP AngleConstraintCreator::create_constraint() const {
	return ConstraintOP( new AngleConstraint( id::AtomID(), id::AtomID(), id::AtomID(), NULL ) ) ;
}

std::string AngleConstraintCreator::keyname() const
{
	return "Angle";
}

AtomPairConstraintCreator::AtomPairConstraintCreator() {}
AtomPairConstraintCreator::~AtomPairConstraintCreator() {}

ConstraintOP AtomPairConstraintCreator::create_constraint() const {
	return ConstraintOP( new AtomPairConstraint( id::AtomID(), id::AtomID(), NULL ) );
}

std::string AtomPairConstraintCreator::keyname() const
{
	return "AtomPair";
}

BigBinConstraintCreator::BigBinConstraintCreator() {}
BigBinConstraintCreator::~BigBinConstraintCreator() {}

ConstraintOP BigBinConstraintCreator::create_constraint() const {
	return ConstraintOP( new BigBinConstraint );
}

std::string BigBinConstraintCreator::keyname() const
{
	return "BigBin";
}

CoordinateConstraintCreator::CoordinateConstraintCreator() {}
CoordinateConstraintCreator::~CoordinateConstraintCreator() {}

ConstraintOP CoordinateConstraintCreator::create_constraint() const {
	return ConstraintOP( new CoordinateConstraint );
}

std::string CoordinateConstraintCreator::keyname() const
{
	return "CoordinateConstraint";
}

DihedralConstraintCreator::DihedralConstraintCreator() {}
DihedralConstraintCreator::~DihedralConstraintCreator() {}

ConstraintOP DihedralConstraintCreator::create_constraint() const {
	return ConstraintOP( new DihedralConstraint( id::AtomID(), id::AtomID(), id::AtomID(), id::AtomID(), NULL ) );
}

std::string DihedralConstraintCreator::keyname() const
{
	return "Dihedral";
}

DihedralPairConstraintCreator::DihedralPairConstraintCreator() {}
DihedralPairConstraintCreator::~DihedralPairConstraintCreator() {}

ConstraintOP DihedralPairConstraintCreator::create_constraint() const {
	return ConstraintOP( new DihedralPairConstraint(
		id::AtomID(), id::AtomID(), id::AtomID(), id::AtomID(),
		id::AtomID(), id::AtomID(), id::AtomID(), id::AtomID(), NULL ) );
}

std::string DihedralPairConstraintCreator::keyname() const
{
	return "DihedralPair";
}

KofNConstraintCreator::KofNConstraintCreator() {}
KofNConstraintCreator::~KofNConstraintCreator() {}

ConstraintOP KofNConstraintCreator::create_constraint() const {
	return ConstraintOP( new KofNConstraint );
}

std::string KofNConstraintCreator::keyname() const
{
	return "KofNConstraint";
}

LocalCoordinateConstraintCreator::LocalCoordinateConstraintCreator() {}
LocalCoordinateConstraintCreator::~LocalCoordinateConstraintCreator() {}

ConstraintOP LocalCoordinateConstraintCreator::create_constraint() const {
	return ConstraintOP( new LocalCoordinateConstraint );
}

std::string LocalCoordinateConstraintCreator::keyname() const
{
	return "LocalCoordinateConstraint";
}

MultiConstraintCreator::MultiConstraintCreator() {}
MultiConstraintCreator::~MultiConstraintCreator() {}

ConstraintOP MultiConstraintCreator::create_constraint() const {
	return ConstraintOP( new MultiConstraint );
}

std::string MultiConstraintCreator::keyname() const
{
	return "MultiConstraint";
}

SiteConstraintCreator::SiteConstraintCreator() {}
SiteConstraintCreator::~SiteConstraintCreator() {}

ConstraintOP SiteConstraintCreator::create_constraint() const {
    return ConstraintOP( new SiteConstraint );
}

std::string SiteConstraintCreator::keyname() const
{
    return "SiteConstraint";
}

SiteConstraintResiduesCreator::SiteConstraintResiduesCreator() {}
SiteConstraintResiduesCreator::~SiteConstraintResiduesCreator() {}

ConstraintOP SiteConstraintResiduesCreator::create_constraint() const {
    return ConstraintOP( new SiteConstraintResidues );
}

std::string SiteConstraintResiduesCreator::keyname() const
{
    return "SiteConstraintResidues";
}

FabConstraintCreator::FabConstraintCreator() {}
FabConstraintCreator::~FabConstraintCreator() {}

ConstraintOP FabConstraintCreator::create_constraint() const {
    return ConstraintOP( new FabConstraint );
}

std::string FabConstraintCreator::keyname() const
{
    return "FabConstraint";
}

NamedAngleConstraintCreator::NamedAngleConstraintCreator() {}
NamedAngleConstraintCreator::~NamedAngleConstraintCreator() {}

ConstraintOP NamedAngleConstraintCreator::create_constraint() const {
	return ConstraintOP( new NamedAngleConstraint( id::NamedAtomID(), id::NamedAtomID(), id::NamedAtomID(), NULL ) ) ;
}

std::string NamedAngleConstraintCreator::keyname() const
{
	return "NamedAngle";
}


} //namespace constraints
} //namespace scoring
} //namespace core
