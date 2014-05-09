// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/silent/SilentStructCreator.hh
/// @brief  Base class for SilentStructCreators for the SilentStruct load-time factory registration scheme
/// @author James Thompson

// Unit Headers
#include <core/io/silent/BasicSilentStructCreators.hh>

// Package Headers

#include <core/io/silent/SilentStruct.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/RigidBodySilentStruct.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/ScoreJumpFileSilentStruct.hh>

//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {


ProteinSilentStruct_SinglePrecCreator::ProteinSilentStruct_SinglePrecCreator() {}
ProteinSilentStruct_SinglePrecCreator::~ProteinSilentStruct_SinglePrecCreator() {}
SilentStructOP ProteinSilentStruct_SinglePrecCreator::create_silent_struct() const {
	return new ProteinSilentStruct_SinglePrec;
}

std::string ProteinSilentStruct_SinglePrecCreator::keyname() const {
	return "protein_float";
}

// class def for ProteinSilentStruct
ProteinSilentStructCreator::ProteinSilentStructCreator() {}
ProteinSilentStructCreator::~ProteinSilentStructCreator() {}
SilentStructOP ProteinSilentStructCreator::create_silent_struct() const {
	return new ProteinSilentStruct;
}

std::string ProteinSilentStructCreator::keyname() const {
	return "protein";
}

// class def for RNA_SilentStruct
RNA_SilentStructCreator::RNA_SilentStructCreator() {}
RNA_SilentStructCreator::~RNA_SilentStructCreator() {}
SilentStructOP RNA_SilentStructCreator::create_silent_struct() const {
	return new RNA_SilentStruct;
}

std::string RNA_SilentStructCreator::keyname() const {
	return "rna";
}

// class def for BinarySilentStruct
BinarySilentStructCreator::BinarySilentStructCreator() {}
BinarySilentStructCreator::~BinarySilentStructCreator() {}
SilentStructOP BinarySilentStructCreator::create_silent_struct() const {
	return new BinarySilentStruct;
}

std::string BinarySilentStructCreator::keyname() const {
	return "binary";
}

// class def for ScoreFileSilentStruct
ScoreFileSilentStructCreator::ScoreFileSilentStructCreator() {}
ScoreFileSilentStructCreator::~ScoreFileSilentStructCreator() {}
SilentStructOP ScoreFileSilentStructCreator::create_silent_struct() const {
	return new ScoreFileSilentStruct;
}

std::string ScoreFileSilentStructCreator::keyname() const {
	return "score";
}

// class def for ScoreJumpFileSilentStruct
ScoreJumpFileSilentStructCreator::ScoreJumpFileSilentStructCreator() {}
ScoreJumpFileSilentStructCreator::~ScoreJumpFileSilentStructCreator() {}
SilentStructOP ScoreJumpFileSilentStructCreator::create_silent_struct() const {
  return new ScoreJumpFileSilentStruct;
}

std::string ScoreJumpFileSilentStructCreator::keyname() const {
  return "score_jump";
}


// class def for ProteinSilentStruct
RigidBodySilentStructCreator::RigidBodySilentStructCreator() {}
RigidBodySilentStructCreator::~RigidBodySilentStructCreator() {}
SilentStructOP RigidBodySilentStructCreator::create_silent_struct() const {
	return new RigidBodySilentStruct;
}

std::string RigidBodySilentStructCreator::keyname() const {
	return "rigid_body";
}


} //namespace silent
} //namespace io
} //namespace core
