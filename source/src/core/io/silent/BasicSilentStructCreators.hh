// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/BasicSilentStructCreator.hh
/// @brief  Base class for BasicSilentStructCreators for the BasicSilentStruct load-time factory registration scheme
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_BasicSilentStructCreators_hh
#define INCLUDED_core_io_silent_BasicSilentStructCreators_hh

// Unit Headers
#include <core/io/silent/SilentStructCreator.hh>

// c++ headers

#include <core/types.hh>
#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {

/// @brief creator for the ProteinSilentStruct_SinglePrec class
class ProteinSilentStruct_SinglePrecCreator : public SilentStructCreator
{
public:
	ProteinSilentStruct_SinglePrecCreator();
	~ProteinSilentStruct_SinglePrecCreator() override;

	SilentStructOP create_silent_struct( SilentFileOptions const & ) const override;
	std::string keyname() const override;
};


/// @brief creator for the ProteinSilentStruct class
class ProteinSilentStructCreator : public SilentStructCreator
{
public:
	ProteinSilentStructCreator();
	~ProteinSilentStructCreator() override;

	SilentStructOP create_silent_struct( SilentFileOptions const & ) const override;
	std::string keyname() const override;
};


/// @brief creator for the RNA_SilentStruct class
class RNA_SilentStructCreator : public SilentStructCreator
{
public:
	RNA_SilentStructCreator();
	~RNA_SilentStructCreator() override;

	SilentStructOP create_silent_struct( SilentFileOptions const & ) const override;
	std::string keyname() const override;
};


/// @brief creator for the BinarySilentStruct class
class BinarySilentStructCreator : public SilentStructCreator
{
public:
	BinarySilentStructCreator();
	~BinarySilentStructCreator() override;

	SilentStructOP create_silent_struct( SilentFileOptions const & ) const override;
	std::string keyname() const override;
};


/// @brief creator for the ScoreFileSilentStruct class
class ScoreFileSilentStructCreator : public SilentStructCreator
{
public:
	ScoreFileSilentStructCreator();
	~ScoreFileSilentStructCreator() override;

	SilentStructOP create_silent_struct( SilentFileOptions const & ) const override;
	std::string keyname() const override;
};


/// @brief creator for the ScoreJumpFileSilentStruct class
class ScoreJumpFileSilentStructCreator : public SilentStructCreator
{
public:
	ScoreJumpFileSilentStructCreator();
	~ScoreJumpFileSilentStructCreator() override;

	SilentStructOP create_silent_struct( SilentFileOptions const & ) const override;
	std::string keyname() const override;
};


/// @brief creator for the RigidBodySilentStruct class
class RigidBodySilentStructCreator : public SilentStructCreator
{
public:
	RigidBodySilentStructCreator();
	~RigidBodySilentStructCreator() override;

	SilentStructOP create_silent_struct( SilentFileOptions const & ) const override;
	std::string keyname() const override;
};


} //namespace silent
} //namespace io
} //namespace core

#endif
