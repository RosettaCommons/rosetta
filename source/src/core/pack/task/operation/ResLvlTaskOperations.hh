// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResLvlTaskOperations.hh
/// @brief  core-level (very general) derived classes that wrap widely-used methods of the ResidueLevelTask interface. These are used by higher-level TaskOperations that allow the user to configure the behavior of PackerTasks that are created by TaskFactory.
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResLvlTaskOperations_hh
#define INCLUDED_core_pack_task_operation_ResLvlTaskOperations_hh

// Unit Headers
#include <core/pack/task/operation/ResLvlTaskOperations.fwd.hh>

#include <core/pack/task/operation/ResLvlTaskOperation.hh>

// Project Headers
#include <core/pack/task/PackerTask.fwd.hh>


#include <utility/vector1.hh>
#include <string>


namespace core {
namespace pack {
namespace task {
namespace operation {

class RestrictToRepackingRLT : public ResLvlTaskOperation
{
public:
	typedef ResLvlTaskOperation parent;
public:
	virtual ~RestrictToRepackingRLT();
	virtual ResLvlTaskOperationOP clone() const;
	virtual void apply( ResidueLevelTask & ) const;
};

class RestrictAbsentCanonicalAASRLT : public ResLvlTaskOperation
{
public:
	typedef ResLvlTaskOperation parent;
public:
	RestrictAbsentCanonicalAASRLT();
	virtual ~RestrictAbsentCanonicalAASRLT();
	virtual ResLvlTaskOperationOP clone() const;
	virtual void apply( ResidueLevelTask & ) const;
	// if an amino acid is not present (false) in the boolean vector, then do not allow it at this position.  The boolean vector is a 20-length vector in alphabetical order by one-letter code.
	void aas_to_keep( utility::vector1< bool > const & );
	void aas_to_keep( std::string const & );
	virtual void parse_tag( TagCOP );
private:
	utility::vector1< bool > canonical_aas_to_keep_;
};

class DisallowIfNonnativeRLT:  public ResLvlTaskOperation
{
public:
	typedef ResLvlTaskOperation parent;
public:
	DisallowIfNonnativeRLT();
	DisallowIfNonnativeRLT( utility::vector1< bool > disallowed_aas );
	virtual ~DisallowIfNonnativeRLT();
	virtual ResLvlTaskOperationOP clone() const;
	virtual void apply(  ResidueLevelTask & rlt ) const;
	//helper functions to define desired AAs
	void clear();
	//define as true which residues are NOT allowed
	void disallow_aas( utility::vector1< bool > const & cannonical_disallowed );
	void disallow_aas( std::string const & aa_string );
	virtual void parse_tag( TagCOP );

private:
	utility::vector1< bool > invert_vector( utility::vector1< bool > disallowed_aas);
	utility::vector1< bool > disallowed_aas_;
	utility::vector1< bool > allowed_aas_;
};

class PreventRepackingRLT : public ResLvlTaskOperation
{
public:
	typedef ResLvlTaskOperation parent;
public:
	virtual ~PreventRepackingRLT();
	virtual ResLvlTaskOperationOP clone() const;
	virtual void apply( ResidueLevelTask & ) const;
};

class AddBehaviorRLT : public ResLvlTaskOperation
{
public:
	typedef ResLvlTaskOperation parent;
public:
	AddBehaviorRLT();
	AddBehaviorRLT( std::string const & behavior );
	virtual ~AddBehaviorRLT();
	virtual ResLvlTaskOperationOP clone() const;
	virtual void apply( ResidueLevelTask & ) const;
	virtual void parse_tag( TagCOP );
private:
	std::string behavior_;
};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
