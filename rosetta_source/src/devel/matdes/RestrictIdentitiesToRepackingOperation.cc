// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/RestrictIdentitiesToRepackingOperation.cc
/// @brief	Restricts specified residue types to only repack, no design.
/// @author Neil King (neilking@uw.edu)

// Unit Headers
#include <devel/matdes/RestrictIdentitiesToRepackingOperation.hh>
#include <devel/matdes/RestrictIdentitiesToRepackingOperationCreator.hh>

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <core/kinematics/Jump.hh>
#include <boost/functional/hash.hpp>



using basic::Error;
using basic::Warning;
static basic::Tracer TR( "devel.matdes.RestrictIdentitiesToRepackingOperation" );

namespace devel {
namespace matdes {
	
	RestrictIdentitiesToRepackingOperation::RestrictIdentitiesToRepackingOperation() {}
	
	RestrictIdentitiesToRepackingOperation::RestrictIdentitiesToRepackingOperation( utility::vector1 < std::string > identities )
	: identities_( identities )
	{
	}
	
	RestrictIdentitiesToRepackingOperation::~RestrictIdentitiesToRepackingOperation() {}
	
	core::pack::task::operation::TaskOperationOP
	RestrictIdentitiesToRepackingOperationCreator::create_task_operation() const
	{
		return new RestrictIdentitiesToRepackingOperation;
	}
	
	core::pack::task::operation::TaskOperationOP RestrictIdentitiesToRepackingOperation::clone() const
	{
		return new RestrictIdentitiesToRepackingOperation( *this );
	}


	void
	RestrictIdentitiesToRepackingOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
	{
		runtime_assert( identities_.size() != 0 );

		for (Size j=1; j<=identities_.size(); j++) {
			for (Size i=1; i<=pose.n_residue(); i++) {
				std::string aa_name = pose.residue(i).name3();
				if (aa_name == identities_[j]) {
	      	task.nonconst_residue_task(i).prevent_repacking();
					TR.Debug << "restrict to repacking: " << aa_name << i << std::endl;	
				}
    	}
		}
	}

	utility::vector1< std::string >
	RestrictIdentitiesToRepackingOperation::get_identities() const
	{
		return identities_;
	}
	
	void
	RestrictIdentitiesToRepackingOperation::set_identities( utility::vector1  < std::string > identities_vec )
	{
		runtime_assert( identities_vec.size() != 0 );
		identities_.clear();
		foreach( std::string const item, identities_vec )
			identities_.push_back( item );
	}

	void
	RestrictIdentitiesToRepackingOperation::parse_tag( TagPtr tag )
	{
		unparsed_identities_ = tag->getOption< std::string >( "identities" ) ;
		if( unparsed_identities_ != "" ){
			
      utility::vector1< std::string > const ids( utility::string_split( unparsed_identities_ , ',' ) );
			identities_.clear();

			TR << "Only repacking, not designing, residues of type(s): ";
      foreach( std::string const id, ids ){
				TR << id << " ";
       	identities_.push_back( id );
      }
			TR.Debug << std::endl;
    }
	}
	
} //namespace matdes
} //namespace devel
