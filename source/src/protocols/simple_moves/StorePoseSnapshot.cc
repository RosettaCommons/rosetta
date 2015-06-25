// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/StorePoseSnapshot.cc
/// @brief Stores current residue indices in the pose as a ReferencePose.  As residues are added or subtracted,
/// this permits the user to set up movers based on the current residue indices rather than the modified indices
/// of a future state.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker Laboratory.

// Unit headers
#include <protocols/simple_moves/StorePoseSnapshot.hh>
#include <protocols/simple_moves/StorePoseSnapshotCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

//parsing
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>


// Utility Headers

// Unit Headers

// C++ headers


namespace protocols {
	namespace simple_moves {

		using namespace core;
		using namespace core::chemical;
		using namespace std;
		using namespace numeric::conversions;

		using core::pose::Pose;
		using core::conformation::Residue;

		static thread_local basic::Tracer TR( "protocols.simple_moves.StorePoseSnapshot" );

		std::string
		StorePoseSnapshotCreator::keyname() const
		{
			return StorePoseSnapshotCreator::mover_name();
		}

		protocols::moves::MoverOP
		StorePoseSnapshotCreator::create_mover() const {
			return protocols::moves::MoverOP( new StorePoseSnapshot );
		}

		std::string
		StorePoseSnapshotCreator::mover_name()
		{
			return "StorePoseSnapshot";
		}

		StorePoseSnapshot::~StorePoseSnapshot() {}

		/// @brief Default constructor
		///
		StorePoseSnapshot::StorePoseSnapshot() :
			parent(),
			reference_pose_name_( "" )
		{}

		/// @brief Copy constructor
		///
		StorePoseSnapshot::StorePoseSnapshot( StorePoseSnapshot const &src ) :
			parent(src),
			reference_pose_name_( src.reference_pose_name_ )
		{}

		/// @brief Apply function -- actually apply this mover to the pose, modifying the pose.
		///
		void StorePoseSnapshot::apply( Pose & pose ) {
			runtime_assert_string_msg(
				reference_pose_name()!="",
				"Error in protocols::simple_moves::StorePoseSnapshot::apply(): The reference pose name is currently blank.  This must be set before the apply() function is called."
			);
			pose.reference_pose_from_current( reference_pose_name() ); //Yes, that's all this mover does.  It's that simple.
			if (TR.visible()) {
				TR << "Stored pose snapshot " << reference_pose_name() << "." << std::endl;
				TR.flush();
			}
			return;
		}

		/// @brief Get the mover name.
		///
		std::string
		StorePoseSnapshot::get_name() const {
			return StorePoseSnapshotCreator::mover_name();
		}

		/// @brief Parse RosettaScripts XML to set up this mover.
		/// @details This is called at script initialization, long before the apply()
		/// function is called.
		void StorePoseSnapshot::parse_my_tag( utility::tag::TagCOP tag,
				basic::datacache::DataMap &,
				protocols::filters::Filters_map const &,
				protocols::moves::Movers_map const &,
				Pose const & //pose
		)
		{
			runtime_assert_string_msg( tag->hasOption("reference_pose_name"), "Error in protocols::simple_moves::StorePoseSnapshot::parse_my_tag():  When parsing options for the StorePoseSnapshot mover, no \"reference_pose_name\" option was found.  This is required." );
			
			set_reference_pose_name( tag->getOption< std::string >( "reference_pose_name", "" ) );
			if (TR.visible()) TR << "Set reference pose name to " << reference_pose_name() << "." << std::endl;

			if (TR.visible()) TR.flush();
			return;
		}
		
		/// @brief Set the name of the reference pose object that will be created and stored in the pose.
		///
		void StorePoseSnapshot::set_reference_pose_name( std::string const &name_in ) {
			runtime_assert_string_msg( name_in!="", "Error in protocols::simple_moves::StorePoseSnapshot::set_reference_pose_name(): The name cannot be an empty string." );
			reference_pose_name_ = name_in;
			return;
		}
		
		/// @brief Return the name of the reference pose object that will be created and stored in the pose.
		///
		std::string StorePoseSnapshot::reference_pose_name( ) const { return reference_pose_name_; }

	} // simple_moves
} // protocols
