// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/StorePoseSnapshot.hh
/// @brief Stores current residue indices in the pose as a ReferencePose.  As residues are added or subtracted,
/// this permits the user to set up movers based on the current residue indices rather than the modified indices
/// of a future state.  Header files for the mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker Laboratory.

#ifndef INCLUDED_protocols_simple_moves_StorePoseSnapshot_hh
#define INCLUDED_protocols_simple_moves_StorePoseSnapshot_hh

#include <protocols/simple_moves/StorePoseSnapshot.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

#include <utility/vector1.hh>


// Utility headers

// C++ headers

// Unit headers

namespace protocols {
	namespace simple_moves {

		/// @brief A mover to store current residue indices in a map that will be updated as residues are added or deleted,
		/// permitting residue indices from the current state to be used to set up movers applied in the future.
		class StorePoseSnapshot : public protocols::moves::Mover
		{
			private:
				typedef protocols::moves::Mover parent;
			public:
				/// @brief Default constructor.
				///
				StorePoseSnapshot();
				
				/// @brief Copy constructor.
				///
				StorePoseSnapshot( StorePoseSnapshot const &src );
				
				/// @brief Destructor.
				///
				virtual ~StorePoseSnapshot();

				/// @brief Apply function -- actually apply this mover to the pose, modifying the pose.
				///
				virtual void apply( core::pose::Pose & pose );

				/// @brief Get the mover name.
				///
				virtual std::string get_name() const;

				/// @brief Makea  copy of this mover, and return an owning pointer to the copy.
				///
				virtual protocols::moves::MoverOP clone() const {
					return (protocols::moves::MoverOP( new protocols::simple_moves::StorePoseSnapshot( *this ) ) );
				}
				
				/// @brief Create a new instance of this mover, initialized to default settings.
				///
				virtual protocols::moves::MoverOP fresh_instance() const {
					return protocols::moves::MoverOP( new StorePoseSnapshot );
				}
				
				/// @brief Parse RosettaScripts XML to set up this mover.
				/// @details This is called at script initialization, long before the apply()
				/// function is called.
				void parse_my_tag( utility::tag::TagCOP tag,
					basic::datacache::DataMap &,
					protocols::filters::Filters_map const &,
					protocols::moves::Movers_map const &,
					core::pose::Pose const & );
					
				/// @brief Set the name of the reference pose object that will be created and stored in the pose.
				///
				void set_reference_pose_name( std::string const &name_in );
				
				/// @brief Return the name of the reference pose object that will be created and stored in the pose.
				///
				std::string reference_pose_name( ) const;

			private:
			
			/// @brief The name of the ReferencePose object that will be created and stored in the pose.
			///
			std::string reference_pose_name_;
			
		};

	} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_StorePoseSnapshot_HH_
