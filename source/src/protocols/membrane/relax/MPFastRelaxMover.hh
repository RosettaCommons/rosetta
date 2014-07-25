// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/relax/MPFastRelaxMover.hh
///
/// @brief      Basic FastRelax Protocol for Membrane Proteins
/// @details    Refinement and minimization of membrane protein structures using an adapted
///				version of the FastRelax protocol. Uses the membrane framework and adapted
///				minimization settings. 			
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/20/14)

#ifndef INCLUDED_protocols_membrane_relax_MPFastRelaxMover_hh
#define INCLUDED_protocols_membrane_relax_MPFastRelaxMover_hh

// Unit Headers
#include <protocls/membrane/relax/MPFastRleaxMover.hh> 

// Project Headers
#include <protocols/moves/Mover.hh> 

namespace protocols {
namespace membrane {
namespace relax {

class MPFastRleaxMover : public Mover {

public: 

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Defailt Constructor
	/// @details Create a default version of the MPFastRelax Protocol
	MPFastRelaxMover();

	/// @brief Custom Constructor - Custom Protocol Setup
	/// @details Create a custom version of the mp fast relax protocol given
	/// a user provided energy function & number of cycles
	MPFastRelaxMover(
		ScoreFunctionOP sfxn, 
		Size cycles
		);

	/// @brief Copy Constructor
	/// @details Create a deep copy of this mover
	MPFastRelaxMover( MPFastRelaxMover const & src ); 

	/// @brief Assignment Operator
	/// @details Create a deep copy of this mover overloading the assignment operator "="
	MPFastRelaxMover & 
	operator=( MPFastRelaxMover const & src ); 

	/// @brief Destructor
	~MPFastRelaxMover(); 

	//////////////////////
	//// Mover Methods ///
	//////////////////////
	
	/// @brief Get the name of this mover (MPFastRelaxMover)
	virtual std::string get_name() const;
	
	/// @brief Perform Membrane Fast Relax Protocol
	virtual void apply( Pose & pose );
	
	////////////////////////////////
	//// Rosetta Scripts Methods ///
	////////////////////////////////
	
	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;
	
	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;
	
	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
	  utility::tag::TagCOP tag,
	  basic::datacache::DataMap &,
	  protocols::filters::Filters_map const &,
	  protocols::moves::Movers_map const &,
	  core::pose::Pose const &
	  );

public: // methods

	/// @Brief Perform Refinement Cycle
	/// @details Perform a single refinement cycle of reweighting, repack, and minimization
	void perform_refinement_cycle( Pose & pose ); 


private: 

	ScoreFunctionOP sfxn_; 
	PackRotamersMover pack_mover_; 
	MinMoverOP min_mover_; 

};

} // relax
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_relax_MPFastRelaxMover_hh
