// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/TransformInfoMembraneMoverCreator.hh
/// @brief      Transform pose into membrane coordinates (Rosetta Scripts Hook)
/// @details	Requires a MembraneInfo object with all of its associated
///				information. This can be done by calling AddMembraneMover
///				beforehand. pose.conformation().is_membrane() should always
///				return true.
///				CAUTION: THIS MOVER ONLY WORKS FOR A FIXED MEMBRANE WHERE THE
///				MEMBRANE VIRTUAL RESIDUE IS AT THE ROOT OF THE FOLDTREE!!!
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TransformIntoMembraneMover_hh
#define INCLUDED_protocols_membrane_TransformIntoMembraneMover_hh

// Unit Headers
#include <protocols/membrane/TransformIntoMembraneMover.fwd.hh>
#include <protocols/membrane/TransformIntoMembraneMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/TranslationRotationMover.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/MembraneInfo.fwd.hh>
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/Span.fwd.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>
#include <protocols/membrane/geometry/EmbeddingDef.fwd.hh>
#include <protocols/membrane/geometry/util.hh>

// Package Headers
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh> 
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;
	  
/// @brief	Takes a pose, calculates the best possible center and normal from
///			the topology object (through spanfile from OCTOPUS) and the
///			coordinates of the pose
///			CAUTION: THIS MOVER ONLY WORKS FOR A FIXED MEMBRANE WHERE THE
///			MEMBRANE VIRTUAL RESIDUE IS AT THE ROOT OF THE FOLDTREE!!!
class TransformIntoMembraneMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details uses membrane coords: center(0,0,0), normal(0,0,15) and membrane jump
	TransformIntoMembraneMover();

	/// @brief Constructor that uses jump number for transformation
	/// @details uses membrane coords: center(0,0,0), normal(0,0,15);
	///				downstream jump will be transformed
	TransformIntoMembraneMover( SSize jump );

	/// @brief Custom Constructor - mainly for PyRosetta
	/// @details user can specify membrane coords into which protein should be
	///          transformed to; uses membrane jump
	TransformIntoMembraneMover(
		Vector mem_center,
		Vector mem_normal,
		std::string spanfile = ""
	);

	/// @brief Custom Constructor using jumpnumber
	/// @details user can specify membrane coords into which protein should be
	///          transformed to; downstream partner will be transformed
	TransformIntoMembraneMover( SSize jump,
							   Vector mem_center,
							   Vector mem_normal,
							   std::string spanfile = ""
							   );

	/// @brief Copy Constructor
	TransformIntoMembraneMover( TransformIntoMembraneMover const & src );
	
	/// @brief Destructor
	virtual ~TransformIntoMembraneMover();
	
	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////
	
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
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////
	
	/// @brief Get the name of this Mover (TransformIntoMembraneMover)
	virtual std::string get_name() const;
		
	/// @brief Move the pose into membrane coordinate frame
	virtual void apply( Pose & pose );
	
private: // methods
	
	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	/// @details Register mover-relevant options with JD2 - includes
	/// mp, setup options: center, normal, spanfiles, jumpnum
	void register_options();
	
	/// @brief Initialize Mover options from the commandline
	/// @details Initialize mover settings from the commandline
	/// mainly in the mp, setup group: center, normal,
	/// spanfiles, jumpnum
	void init_from_cmd();
	
	
private: // data

	// Pose residue typeset, and visualization
	bool fullatom_;

	// jump number: downstream partner will be transformed
	SSize jump_;
	
	// center and normal coordinates with respect to fixed membrane
	Vector mem_center_;
	Vector mem_normal_;

	// spanfile name
	std::string spanfile_;

	// SpanningTopology
	SpanningTopologyOP topology_;
	
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TransformIntoMembraneMover_hh
