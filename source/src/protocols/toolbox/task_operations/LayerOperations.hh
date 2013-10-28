// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/LayerOperations.hh
/// @brief select residues depending on layer: core, boundary, and surface
/// and the layer are defined by accessible surface area of each residue, like nobu's just more general
/// @author Eva

#ifndef INCLUDED_protocols_toolbox_LayerOperations_hh
#define INCLUDED_protocols_toolbox_LayerOperations_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/AtomType.hh>

#include <core/scoring/sasa.hh>
//#include <protocols/toolbox/LayerOperations.fwd.hh>
#include <map>
#include <utility/vector1.hh>

//package headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/types.hh>


// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace toolbox {
namespace task_operations{

class LayerOperations : public core::pack::task::operation::TaskOperation {
public:
	typedef core::Size Size;
	//typedef core::Real Real;
	typedef std::string String;
	//typedef core::chemical::AA AA;
	typedef core::pose::Pose Pose;
	typedef core::Real Real;
	//typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;

public: // constructor/destructor

	/// @brief default constructor
	LayerOperations();

	/// @brief value constructor
	LayerOperations( bool const pick_core, bool const pick_boundary, bool const pick_surface );

	/// @brief value constructor
	/// @detail selected layer can be given by string, for example, core_surface or core_boundary_surface
	LayerOperations( String const pick );

	/// @brief destructor
	~LayerOperations();
    
    /// @brief make clone
	virtual TaskOperationOP clone() const;

public:

	void initialize( Real burial, Real surface );

public:  // mutator

	/// @brief
	void
	set_design_layer( bool const pick_core, bool const pick_boundary, bool const pick_surface )
	{
		pick_core_ = pick_core;
		pick_boundary_ = pick_boundary;
		pick_surface_ = pick_surface;
	}

	/// @brief set pore radius for colculating asa
	void pore_radius( Real const ps ){
		pore_radius_ = ps;
	}

	/// @brief accessible surface for evaluating residues are in surface or not
	void sasa_surface( Real const r, String const ss="" );

	/// @brief accessible surface for evaluating residues are in core or not
	void sasa_core( Real const r, String const ss="" );

	/// @brief accessible surface for evaluating residues are in core or not
	void make_rasmol_format_file( bool const b ){
		make_rasmol_format_file_ = b;
	}

	/// @brief
	//void exclude_aatypes_for_selection( utility::vector1< AA > const & aas );

	/// @brief
	//void restrict_aatypes_for_selection( utility::vector1< AA > const & aas );

public:  // accessor

	/// @brief accessbile surface are of each residue
	Real rsd_sasa( Size const i ) const;

	/// @brief return defined layer for each residue
	String layer( Size const i ) const;

	/// @brief selected residues on boundary
	utility::vector1< Size > const & selected_boundary_residues() const;

	/// @brief selected residues in core
	utility::vector1< Size > const & selected_core_residues() const;

	/// @brief selected residues on surface
	utility::vector1< Size > const & selected_surface_residues() const;


public: // main operation


	/// @brief apply
	utility::vector1< Size > const
	compute( Pose const & pose );

private: // helper functions

	/// @brief
	LayerOperations(LayerOperations const &);

	/// @brief
	utility::vector1< Real > const
	calc_rsd_sasa( Pose const & pose ) const;

protected:

	/// @brief design/repack core ?
	bool pick_core_;

	/// @brief design/repack boundary ?
	bool pick_boundary_;

	/// @brief design surface ?
	bool pick_surface_;

	/// @brief pore radius for calculating asa( accessible surface area )
	Real pore_radius_;

	/// @brief asa pore radius for calculating asa
	std::map< char, Real > burial_;

	/// @brief pore radius for calculating asa
	std::map< char, Real > surface_;

	/// @brief amino acid types excluded for selection
	//utility::vector1< AA > excluded_aatypes_for_selection_;

	/// @brief amino acid types restricted for selection
	//utility::vector1< AA > restricted_aatypes_for_selection_;

	/// @brief selected residues in core
	utility::vector1< Size > selected_core_residues_;

  /// @brief selected residues at boundary
	utility::vector1< Size > selected_boundary_residues_;

  /// @brief selected residues on surface
	utility::vector1< Size > selected_surface_residues_;

	/// @brief create rasmol format script for selected residues; red: surface, blue: core, green: boundary
	bool make_rasmol_format_file_;

	/// @brief asa for each residue
	utility::vector1< Real > rsd_sasa_;

	/// @brief
	utility::vector1< String > rsd_layer_;

	/// @brief design selection?
	bool design_;

	/// @brief repack the rest of the protein?
	bool repack_rest_;

	/// @brief which chain to pick
	core::Size chain_;
	
	/// @brief take chains apart, as defined by user specified jump number
	core::Size jump_;

};

}
} 
} // protocols

#endif
