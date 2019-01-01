// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSLigandTransformMover.hh
/// @brief   Finds a ligand position and orientation that minimizes experimental PCS data
///          through a grid search. The pose found by PCSLigandTransformMover can be used
///          as starting position for subsequent ligand docking.
/// @details last Modified: 06/02/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pcs_PCSLigandTransformMover_HH
#define INCLUDED_protocols_nmr_pcs_PCSLigandTransformMover_HH

// Unit headers
#include <protocols/nmr/pcs/PCSLigandTransformMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/scoring/nmr/pcs/PCSData.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Numeric headers
#include <numeric/VoxelGrid.impl.hh>
#include <numeric/geometry/BoundingBox.hh>

// C++ headers
#include <iostream>
#include <string>
#include <map>

namespace protocols {
namespace nmr {
namespace pcs {

/// Two utility classes used by the PCSLigandTranformMover
class AtomGridPoint : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor
	AtomGridPoint();

	/// @brief Construct from atom name, ID and 3D coordinates
	AtomGridPoint(
		std::string const & atom,
		core::Size const id,
		core::Vector const & coords,
		core::conformation::Residue const & residue
	);

	/// @brief Destructor
	~AtomGridPoint();

	/// @brief Is relevant for neighbor search
	bool is_relevant_neighbor() const { return true; }

	// Getters
	core::Vector const & get_coordinates() const { return coords_; }
	std::string atom_name() const { return atom_name_; }
	core::Size id() const { return id_; }
	core::conformation::Residue const * residue() const { return residue_; }
	void set_coordinates(core::Vector const & v) { coords_ = v; }

private: // Data

	std::string atom_name_;
	core::Size id_;
	core::Vector coords_;
	core::conformation::Residue const * residue_;

};

class AtomGrid : public numeric::VoxelGrid< AtomGridPoint > {

public: // Methods

	/// @brief default constructor
	AtomGrid(
		Real const & resolution,
		bool const & cache_edges = false
	);

	/// @brief Construct from a vector of points
	AtomGrid(
		Real const & resolution,
		utility::vector1< AtomGridPoint > const & points,
		bool const & cache_edges = false
	);

	/// @brief Construct from a vector of points
	AtomGrid(
		Real const & resolution,
		utility::vector1< AtomGridPoint const * > const & points,
		bool const & cache_edges = false
	);

	/// @brief Destructor
	~AtomGrid() override;

	/// @brief Extract the 3D coordinate of a given object of type VoxelGridPoint.
	Vector const *
	ExtractPosition(
		AtomGridPoint const & point
	) const override;

	/// @brief Check if two grid points are the same.
	bool
	IsSameItem(
		AtomGridPoint const & point1,
		AtomGridPoint const & point2
	) const override;

	/// @brief Check if this item is relevant for neighbor search
	bool
	IsRelevantItem(
		AtomGridPoint const & point
	) const override;

};

class PCSLigandTransformMover : public protocols::moves::Mover {

public: // Typedefs

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::scoring::nmr::Matrix Matrix;
	typedef core::scoring::nmr::pcs::PCSData PCSData;
	typedef core::scoring::nmr::pcs::PCSDataOP PCSDataOP;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::conformation::Residue Residue;
	typedef numeric::geometry::BoundingBox< core::Vector > BoundingBox;
	typedef core::scoring::nmr::Vec5 Vec5;

private:

	/// @brief Utility class of PCSLigandTransformMover which
	///        holds references to PCSData used in lmmin
	///        (Levenberg Marquardt minimization) function
	class LMMinPCSDataRef {

	public:

		PCSData const * pcs_data_;
		std::map<Size,Vector> const * atom_xyz_map_;
		Vector const * xyz_start_;
		Real xyz_delta_;

	public:

		LMMinPCSDataRef() :
			pcs_data_(nullptr),
			atom_xyz_map_(nullptr),
			xyz_start_(nullptr),
			xyz_delta_()
		{}

		LMMinPCSDataRef(
			PCSData const * data,
			std::map<Size,Vector> const * map,
			Vector const * start,
			Real delta

		) :
			pcs_data_(data),
			atom_xyz_map_(map),
			xyz_start_(start),
			xyz_delta_(delta)
		{}

	};

public: // Methods

	/// @brief Default constructor
	PCSLigandTransformMover();

	// @brief Construct from PCSData
	PCSLigandTransformMover(
		PCSDataOP data,
		ScoreFunctionOP sfxn
	);

	/// @brief Copy constructor
	PCSLigandTransformMover(PCSLigandTransformMover const & other);

	/// @brief Assignment operator
	PCSLigandTransformMover & operator=(PCSLigandTransformMover const & rhs);

	/// @brief Destructor
	~PCSLigandTransformMover() override;

	/// @brief Get the name of this mover
	std::string get_name() const override;

	static std::string mover_name();

	/// @brief Perform grid search of ligand position by minimizing the PCS score
	void apply(Pose & pose) override;

	/// @brief Return a clone of the Mover object.
	protocols::moves::MoverOP clone() const override;

	/// @brief Generates a new Mover object freshly created with the default ctor.
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Parse tags of XML script
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose
	) override;

	/// @brief Create XML schema definition for PREMover
	static void provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd);

	// Getters & Setters
	PCSData const & get_pcs_data() const { return *pcs_data_; }
	ScoreFunctionOP get_scorefunction() { return sfxn_; }
	BoundingBox const & get_gridsearch_range() const { return gs_box_; }
	Real get_trans_step() const { return trans_step_; }
	Real get_rot_step() const { return rot_step_; }
	Real get_resolution_damping() const { return rsol_damping_; }
	bool optimized_transform() const { return optimize_transform_; }
	char get_chain() const { return chain_; }

	void set_scorefunction(ScoreFunctionOP sfxn) { sfxn_ = sfxn; }
	void set_gridsearch_range(BoundingBox const & box) { gs_box_ = box; }
	void set_gridsearch_range(Vector const & lower, Vector const & upper) {
		runtime_assert(lower.x() <= upper.x());
		runtime_assert(lower.y() <= upper.y());
		runtime_assert(lower.z() <= upper.z());
		gs_box_.set_lower(lower);
		gs_box_.set_upper(upper);
	}
	void set_trans_step(Real step) { trans_step_ = step; }
	void set_rot_step(Real step) { rot_step_ = step; }
	void set_resolution_damping(Real damping) {
		runtime_assert(damping >= 0.0);
		runtime_assert(damping <= 1.0);
		rsol_damping_ = damping;
	}
	void do_optimized_transform() { optimize_transform_ = true; }
	void undo_optimized_transform() { optimize_transform_ = false; }
	void set_chain(char const c) { chain_ = c; }

private: // Methods

	/// @brief Build grid around protein and fill it with non-ligand atoms
	///        Update the bounding box of the xyz grid search to span the
	///        range: [ min_xyz atom_grid - nbr_radius, max_xyz atom_grid + nbr_radius ]
	void reset_grid_and_bounding_box(Pose const & pose, Size const & ligand_resid);

	/// @brief The local frame of the ligand is located in the ligand's
	///        center of mass and has unit basis
	Vector define_ligand_frame_origin(Residue const & ligand);

	/// @brief Traverse ligand center to next xyz point (in Angstrom)
	bool next_ligand_position(Vector & xyz);

	/// @brief Change ligand orientation to new alpha, beta, gamma angles (in degree)
	bool next_ligand_orientation(Vector & abc);

	/// @brief Perform rotation of vector v by matrix R in local frame centered
	///        at point ori and with unit basis
	Vector local_rotation(Matrix const & R, Vector const & ori, Vector const & v);

	/// @brief Build map of atom index and xyz coordinates for ligand residue
	std::map<Size,Vector> build_ligand_atom_xyz_table(Residue const & ligand);

	/// @brief Find the best position and orientation for this residue according
	///        to the PCS data using a grid search.
	void find_best_ligand_pose_with_grid_search(Size const & resid, Residue const & ligand, Vector & position, Vector & orientation);

	/// @brief Optimize the ligand position and orientation by NLS
	void optimize_ligand_pose_with_nls(Residue const & ligand, Vector & position, Vector & orientation);

	/// @brief pcs error function used in the lmmin function
	///        * par is an array of the positional parameters to be fitted [x, y, z, alpha, beta, gamma]
	///        * data is a pointer to the LMMinPCSDataRef container reference
	///        * fvc is an array holding the residuals of the fit calculation
	friend void pcs_lmmin(Real const *par, int /*m_dat*/, void const *data, Real *fvec, int * /*info*/);

	/// @brief Move ligand residue closer to protein surface. To this end, the ligand
	///        is simply moved along the jump with the protein in 0.1 Ang steps until
	///        an increase in the fa_rep score is monitored.
	void move_ligand_close_to_surface(Pose & pose, Size const & ligand_resid);

	/// @brief Calculate the damped resolution in the range [rmax, rmin] using the damping parameter
	Real damped_resolution(Real rmax, Real rmin);

private: // Data

	/// @brief experimental PCS data for all lanthanides and tagging sites
	PCSDataOP pcs_data_;

	/// @brief a scorefunction to monitor if ligand touches the protein
	ScoreFunctionOP sfxn_;

	/// @brief voxelgrid that spans around the protein and which is used for clash check
	AtomGrid grid_;

	/// @brief Bounding box of the ligand xyz position grid search
	BoundingBox gs_box_;

	/// @brief translational stepsize (in Ang.) of grid search
	Real trans_step_;

	/// @brief rotational stepsize (in degree) of grid search
	Real rot_step_;

	/// @brief smoothening parameter that can be used to damp the resolution
	///        of the atom clash check with the voxel grid. Takes values between 0.0 and 1.0
	Real rsol_damping_;

	/// @brief do further optimization of position and orientation as by grid search
	bool optimize_transform_;

	/// @brief ligand chain
	char chain_;

};

} // namespace pcs
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_pcs_PCSLigandTransformMover_HH
