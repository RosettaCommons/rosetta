// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSLigandTransformMover.cc
/// @brief   Implementation of PCSLigandTransformMover
/// @details last Modified: 06/02/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/nmr/pcs/PCSLigandTransformMover.hh>
#include <protocols/nmr/pcs/PCSLigandTransformMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/scoring/nmr/pcs/PCSData.hh>
#include <core/scoring/nmr/pcs/PCSMultiSet.hh>
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/util.hh>

// Core headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/mover_schemas.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <ObjexxFCL/format.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/nls/lmmin.hh>
#include <numeric/random/random.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>
#include <limits>

namespace protocols {
namespace nmr {
namespace pcs {

static basic::Tracer TR( "protocols.nmr.pcs.PCSLigandTransformMover" );

AtomGridPoint::AtomGridPoint() :
	atom_name_(""),
	id_(),
	coords_(),
	residue_(nullptr)
{}

AtomGridPoint::AtomGridPoint(
	std::string const & atom,
	core::Size const id,
	core::Vector const & coords,
	core::conformation::Residue const & residue
) :
	atom_name_(atom),
	id_(id),
	coords_(coords),
	residue_(&residue)
{}

AtomGridPoint::~AtomGridPoint() {}


AtomGrid::AtomGrid(
	Real const & resolution,
	bool const & cache_edges
) :
	numeric::VoxelGrid< AtomGridPoint >(resolution, cache_edges)
{}

AtomGrid::AtomGrid(
	Real const & resolution,
	utility::vector1< AtomGridPoint > const & points,
	bool const & cache_edges
) :
	numeric::VoxelGrid< AtomGridPoint >(resolution, cache_edges)
{
	utility::vector1< AtomGridPoint const * > points_ptr;
	points_ptr.reserve(points.size());
	for ( auto & p : points ) {
		points_ptr.push_back(&p);
	}
	SetObjects(points_ptr);
}

AtomGrid::AtomGrid(
	Real const & resolution,
	utility::vector1< AtomGridPoint const * > const & points,
	bool const & cache_edges
) :
	numeric::VoxelGrid< AtomGridPoint >(resolution, cache_edges)
{
	SetObjects(points);
}

AtomGrid::~AtomGrid() {}

core::Vector const *
AtomGrid::ExtractPosition(
	AtomGridPoint const & point
) const
{
	return &(point.get_coordinates());
}

bool
AtomGrid::IsSameItem(
	AtomGridPoint const & point1,
	AtomGridPoint const & point2
) const
{
	return point1.id() == point2.id()
		&& point1.atom_name() == point2.atom_name()
		&& point1.residue() == point2.residue();
}

bool
AtomGrid::IsRelevantItem(
	AtomGridPoint const & point
) const
{
	return point.is_relevant_neighbor();
}

PCSLigandTransformMover::PCSLigandTransformMover() :
	protocols::moves::Mover(),
	pcs_data_(nullptr),
	sfxn_(nullptr),
	grid_(10.0, false),
	gs_box_(),
	trans_step_(5.0),
	rot_step_(20.0),
	rsol_damping_(0.0),
	optimize_transform_(false),
	chain_('X')
{}

PCSLigandTransformMover::PCSLigandTransformMover(
	PCSDataOP data,
	ScoreFunctionOP sfxn
) :
	protocols::moves::Mover(),
	pcs_data_(data),
	sfxn_(sfxn),
	grid_(10.0, false),
	gs_box_(),
	trans_step_(5.0),
	rot_step_(20.0),
	rsol_damping_(0.0),
	optimize_transform_(false),
	chain_('X')
{}

PCSLigandTransformMover::PCSLigandTransformMover(PCSLigandTransformMover const & other) :
	protocols::moves::Mover(other),
	pcs_data_( new PCSData( *(other.pcs_data_) ) ),
	sfxn_( other.sfxn_->clone() ),
	grid_(other.grid_),
	gs_box_(other.gs_box_),
	trans_step_(other.trans_step_),
	rot_step_(other.rot_step_),
	rsol_damping_(other.rsol_damping_),
	optimize_transform_(other.optimize_transform_),
	chain_(other.chain_)
{}

PCSLigandTransformMover &
PCSLigandTransformMover::operator=(PCSLigandTransformMover const & rhs) {
	if ( this == &rhs ) {
		return *this;
	}
	return *( new PCSLigandTransformMover( *this ) );
}

PCSLigandTransformMover::~PCSLigandTransformMover() {}

std::string
PCSLigandTransformMover::mover_name() {
	return "PCSLigandTransformMover";
}

std::string
PCSLigandTransformMover::get_name() const {
	return mover_name();
}

std::string
PCSLigandTransformMoverCreator::keyname() const {
	return PCSLigandTransformMover::mover_name();
}

protocols::moves::MoverOP
PCSLigandTransformMover::clone() const {
	return protocols::moves::MoverOP(new PCSLigandTransformMover(*this));
}

protocols::moves::MoverOP
PCSLigandTransformMover::fresh_instance() const {
	return protocols::moves::MoverOP(new PCSLigandTransformMover);
}

protocols::moves::MoverOP
PCSLigandTransformMoverCreator::create_mover() const {
	return protocols::moves::MoverOP(new PCSLigandTransformMover);
}

void
PCSLigandTransformMover::apply(Pose & pose) {

	using namespace core::pose;
	using namespace core::conformation;

	Size lig_resid = get_resnums_for_chain(pose, chain_)[1];
	Residue const & lig = pose.residue(lig_resid);

	// Create grid around protein and fill it with atoms
	// Grid resolution should be the ligand's neighbor radius
	// Additionally, set up the dimensions of the grid search
	TR.Info << "Setting up AtomGrid around protein" << std::endl;
	reset_grid_and_bounding_box(pose, lig_resid);

	// Find the best ligand position and orientation given the
	// PCS data and using a grid search
	Vector lig_frame_ori = define_ligand_frame_origin(lig);
	Vector best_gs_xyz, best_gs_abc;
	TR.Info << "Running PCS grid search with " << trans_step_ << " (Ang) and " << rot_step_ << " (deg) steps" << std::endl;
	find_best_ligand_pose_with_grid_search(lig_resid, lig, best_gs_xyz, best_gs_abc);

	// Optionally, do further optimization
	if ( optimize_transform_ ) {
		TR.Info << "Optimizing ligand transformation" << std::endl;
		optimize_ligand_pose_with_nls(lig, best_gs_xyz, best_gs_abc);
	}

	// Translation to best grid point position
	Vector vt = best_gs_xyz - lig_frame_ori;
	// Rotation to best orientation
	Matrix Rot = numeric::rotation_matrix_from_euler_angles_ZXZ(best_gs_abc);

	// Update ligand residue coordinates in pose
	for ( Size i(1); i <= lig.natoms(); ++i ) {
		Vector atom_xyz_transformed = local_rotation(Rot, lig_frame_ori, lig.xyz(i));
		atom_xyz_transformed += vt;
		pose.set_xyz(core::id::AtomID(i,lig_resid), atom_xyz_transformed);
	}

	// Mover ligand closer to protein surface until it touches
	move_ligand_close_to_surface(pose, lig_resid);
}

void
PCSLigandTransformMover::reset_grid_and_bounding_box(
	Pose const & pose,
	Size const & ligand_resid
)
{
	using namespace core::pose;
	using namespace core::conformation;

	Real resolution = pose.residue(ligand_resid).type().nbr_radius();
	Real damped_rsl = damped_resolution(resolution, 0.2*resolution);
	grid_.Clear();

	utility::vector1< AtomGridPoint > grid_points;
	grid_points.reserve(pose.total_atoms());

	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		if ( pose.residue(i).is_virtual_residue() || i == ligand_resid ) {
			continue;
		}
		for ( Size j(1); j <= pose.residue(i).natoms(); ++j ) {
			if ( !pose.residue(i).atom_is_hydrogen(j) ) {
				AtomGridPoint pt(pose.residue(i).atom_name(j), i, pose.residue(i).xyz(j), pose.residue(i));
				grid_points.push_back(pt);
			}
		}
	}
	grid_ = AtomGrid(damped_rsl, grid_points, false);
	gs_box_.set_lower(grid_.GetBoundingBox().lower() - pose.residue(ligand_resid).type().nbr_radius());
	gs_box_.set_upper(grid_.GetBoundingBox().upper() + pose.residue(ligand_resid).type().nbr_radius());
}

core::Vector
PCSLigandTransformMover::define_ligand_frame_origin(Residue const & ligand) {
	using namespace core::conformation;

	utility::vector1<Vector> lig_coords;
	lig_coords.reserve(ligand.natoms());
	for ( utility::vector1<Atom>::const_iterator it = ligand.atom_begin(), it_end = ligand.atom_end(); it != it_end; ++it ) {
		lig_coords.push_back(it->xyz());
	}
	runtime_assert_msg(lig_coords.size() > 0, "Number of ligand atom coordinates must be > 0.");

	return numeric::center_of_mass(lig_coords);
}

bool
PCSLigandTransformMover::next_ligand_position(Vector & xyz)
{
	xyz.x() += trans_step_;

	if ( xyz.x() - gs_box_.upper().x() > trans_step_/2.0 ) {
		xyz.x() = gs_box_.lower().x();
		xyz.y() += trans_step_;
	}
	if ( xyz.y() - gs_box_.upper().y() > trans_step_/2.0 ) {
		xyz.y() = gs_box_.lower().y();
		xyz.z() += trans_step_;
	}
	if ( xyz.z() - gs_box_.upper().z() > trans_step_/2.0 ) {
		xyz.z() = gs_box_.lower().z();
		return false;
	}
	return true;
}

bool
PCSLigandTransformMover::next_ligand_orientation(Vector & abc)
{
	abc.x() += rot_step_;

	if ( abc.x() > 360.0 ) {
		abc.x() = 0.0;
		abc.y() += rot_step_;
	}
	if ( abc.y() > 180.0 ) {
		abc.y() = 0.0;
		abc.z() += rot_step_;
	}
	if ( abc.z() > 360.0 ) {
		abc.z() = 0.0;
		return false;
	}
	return true;
}

core::Vector
PCSLigandTransformMover::local_rotation(
	Matrix const & R,
	Vector const & ori,
	Vector const & v)
{
	// 1) Transform coordinates into local frame
	Vector v_t = v - ori;
	// 2) Apply rotation in local frame
	v_t = R * v_t;
	// 3) Transform coordinates back into global frame
	v_t = v_t + ori;
	return v_t;
}

std::map<core::Size,core::Vector>
PCSLigandTransformMover::build_ligand_atom_xyz_table(Residue const & ligand)
{
	std::map<Size,Vector> table;
	for ( Size i(1); i <= ligand.natoms(); ++i ) {
		table[i] = ligand.xyz(i);
	}
	return table;
}

core::Real
PCSLigandTransformMover::damped_resolution(
	Real rmax,
	Real rmin
)
{
	return (1.0-rsol_damping_)*(1.0-rsol_damping_)*(rmax-rmin)+rmin;
}

void
PCSLigandTransformMover::find_best_ligand_pose_with_grid_search(
	Size const & resid,
	Residue const & ligand,
	Vector & position,
	Vector & orientation
)
{
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pcs;

	// Define local ligand frame centered in its centroid
	Vector lig_frame_ori = define_ligand_frame_origin(ligand);

	// Create a test grid point for the ligand residue for neighbor search
	AtomGridPoint lig_grid_pt("ORI", resid, lig_frame_ori, ligand);

	// Closure to calculate translation vector between ligand center and grid point p
	// In addition, it sets our test grid point to p
	auto translation_vector = [&] (Vector const & p) {
		lig_grid_pt.set_coordinates(p);
		return p - lig_frame_ori;
	};

	// Shift ligand to the lower left front corner of the grid search bounding box
	Vector vt = translation_vector(gs_box_.lower());

	// Starting position and orientation of ligand
	Vector current_xyz(gs_box_.lower());
	Vector current_abc(0.0,0.0,0.0);
	TR.Trace << "Ligand start position before PCS grid search " << current_xyz.to_string() << " (Ang, Ang, Ang)" << std::endl;
	TR.Trace << "Ligand start orientation before PCS grid search " << current_abc.to_string() << " (deg, deg, deg)" << std::endl;
	Matrix Rot;
	Real search_radius = damped_resolution(ligand.type().nbr_radius(), 0.2*ligand.type().nbr_radius());

	// PCS data
	utility::vector1<PCSMultiSetOP> const & pcs_multiset_vec = pcs_data_->get_pcs_multiset_vec();
	Size n_pcs_tags = pcs_data_->get_number_tags();

	Real best_score = std::numeric_limits< Real >::max(), current_score = std::numeric_limits< Real >::max();
	Vector best_xyz = current_xyz, best_abc = current_abc;
	bool skip_data = false;

	while ( next_ligand_position(current_xyz) ) {

		// Shift ligand center to new grid point
		vt = translation_vector(current_xyz);

		// Continue only if ligand center has no neighboring protein atoms with search_radius
		if ( !grid_.HasNeighbors(lig_grid_pt,search_radius) ) {

			while ( next_ligand_orientation(current_abc) ) {

				// Rotation matrix for current set of euler angles
				Rot = numeric::rotation_matrix_from_euler_angles_ZXZ(current_abc);

				// Calculate PCS score for this ligand position and orientation
				current_score = 0.0;
				skip_data = false;

				// Number of tags
				for ( Size i(1); i <= n_pcs_tags; ++i ) {

					if ( skip_data ) { // If current score is already higher than best score, we skip the next datasets
						break;
					}

					Size n_metals = pcs_multiset_vec[i]->get_number_metal_ions();
					utility::vector1<PCSSingleSetOP> const & pcs_singleset_vec = pcs_multiset_vec[i]->get_pcs_singleset_vec();

					// Number of metals (i.e. PCS experiments) for this tag
					for ( Size j(1); j <= n_metals; ++j ) {

						Vector tensor_angles( pcs_singleset_vec[j]->get_tensor_const()->get_alpha(),
							pcs_singleset_vec[j]->get_tensor_const()->get_beta(),
							pcs_singleset_vec[j]->get_tensor_const()->get_gamma() );
						Matrix tensor_rotation = rotation_matrix_from_euler_angles(tensor_angles,
							pcs_singleset_vec[j]->get_tensor_const()->get_euler_convention());
						Vec5 tensor_params = { pcs_singleset_vec[j]->get_tensor_const()->get_metal_center().x(),
							pcs_singleset_vec[j]->get_tensor_const()->get_metal_center().y(),
							pcs_singleset_vec[j]->get_tensor_const()->get_metal_center().z(),
							pcs_singleset_vec[j]->get_tensor_const()->get_ax(),
							pcs_singleset_vec[j]->get_tensor_const()->get_rh() };
						Real singleset_score(0.0);
						utility::vector1<PCSSingle> const & pcs_single_vec = pcs_singleset_vec[j]->get_single_pcs_vec();
						// PCS values for this dataset
						for ( utility::vector1<PCSSingle>::const_iterator pcs_it = pcs_single_vec.begin();
								pcs_it != pcs_single_vec.end(); ++pcs_it ) {

							Real pcs_calc(0.0);
							// Equivalent spins that contribute to this PCS
							for ( utility::vector1<core::id::AtomID>::const_iterator atom_it = pcs_it->get_protein_spins().begin();
									atom_it != pcs_it->get_protein_spins().end(); ++atom_it ) {

								if ( atom_it->rsd() != resid || atom_it->atomno() > ligand.natoms() ) {
									utility_exit_with_message("Could not find AtomID in ligand residue.");
								}
								// Calculate transformed coordinates
								Vector atom_xyz_transformed = local_rotation(Rot, lig_frame_ori, ligand.xyz(atom_it->atomno()));
								atom_xyz_transformed += vt;

								pcs_calc += pcs_func(tensor_params, tensor_rotation, atom_xyz_transformed) ;
							} // One single PCS value
							pcs_calc /= pcs_it->get_protein_spins().size();
							Real dev = pcs_calc - pcs_it->get_pcs_exp();
							singleset_score += dev * dev * pcs_it->get_weight();

						} // PCS values for this dataset
						current_score += std::sqrt(singleset_score) * pcs_singleset_vec[j]->get_weight();
						if ( current_score > best_score ) {
							skip_data = true;
							continue;
						}

					} // Number of metals (i.e. PCS experiments) for this tag

				} // Number of tags
				if ( current_score < best_score ) {
					best_score = current_score;
					best_xyz = current_xyz;
					best_abc = current_abc;
				}
			}
		}
	}

	TR.Trace << "Best ligand position after PCS grid search " << best_xyz.to_string() << " (Ang, Ang, Ang)" << std::endl;
	TR.Trace << "Ligand start orientation after PCS grid search " << best_abc.to_string() << " (deg, deg, deg)" << std::endl;
	position = best_xyz;
	orientation = best_abc;
}

void
PCSLigandTransformMover::move_ligand_close_to_surface(
	Pose & pose,
	Size const & ligand_resid
)
{
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::scoring;

	utility::vector1<bool> chain_X_resid(pose.total_residue(), false);
	runtime_assert_msg(ligand_resid <= chain_X_resid.size(), "Ligand residue ID outside of pose.total_residue()");
	chain_X_resid[ligand_resid] = true;

	Size jump_id = jump_which_partitions(pose.fold_tree(), chain_X_resid);
	if ( jump_id == 0 ) {
		TR.Warning << "No jump ID for ligand residue found while trying to move ligand closer to protein surface. Trying to reset fold tree." << std::endl;
		FoldTree new_fold_tree(get_foldtree_which_partitions(pose.fold_tree(), chain_X_resid));
		pose.fold_tree(new_fold_tree);
		jump_id = jump_which_partitions(pose.fold_tree(),chain_X_resid);
		runtime_assert(jump_id != 0);
	}

	(*sfxn_)(pose);
	Real initial_fa_rep = pose.energies().total_energies()[fa_rep];

	// Mover towards protein to detect of ligand is touching protein surface
	protocols::rigid::RigidBodyTransMover trans_mover(pose, jump_id);
	trans_mover.trans_axis(trans_mover.trans_axis().negate());
	trans_mover.step_size(1);
	trans_mover.apply(pose);
	(*sfxn_)(pose);
	Real new_fa_rep = pose.energies().total_energies()[fa_rep];
	Real fa_rep_delta = initial_fa_rep - new_fa_rep;
	if ( std::abs(fa_rep_delta) < 1.0e-4 ) { // Move closer to surface
		TR.Warning << "Ligand is not touching the protein. Move ligand closer to protein." << std::endl;
		trans_mover.step_size(0.1);
		do {
			trans_mover.apply(pose);
			(*sfxn_)(pose);
			new_fa_rep = pose.energies().total_energies()[fa_rep];
			fa_rep_delta = initial_fa_rep - new_fa_rep;
		} while ( std::abs(fa_rep_delta) < 1.0e-4 );
	} else {        // Move back
		trans_mover.trans_axis(trans_mover.trans_axis().negate());
		trans_mover.apply(pose);
	}
}

void
pcs_lmmin(
	core::Real const *par,
	int /*m_dat*/,
	void const *data,
	core::Real *fvec,
	int * /*info*/
)
{
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pcs;

	PCSLigandTransformMover::LMMinPCSDataRef * min_data = static_cast< PCSLigandTransformMover::LMMinPCSDataRef * >(const_cast< void* >(data));

	core::Real * nonconst_par = const_cast< core::Real* >(par);

	// Constrain the fit parameter to reasonable ranges.
	if ( !( (min_data->xyz_start_->x() - min_data->xyz_delta_) <= par[0] ) ||
			!( (min_data->xyz_start_->x() + min_data->xyz_delta_) >= par[0] ) ) {
		nonconst_par[0] = min_data->xyz_delta_*std::tanh(par[0]) + min_data->xyz_start_->x();
	}
	if ( !( (min_data->xyz_start_->y() - min_data->xyz_delta_) <= par[1] ) ||
			!( (min_data->xyz_start_->y() + min_data->xyz_delta_) >= par[1] ) ) {
		nonconst_par[1] = min_data->xyz_delta_*std::tanh(par[1]) + min_data->xyz_start_->y();
	}
	if ( !( (min_data->xyz_start_->z() - min_data->xyz_delta_) <= par[2] ) ||
			!( (min_data->xyz_start_->z() + min_data->xyz_delta_) >= par[2] ) ) {
		nonconst_par[2] = min_data->xyz_delta_*std::tanh(par[2]) + min_data->xyz_start_->z();
	}
	if ( !( 0.0 <= par[3] ) || !( par[3] < 360.0 ) ) {
		nonconst_par[3] = 180.0*std::tanh(par[3])+180.0;
	}
	if ( !( 0.0 <= par[4] ) || !( par[4] < 180.0 ) ) {
		nonconst_par[4] = 90.0*std::tanh(par[4])+90.0;
	}
	if ( !( 0.0 <= par[5] ) || !( par[5] < 360.0 ) ) {
		nonconst_par[5] = 180.0*std::tanh(par[5])+180.0;
	}

	core::Vector local_ori(0.0,0.0,0.0);
	for ( auto const & atm : *(min_data->atom_xyz_map_) ) {
		local_ori += atm.second;
	}
	local_ori /= min_data->atom_xyz_map_->size();

	// Current set of transformations
	core::Vector current_position(nonconst_par[0], nonconst_par[1], nonconst_par[2]);
	core::Vector current_orientation(nonconst_par[3], nonconst_par[4], nonconst_par[5]);
	Matrix Rot(numeric::rotation_matrix_from_euler_angles_ZXZ(current_orientation));
	core::Vector vt(current_position - local_ori);

	int ii(0);
	utility::vector1<PCSMultiSetOP> const & pcs_multiset_vec = min_data->pcs_data_->get_pcs_multiset_vec();
	for ( utility::vector1<PCSMultiSetOP>::const_iterator multiset_it = pcs_multiset_vec.begin();
			multiset_it != pcs_multiset_vec.end(); ++multiset_it ) {

		utility::vector1<PCSSingleSetOP> const & pcs_singleset_vec = (*multiset_it)->get_pcs_singleset_vec();
		for ( utility::vector1<PCSSingleSetOP>::const_iterator singleset_it = pcs_singleset_vec.begin();
				singleset_it != pcs_singleset_vec.end(); ++singleset_it ) {

			core::Vector tensor_angles( (*singleset_it)->get_tensor_const()->get_alpha(),
				(*singleset_it)->get_tensor_const()->get_beta(),
				(*singleset_it)->get_tensor_const()->get_gamma() );
			Matrix tensor_rotation = rotation_matrix_from_euler_angles(tensor_angles,
				(*singleset_it)->get_tensor_const()->get_euler_convention());
			Vec5 tensor_params = { (*singleset_it)->get_tensor_const()->get_metal_center().x(),
				(*singleset_it)->get_tensor_const()->get_metal_center().y(),
				(*singleset_it)->get_tensor_const()->get_metal_center().z(),
				(*singleset_it)->get_tensor_const()->get_ax(),
				(*singleset_it)->get_tensor_const()->get_rh() };

			utility::vector1<PCSSingle> const & pcs_single_vec = (*singleset_it)->get_single_pcs_vec();
			for ( utility::vector1<PCSSingle>::const_iterator pcs_it = pcs_single_vec.begin();
					pcs_it != pcs_single_vec.end(); ++pcs_it ) {

				core::Real pcs_calc(0.0);
				for ( utility::vector1<core::id::AtomID>::const_iterator atom_it = pcs_it->get_protein_spins().begin();
						atom_it != pcs_it->get_protein_spins().end(); ++atom_it ) {

					// Calculate transformed coordinates
					core::Vector atom_xyz_transformed = min_data->atom_xyz_map_->at(atom_it->atomno()) - local_ori;
					atom_xyz_transformed = Rot * atom_xyz_transformed;
					atom_xyz_transformed = atom_xyz_transformed + local_ori + vt;

					pcs_calc += pcs_func(tensor_params, tensor_rotation, atom_xyz_transformed) ;
				}
				pcs_calc /= pcs_it->get_protein_spins().size();
				core::Real dev = pcs_it->get_pcs_exp() - pcs_calc;

				fvec[ii++] = dev * pcs_it->get_weight();

			} // one single PCS

		} // one PCS experiment

	} // all PCS tags
}

void
PCSLigandTransformMover::optimize_ligand_pose_with_nls(
	Residue const & ligand,
	Vector & position,
	Vector & orientation
)
{
	TR.Trace << "Ligand position before optimization " << position.to_string() << " (Ang, Ang, Ang)" << std::endl;
	TR.Trace << "Ligand orientation before optimization " << orientation.to_string() << " (deg, deg, deg)" << std::endl;

	utility::fixedsizearray1<Real,6> fit_params = { position.x(), position.y(), position.z(), orientation.x(), orientation.y(), orientation.z() };
	utility::fixedsizearray1<Real,6> best_params(fit_params);
	numeric::nls::lm_status_struct status;
	Real bestnorm = std::numeric_limits< Real >::infinity();
	int no_pcs = static_cast<int>(pcs_data_->get_total_number_pcs());
	Size no_repeats(10);
	std::map<Size,Vector> atom_xyz_table = build_ligand_atom_xyz_table(ligand);
	LMMinPCSDataRef mindata(pcs_data_.get(), &atom_xyz_table, &position, trans_step_);

	for ( Size i(1); i <= no_repeats; ++i ) {
		numeric::nls::lmmin( 6, &fit_params[1], no_pcs, (const void*) &mindata, pcs_lmmin, &status, numeric::nls::lm_printout_std);
		if ( status.fnorm < bestnorm ) {
			bestnorm=status.fnorm;
			best_params = fit_params;
		}
		fit_params[4] = 360.0 * numeric::random::uniform();
		fit_params[5] = 180.0 * numeric::random::uniform();
		fit_params[6] = 360.0 * numeric::random::uniform();
	}

	position.x() = best_params[1];
	position.y() = best_params[2];
	position.z() = best_params[3];
	orientation.x() = best_params[4];
	orientation.y() = best_params[5];
	orientation.z() = best_params[6];

	TR.Trace << "Ligand position after optimization " << position.to_string() << " (Ang, Ang, Ang)" << std::endl;
	TR.Trace << "Ligand orientation after optimization " << orientation.to_string() << " (deg, deg, deg)" << std::endl;
}

/// @brief Parse tags of XML script
void
PCSLigandTransformMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose
)
{
	try {
		// PCSData
		if ( tag->hasOption("pcs_input_file") ) {
			std::string filename;
			filename = tag->getOption< std::string >( "pcs_input_file", "" );
			if ( filename == "" ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'PCSLigandTransformMover' requires 'pcs_input_file' tag.");
			} else {
				pcs_data_ = PCSDataOP(new PCSData(filename, pose));
			}
		}
		// ScoreFunction
		if ( tag->hasOption("scorefxn") ) {
			sfxn_ = protocols::rosetta_scripts::parse_score_function(tag, "scorefxn", datamap);
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'PCSLigandTransformMover' requires 'scorefxn' tag");
		}
		// Chain
		if ( tag->hasOption("chain") ) {
			std::string chain = tag->getOption< std::string >("chain", "X");
			if ( chain == "" ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'PCSLigandTransformMover' requires 'chain' tag.");
			} else {
				chain_ = char(chain[0]);
			}
		}
	} catch ( utility::excn::RosettaScriptsOptionError & excn ) {
		TR << "caught exception " << excn.msg() << std::endl;
	}
}

/// @brief Create XML schema definition for PCSLigandTransformMover
void
PCSLigandTransformMover::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd) {
	using namespace utility::tag;

	// Basic attributes
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "pcs_input_file", xs_string, "PCS data input file" )
		+ XMLSchemaAttribute::required_attribute( "scorefxn", xs_string, "Scorefunction to be used by PCSLigandTransformMover" )
		+ XMLSchemaAttribute::attribute_w_default( "chain", xsct_char, "Ligand chain", "X" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"This Mover performs a PCS data driven grid search of the ligand position and orientation", attlist);
}

void
PCSLigandTransformMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PCSLigandTransformMover::provide_xml_schema( xsd );
}

} // namespace pcs
} // namespace nmr
} // namespace protocols
