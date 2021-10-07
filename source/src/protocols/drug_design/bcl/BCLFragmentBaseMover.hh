// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

////////////////////////////////////////////////////////////////////////////////
/// @file
///
/// @brief
/// A base mover class for perturbing small molecule fragments with the BCL
///
/// @details
/// This is a base class for perturbing small molecules with the BCL. It
/// primarily contains options relating to generating 3D conformers of fragments,
/// which is very common for most perturbation types. It also contains a
/// BCLFragmentHandler object as a member, which is what we use to interconvert
/// between BCL FragmentComplete objects and Rosetta MutableResidueType objects.
///
/// @author
/// Benjamin P. Brown (benjamin.p.brown17@gmail.com)
////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_protocols_drug_design_bcl_BCLFragmentBaseMover_HH
#define INCLUDED_protocols_drug_design_bcl_BCLFragmentBaseMover_HH

// unit forward includes
#include <protocols/drug_design/bcl/BCLFragmentBaseMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// core includes
#include <core/chemical/bcl/BCLFragmentHandler.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/AtomRefMapping.hh>

// package includes
#include <core/chemical/AtomRefMapping.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// utility includes
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// citation manager
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// BCL includes
#ifdef USEBCL
#include <bcl/include/chemistry/bcl_chemistry_atoms_complete_standardizer.h>
#include <bcl/include/chemistry/bcl_chemistry_bond_isometry_handler.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_complete.h>
#include <bcl/include/chemistry/bcl_chemistry_stereocenters_handler.h>
#include <bcl/include/math/bcl_math_mutate_interface.h>
#include <bcl/include/sdf/bcl_sdf.h>
#include <bcl/include/sdf/bcl_sdf_atom_info.h>
#include <bcl/include/sdf/bcl_sdf_bond_info.h>
#include <bcl/include/util/bcl_util_implementation.h>
#include <bcl/include/util/bcl_util_loggers.h>
#endif

namespace protocols
{
namespace drug_design
{
namespace bcl
{
class BCLFragmentBaseMover :
	public protocols::moves::Mover
{
public:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	//! @brief default constructor
	BCLFragmentBaseMover();

#ifdef USEBCL
	//! @brief full constructor
	BCLFragmentBaseMover
	(
		core::chemical::bcl::BCLFragmentHandler const &handler,
		std::string const &ligand_chain_id,
		std::string const &conformation_comparer,
		double const tolerance,
		size_t const &n_confs,
		size_t const &n_max_iterations,
		bool const change_chirality,
		double const random_dihedral_change_prob,
		bool const generate_3d,
		double const clash_tolerance,
		bool const cluster,
		double const clash_resolution,
		bool const sample_dihed,
		bool const sample_rings,
		bool const sample_bonds
	);
#endif

	//! @brief create copy constructor
	protocols::moves::MoverOP clone() const override;

	//! @brief create this type of objectt
	protocols::moves::MoverOP fresh_instance() const override;

	/////////////////
	// data access //
	/////////////////

	static std::string mover_name();

	std::string get_name() const override;

	// data access for private data
protected:
	core::chemical::bcl::BCLFragmentHandler &get_base_handler();
	std::string &get_base_ligand_chain();
	std::string &get_base_conf_comparer();
	double &get_base_conf_tolerance();
	core::Size &get_base_n_max_confs();
	core::Size &get_base_n_max_conf_iterations();
	bool &get_base_conf_sample_chirality();
	double &get_base_conf_rand_dihedral_prob();
	bool &get_base_conf_generate_3d();
	double &get_base_conf_clash_tolerance();
	bool &get_base_conf_cluster();
	double &get_base_conf_clash_resolution();

public:

	////////////////
	// operations //
	////////////////

	//! @brief apply BCLFragmentBaseMover
	void apply( core::pose::Pose & POSE ) override;

	//! @brief set the residue ID of the ligand to be mutated
	void set_ligand( std::string const &ligand_chain_id);

	//! @brief set the fragment conformation comparer type
	void set_conf_comparer( std::string const &conformation_comparer);

	//! @brief set the tolerance for different conformers
	void set_conf_tolerance( double const tolerance);

	//! @brief set the max number of desired fragment conformers
	void set_n_max_confs( size_t const n_confs);

	//! @brief set the max number of iterations for conformer sampling
	void set_n_max_conf_iterations( size_t const n_max_iterations);

	//! @brief allow R/S sampling at chiral centers
	void set_sample_chirality( bool const change_chirality);

	//! @brief set the probability of accepting a random dihedral
	void set_rand_dihed_prob( double const random_dihedral_change_prob);

	//! @brief set de novo conformer generation
	void set_conf_generate_3d( bool const generate_3d);

	//! @brief set the conformer clash tolerance
	void set_conf_clash_tolerance( double const clash_tolerance);

	//! @brief set leader-follower clustering during conf generation
	void set_conf_cluster( bool const cluster);

	//! @brief set the conformer clash resolution
	void set_conf_clash_resolution( double const clash_resolution);

	//! @brief set dihedral bin sampling preference
	void set_conf_sample_dihed( bool const sample_dihed);

	//! @brief set ring sampling preference
	void set_conf_sample_rings( bool const sample_rings);

	//! @brief set bond angle/length sampling preference
	void set_conf_sample_bond_angles( bool const sample_bonds);

	//////////////////////
	// helper functions //
	//////////////////////

#ifdef USEBCL
	//! @brief Convert a Pose Residue to a BCL Fragment.
	//! @param POSE the pose containing the ligand
	//! @param RESNUM number of the residue in the pose
	::bcl::chemistry::FragmentComplete pose_residue_to_fragment
	(
		core::pose::Pose const &pose,
		core::Size const resnum
	);
#endif

	//! @brief initial setup of conformer sampler upon object construction
	void init_fragment_sample_confs();

	//! @brief parse xml file
	void parse_my_tag( TagCOP TAG, basic::datacache::DataMap & data ) override;

	static void base_attributes( utility::tag::AttributeList &attlist);

	//! @brief Provide citations to the passed CitationCollectionList.
	//! @details This mover is unpublished. Hopefully we can publish it soon and I can change
	//! this to a citation for an actual manuscript
	void provide_citation_info( basic::citation_manager::CitationCollectionList &citations) const override;

private:

	//////////
	// data //
	//////////

	//! the handler to convert between BCL and Rosetta small molecule objects
	core::chemical::bcl::BCLFragmentHandler handler_;

	//! the pose residue chain ID of the ligand to be mutated
	std::string ligand_chain_;

	// what follows are a bunch of BCL::Conf options;
	// it is important that we have an internal SampleConformations object
	// otherwise we will not be able to build a new 3D Pose;

	//! @brief conformation comparer tells us how to distinguish conformers
	std::string conf_comparer_; // recommend SymmetryRMSD

	//! @brief tolerance tells us how far away two potential conformers must be in conf_comparer space to be different
	double conf_tolerance_; // recommend 0.25 with SymmetryRMSD, can consider smaller e.g. 0.125 for local sampling

	//! @brief max number of conformers to generate
	core::Size n_max_confs_; // general use, recommend 10 - 250

	//! @brief max number of conformer generation iterations
	core::Size n_max_conf_iterations_; // recommend 4-5x n_max_confs

	//! @brief whether to sample chirality or not
	bool conf_sample_chirality_; // false in most cases

	//! @brief probability of including a randomly selected dihedral for inclusion in conformer
	double conf_rand_dihedral_prob_; // recommend keeping this turned off (prob = 0.0)

	//! @brief create a 3D conformer from scratch; no input coordinate info is used
	bool conf_generate_3d_; // be cafeful setting this to true for complex ring systems

	//! @brief clash tolerance in Angstroms
	double conf_clash_tolerance_; // recommend 0.1

	//! @brief leader-follower clustering during conf sampling
	bool conf_cluster_; // recommend true; improves native conf recovery

	//! @brief clash resolution in Angstroms
	double conf_clash_resolution_ ; // recommend 0.1

	//! @brief skip dihedral sampling
	bool conf_sample_dihed_; // for local set false

	//! @brief skip ring sampling
	bool conf_sample_rings_; // for local set false

	//! @brief skip bond angle/length sampling
	bool conf_sample_bond_angles_; // generally keep true

}; // BCLFragmentBaseMover
} // bcl
} // protocols
} // drug_design

#endif // INCLUDED_protocols_drug_design_bcl_BCLFragmentBaseMover_HH
