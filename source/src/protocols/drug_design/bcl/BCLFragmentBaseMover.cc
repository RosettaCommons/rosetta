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

// unit include
#include <protocols/drug_design/bcl/BCLFragmentBaseMover.hh>

// protocols includes
#include <protocols/drug_design/bcl/BCLSampleConfsManager.hh>
#include <protocols/drug_design/bcl/BCLReferenceSDFilesManager.hh>

// rosetta scripts
#include <protocols/rosetta_scripts/util.hh>

// projects includes
#include <core/chemical/Atom.hh>
#include <core/chemical/Atom.fwd.hh>
#include <core/chemical/AtomProperty.hh>
#include <core/chemical/AtomProperties.fwd.hh>
#include <core/chemical/bcl/util.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/Bond.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>

// utility includes
#include <utility/numbers.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ includes
#include <string>

// debug
#include <utility/io/ozstream.hh>
#include <numeric/xyzVector.hh>

// BCL includes
#ifdef USEBCL
#include <bcl/include/chemistry/bcl_chemistry_atom_complete.h>
#include <bcl/include/chemistry/bcl_chemistry_atom_vector.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_track_mutable_atoms.h>
#include <bcl/include/chemistry/bcl_chemistry_rotamer_library_file.h>
#include <bcl/include/io/bcl_io_ofstream.h>
#include <bcl/include/io/bcl_io_file.h>
#include <bcl/include/linal/bcl_linal_vector_3d.h>
#include <bcl/include/storage/bcl_storage_vector.h>
#include <bcl/include/util/bcl_util_format.h>
#endif

namespace protocols
{
namespace drug_design
{
namespace bcl
{
// tracer
static basic::Tracer TR("protocols.drug_design.bcl.BCLFragmentBaseMover");

//////////
// data //
//////////

//////////////////////////////////
// construction and destruction //
//////////////////////////////////

//! @brief default constructor
BCLFragmentBaseMover::BCLFragmentBaseMover() :
	Mover(),
	handler_( core::chemical::bcl::BCLFragmentHandler()),
	ligand_chain_("X"),
	conf_comparer_( "SymmetryRMSD"),
	conf_tolerance_( 0.25),
	n_max_confs_( 100),
	n_max_conf_iterations_( 500),
	conf_sample_chirality_( false),
	conf_rand_dihedral_prob_( 0.0),
	conf_generate_3d_( false),
	conf_clash_tolerance_( 0.1),
	conf_cluster_( true),
	conf_clash_resolution_( 0.1),
	conf_sample_dihed_( true),
	conf_sample_rings_( true),
	conf_sample_bond_angles_( true)
{
	// now setup the conformer sampler with the member data
	init_fragment_sample_confs();
}

//! @brief default full constructor
#ifdef USEBCL
BCLFragmentBaseMover::BCLFragmentBaseMover
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
) :
	Mover(),
	handler_( handler),
	ligand_chain_( ligand_chain_id),
	conf_comparer_( conformation_comparer),
	conf_tolerance_( tolerance),
	n_max_confs_( n_confs),
	n_max_conf_iterations_( n_max_iterations),
	conf_sample_chirality_( change_chirality),
	conf_rand_dihedral_prob_( random_dihedral_change_prob),
	conf_generate_3d_( generate_3d),
	conf_clash_tolerance_( clash_tolerance),
	conf_cluster_( cluster),
	conf_clash_resolution_( clash_resolution),
	conf_sample_dihed_( sample_dihed),
	conf_sample_rings_( sample_rings),
	conf_sample_bond_angles_( sample_bonds)
{
	// obligatory redundant require_bcl
	core::chemical::bcl::require_bcl();

	// now setup the conformer sampler with the member data
	init_fragment_sample_confs();
}
#endif


//! @brief clone this object
protocols::moves::MoverOP
BCLFragmentBaseMover::clone() const
{
	return utility::pointer::make_shared< BCLFragmentBaseMover >( *this );
}

//! @brief create this type of object
protocols::moves::MoverOP
BCLFragmentBaseMover::fresh_instance() const
{
	return utility::pointer::make_shared< BCLFragmentBaseMover >();
}

/////////////////
// data access //
/////////////////

std::string BCLFragmentBaseMover::mover_name()
{
	return "BCLFragmentBaseMover";
}

std::string BCLFragmentBaseMover::get_name() const
{
	return mover_name();
}

// the following functions retrieve base class member data in lieu of making them protected
core::chemical::bcl::BCLFragmentHandler &BCLFragmentBaseMover::get_base_handler()
{
	return handler_;
}
std::string &BCLFragmentBaseMover::get_base_ligand_chain()
{
	return ligand_chain_;
}

std::string &BCLFragmentBaseMover::get_base_conf_comparer()
{
	return conf_comparer_;
}
double &BCLFragmentBaseMover::get_base_conf_tolerance()
{
	return conf_tolerance_;
}
core::Size &BCLFragmentBaseMover::get_base_n_max_confs()
{
	return n_max_confs_;
}
core::Size &BCLFragmentBaseMover::get_base_n_max_conf_iterations()
{
	return n_max_conf_iterations_;
}
bool &BCLFragmentBaseMover::get_base_conf_sample_chirality()
{
	return conf_sample_chirality_;
}
double &BCLFragmentBaseMover::get_base_conf_rand_dihedral_prob()
{
	return conf_rand_dihedral_prob_;
}
bool &BCLFragmentBaseMover::get_base_conf_generate_3d()
{
	return conf_generate_3d_;
}
double &BCLFragmentBaseMover::get_base_conf_clash_tolerance()
{
	return conf_clash_tolerance_;
}
bool &BCLFragmentBaseMover::get_base_conf_cluster()
{
	return conf_cluster_;
}
double &BCLFragmentBaseMover::get_base_conf_clash_resolution()
{
	return conf_clash_resolution_;
}

////////////////
// operations //
////////////////

//! @brief apply BCLFragmentBaseMover
void BCLFragmentBaseMover::apply( core::pose::Pose & /* POSE */ ) {
	// extras=bcl required for construction of this class object
	core::chemical::bcl::require_bcl();
}

//! @brief set the residue ID of the ligand to be mutated
void BCLFragmentBaseMover::set_ligand( std::string const &ligand_chain_id)
{
	ligand_chain_ = ligand_chain_id;
}

//! @brief set the fragment conformation comparer type
void BCLFragmentBaseMover::set_conf_comparer( std::string const &conformation_comparer)
{
	conf_comparer_ = conformation_comparer;
}

//! @brief set the tolerance for different conformers
void BCLFragmentBaseMover::set_conf_tolerance( double const tolerance)
{
	conf_tolerance_ = tolerance;
}

//! @brief set the max number of desired fragment conformers
void BCLFragmentBaseMover::set_n_max_confs( size_t const n_confs)
{
	n_max_confs_ = n_confs;
}

//! @brief set the max number of iterations for conformer sampling
void BCLFragmentBaseMover::set_n_max_conf_iterations( size_t const n_max_iterations)
{
	n_max_conf_iterations_ = n_max_iterations;
}

//! @brief allow R/S sampling at chiral centers
void BCLFragmentBaseMover::set_sample_chirality( bool const change_chirality)
{
	conf_sample_chirality_ = change_chirality;
}

//! @brief set the probability of accepting a random dihedral
void BCLFragmentBaseMover::set_rand_dihed_prob( double const random_dihedral_change_prob)
{
	conf_rand_dihedral_prob_ = random_dihedral_change_prob;
}

//! @brief set de novo conformer generation
void BCLFragmentBaseMover::set_conf_generate_3d( bool const generate_3d)
{
	conf_generate_3d_ = generate_3d;
}

//! @brief set the conformer clash tolerance
void BCLFragmentBaseMover::set_conf_clash_tolerance( double const clash_tolerance)
{
	conf_clash_tolerance_ = clash_tolerance;
}

//! @brief set leader-follower clustering during conf generation
void BCLFragmentBaseMover::set_conf_cluster( bool const cluster)
{
	conf_cluster_ = cluster;
}

//! @brief set the conformer clash resolution
void BCLFragmentBaseMover::set_conf_clash_resolution( double const clash_resolution)
{
	conf_clash_resolution_ = clash_resolution;
}

//! @brief set dihedral bin sampling preference
void BCLFragmentBaseMover::set_conf_sample_dihed( bool const sample_dihed)
{
	conf_sample_dihed_ = sample_dihed;
}

//! @brief set ring sampling preference
void BCLFragmentBaseMover::set_conf_sample_rings( bool const sample_rings)
{
	conf_sample_rings_ = sample_rings;
}

//! @brief set bond angle/length sampling preference
void BCLFragmentBaseMover::set_conf_sample_bond_angles( bool const sample_bonds)
{
	conf_sample_bond_angles_ = sample_bonds;
}

//////////////////////
// helper functions //
//////////////////////

#ifdef USEBCL
//! @brief Convert a Pose Residue to a BCL Fragment.
//! @param POSE the pose containing the ligand
//! @param RESNUM number of the residue in the pose
::bcl::chemistry::FragmentComplete BCLFragmentBaseMover::pose_residue_to_fragment
(
	core::pose::Pose const &pose,
	core::Size const resnum
)
{
	// get pose residue, which contains xyz info
	core::conformation::Residue raw_res( pose.residue( resnum));

	// get ligand restype, which does not contain xyz info
	core::chemical::MutableResidueTypeOP restype
		(
		utility::pointer::make_shared< core::chemical::MutableResidueType >( pose.residue_type( resnum) )
	);

	// iterate over pose atoms in pose residue and save the xyz coords
	::bcl::storage::Vector< ::bcl::linal::Vector3D> orig_xyz;
	for
		(
				auto atoms_itr( raw_res.atom_begin()), atoms_itr_end( raw_res.atom_end());
				atoms_itr != atoms_itr_end;
				++atoms_itr
				) {
		// let's save some xyz coordinates 0-indexed by atom
		::bcl::linal::Vector3D coords
			(
			atoms_itr->xyz().x(),
			atoms_itr->xyz().y(),
			atoms_itr->xyz().z()
		);
		orig_xyz.PushBack( coords);
	}

	// convert to BCL fragment
	::bcl::chemistry::FragmentComplete fragment( handler_.restype_to_fragment( (*restype), false));

	// the change to restype removed atom coords; reassign coords from the pose residue
	auto coord_itr( orig_xyz.Begin());
	::bcl::storage::Vector< ::bcl::sdf::AtomInfo> frag_atom_info( fragment.GetAtomInfo());
	::bcl::storage::Vector< ::bcl::sdf::BondInfo> frag_bond_info( fragment.GetBondInfo());
	for
		(
				auto atom_itr( frag_atom_info.Begin()), atom_itr_end( frag_atom_info.End());
				atom_itr != atom_itr_end;
				++atom_itr, ++coord_itr
				) {
		atom_itr->SetCoordinates( *coord_itr);
	}

	// make new fragment with correct coordinates
	::bcl::chemistry::AtomVector< ::bcl::chemistry::AtomComplete> atom_v( frag_atom_info, frag_bond_info);
	fragment = ::bcl::chemistry::FragmentComplete( atom_v, raw_res.name());

	// return the BCL fragment
	return fragment;
}
#endif


void BCLFragmentBaseMover::init_fragment_sample_confs()
{
#ifdef USEBCL
	// create our sample confs object
	BCLSampleConfsManager::get_instance()->set_sample_confs(
	::bcl::chemistry::SampleConformations
	(
			BCLSampleConfsManager::get_instance()->get_rotamer_library_file(),
			conf_comparer_,
			conf_tolerance_,
			n_max_confs_,
			n_max_conf_iterations_,
			conf_sample_chirality_,
			conf_rand_dihedral_prob_,
			conf_generate_3d_,
			conf_clash_tolerance_,
			conf_cluster_,
			conf_clash_resolution_
	));
	// sampling bools: dihedrals, rings, bond angles, chirality
	BCLSampleConfsManager::get_instance()->get_sample_confs().SetSamplingPreferences
	(
		conf_sample_dihed_, conf_sample_rings_, conf_sample_bond_angles_, conf_sample_chirality_
	);

#endif
}

//! @brief parse xml file
void BCLFragmentBaseMover::parse_my_tag( TagCOP tag, basic::datacache::DataMap & /*data*/ )
{
	// extras=bcl required for construction of this class object
	core::chemical::bcl::require_bcl();

	// set the ligand pose residue index
	this->set_ligand( tag->getOption< std::string>( "ligand_chain", ligand_chain_));

	// set the conformation comparer type
	this->set_conf_comparer( tag->getOption< std::string>( "conformation_comparer", conf_comparer_));

	// set the conformer tolerance distinguishing conformers
	this->set_conf_tolerance( tag->getOption< double>( "conformer_tolerance", conf_tolerance_));

	// set the maximum number of conformers to generate
	this->set_n_max_confs( tag->getOption< size_t>( "n_max_confs", n_max_confs_));

	// set max iterations
	this->set_n_max_conf_iterations( tag->getOption< size_t>( "n_max_conf_iterations", n_max_conf_iterations_));

	// set sample chirality
	this->set_sample_chirality( tag->getOption< bool>( "sample_chirality", conf_sample_chirality_));

	// set rand dihed mutate prob
	this->set_rand_dihed_prob( tag->getOption< double>( "conf_random_dihed_mutate_prob", conf_rand_dihedral_prob_));

	// set generate 3d
	this->set_conf_generate_3d( tag->getOption< bool>( "conf_generate_3d", conf_generate_3d_));

	// set clash tolerance
	this->set_conf_clash_tolerance( tag->getOption< double>( "conf_clash_tolerance", conf_clash_tolerance_));

	// set cluster
	this->set_conf_cluster( tag->getOption< bool>( "conf_cluster", conf_cluster_));

	// set clash resolution
	this->set_conf_clash_resolution( tag->getOption< double>( "conf_clash_resolution", conf_clash_resolution_));

	// set sample dihed preference
	this->set_conf_sample_dihed( tag->getOption< bool>( "conf_sample_dihed", conf_sample_dihed_));

	// set sample rings preference
	this->set_conf_sample_rings( tag->getOption< bool>( "conf_sample_rings", conf_sample_rings_));

	// set sample bond angles preference
	this->set_conf_sample_bond_angles( tag->getOption< bool>( "conf_sample_bond_angles", conf_sample_bond_angles_));
}

//! @brief set xml options
void BCLFragmentBaseMover::base_attributes( utility::tag::AttributeList &attlist)
{
	// initialize attribute list
	using namespace utility::tag;

	// load attributes
	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"ligand_chain", xs_string,
		"The chain ID of the ligand residue to be mutated. "
		"The specified chain may not contain more than one residue/molecule.",
		"X"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conformation_comparer", xs_string,
		"method with which to compare conformers",
		"SymmetryRMSD"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conformer_tolerance", xsct_real,
		"The amount of tolerance distinguishing conformers",
		"0.25"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"n_max_confs", xsct_non_negative_integer,
		"The maximum number of rconformers to generate",
		"100"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"n_max_conf_iterations", xsct_non_negative_integer,
		"The maximum number of conformer generation iterations to perform; "
		"usually 4-5x the maximum number of desired conformers is sufficient",
		"500"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"sample_chirality", xsct_rosetta_bool,
		"If true, also sample R/S configurations; if false, chirality "
		"will not be changed by BCL::Conf",
		"false"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conf_random_dihed_mutate_prob", xsct_real,
		"The probability of performing a random dihedral perturbation "
		"(ignoring the statistical potential) when generating conformers; "
		"generally not recommended outside of benchmarking",
		"0.0"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conf_generate_3d", xsct_rosetta_bool,
		"If true, ignore all input coordinate information when generating conformers; "
		"if false, retain input coordinate information to initiate sampling; note that "
		"BCL::Conf may have difficulty performing de novo 3D conformer generation for "
		"large macrocycles and complex ring systems",
		"false"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conf_clash_tolerance", xsct_real,
		"The maximum average angstroms clash present across all atoms in the molecule. "
		"A clash is defined as atoms more than three atoms separated that are closer than "
		"the sum of their vdw radii - 0.7 Angstroms ",
		"0.1"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conf_cluster", xsct_rosetta_bool,
		"Perform leader-follower clustering when generating small molecule conformers; "
		"generally recommended as it improves native conformer recovery",
		"true"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conf_clash_resolution", xsct_real,
		"maximum number of tries  (times number of dihedral angles plus number of bond angles) "
		"to resolve clashes permolecule, before clustering. Larger numbers slow the calculation "
		"and make it more likely to find kinetically inaccessible conformational space, but a value of "
		"0 makes it hard to sample conformations of strained molecules, so a value around "
		"0.05-0.5 is usually best",
		"0.1"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conf_sample_dihed", xsct_rosetta_bool,
		"If true, sample across dihedral bins; if false, restrict conformational sampling to current dihedral bins",
		"true"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conf_sample_rings", xsct_rosetta_bool,
		"If true, sample rings; if false, restrict conformational sampling to current ring conformer",
		"true"
	);

	attlist + XMLSchemaAttribute::attribute_w_default
		(
		"conf_sample_bond_angles", xsct_rosetta_bool,
		"If true, sample bond angles and lengths; if false, restrict conformational sampling to current bond angles and lengths",
		"true"
	);
}

//! @brief Provide citations to the passed CitationCollectionList.
//! @details This mover is unpublished. Hopefully we can publish it soon and I can change
//! this to a citation for an actual manuscript
void
BCLFragmentBaseMover::provide_citation_info
(
	basic::citation_manager::CitationCollectionList &citations
) const
{
	using namespace basic::citation_manager;
	UnpublishedModuleInfoOP mover_info( utility::pointer::make_shared< UnpublishedModuleInfo>
		(
		"BCLFragmentBaseMover",
		CitedModuleType::Mover
		));

	// add me
	mover_info->add_author
		(
		"Benjamin P. Brown", "Vanderbilt University", "benjamin.p.brown17@gmail.com"
	);

	// add Jeff
	mover_info->add_author
		(
		"Jeffrey Mendenhall", "Vanderbilt University", "jeffreymendenhall@gmail.com"
	);

	// add Jens the bossman
	mover_info->add_author
		(
		"Jens Meiler", "Vanderbilt University", "jens@meilerlab.org"
	);

	citations.add( mover_info);
}

} // bcl
} // drug_design
} // protocols
