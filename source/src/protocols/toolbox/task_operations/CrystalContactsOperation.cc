// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/CrystalContactsOperation.cc
/// @brief  Exclude crystal contacts from design
/// @author Patrick Conway (ptconway@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/CrystalContactsOperation.hh>
#include <protocols/toolbox/task_operations/CrystalContactsOperationCreator.hh>

// Project Headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/sasa.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C++ Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.task_operations.CrystalContactsOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

core::pack::task::operation::TaskOperationOP
CrystalContactsOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new CrystalContactsOperation );
}


CrystalContactsOperation::CrystalContactsOperation( core::Real all_gap, core::Real polar_gap, core::Real max_buried_sasa, bool invert, bool nbr_radius_to_nbr_radius, bool nbr_radius_to_atoms, bool atoms_to_atoms ):
	all_gap_(all_gap),     // add this to all calculated distances
	polar_gap_(polar_gap),    // if either residue is polar - add this to calculated distances
	max_buried_sasa_(max_buried_sasa),  // ignore buried residues as defined by maximum allowed sasa
	invert_(invert),      // design residues in contact
	nbr_radius_to_nbr_radius_(nbr_radius_to_nbr_radius), // contact determined by nbr radius overlap (CBeta to CBeta)
	nbr_radius_to_atoms_(nbr_radius_to_atoms),      // contact determined by nbr radius to atom distance (CBeta to any atom on symmetric partner)
	atoms_to_atoms_(atoms_to_atoms)           // contact determined by atom to atom distances (any atom to any atom)
{}

CrystalContactsOperation::~CrystalContactsOperation() {}

core::pack::task::operation::TaskOperationOP CrystalContactsOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new CrystalContactsOperation( *this ) );
}

void
CrystalContactsOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core;
	using namespace basic;
	using namespace pose;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace utility;

	if ( !is_symmetric(pose) ) {
		utility_exit_with_message( "Cannot evaluate crystal contacts on an asymmetric pose" );
	}

	// calc sc sasa (mc can be exposed as long as sc is buried)
	Pose asymm_pose;
	extract_asymmetric_unit( pose, asymm_pose, true );
	utility::vector1< core::Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa_sc( asymm_pose, rsd_sasa, true /*normalize*/);

	// find crystal contacts - iterate through all residue pairs between subunit 1 and the other subunits
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	std::set<Size> contacts;
	for ( Size ir=1; ir<=sym_info->num_independent_residues(); ir++ ) {
		TR.Debug << "sasa: " << ir << " " << rsd_sasa[ir] << std::endl;
		if ( rsd_sasa[ir] < max_buried_sasa_ ) continue;      // if residue has no sc SASA, ignore
		if ( !pose.residue(ir).is_protein() ) continue;

		for ( Size jr=sym_info->num_independent_residues()+1; jr<=sym_info->num_total_residues_without_pseudo(); jr++ ) {
			if ( !pose.residue(jr).is_protein() ) continue;

			if ( is_crystal_contact( pose.residue(ir), pose.residue(jr) ) ) {
				contacts.insert(ir);
				break;
			}
		}
	}

	std::string design_select = "";
	for ( Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++ ) {
		// if find residue in contacts or greater than asymm.nres, exclude
		// use != invert_ to flip Boolean, exclude residues not in contact or greater than asymm nres
		if ( ( ((contacts.find(ir) != contacts.end()) != invert_ )) || ir > sym_info->num_independent_residues() ) {
			task.nonconst_residue_task(ir).prevent_repacking();
		} else {
			design_select += ObjexxFCL::string_of(ir)+"+";
		}
	}
	TR << "sele resi " << design_select << std::endl;
}

bool
CrystalContactsOperation::is_crystal_contact( core::conformation::Residue const & asymm_residue, core::conformation::Residue const & symm_residue ) const {
	bool is_contact = false;
	std::string atom_asymm = (asymm_residue.name3() == "GLY") ? "CA" : "CB";
	std::string atom_symm = (symm_residue.name3() == "GLY") ? "CA" : "CB";
	core::Real contact_distance;

	// is CBeta to CBeta distance less than sum of nbr_radii plus gaps?
	if ( nbr_radius_to_nbr_radius_ ) {
		contact_distance = asymm_residue.nbr_radius() + symm_residue.nbr_radius() + all_gap_;
		if ( asymm_residue.is_polar() /*&& asymm_residue.is_polar()*/ ) {
			contact_distance += polar_gap_;
		}
		if ( asymm_residue.xyz(atom_asymm).distance_squared(symm_residue.xyz(atom_symm)) <= contact_distance*contact_distance ) {
			is_contact = true;
			TR.Debug << "contact: " << asymm_residue.seqpos() << " " << symm_residue.seqpos() << " " << contact_distance << std::endl;
		}
	} else if ( nbr_radius_to_atoms_ ) {
		// is CBeta to atom distance less than asymm nbr_radius plus gaps?
		for ( Size atom_symm = symm_residue.natoms(); atom_symm > 0; --atom_symm ) {
			contact_distance = asymm_residue.nbr_radius() + all_gap_;
			if ( asymm_residue.is_polar() && (symm_residue.atom_type(atom_symm).is_acceptor() || symm_residue.atom_type(atom_symm).is_polar_hydrogen() || symm_residue.atom_type(atom_symm).is_donor() ) ) {
				contact_distance += polar_gap_;
			}
			if ( asymm_residue.xyz(atom_asymm).distance_squared(symm_residue.xyz(atom_symm)) <= contact_distance*contact_distance ) {
				is_contact = true;
				TR.Debug << "contact: " << asymm_residue.seqpos() << " " << symm_residue.seqpos() << " " << contact_distance << std::endl;
				TR.Debug << "contact: " << asymm_residue.name3() << " " << symm_residue.name3() << " " << contact_distance << std::endl;
				TR.Debug << "contact: " << atom_asymm << " " << symm_residue.atom_name(atom_symm) << " " << contact_distance << std::endl;
				TR.Debug << "actual distance squared: " << asymm_residue.xyz(atom_asymm).distance_squared(symm_residue.xyz(atom_symm)) << " contact squared: " << contact_distance*contact_distance << std::endl;
				//TR.Debug << "contact: " << asymm_residue.xyz(atom_asymm)[0] << " " << asymm_residue.xyz(atom_asymm)[1] << " " << asymm_residue.xyz(atom_asymm)[2] << " " << symm_residue.xyz(atom_symm)[0] << " " << symm_residue.xyz(atom_symm)[1] << " " << symm_residue.xyz(atom_symm)[2] << " " << contact_distance << std::endl;
				break;
			}
		}
	} else if ( atoms_to_atoms_ ) {
		// are atom to atom distances less than gaps?
		for ( Size atom_asymm = asymm_residue.natoms(); atom_asymm > 0; --atom_asymm ) {
			for ( Size atom_symm = symm_residue.natoms(); atom_symm > 0; --atom_symm ) {
				contact_distance = all_gap_;
				if ( (asymm_residue.atom_type(atom_symm).is_acceptor() && symm_residue.atom_type(atom_symm).is_polar_hydrogen()) ||
						(asymm_residue.atom_type(atom_symm).is_polar_hydrogen() && symm_residue.atom_type(atom_symm).is_acceptor()) ) {
					contact_distance += polar_gap_;
				}
				if ( asymm_residue.xyz(atom_asymm).distance_squared(symm_residue.xyz(atom_symm)) <= contact_distance*contact_distance ) {
					is_contact = true;
					TR.Debug << "contact: " << asymm_residue.seqpos() << " " << symm_residue.seqpos() << " " << contact_distance << std::endl;
					break;
				}
			}
			if ( is_contact ) {
				break;
			}
		}
	}

	return is_contact;
}



void
CrystalContactsOperation::parse_tag( TagCOP tag, DataMap & )
{
	all_gap_ = tag->getOption<core::Real>("all_gap", 0.5);
	polar_gap_ = tag->getOption<core::Real>("polar_gap", 2.5);
	max_buried_sasa_ = tag->getOption<core::Real>("max_buried_sasa", 0.01);
	invert_ = tag->getOption< bool >("invert",0);

	nbr_radius_to_nbr_radius_ = tag->getOption<bool>("nbr_radius_to_nbr_radius", 0);
	nbr_radius_to_atoms_ = tag->getOption<bool>("nbr_radius_to_atoms", 1);
	atoms_to_atoms_ = tag->getOption<bool>("atoms_to_atoms", 0);
}

} //namespace task_operations
} //namespace toolbox
} //namespace protocols
