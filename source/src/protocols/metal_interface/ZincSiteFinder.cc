// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/metal_interface/bder/ZincSiteFinder.cc
/// @brief  Searches pose for a zinc residue, then fills a vector of MetalSiteResidue objects with info including sequence position of coordinating sidechains, ligand atom xyz, ligand atom name, and atom ids to provide a convenient way for protocols to add metalsite constraints (ligand refers to protein sidechains)
/// @author Bryan Der

#include <protocols/metal_interface/ZincSiteFinder.hh>
#include <protocols/metal_interface/MetalSiteResidue.hh>
#include <protocols/metal_interface/FindClosestAtom.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.metal_interface.ZincSiteFinder" );

typedef numeric::xyzVector<core::Real> point;
using namespace core;

namespace protocols {
namespace metal_interface {


ZincSiteFinder::ZincSiteFinder() // default constructor
: n_ligands_ (0), zinc_res_(0), parse_error_(false) // 4 ligands to zinc are expected by default - if not 4, then the value must be set in 'set_expecting_n_ligands'
{
	msr_.clear();
}
ZincSiteFinder::ZincSiteFinder( core::Size zinc_res ) // constructor for known zinc sites, particularly useful if there are multiple zinc sites because this protocol would stop after finding the first one
: n_ligands_ (0), zinc_res_(zinc_res), parse_error_(false)
{
	msr_.clear();
}

ZincSiteFinder::~ZincSiteFinder() // destructor
{
}


void
ZincSiteFinder::set_expecting_n_ligands( Size n ) {
	//To check if we've successfully identified the entire metal site.
	//Sometimes large coordination distances can result in not finding all suppos'ed coordinating residues.
	//If set to 0, the program will find however many it finds based on the zn-distance cutoff of 3.0 Angstroms
	n_ligands_ = n;
}

bool
ZincSiteFinder::check_for_parse_error() { //allows pose to be skipped if metal site was not properly found (often due tocoordination distances that are too large)
	return parse_error_;
}


/// @details First finds zinc, then iterates through protein residues until a Cys/His/Asp/Glu sidechain atom (S, N, O) is within 3 Angstroms of the zinc.  Upon finding this residue, it appends the vector of MetalSiteResidue objects.

utility::vector1< protocols::metal_interface::MetalSiteResidueOP >
ZincSiteFinder::find_zinc_site( pose::Pose const & pose )
{
	point zinc;
	Size pose_length = pose.n_residue();
	Size index( 0 ); // used to fill the various metalsite vectors

	//if the zinc residue position was given in the constructor...
	if ( zinc_res_ > 0 && zinc_res_ <= pose_length ) {
		std::string name3 = pose.residue(zinc_res_).name3();
		assert( name3 == " ZN" || name3 == "ZN " || name3 == "ZN" || name3 == "ZNX" || name3 == "HIZ");
		msr_.push_back( protocols::metal_interface::MetalSiteResidueOP( new protocols::metal_interface::MetalSiteResidue ) );
		index++;
		zinc = pose.residue(zinc_res_).atom(1).xyz();
		msr_[1]->set_seqpos( zinc_res_ );
		msr_[1]->set_ligand_atom_xyz( zinc );
		msr_[1]->set_ligand_atom_name("ZN");
		msr_[1]->set_ligand_atom_id( core::id::AtomID( 1 /*zinc is only atom*/, zinc_res_ ) );
	} else {
		//The zinc atom must be found FIRST in order to subsequently find the liganding residues
		for ( Size i(1); i <= pose_length; ++i ) {
			std::string name3 = pose.residue(i).name3();
			if ( name3 == " ZN" || name3 == "ZN " || name3 == "ZN" || name3 == "ZNX" ) {
				msr_.push_back( protocols::metal_interface::MetalSiteResidueOP( new protocols::metal_interface::MetalSiteResidue ) );
				index++;
				zinc = pose.residue(i).atom(1).xyz();
				msr_[1]->set_seqpos( i );
				msr_[1]->set_ligand_atom_xyz( zinc );
				msr_[1]->set_ligand_atom_name("ZN");
				msr_[1]->set_ligand_atom_id( core::id::AtomID( 1 /*zinc is only atom*/, i ) );
				TR << "Found zinc: res " << i << " name " << name3 << std::endl;
				break; //found zinc
			} else if ( pose.residue(i).name3() == "HIZ" ) {
				//HIZ contains a histidine + zinc, I used this as the transition state in RosettaMatch
				msr_.push_back( protocols::metal_interface::MetalSiteResidueOP( new protocols::metal_interface::MetalSiteResidue ) );
				index++;
				zinc = pose.residue(i).atom(5).xyz(); //ZN1 is atom 5 of the HIZ residue
				msr_[1]->set_seqpos( i );
				msr_[1]->set_ligand_atom_xyz( zinc );
				msr_[1]->set_ligand_atom_name("ZN1");
				msr_[1]->set_ligand_atom_id( core::id::AtomID( 5 /*ZN1 is 5th atom*/, i ) );
				TR << "Found zinc: res " << i << " name HIZ" << std::endl;
				break; //found zinc
			}
		}
	}

	//Now that we have the zinc, iterate through Cys/His/Asp/Glu residues, if coordinating atom is within 3.0 Angstroms, store the residue number, residue name, atom ids

	for ( Size i(1); i <= pose_length; ++i ) {
		// AMW: assignment based on residue identity and all the other stuff are separate tasks.
		std::string lig_atom, pre_lig_atom, pre_pre_lig_atom;
		if ( pose.residue(i).name3() == "CYS" && !pose.residue(i).has_variant_type( chemical::DISULFIDE ) ) {
			lig_atom = " SG ";
			pre_lig_atom = "CB";
			pre_pre_lig_atom = "CA";
		} else if ( pose.residue(i).name3() == "HIS" ) {
			lig_atom = protocols::metal_interface::find_closest_atom( pose.residue(i), zinc );
			if ( lig_atom == " ND1" ) {
				pre_lig_atom = "CG";
				pre_pre_lig_atom = "CB";
			} else {
				pre_lig_atom = "CD2";
				pre_pre_lig_atom = "CG";
			}
		} else if ( pose.residue(i).name3() == "ASP" ) {
			lig_atom = protocols::metal_interface::find_closest_atom( pose.residue(i), zinc );
			pre_lig_atom = "CG";
			pre_pre_lig_atom = "CB";
		} else if ( pose.residue(i).name3() == "GLU" ) {
			lig_atom = protocols::metal_interface::find_closest_atom( pose.residue(i), zinc );
			pre_lig_atom = "CD";
			pre_pre_lig_atom = "CG";
		} else {
			continue;
		}

		point p = pose.residue(i).atom(lig_atom).xyz();
		Real dist = zinc.distance( p );
		if ( dist * dist < 9.0 ) {
			msr_.push_back( protocols::metal_interface::MetalSiteResidueOP( new protocols::metal_interface::MetalSiteResidue ) );
			++index;
			msr_[index]->set_resname(pose.residue(i).name3());
			msr_[index]->set_seqpos( i );
			msr_[index]->set_ligand_atom_xyz( p );
			msr_[index]->set_ligand_atom_name(lig_atom.c_str());
			//atom_ids are required for metalsite constraints
			msr_[index]->set_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(lig_atom.c_str()), i ) );
			msr_[index]->set_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(pre_lig_atom.c_str()), i ) );
			msr_[index]->set_pre_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(pre_pre_lig_atom.c_str()), i ) );
		}
	}

	for ( Size ii(1); ii <=index; ii++ ) {
		TR << "msr_" << ii << " " << msr_[ii]->get_resname() << " " << msr_[ii]->get_seqpos() << " " << msr_[ii]->get_ligand_atom_xyz() << " " << msr_[ii]->get_ligand_atom_name() << std::endl;
		TR << msr_[ii]->get_ligand_atom_id() << "  " << msr_[ii]->get_pre_ligand_atom_id() << msr_[ii]->get_pre_pre_ligand_atom_id() << std::endl;
	}


	// if zinc site was not fully found, use pare_error_ to skip this pose but not exit the job
	if ( msr_[1]->get_ligand_atom_name() != "ZN" && msr_[1]->get_ligand_atom_name() != "ZN1" && msr_[1]->get_ligand_atom_name() != "ZNX" ) {
		TR << "Zinc not found after parsing - parse_error=true" << std::endl;
		parse_error_ = true;
	}
	if ( n_ligands_ == 0 ) {
		TR << "Didn't know how many ligands to expect, found " << index - 1 << " after parsing." << std::endl;
	} else if ( index < n_ligands_ + 1 /*add 1 for the metal*/ ) {
		TR << "Metalsite incomplete after parsing (probably due to very bad geometry or incorrect number of expected ligands) - parse_error=true" << std::endl;
		parse_error_ = true;
	} else if ( index > n_ligands_ + 1 /*add 1 for the metal*/ ) {
		TR << "Found too many ligands - parse_error=true" << std::endl;
		parse_error_ = true;
	}

	return msr_;
}//find_zinc_site


}//namespace metal_interface
}//namespace protocols
