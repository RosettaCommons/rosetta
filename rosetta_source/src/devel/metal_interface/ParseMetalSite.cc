// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta% Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   devel/metal_interface/bder/ParseMetalSite.hh
/// @brief  Searches pose for a zinc residue, then fills a vector of MetalSiteResidue objects with info including sequence position of coordinating sidechains, ligand atom xyz, ligand atom name, and atom ids to provide a convenient way for protocols to add metalsite constraints (ligand refers to protein sidechains)
/// @author Bryan Der

#include <devel/metal_interface/ParseMetalSite.hh>
#include <devel/metal_interface/MetalSiteResidue.hh>
#include <devel/metal_interface/FindClosestAtom.hh>

#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh> // to print a point
#include <utility/exit.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.io.hh>



static basic::Tracer TR("devel.metal_interface.ParseMetalSite");

typedef numeric::xyzVector<core::Real> point;
using namespace core;

namespace devel{
namespace metal_interface{


ParseMetalSite::ParseMetalSite() // default constructor
	: n_ligands_ (4), zinc_res_(0), parse_error_(false) // 4 ligands to zinc are expected by default - if not 4, then the value must be set in 'set_expecting_n_ligands'
{
	msr_.clear();
}
ParseMetalSite::ParseMetalSite( core::Size zinc_res ) // constructor for known zinc sites, particularly useful if there are multiple zinc sites because this protocol would stop after finding the first one
	: n_ligands_ (4), zinc_res_(zinc_res), parse_error_(false)
{
	msr_.clear();
}

ParseMetalSite::~ParseMetalSite() // destructor
{
}


void
ParseMetalSite::set_expecting_n_ligands( Size n ) { // will tell us if we've successfully identified the entire metal site.  Sometimes large coordination distances can result in not finding all suppos'ed coordinating residues.
	n_ligands_ = n;
}

bool
ParseMetalSite::check_for_parse_error() { //allows pose to be skipped if metal site was not properly found (often due tocoordination distances that are too large)
	return parse_error_;
}


///@details First finds zinc, then iterates through protein residues until a Cys/His/Asp/Glu sidechain atom (S, N, O) is within 3 Angstroms of the zinc.  Upon finding this residue, it appends the vector of MetalSiteResidue objects.

utility::vector1< devel::metal_interface::MetalSiteResidueOP >
ParseMetalSite::parse_metalsite( pose::Pose const & pose )
{
	point zinc;
	Size pose_length = pose.n_residue();
	Size index( 0 ); // used to fill the various metalsite vectors

	//if the zinc residue position was given in the constructor...
	if(zinc_res_ > 0 && zinc_res_ <= pose_length) {
		std::string name3 = pose.residue(zinc_res_).name3();
		assert( name3 == " ZN" || name3 == "ZN " || name3 == "ZN" || name3 == "ZNX");
		msr_.push_back( new devel::metal_interface::MetalSiteResidue );
		index++;
		zinc = pose.residue(zinc_res_).atom(1).xyz();
		msr_[1]->set_seqpos( zinc_res_ );
		msr_[1]->set_ligand_atom_xyz( zinc );
		msr_[1]->set_ligand_atom_name("ZN");
		msr_[1]->set_ligand_atom_id( core::id::AtomID( 1 /*zinc is only atom*/, zinc_res_ ) );
	}

	else {
		//The zinc atom must be found FIRST in order to subsequently find the liganding residues
		for ( Size i(1); i <= pose_length; ++i ) {

			if ( pose.residue(i).name3() == " ZN" || pose.residue(i).name3() == "ZN " || pose.residue(i).name3() == "ZN" ) {
				msr_.push_back( new devel::metal_interface::MetalSiteResidue );
				index++;
				zinc = pose.residue(i).atom(1).xyz();
				msr_[1]->set_seqpos( i );
				msr_[1]->set_ligand_atom_xyz( zinc );
				msr_[1]->set_ligand_atom_name("ZN");
				msr_[1]->set_ligand_atom_id( core::id::AtomID( 1 /*zinc is only atom*/, i ) );
				break; //found zinc
			}
		}
	}

	//Now that we have the zinc, iterate through Cys/His/Asp/Glu residues, if coordinating atom is within 3.0 Angstroms, store the residue number, residue name, atom ids
	for ( Size i(1); i <= pose_length; ++i ) {

		if ( pose.residue(i).name3() == "CYS" ) {
			point p = pose.residue(i).atom(" SG ").xyz();
			Real dist = zinc.distance( p );
			if ( dist * dist < 9.0 ) {
				msr_.push_back( new devel::metal_interface::MetalSiteResidue );
				++index;
				msr_[index]->set_resname("CYS");
				msr_[index]->set_seqpos( i );
				msr_[index]->set_ligand_atom_xyz( p );
				msr_[index]->set_ligand_atom_name(" SG ");
				//atom_ids are required for metalsite constraints
				msr_[index]->set_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(" SG "), i ) );
				msr_[index]->set_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CB"), i ) );
				msr_[index]->set_pre_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CA"), i ) );

			}
		}

		else if ( pose.residue(i).name3() == "HIS" ) {
			std::string atom_n = devel::metal_interface::find_closest_atom( pose.residue(i), zinc );

			if ( atom_n == " ND1" ) {
				point p = pose.residue(i).atom(" ND1").xyz();
				Real dist = zinc.distance( p );
				if ( dist * dist < 9.0 ) {
					msr_.push_back( new devel::metal_interface::MetalSiteResidue );
					++index;
					msr_[index]->set_resname("HIS");
					msr_[index]->set_seqpos( i );
					msr_[index]->set_ligand_atom_xyz( p );
					msr_[index]->set_ligand_atom_name(" ND1");
					//atom_ids are required for metalsite constraints
					msr_[index]->set_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(" ND1"), i ) );
					msr_[index]->set_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CG"), i ) );
					msr_[index]->set_pre_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CB"), i ) );
				}

			}
			else if ( atom_n == " NE2" ) {
				point p = pose.residue(i).atom(" NE2").xyz();
				Real dist = zinc.distance( p );
				if ( dist * dist < 9.0 ) {
					msr_.push_back( new devel::metal_interface::MetalSiteResidue );
					++index;
					msr_[index]->set_resname("HIS");
					msr_[index]->set_seqpos( i );
					msr_[index]->set_ligand_atom_xyz( p );
					msr_[index]->set_ligand_atom_name(" NE2");
					//atom_ids are required for metalsite constraints
					msr_[index]->set_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(" NE2"), i ) );
					msr_[index]->set_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CD2"), i ) );
					msr_[index]->set_pre_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CG"), i ) );
				}
			}
		}


		else if ( pose.residue(i).name3() == "ASP" ) {
			std::string atom_n = devel::metal_interface::find_closest_atom( pose.residue(i), zinc );

			if ( atom_n == " OD1" ) {
				point p = pose.residue(i).atom(" OD1").xyz();
				Real dist = zinc.distance( p );
				if ( dist * dist < 9.0 ) {
					msr_.push_back( new devel::metal_interface::MetalSiteResidue );
					++index;
					msr_[index]->set_resname("ASP");
					msr_[index]->set_seqpos( i );
					msr_[index]->set_ligand_atom_xyz( p );
					msr_[index]->set_ligand_atom_name(" OD1");
					//atom_ids are required for metalsite constraints
					msr_[index]->set_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(" OD1"), i ) );
					msr_[index]->set_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CG"), i ) );
					msr_[index]->set_pre_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CB"), i ) );
				}
			}

			else if ( atom_n == " OD2" ) {
				point p = pose.residue(i).atom(" OD2").xyz();
				Real dist = zinc.distance( p );
				if ( dist * dist < 9.0 ) {
					msr_.push_back( new devel::metal_interface::MetalSiteResidue );
					++index;
					msr_[index]->set_resname("ASP");
					msr_[index]->set_seqpos( i );
					msr_[index]->set_ligand_atom_xyz( p );
					msr_[index]->set_ligand_atom_name(" OD2");
					//atom_ids are required for metalsite constraints
					msr_[index]->set_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(" OD2"), i ) );
					msr_[index]->set_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CG"), i ) );
					msr_[index]->set_pre_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CB"), i ) );
				}
			}
		}

		else if ( pose.residue(i).name3() == "GLU" ) {
			std::string atom_n = devel::metal_interface::find_closest_atom( pose.residue(i), zinc );
			if ( atom_n == " OE1" ) {
				point p = pose.residue(i).atom(" OE1").xyz();
				Real dist = zinc.distance( p );
				if ( dist * dist < 9.0 ) {
					msr_.push_back( new devel::metal_interface::MetalSiteResidue );
					++index;
					msr_[index]->set_resname("GLU");
					msr_[index]->set_seqpos( i );
					msr_[index]->set_ligand_atom_xyz( p );
					msr_[index]->set_ligand_atom_name(" OE1");
					//atom_ids are required for metalsite constraints
					msr_[index]->set_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(" OE1"), i ) );
					msr_[index]->set_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CD"), i ) );
					msr_[index]->set_pre_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CG"), i ) );
				}
			}

			else if ( atom_n == " OE2" ) {
				point p = pose.residue(i).atom(" OE2").xyz();
				Real dist = zinc.distance( p );
				if ( dist * dist < 9.0 ) {
					msr_.push_back( new devel::metal_interface::MetalSiteResidue );
					++index;
					msr_[index]->set_resname("GLU");
					msr_[index]->set_seqpos( i );
					msr_[index]->set_ligand_atom_xyz( p );
					msr_[index]->set_ligand_atom_name(" OE2");
					//atom_ids are required for metalsite constraints
					msr_[index]->set_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index(" OE2"), i ) );
					msr_[index]->set_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CD"), i ) );
					msr_[index]->set_pre_pre_ligand_atom_id( core::id::AtomID( pose.residue(i).atom_index("CG"), i ) );
				}
			}
		}


	}

	for (Size ii(1); ii <=index; ii++) {
		TR << "msr_" << ii << " " << msr_[ii]->get_resname() << " " << msr_[ii]->get_seqpos() << " " << msr_[ii]->get_ligand_atom_xyz() << " " << msr_[ii]->get_ligand_atom_name() << std::endl;
		TR << msr_[ii]->get_ligand_atom_id() << "  " << msr_[ii]->get_pre_ligand_atom_id() << msr_[ii]->get_pre_pre_ligand_atom_id() << std::endl;
	}


	// ensure proper parsing of the metalsite in such a way that will skip this pose but not exit the job
	if (msr_[1]->get_ligand_atom_name() != "ZN") {
		//utility_exit_with_message("Metalsite ZINC was not properly identified.");
		TR << "Zinc not found after parsing - parse_error=true" << std::endl;
		parse_error_ = true;
	}
	if (index != n_ligands_ + 1 /*add 1 for the metal*/) {
		TR << "Metalsite incomplete after parsing (probably due to very bad geometry or incorrect number of expected ligands) - parse_error=true" << std::endl;
		parse_error_ = true;
	}


  return msr_;
}//parse_metal_site



}//namespace metal_interface
}//namespace devel
