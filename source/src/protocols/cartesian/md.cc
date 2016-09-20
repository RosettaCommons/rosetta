// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief


/// Note: This is an extremely crude and slow implementation of cartesian MD in rosetta.
/// I make no promises as to this codes correctness. use at own risk.
/// Nor is this code efficient. Nor is it intended for production simulations.


// Unit headers
#include <protocols/cartesian/md.hh>

// Package headers
#include <core/id/AtomID.hh>
#include <core/optimization/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

// // ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// // Numeric headers

#include <iostream>
#include <fstream>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/random/random.fwd.hh>


using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace cartesian {

static THREAD_LOCAL basic::Tracer TR( "core.optimization.md" );

float sqr(float t){ return t*t; }

MolecularDynamics::MolecularDynamics(
	core::pose::PoseOP & inputpose,
	core::scoring::ScoreFunction const & scorefxn
):
	pose( inputpose )
{

	scorefxn( *pose );
	mm.set_bb (  true );
	mm.set_chi(  true );
	min_map.setup( *pose, mm );
	core::kinematics::DomainMap const & dmap = min_map.domain_map();

	pose->energies().set_use_nblist( *pose, dmap, true );
	scorefxn.setup_for_minimizing( *pose, min_map );

	createCartesianArray( );
	createBondList(  );
	createAngleList(  );
	createDihedralList(  );
	setDihedralDerivatives();

}

void MolecularDynamics::createCartesianArray( )
{
	using namespace core;

	Size const nres( pose->size() );
	cartom.clear();
	for ( Size ir = 1; ir <= nres; ++ir ) {
		conformation::Residue const & rsd( pose->residue( ir ) );
		for ( Size i = 1; i <= rsd.natoms(); i++ ) {
			id::AtomID atom_id(  i, ir );
			CartesianAtom atom;
			atom.index = i;
			atom.res = ir;
			atom.atom_id = atom_id;

			atom.force = core::Vector(0,0,0);
			atom.velocity = core::Vector(0,0,0);
			atom.position = pose->xyz( atom_id );

			atom.old_force = core::Vector(0,0,0);
			atom.old_velocity = core::Vector(0,0,0);
			atom.old_position = core::Vector(0,0,0);

			atom.mass = 12.001;  // assume carbon for now, we can be more accurate later.
			cartom.push_back(atom);
		}
	}

}

void MolecularDynamics::setCartesianPositionsFromPose( )
{
	using namespace core;

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		cartom[i].position = pose->xyz( cartom[i].atom_id );
	}
}


void MolecularDynamics::setPosePositionsFromCartesian( )
{
	using namespace core;

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		pose->set_xyz( cartom[i].atom_id,  cartom[i].position );
	}
}


void MolecularDynamics::zeroForces( )
{
	using namespace core;

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		cartom[i].force = core::Vector( 0,0,0 );
	}
}

int MolecularDynamics::findCartomAtom(
	const core::id::AtomID  &id1
)
{
	using namespace core;

	int number=-1; int i;
	for ( i = 1; i <= (int)cartom.size(); i ++ ) {
		if ( (cartom[i].atom_id.rsd() == id1.rsd() ) &&
				(cartom[i].atom_id.atomno() == id1.atomno() )
				) { number = i; break; }
	}
	if ( i > (int)cartom.size() ) number = -1;
	return number;
}

void MolecularDynamics::getCartesianDerivatives(
	core::scoring::ScoreFunction const & scorefxn
)
{
	using namespace core;

	using namespace core::optimization;
	Size const nres( pose->size() );
	scorefxn.setup_for_derivatives( *pose);
	int count = 1;
	for ( Size ir = 1; ir <= nres; ++ir ) {
		conformation::Residue const & rsd( pose->residue( ir ) );
		for ( Size i = 1; i <= rsd.natoms(); i++ ) {
			id::AtomID atom_id(  i, ir );
			core::Vector F1(0,0,0),F2(0,0,0);
			scorefxn.eval_npd_atom_derivative( atom_id, *pose, min_map.domain_map(), F1, F2 );
			cartom[count].force = F2;
			count += 1;
		}
	}

	//std::cout << "Get cartesian derivates \n";
	// Now handle torsonal derivatives aka dEdphi


	// now loop over the torsions in the map (the map MUST be a map of everything!)

	int imap( 1 ); // for indexing into de_dvars( imap )
	for ( auto it=min_map.begin(), ite=min_map.end();
			it != ite; ++it, ++imap ) {
		using namespace id;

		DOF_Node const & dof_node( **it );

		// NOTE: deriv is in the units of the degree of freedom as
		// represented internally, without any scale factors applied
		// ie the units returned by pose->get_atom_tree_torsion(...)
		//
		//
		// type  -- units
		// -------------------
		// PHI   -- radians
		// THETA -- radians
		// D     -- angstroms
		// RB1-3 -- angstroms
		// RB4-6 -- degrees    (!)

		//Real deriv = 0.0;

		/*kinematics::tree::Atom const & atom
		( pose->atom_tree().atom( dof_node.atom_id() ) );*/

		// eg rama,Paa,dunbrack,and torsional constraints
		Real deriv = scorefxn.eval_dof_derivative
			( dof_node.dof_id(), dof_node.torsion_id(), *pose );
		// deriv /= min_map.torsion_scale_factor( dof_node );


		if ( deriv == 0 ) continue;
		// work out which atoms are involved !
		AtomID id1;
		AtomID id2;
		AtomID id3;
		AtomID id4;

		pose->conformation().get_torsion_angle_atom_ids(
			dof_node.torsion_id(),
			id1,id2,id3,id4 );

		// std::cout << " " << id1.atomno() << "  " << id2.atomno() << "  " << id3.atomno() << "  " << id4.atomno() << "  -- ";

		// work out the inices in the cartom array !
		//Size i;
		int no1= findCartomAtom( id1 ); if ( no1 < 0 ) std::cerr << "ERROR ERROR ERROR id1 \n";
		int no2= findCartomAtom( id2 ); if ( no2 < 0 ) std::cerr << "ERROR ERROR ERROR id2 \n";
		int no3= findCartomAtom( id3 ); if ( no3 < 0 ) std::cerr << "ERROR ERROR ERROR id3 \n";
		int no4= findCartomAtom( id4 ); if ( no4 < 0 ) std::cerr << "ERROR ERROR ERROR id4 \n";
		//  int no2=-1; for ( i = 1; i <= cartom.size(); i ++ ) if( cartom[i].atom_id == id2 ) { no2 = i; break; } if( i > cartom.size()) std::cerr << "ERROR ERROR ERROR no2 \n";
		//  int no3=-1; for ( i = 1; i <= cartom.size(); i ++ ) if( cartom[i].atom_id == id3 ) { no3 = i; break; } if( i > cartom.size()) std::cerr << "ERROR ERROR ERROR no3 \n";
		//  int no4=-1; for ( i = 1; i <= cartom.size(); i ++ ) if( cartom[i].atom_id == id4 ) { no4 = i; break; } if( i > cartom.size()) std::cerr << "ERROR ERROR ERROR no4 \n";


		//std::cout << " " << cartom[no1].atom_id.rsd() << " " << cartom[no2].atom_id.rsd()
		//         << " " << cartom[no3].atom_id.rsd() << " " << cartom[no4].atom_id.rsd();
		//std::cout << " " << deriv << "  ";


		// convert to cartesian derivatives on those atoms using the usual MD code - lift from PD


		//float  epot_dihedral_this=0;
		float  /*phi,*/ sin_phi, cos_phi;

		core::Vector fi(0), fab(0), fj(0);

		core::Vector vti_vta = cartom[no1].position - cartom[no2].position;
		core::Vector vta_vtb = cartom[no2].position - cartom[no3].position;
		core::Vector vtb_vtj = cartom[no3].position - cartom[no4].position;

		core::Vector nrml1 = cross_product( vti_vta, vta_vtb );
		core::Vector nrml2 = cross_product( vta_vtb, vtb_vtj );
		core::Vector nrml3 = cross_product( vta_vtb, nrml1   );

		float inv_nrml1_mag = 1.0 / nrml1.magnitude();
		float inv_nrml2_mag = 1.0 / nrml2.magnitude();
		float inv_nrml3_mag = 1.0 / nrml3.magnitude();

		cos_phi = dot( nrml1, nrml2 ) * inv_nrml1_mag * inv_nrml2_mag;
		sin_phi = dot( nrml3, nrml2 ) * inv_nrml3_mag * inv_nrml2_mag;

		nrml2 *= inv_nrml2_mag;

		if ( fabs( sin_phi ) > 0.1 ) {
			nrml1 *= inv_nrml1_mag;
			core::Vector dcosdnrml1 = inv_nrml1_mag * ( nrml1 * cos_phi - nrml2 );
			core::Vector dcosdnrml2 = inv_nrml2_mag * ( nrml2 * cos_phi - nrml1 );

			deriv *= -1;

			deriv /= sin_phi;
			fi += deriv * cross_product( vta_vtb, dcosdnrml1 );
			fj += deriv * cross_product( vta_vtb, dcosdnrml2 );
			fab += deriv * ( cross_product( vti_vta, dcosdnrml1 ) + cross( vtb_vtj, dcosdnrml2 ) );
		} else {
			nrml3 *= inv_nrml3_mag;
			core::Vector dsindnrml3 = inv_nrml3_mag * ( nrml3 * sin_phi - nrml2 );
			core::Vector dsindnrml2 = inv_nrml2_mag * ( nrml2 * sin_phi - nrml3 );

			deriv *= -1;

			deriv /= -cos_phi;
			fi += deriv * update_operation( vta_vtb, dsindnrml3 );
			fj += deriv * ( cross( dsindnrml2, vta_vtb ) );
			fab += deriv * update_5way_operation( vta_vtb, vti_vta, dsindnrml3, dsindnrml2, vtb_vtj );
		}

		//  fi.mul(PhysicsConst::invAngstrom);
		//  fab.mul(PhysicsConst::invAngstrom);
		//  fj.mul(PhysicsConst::invAngstrom);

		// add to cartesian derivatives

		//    std::cout << fi.x() << "   " ;
		//    std::cout << fi.y() << "   " ;
		//    std::cout << fi.z() << "   " ;
		//    std::cout << "   |  ";
		//    std::cout << fab.x() - fi.x() << "   " ;
		//    std::cout << fab.y() - fi.y() << "   " ;
		//    std::cout << fab.z() - fi.z() << "   " ;
		//    std::cout << "   |  ";
		//    std::cout << fj.x() - fab.x() << "   " ;
		//    std::cout << fj.y() - fab.y() << "   " ;
		//    std::cout << fj.z() - fab.z() << "   " ;
		//    std::cout << "   |  ";
		//    std::cout << fj.x() << "   " ;
		//    std::cout << fj.y() << "   " ;
		//    std::cout << fj.z() << "   " ;
		//    std::cout << std::endl;

		cartom[no1].force += fi;
		cartom[no2].force += fab - fi;
		cartom[no3].force += fj - fab;
		cartom[no4].force -= fj;

		//epot_dihedral += epot_dihedral_this;

		// DONE

	} // loop over map


}


void MolecularDynamics::createBondList( )
{
	using namespace core;

	Size const nres( pose->size() );
	bondlist.clear();

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		// for each atom get neighbours within the residue;
		conformation::Residue const & rsd( pose->residue( cartom[i].atom_id.rsd() ) );
		core::chemical::AtomIndices indices = rsd.bonded_neighbor( cartom[i].atom_id.atomno() );
		for ( Size j = 1; j <= indices.size(); j++ ) {
			if ( int(indices[j]) < int(cartom[i].atom_id.atomno()) ) continue;
			//std::cout << indices[j] << "  ";
			MD_Bond newbond;
			newbond.atom_id_1 = cartom[i].atom_id;
			newbond.atom_id_2 = id::AtomID( indices[j], cartom[i].atom_id.rsd()  );
			newbond.index1 = -1;
			newbond.index2 = -1;
			newbond.length = 0.0;

			bondlist.push_back( newbond );
		}
	}

	// add residue bonds!
	for ( Size ir = 1; ir < nres; ++ir ) {
		conformation::Residue const & rsd1( pose->residue( ir ) );
		conformation::Residue const & rsd2( pose->residue( ir+1 ) );

		if ( !rsd1.is_bonded( rsd2 ) ) continue;

		//`nstd::cout << "RES " << ir << std::endl;

		int atom1 = rsd1.connect_atom( rsd2 );
		int atom2 = rsd2.connect_atom( rsd1 );

		MD_Bond newbond;
		newbond.atom_id_1 = id::AtomID( atom1, ir  );
		newbond.atom_id_2 = id::AtomID( atom2, ir+1  );
		newbond.index1 = atom1;
		newbond.index2 = atom2;
		newbond.length = 0.0;

		bondlist.push_back( newbond );
	}

	// post process the list
	// find all the indices
	for ( Size i = 1; i <= bondlist.size(); i ++ ) {
		bondlist[i].index1 = -1;
		bondlist[i].index2 = -1;
		for ( Size j = 1; j <= cartom.size(); j ++ ) {
			if ( bondlist[i].atom_id_1 == cartom[j].atom_id ) bondlist[i].index1 = j;
			if ( bondlist[i].atom_id_2 == cartom[j].atom_id ) bondlist[i].index2 = j;
		}
		if ( bondlist[i].index1 <= 0 ) std::cout << "CODE ERROR!\n";
		if ( bondlist[i].index2 <= 0 ) std::cout << "CODE ERROR!\n";

		// measure length

		bondlist[i].length = pose->xyz( bondlist[i].atom_id_1 ).distance( pose->xyz( bondlist[i].atom_id_2 ) );
		/*
		std::cout <<
		bondlist[i].atom_id_1.rsd()       << "  " <<
		bondlist[i].atom_id_1.atomno()    << "  " <<
		bondlist[i].index1                << "  " <<
		bondlist[i].atom_id_2.rsd()       << "  " <<
		bondlist[i].atom_id_2.atomno()    << "  " <<
		bondlist[i].index2                << "  " <<
		bondlist[i].length                << "  " <<
		std::endl;
		*/
	}
}


void MolecularDynamics::createAngleList( )
{
	using namespace core;

	//Size const nres( pose->size() );
	anglelist.clear();

	for ( Size i = 1; i <= bondlist.size(); i ++ ) {
		for ( Size j = 1; j <= bondlist.size(); j ++ ) {
			if ( i == j ) continue;

			MD_Angle newangle;
			if (  bondlist[i].index1 ==  bondlist[j].index1 ) {
				newangle.index1 = bondlist[i].index2;
				newangle.index2 = bondlist[i].index1;
				newangle.index3 = bondlist[j].index2;
			} else if (  bondlist[i].index1 ==  bondlist[j].index2 ) {
				newangle.index1 = bondlist[i].index2;
				newangle.index2 = bondlist[i].index1;
				newangle.index3 = bondlist[j].index1;
			} else if (  bondlist[i].index2 ==  bondlist[j].index1 ) {
				newangle.index1 = bondlist[i].index1;
				newangle.index2 = bondlist[i].index2;
				newangle.index3 = bondlist[j].index2;
			} else if (  bondlist[i].index2 ==  bondlist[j].index2 ) {
				newangle.index1 = bondlist[i].index1;
				newangle.index2 = bondlist[i].index2;
				newangle.index3 = bondlist[j].index1;
			} else {
				continue;
			}

			if ( newangle.index1 > newangle.index3 ) continue;

			newangle.atom_id_1 = cartom[ newangle.index1 ].atom_id;
			newangle.atom_id_2 = cartom[ newangle.index2 ].atom_id;
			newangle.atom_id_3 = cartom[ newangle.index3 ].atom_id;
			newangle.length = pose->xyz( newangle.atom_id_1 ).distance( pose->xyz( newangle.atom_id_3 ) );
			/*
			std::cout <<
			newangle.atom_id_1.rsd()       << "  " <<
			newangle.atom_id_1.atomno()    << "  " <<
			newangle.index1                << "  " <<
			newangle.atom_id_2.rsd()       << "  " <<
			newangle.atom_id_2.atomno()    << "  " <<
			newangle.index2                << "  " <<
			newangle.atom_id_3.rsd()       << "  " <<
			newangle.atom_id_3.atomno()    << "  " <<
			newangle.index3                << "  " <<
			newangle.length                << "  " <<
			std::endl;
			*/
			anglelist.push_back( newangle );
		}
	}
}

MD_HarmonicDihedral MolecularDynamics::createDihedral(
	const core::conformation::Residue &rsd1,
	const core::conformation::Residue &rsd2,
	const core::conformation::Residue &rsd3,
	const core::conformation::Residue &rsd4,
	std::string const & name1,
	std::string const & name2,
	std::string const & name3,
	std::string const & name4
){
	using namespace core;

	MD_HarmonicDihedral newdihed;
	newdihed.atom_id_1 = id::AtomID( rsd1.atom_index( name1 ) ,rsd1.seqpos() );
	newdihed.atom_id_2 = id::AtomID( rsd2.atom_index( name2 ) ,rsd2.seqpos() );
	newdihed.atom_id_3 = id::AtomID( rsd3.atom_index( name3 ) ,rsd3.seqpos() );
	newdihed.atom_id_4 = id::AtomID( rsd4.atom_index( name4 ) ,rsd4.seqpos() );
	newdihed.index1  = findCartomAtom(  newdihed.atom_id_1 );
	newdihed.index2  = findCartomAtom(  newdihed.atom_id_2 );
	newdihed.index3  = findCartomAtom(  newdihed.atom_id_3 );
	newdihed.index4  = findCartomAtom(  newdihed.atom_id_4 );

	if ( (newdihed.index1 < 0) ) std::cout << "Can't find" << name1 << std::endl;
	if ( (newdihed.index2 < 0) ) std::cout << "Can't find" << name2 << std::endl;
	if ( (newdihed.index3 < 0) ) std::cout << "Can't find" << name3 << std::endl;
	if ( (newdihed.index4 < 0) ) std::cout << "Can't find" << name4 << std::endl;
	if (
			(newdihed.index1 < 0) ||
			(newdihed.index2 < 0) ||
			(newdihed.index3 < 0) ||
			(newdihed.index4 < 0)
			) {
		std::cout
			<< name1 << " "
			<< name2 << " "
			<< name3 << " "
			<< name4 << " "
			<< newdihed.atom_id_1.rsd() << " "
			<< newdihed.atom_id_2.rsd()  << " "
			<< newdihed.atom_id_3.rsd()  << " "
			<< newdihed.atom_id_4.rsd()  << " "
			<< newdihed.atom_id_1.atomno() << " "
			<< newdihed.atom_id_2.atomno() << " "
			<< newdihed.atom_id_3.atomno() << " "
			<< newdihed.atom_id_4.atomno() << " "
			<< std::endl;
		exit(0);
	}

	newdihed.angle = 0;
	return newdihed;
}


MD_HarmonicDihedral MolecularDynamics::createDihedral(
	const core::conformation::Residue &rsd,
	std::string const & name1,
	std::string const & name2,
	std::string const & name3,
	std::string const & name4
){
	return createDihedral( rsd, rsd, rsd, rsd, name1, name2, name3, name4 );
}

/// Yes i know this is hard coded stuff. Read warning at the top of this file.
void MolecularDynamics::createDihedralList(  )
{
	using namespace core;
	int ir;
	for ( ir = 1; ir < (int)pose->size(); ir ++ ) {
		const  conformation::Residue &rsd = pose->residue( ir );


		if ( ( !(rsd.aa() == chemical::aa_pro) )  &&
				( ir != 1 ) &&
				( ir <  ( (int)pose->size() - 1 ) )
				) {

			dihedrallist.push_back( createDihedral(
				pose->residue( ir - 1 ),
				pose->residue( ir ),
				pose->residue( ir ),
				pose->residue( ir ),
				"C",  "CA",   "N",   "H"  ) );

			dihedrallist.push_back( createDihedral(
				pose->residue( ir ),
				pose->residue( ir + 1 ),
				pose->residue( ir ),
				pose->residue( ir ),
				"CA",  "N",   "C",   "O" ) );

		}

		//if( rsd.aa() == chemical::aa_ala )
		//if( rsd.aa() == chemical::aa_cys )
		if ( rsd.aa() == chemical::aa_asp ) {
			//std::cout << "Preparring ASP\n";
			dihedrallist.push_back( createDihedral( rsd,   "CB",  "OD1",  "CG",  "OD2" ) );
		}
		if ( rsd.aa() == chemical::aa_glu ) {
			//std::cout << "Preparring ASP\n";
			dihedrallist.push_back( createDihedral( rsd,   "CG",  "OE1",  "CD",  "OE2" ) );
		}
		if ( rsd.aa() == chemical::aa_phe ) {
			//std::cout << "Preparring PHE\n";

			// Ring Dihedrals
			dihedrallist.push_back( createDihedral( rsd,    "CB",  "CG", "CD2", "HD2"   ) );

			dihedrallist.push_back( createDihedral( rsd,    "HE2", "CE2", "CZ", "HZ"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "HZ",  "CZ", "CE1","HE1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HE1",  "CE1","CD1", "HD1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HD1", "CD1", "CG", "CB"   ) );

			dihedrallist.push_back( createDihedral( rsd,    "CD1", "CG", "CD2", "HD2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CG",  "CD2","CE2", "HE2"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CD2", "CE2", "CZ", "HZ"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CE2", "CZ", "CE1","HE1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CZ",  "CE1","CD1", "HD1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CE1", "CD1", "CG", "CB"   ) );

			dihedrallist.push_back( createDihedral( rsd,    "CB",  "CG", "CD2", "CE2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "HD2",  "CD2","CE2", "CZ"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HE2", "CE2", "CZ", "CE1"  ) );
			dihedrallist.push_back( createDihedral( rsd,    "HZ",  "CZ", "CE1","CD1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HE1",  "CE1","CD1", "CG"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HD1", "CD1", "CG", "CD2"  ) );

			// Impropers
			dihedrallist.push_back( createDihedral( rsd,    "CG",  "CE2", "CD2", "HD2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CD2",  "CZ", "CE2", "HE2"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CZ",  "CD2", "CE2", "HE2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CE1",  "CE2", "CZ", "HZ"     ) );
			dihedrallist.push_back( createDihedral( rsd,    "CD1",  "CZ", "CE1", "HE1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CG",  "CE1", "CD1", "HD1"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CD1",  "CD2", "CG", "CB"     ) );
		}

		//if( rsd.aa() == chemical::aa_gly )
		if ( rsd.aa() == chemical::aa_his ) {

			if ( rsd.has( "HE2" ) ) {
				dihedrallist.push_back( createDihedral( rsd,    "CE1",  "CD2", "NE2","HE2"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "CG",  "NE2","CD2", "HD2"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "ND1", "NE2", "CE1", "HE1"   ) );
				dihedrallist.push_back( createDihedral( rsd,    "ND1",  "CD2", "CG", "CB"     ) );

				dihedrallist.push_back( createDihedral( rsd,    "CB",  "CG", "CD2","HD2"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "HD2","CD2",  "NE2","HE2"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "HE2", "NE2", "CE1", "HE1"   ) );
				dihedrallist.push_back( createDihedral( rsd,    "HE1",  "CE1", "ND1","CG"     ) );

				dihedrallist.push_back( createDihedral( rsd,    "CB",  "CG", "CD2","NE2"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "HD2","CD2",  "NE2","CE1"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "HE2", "NE2", "CE1", "ND1"   ) );
			} else {
				// HD1 torsions misssing - fix this - mike
				dihedrallist.push_back( createDihedral( rsd,    "CG",  "NE2","CD2", "HD2"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "ND1", "NE2", "CE1", "HE1"   ) );
				dihedrallist.push_back( createDihedral( rsd,    "ND1",  "CD2", "CG", "CB"     ) );

				dihedrallist.push_back( createDihedral( rsd,    "CB",  "CG", "CD2","HD2"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "HE1",  "CE1", "ND1","CG"     ) );

				dihedrallist.push_back( createDihedral( rsd,    "CB",  "CG", "CD2","NE2"    ) );
				dihedrallist.push_back( createDihedral( rsd,    "HD2","CD2",  "NE2","CE1"    ) );


			}
		}
		//if( rsd.aa() == chemical::aa_ile )
		//if( rsd.aa() == chemical::aa_lys )
		//if( rsd.aa() == chemical::aa_leu )
		//if( rsd.aa() == chemical::aa_met )
		if ( rsd.aa() == chemical::aa_asn ) {
			//std::cout << "Preparring ASN\n";
			dihedrallist.push_back( createDihedral( rsd,    "CB",  "ND2",  "CG",  "OD1"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CG",  "1HD2", "ND2", "2HD2"  ) );

			dihedrallist.push_back( createDihedral( rsd,    "OD1", "CG",  "ND2", "1HD2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "OD1", "CG",  "ND2", "2HD2"  ) );
			dihedrallist.push_back( createDihedral( rsd,    "CB", "CG",  "ND2", "1HD2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CB", "CG",  "ND2", "2HD2"  ) );

		}
		//if( rsd.aa() == chemical::aa_pro )
		if ( rsd.aa() == chemical::aa_gln ) {
			//std::cout << "Preparring GLN\n";
			dihedrallist.push_back( createDihedral( rsd,  "CG", "NE2",  "CD",  "OE1" ) );
			dihedrallist.push_back( createDihedral( rsd,  "CD", "1HE2", "NE2", "2HE2" ) );

			dihedrallist.push_back( createDihedral( rsd,    "OE1", "CG",  "NE2", "1HE2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "OE1", "CG",  "NE2", "2HE2"  ) );
			dihedrallist.push_back( createDihedral( rsd,    "CG", "CD",  "NE2", "1HE2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CG", "CD",  "NE2", "2HE2"  ) );
		}
		if ( rsd.aa() == chemical::aa_arg ) {
			//std::cout << "Preparring ARG\n";
			dihedrallist.push_back( createDihedral( rsd,  "NE",   "NH1", "CZ",  "NH2"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "CD",   "CZ", "NE",  "HE"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "CZ",   "1HH1", "NH1",  "2HH1"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "CZ",   "1HH2", "NH2",  "2HH2"   ) );

			dihedrallist.push_back( createDihedral( rsd,  "NH1",   "CZ", "NH2",  "1HH2"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NH1",   "CZ", "NH2",  "2HH2"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NE",    "CZ", "NH2",  "1HH2"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NE",    "CZ", "NH2",  "2HH2"   ) );

			dihedrallist.push_back( createDihedral( rsd,  "NH2",   "CZ", "NH1",  "1HH1"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NH2",   "CZ", "NH1",  "2HH1"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NE",    "CZ", "NH1",  "1HH1"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NE",    "CZ", "NH1",  "2HH1"   ) );

			dihedrallist.push_back( createDihedral( rsd,  "NH2",   "CZ", "NE",  "CD"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NH2",   "CZ", "NE",  "HE"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NH1",    "CZ", "NE",  "CD"   ) );
			dihedrallist.push_back( createDihedral( rsd,  "NH2",    "CZ", "NE",  "HE"   ) );


		};
		//if( rsd.aa() == chemical::aa_ser )
		//if( rsd.aa() == chemical::aa_thr
		//if( rsd.aa() == chemical::aa_val )
		if ( rsd.aa() == chemical::aa_trp ) {

			//std::cout << "Preparring TRP\n";

			dihedrallist.push_back( createDihedral( rsd,  "CD1",  "CE2",     "NE1",     "HE1"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CE2",  "CH2",     "CZ2",     "HZ2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CZ2",  "CZ3",     "CH2",     "HH2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CH2",  "CE3",     "CZ3",     "HZ3"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CZ3",  "CD2",     "CE3",     "HE3"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CG",   "NE1",     "CD1",     "HD1"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CD1",  "CD2",     "CG",      "CB"   )  );

			dihedrallist.push_back( createDihedral( rsd,  "CB",   "CG",      "CD1",     "HD1"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HD1",  "CD1",     "NE1",     "HE1"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HE1",  "NE1",     "CE2",     "CZ2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "NE1",  "CE2",     "CZ2",     "HZ2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CE2",  "CZ2",     "CH2",     "HH2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HH2",  "CH2",     "CZ3",     "HZ3"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HZ3",  "CZ3",     "CE3",     "HE3"  )  );
			dihedrallist.push_back( createDihedral( rsd,  "HE3",  "CE3",     "CD2",     "CG"  )  );
			dihedrallist.push_back( createDihedral( rsd,  "CE3",  "CD2",     "CG",      "CB"  )  );

			dihedrallist.push_back( createDihedral( rsd,  "CB",   "CG",      "CD1",     "NE1"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HD1",  "CD1",     "NE1",     "CE2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HE1",  "NE1",     "CE2",     "CD2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "NE1",  "CE2",     "CZ2",     "CH2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CE2",  "CZ2",     "CH2",     "CZ3"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HH2",  "CH2",     "CZ3",     "CE3"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HZ3",  "CZ3",     "CE3",     "CD2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "HE3",  "CE3",     "CD2",     "CE2"   )  );
			dihedrallist.push_back( createDihedral( rsd,  "CE3",  "CD2",     "CG",      "CD1"   )  );


			dihedrallist.push_back( createDihedral( rsd,  "CD2", "CE2", "CZ2", "CH2" )  );
			dihedrallist.push_back( createDihedral( rsd,  "CE2", "CZ2", "CH2", "CZ3" )  );
			dihedrallist.push_back( createDihedral( rsd,  "CZ2", "CH2", "CZ3", "CE3" )  );
			dihedrallist.push_back( createDihedral( rsd,  "CH2", "CZ3", "CE3", "CD2" )  );
			dihedrallist.push_back( createDihedral( rsd,  "CZ3", "CE3", "CD2", "CE2" )  );
			dihedrallist.push_back( createDihedral( rsd,  "CE3", "CD2", "CE2", "CZ2" )  );


		}

		if ( rsd.aa() == chemical::aa_tyr ) {
			//std::cout << "Preparring TYR\n";

			// Ring Dihedrals
			dihedrallist.push_back( createDihedral( rsd,    "CB",  "CG", "CD2", "HD2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "HD2",  "CD2","CE2", "HE2"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HE2", "CE2", "CZ", "OH"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "OH",  "CZ", "CE1","HE1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HE1",  "CE1","CD1", "HD1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HD1", "CD1", "CG", "CB"   ) );

			dihedrallist.push_back( createDihedral( rsd,    "CD1", "CG", "CD2", "HD2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CG",  "CD2","CE2", "HE2"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CD2", "CE2", "CZ", "OH"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CE2", "CZ", "CE1","HE1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CZ",  "CE1","CD1", "HD1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CE1", "CD1", "CG", "CB"   ) );

			dihedrallist.push_back( createDihedral( rsd,    "CB",  "CG", "CD2", "CE2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "HD2",  "CD2","CE2", "CZ"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HE2", "CE2", "CZ", "CE1"  ) );
			dihedrallist.push_back( createDihedral( rsd,    "OH",  "CZ", "CE1","CD1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HE1",  "CE1","CD1", "CG"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "HD1", "CD1", "CG", "CD2"  ) );

			// Impropers
			dihedrallist.push_back( createDihedral( rsd,    "CG",  "CE2", "CD2", "HD2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CD2",  "CZ", "CE2", "HE2"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CZ",  "CD2", "CE2", "HE2"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CE1",  "CE2", "CZ", "OH"     ) );
			dihedrallist.push_back( createDihedral( rsd,    "CD1",  "CZ", "CE1", "HE1"    ) );
			dihedrallist.push_back( createDihedral( rsd,    "CG",  "CE1", "CD1", "HD1"   ) );
			dihedrallist.push_back( createDihedral( rsd,    "CD1",  "CD2", "CG", "CB"     ) );
		}
	}
}


void MolecularDynamics::setDihedralDerivatives( ){
	using namespace core;
	for ( Size i = 1; i <= dihedrallist.size(); i ++ ) {

		// work out which atoms are involved !
		int no1 = dihedrallist[i].index1;
		int no2 = dihedrallist[i].index2;
		int no3 = dihedrallist[i].index3;
		int no4 = dihedrallist[i].index4;

		//std::cout << no1 << "  " << no2 << "  " << no3 << "  " << no4 << "  "
		//          << std::endl;

		//float  epot_dihedral_this=0;
		float  phi, sin_phi, cos_phi;
		core::Vector vti_vta, vta_vtb, vtb_vtj;
		core::Vector nrml1, nrml2, nrml3;
		float  inv_nrml1_mag, inv_nrml2_mag, inv_nrml3_mag;

		core::Vector dcosdnrml1;
		core::Vector dcosdnrml2;
		core::Vector dsindnrml3;
		core::Vector dsindnrml2;

		vti_vta.x() = cartom[no1].position.x() - cartom[no2].position.x();
		vti_vta.y() = cartom[no1].position.y() - cartom[no2].position.y();
		vti_vta.z() = cartom[no1].position.z() - cartom[no2].position.z();

		vta_vtb.x() = cartom[no2].position.x() - cartom[no3].position.x();
		vta_vtb.y() = cartom[no2].position.y() - cartom[no3].position.y();
		vta_vtb.z() = cartom[no2].position.z() - cartom[no3].position.z();

		vtb_vtj.x() = cartom[no3].position.x() - cartom[no4].position.x();
		vtb_vtj.y() = cartom[no3].position.y() - cartom[no4].position.y();
		vtb_vtj.z() = cartom[no3].position.z() - cartom[no4].position.z();

		nrml1.x() = (vti_vta.y() * vta_vtb.z() - vti_vta.z() * vta_vtb.y());
		nrml1.y() = (vti_vta.z() * vta_vtb.x() - vti_vta.x() * vta_vtb.z());
		nrml1.z() = (vti_vta.x() * vta_vtb.y() - vti_vta.y() * vta_vtb.x());

		nrml2.x() = (vta_vtb.y() * vtb_vtj.z() - vta_vtb.z() * vtb_vtj.y());
		nrml2.y() = (vta_vtb.z() * vtb_vtj.x() - vta_vtb.x() * vtb_vtj.z());
		nrml2.z() = (vta_vtb.x() * vtb_vtj.y() - vta_vtb.y() * vtb_vtj.x());

		nrml3.x() = (vta_vtb.y() * nrml1.z() - vta_vtb.z() * nrml1.y());
		nrml3.y() = (vta_vtb.z() * nrml1.x() - vta_vtb.x() * nrml1.z());
		nrml3.z() = (vta_vtb.x() * nrml1.y() - vta_vtb.y() * nrml1.x());

		inv_nrml1_mag = 1.0 / sqrt(sqr(nrml1.x()) + sqr(nrml1.y()) + sqr(nrml1.z()));
		inv_nrml2_mag = 1.0 / sqrt(sqr(nrml2.x()) + sqr(nrml2.y()) + sqr(nrml2.z()));
		inv_nrml3_mag = 1.0 / sqrt(sqr(nrml3.x()) + sqr(nrml3.y()) + sqr(nrml3.z()));

		cos_phi = (nrml1.x() * nrml2.x() + nrml1.y() * nrml2.y() + nrml1.z() * nrml2.z()) * inv_nrml1_mag * inv_nrml2_mag;
		sin_phi = (nrml3.x() * nrml2.x() + nrml3.y() * nrml2.y() + nrml3.z() * nrml2.z()) * inv_nrml3_mag * inv_nrml2_mag;

		nrml2.x() *= inv_nrml2_mag;
		nrml2.y() *= inv_nrml2_mag;
		nrml2.z() *= inv_nrml2_mag;

		phi = -atan2(sin_phi, cos_phi);
		dihedrallist[i].angle = phi;
	}
}


void MolecularDynamics::doBondDerivatives( float &totalepot )
{
	using namespace core;

	float k = 600;
	for ( Size i = 1; i <= bondlist.size(); i ++ ) {

		core::Vector f2 = pose->xyz( bondlist[i].atom_id_1 )   -    pose->xyz( bondlist[i].atom_id_2 );
		float length =  pose->xyz( bondlist[i].atom_id_1 ).distance( pose->xyz( bondlist[i].atom_id_2 ) );

		totalepot +=  0.5 * k * sqr(length - bondlist[i].length);

		float deriv =  k * (length - bondlist[i].length);

		f2 *= deriv / length;

		cartom[ bondlist[i].index1].force += f2;
		cartom[ bondlist[i].index2].force -= f2;
	}
}


void MolecularDynamics::doAngleDerivatives( float &totalepot ) {
	using namespace core;

	float k = 600;
	for ( Size i = 1; i <= anglelist.size(); i ++ ) {

		core::Vector f2 = pose->xyz( anglelist[i].atom_id_1 )   -    pose->xyz( anglelist[i].atom_id_3 );
		float length =  pose->xyz( anglelist[i].atom_id_1 ).distance( pose->xyz( anglelist[i].atom_id_3 ) );

		totalepot +=  0.5 * k * sqr(length - anglelist[i].length);

		float deriv =  k * (length - anglelist[i].length);

		f2 *= deriv / length;

		cartom[ anglelist[i].index1 ].force += f2;
		cartom[ anglelist[i].index3 ].force -= f2;
	}
}


void MolecularDynamics::doDihedralDerivatives( float & totalepot ) {
	using namespace core;

	//  float totene_b4 = totalepot;

	float k = 200;
	for ( Size i = 1; i <= dihedrallist.size(); i ++ ) {

		// work out which atoms are involved !

		int no1 = dihedrallist[i].index1;
		int no2 = dihedrallist[i].index2;
		int no3 = dihedrallist[i].index3;
		int no4 = dihedrallist[i].index4;

		//float  epot_dihedral_this=0;
		float  phi, sin_phi, cos_phi;
		core::Vector vti_vta, vta_vtb, vtb_vtj;
		core::Vector nrml1, nrml2, nrml3;
		float  inv_nrml1_mag, inv_nrml2_mag, inv_nrml3_mag;

		core::Vector dcosdnrml1(0,0,0);
		core::Vector dcosdnrml2(0,0,0);
		core::Vector dsindnrml3(0,0,0);
		core::Vector dsindnrml2(0,0,0);
		core::Vector f, fi, fab, fj;


		vti_vta = cartom[no1].position - cartom[no2].position;
		vta_vtb = cartom[no2].position - cartom[no3].position;
		vtb_vtj = cartom[no3].position - cartom[no4].position;

		fi = fab = fj = 0;

		nrml1 = cross( vti_vta, vta_vtb );
		nrml2 = cross( vta_vtb, vtb_vtj );
		nrml3 = cross( vta_vtb, nrml1   );

		inv_nrml1_mag = 1.0 / nrml1.magnitude();
		inv_nrml2_mag = 1.0 / nrml2.magnitude();
		inv_nrml3_mag = 1.0 / nrml3.magnitude();

		cos_phi = dot( nrml1, nrml2 ) * inv_nrml1_mag * inv_nrml2_mag;
		sin_phi = dot( nrml3, nrml2 ) * inv_nrml3_mag * inv_nrml2_mag;

		nrml2 *= inv_nrml2_mag;
		phi = -atan2(sin_phi, cos_phi);

		if ( fabs( sin_phi ) > 0.1 ) {
			nrml1 *= inv_nrml1_mag;
			dcosdnrml1 = inv_nrml1_mag * ( nrml1 * cos_phi - nrml2 );
			dcosdnrml2 = inv_nrml2_mag * ( nrml2 * cos_phi - nrml1 );
		} else {
			nrml3 *= inv_nrml3_mag;
			dsindnrml3 = inv_nrml3_mag * ( nrml3 * sin_phi - nrml2 );
			dsindnrml2 = inv_nrml2_mag * ( nrml2 * sin_phi - nrml3 );
		}

		double diff = phi - dihedrallist[i].angle;
		if ( diff < -numeric::NumericTraits< double >::pi() )     diff += 2*numeric::NumericTraits< double >::pi();
		if ( diff < -numeric::NumericTraits< double >::pi() )     diff += 2*numeric::NumericTraits< double >::pi();
		if ( diff >  numeric::NumericTraits< double >::pi() )     diff -= 2*numeric::NumericTraits< double >::pi();
		if ( diff >  numeric::NumericTraits< double >::pi() )     diff -= 2*numeric::NumericTraits< double >::pi();
		totalepot += k * sqr(diff);
		double deriv = -2.0 * k * diff;
		//std::cout << "BAH! " << phi << "  " <<  dihedrallist[i].angle << "  " << diff << "  " << deriv << "  " << numeric::NumericTraits< double >::pi() << std::endl;
		// forces
		if ( fabs( sin_phi ) > 0.1 ) {
			deriv /= sin_phi;
			fi += deriv * cross( vta_vtb, dcosdnrml1 );
			fj += deriv * cross( vta_vtb, dcosdnrml2 );
			fab += deriv * ( cross( vti_vta, dcosdnrml1 ) + cross( vtb_vtj, dcosdnrml2 ) );
		} else {
			deriv /= -cos_phi;
			fi += deriv * update_operation( vta_vtb, dsindnrml3 );
			fj += deriv * cross( dsindnrml2, vta_vtb );
			fab += deriv * update_5way_operation( vta_vtb, vti_vta, dsindnrml3, dsindnrml2, vtb_vtj );
		}

		//  fi.mul(PhysicsConst::invAngstrom);
		//  fab.mul(PhysicsConst::invAngstrom);
		//  fj.mul(PhysicsConst::invAngstrom);

		// add to cartesian derivatives

		cartom[ no1 ].force += fi;
		cartom[ no2 ].force += fab - fi;
		cartom[ no3 ].force += fj - fab;
		cartom[ no4 ].force -= fj;
	}
}

void MolecularDynamics::createCartesianDerivatives( core::scoring::ScoreFunction const & scorefxn )
{
	using namespace core;
	//std::cout << "setup_for_scoring" << std::endl;
	Size const nres( pose->size() );

	scorefxn.setup_for_derivatives( *pose);
	cartom.clear();
	for ( Size ir = 1; ir <= nres; ++ir ) {
		// get the appropriate residue from the pose->
		conformation::Residue const & rsd( pose->residue( ir ) );

		//std::cout << "Res: " << ir << "  " <<  rsd.natoms() << std::endl;

		// iterate across neighbors within 12 angstroms
		for ( Size i = 1; i <= rsd.natoms(); i++ ) {
			id::AtomID atom_id(  i, ir );

			CartesianAtom atom;
			core::Vector F1(0,0,0),F2(0,0,0);
			scorefxn.eval_npd_atom_derivative( atom_id, *pose, min_map.domain_map(), F1, F2 );
			atom.index = i;
			atom.res = ir;
			atom.atom_id = atom_id;
			atom.position = core::Vector(0,0,0); // RM: Added to avoid uninitialized variable error
			atom.velocity = core::Vector(0,0,0); // RM: Added to avoid uninitialized variable error
			atom.force = F2;
			atom.old_force = core::Vector(0,0,0);
			atom.old_velocity = core::Vector(0,0,0);
			atom.old_position = core::Vector(0,0,0);
			atom.mass = 12.001;  // assume carbon for now, we can be more accurate later. RN: Added to avoid uninitialized variable error
			cartom.push_back(atom);

		}
	}
}


void MolecularDynamics::setInitialSpeeds(double tgtTemp )
{
	using namespace core;
	//double sigma;
	//int natom = cartom.size();

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		double sigma = sqrt( 1.38065E-23 * tgtTemp / (cartom[i].mass * 1.66053878E-27) );
		cartom[i].velocity.x() = numeric::random::gaussian() * sigma;
		cartom[i].velocity.y() = numeric::random::gaussian() * sigma;
		cartom[i].velocity.z() = numeric::random::gaussian() * sigma;
	}
}


void MolecularDynamics::calcKineticEnergy(
	float &ekin,
	float &Temp
) {
	using namespace core;
	ekin = 0;

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		ekin += (0.5 * cartom[i].mass *  1.66053878E-27 *  inner_product(cartom[i].velocity,cartom[i].velocity));
	}

	// Calculate Temperature from K = 2NkBT/2
	Temp = 2.0 * ekin / (3 * cartom.size() * 8.31 / 6.022E23 );  //

}


void MolecularDynamics::applyForces_BeeMan(
	float &kin,
	float &temp
) {
	using namespace core;

	double t;
	double tinvmass;
	core::Vector dr, dv;
	core::Vector t2a;

	t = 1E-15;

	float forcemul =  -( 1E10 / 6.0221418E23 * 4184.0);

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		tinvmass = t / (6.0 * (cartom[i].mass *  1.66053878E-27   ) ); // t div 6m

		// finish change in velocity using a(n+1)
		dv = ((2.0) * cartom[i].force * forcemul ) * tinvmass;
		cartom[i].velocity += dv;
	}

	calcKineticEnergy( kin, temp ); //calculate instanteneous temperature & pressure

	kin *= 6.0221418E23 /  4184.0;

	//applyThermostat();   //adjust velocities according to Thermostat

	// predict new positions
	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		tinvmass = t / (6.0 * cartom[i].mass *  1.66053878E-27 ); // t div 6m

		// change in position
		dr = t * (cartom[i].velocity + (4.0 * cartom[i].force - cartom[i].old_force )  * forcemul* tinvmass);
		dr /= 1E-10;
		dv = (5.0 * cartom[i].force - cartom[i].old_force )  * forcemul* tinvmass;

		//update positions and speeds.
		//std::cout << dr.x() << "  " << dr.z() << "  "<< dr.z() << "  " << std::endl;
		cartom[i].position += dr;
		cartom[i].velocity += dv;

		// save current force as next Step's old force
		cartom[i].old_force = cartom[i].force;
	}
}


void MolecularDynamics::applyForces_LangevinIntegration(
	double T,
	float &kin,
	float &temp
) {
	using namespace core;
	core::Vector dr, dv;
	double t = 1E-15;
	double tinvmass;
	float forcemul =  -( 1E10 / 6.0221418E23 * 4184.0);

	double gamma = 10E12;
	double gammat = gamma * t;
	double gt = gammat; // these are the consecutive powers for the series expnsions
	double gt2 = gt * gt;
	double gt3 = gt * gt2;
	double gt4 = gt2 * gt2;
	double gt5 = gt2 * gt3;
	double gt6 = gt3 * gt3;
	double gt7 = gt3 * gt4;
	double gt8 = gt4 * gt4;
	double gt9 = gt4 * gt5;

	double c0; // langevin coefficients
	double c1;
	double c2;

	//double egammat;
	double varv;
	double varr;
	double crv;

	core::Vector deltav;
	core::Vector deltar;

	if ( gammat > 0.05 ) {
		// when the friction coeffcient is relatively large, run langevin dynamics
		// the parameters must be calculated using exponential functions
		c0 = exp(-gammat); // langevin coefficients
		c1 = (1 - c0) / gammat;
		c2 = (1 - c1) / gammat;

		double egammat = exp(-gammat);
		varv = (1.0 - sqr(egammat));
		varr = (2.0 * gammat - 3.0 + (4.0 - egammat) * egammat);
		crv = sqr(1.0 - egammat) / (sqrt(varv * varr));
	} else {
		// when the friction coeffcient is in the midrange we can approximate the paramters
		// by taylor expansion which is a little faster
		c0 =
			1 - gt + 0.5 * gt2 - gt3 / 6.0 + gt4 / 24.0 - gt5 / 120.0 + gt6 / 720.0 - gt7 / 5040.0 + gt8 / 40320.0 -
			gt9 / 362880.0;
		c1 =
			1 - (0.5 * gt - gt2 / 6.0 + gt3 / 24.0 - gt4 / 120.0 + gt5 / 720.0 - gt6 / 5040.0 + gt7 / 40320.0 -
			gt8 / 362880.0);
		c2 =
			0.5 - gt / 6.0 + gt2 / 24.0 - gt3 / 120.0 + gt4 / 720.0 - gt5 / 5040.0 + gt6 / 40320.0 - gt7 / 362880.0;

		varr = 2.0 * gt3 / 3.0 - gt4 / 2.0 + 7.0 * gt5 / 30.0 - gt6 / 12.0
			+ 31.0 * gt7 / 1260.0 - gt8 / 160.0 + 127.0 * gt9 / 90720.0;
		varv = 2.0 * gt - 2.0 * gt2 + 4.0 * gt3 / 3.0 - 2.0 * gt4 / 3.0 + 4.0 * gt5 / 15.0
			- 4.0 * gt6 / 45.0 + 8.0 * gt7 / 315.0 - 2.0 * gt8 / 315.0 + 4.0 * gt9 / 2835.0;
		crv = sqrt(3.0) * (0.5 - 3.0 * gt / 16.0 - 17.0 * gt2 / 1280.0 + 17.0 * gt3 / 6144.0
			+ 40967.0 * gt4 / 34406400.0 - 57203.0 * gt5 / 275251200.0 -
			1429487.0 * gt6 / 13212057600.0);
	}

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		tinvmass = t / (cartom[i].mass  *  1.66053878E-27); // t div m

		// finish change in velocity using a(n+1)
		dv = (1.0 - c0 / c1) * cartom[i].force * forcemul * tinvmass / gammat;
		cartom[i].velocity += dv;
	}

	// The velocities are now "real", i.e. finished and thus now is the time
	// to calculate the kinetic energy and apply any barostat
	// note no other thermostats are called here because
	// langevin dynamics regultes the temperature itself - thermostat
	// is built-in so to speak
	calcKineticEnergy( kin, temp);
	kin *= 6.0221418E23 /  4184.0;

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		tinvmass = t / (cartom[i].mass  *  1.66053878E-27); // t div m

		double kTdivm =  1.38065E-23 * T / (cartom[i].mass  * 1.66053878E-27);
		double sigmav = sqrt(kTdivm * varv);
		double sigmar = sqrt(kTdivm * varr) / gamma;

		double n1x = numeric::random::gaussian();
		double n2y = numeric::random::gaussian();
		double n2x = numeric::random::gaussian();
		double n1z = numeric::random::gaussian();
		double n1y = numeric::random::gaussian();
		double n2z = numeric::random::gaussian();

		deltav = core::Vector(sigmav * (crv * n1x + sqrt(1.0 - sqr(crv)) * n2x),
			sigmav * (crv * n1y + sqrt(1.0 - sqr(crv)) * n2y),
			sigmav * (crv * n1z + sqrt(1.0 - sqr(crv)) * n2z));
		deltar = core::Vector(sigmar * n1x, sigmar * n1y, sigmar * n1z);

		// change in position
		dr = t * (c1 * cartom[i].velocity + c2 * cartom[i].force * forcemul * tinvmass) + deltar;
		dr *= 1E10;

		// finish change in velocity using a(n+1)
		dv = (c0 * c2) * cartom[i].force * forcemul * tinvmass / c1 + deltav;

		cartom[i].position += dr;
		cartom[i].velocity *= c0;
		cartom[i].velocity += dv;
	}
}


void MolecularDynamics::applyForces_ConjugateGradient(
	int Step,
	float &current_energy,
	float &m_OldEnergy
) {
	using namespace core;

	core::Vector f;
	static double m_StepMultiplier = 1.0;
	float StepSize = 0.000001;
	float forcemul = 1; // ( 1E10 / 6.0221418E23 * 4184.0);

	//std::cout << "-----> " << m_StepMultiplier << std::endl;
	if ( Step != 0 ) {
		if ( (current_energy < m_OldEnergy) ) { // energy did go down on last Step

			//if( !backstep ){
			double  fmagold, fmagnew, gamma;
			m_StepMultiplier *= 1.2;
			//std::cout << "ENERGY_IMPROVEMENT\n";
			fmagold = 0;
			fmagnew = 0;
			for ( Size i = 1; i <= cartom.size(); i ++ ) {   // calculate the inner dots of the current force
				fmagold += cartom[i].old_force.dot( cartom[i].old_force );
				core::Vector ftemp =  cartom[i].old_force - cartom[i].force;
				fmagnew += ftemp.dot( cartom[i].force );
			}

			//b = ( (gnew - gprev)*(new)   /   prev^2

			if ( fmagold != 0 ) {
				gamma = fmagnew / fmagold;
			} else {
				gamma = 0;
			}

			if ( gamma < 0 ) gamma = 0;
			std::cout << "gamma: " << gamma << std::endl;

			//gamma = 0;
			for ( Size i = 1; i <= cartom.size(); i ++ ) {   // save positions & forces
				cartom[i].velocity    = cartom[i].old_velocity; // set to old direction
				cartom[i].velocity   *= gamma;                 // multiply the old direction with gamma
				cartom[i].velocity   -= cartom[i].force * forcemul;      // and add current force to make current direction
				cartom[i].velocity   *= 1;
				cartom[i].old_position = cartom[i].position;  // save position
				cartom[i].old_force    = cartom[i].force * forcemul;     // save old forces
				cartom[i].old_velocity = cartom[i].velocity;; // save old directions
			}
			m_OldEnergy = current_energy;
			//backstep = false;  // set but never used ~Labonte
			//} else {
			//  std::cout << "BACKSTEP\n";

			//}
		} else { // we went up the other side
			//std::cout << "ENERGY_UP\n" << current_energy << "   " << m_OldEnergy << "  ";

			m_StepMultiplier *= 0.5 / 1.2;
			for ( Size i = 1; i <= cartom.size(); i ++ ) {
				cartom[i].position = cartom[i].old_position; // restore old position
				//cartom[i].force    = cartom[i].old_force / forcemul;   // restore old forces
			}
			//m_OldEnergy += 1; // no idea what this is !?

			//backstep = true;  // set but never used ~Labonte
			//return;

		}
	} else { // this block is for Step == 0 - its just a standard SD Step
		for ( Size i = 1; i < cartom.size(); i++ ) {
			cartom[i].old_position    =  cartom[i].position; // save position (old position = current position)
			cartom[i].old_force       =  cartom[i].force * forcemul; // save old forces
			cartom[i].old_velocity    =  -cartom[i].force * forcemul; // save old directions, equal to old force
			cartom[i].velocity        =  -cartom[i].force * forcemul; // save old directions, equal to old force
		}
		m_OldEnergy = current_energy;
	}

	for ( Size i = 1; i <= cartom.size(); i ++ ) {
		f = cartom[i].velocity;
		f *= StepSize * m_StepMultiplier;
		//if(  sqr(f.x()) + sqr(f.y()) + sqr (f.z())  > 0.2 )
		//   std::cout << "DISP:  " << i << f.x() << "  " << f.y() << "  " << f.z() << std::endl;
		cartom[i].position += f;
	}
}


void MolecularDynamics::doMinimising(
	core::scoring::ScoreFunction const & scorefxn
) {
	using namespace core;

	const int Steps = 400;
	int Step = 0;
	float kin = 0;
	int mcount = 1;

	std::ofstream pdbfile;
	std::string filename ( "min.pdb" );
	pdbfile.open( filename.c_str() , std::ios::out );

	float m_OldEnergy(0.0);

	for ( Step = 0; Step < Steps; Step ++ ) {
		pose::Pose pose2( *pose);
		float current_energy =  scorefxn( pose2 );
		float pot = current_energy;
		float cov_epot = 0;

		zeroForces( );
		getCartesianDerivatives( scorefxn );
		doBondDerivatives( cov_epot );
		doAngleDerivatives( cov_epot );
		doDihedralDerivatives( cov_epot );
		current_energy += cov_epot;

		std::cout << Step << "  " << pot << "  " << cov_epot << "  " << kin << "  " << pot+kin << std::endl;

		if ( Step == 0 ) m_OldEnergy = current_energy;
		applyForces_ConjugateGradient( Step, current_energy, m_OldEnergy );

		setPosePositionsFromCartesian( );

		//if( ( Step % 50 ) == 0 ){
		pdbfile << "MODEL    " << I(5, mcount ) << std::endl;
		pose->dump_pdb( pdbfile );
		pdbfile << "ENDMDL" << std::endl;
		mcount++;
		//}
		pdbfile.flush();
	}

	pdbfile.close();
}


void MolecularDynamics::doMD( core::scoring::ScoreFunction const & scorefxn,
	int Steps,
	float startTemp,
	float endTemp
) {
	using namespace core;
	setInitialSpeeds(startTemp);

	int Step = 0;
	float kin = 0;
	float temp = 0;
	float pot = 0;
	int mcount = 1;

	std::ofstream pdbfile;
	std::string filename ( "md." + right_string_of(startTemp,4,'0') + ".pdb" );
	pdbfile.open( filename.c_str() , std::ios::out );

	float cov_epot=0;
	float TargetTemp=startTemp;
	for ( Step = 0; Step < Steps; Step ++ ) {
		pose::Pose pose2(*pose);
		pot +=  scorefxn( pose2 );
		if ( ( Step % 20 ) == 0 ) {
			std::cout << Step << "  " << pot << "  " << cov_epot << "  " << kin << "  " << pot+kin+cov_epot << "   " << temp << "(" << TargetTemp << ")" << std::endl;
		}
		getCartesianDerivatives( scorefxn );
		pot = 0;
		cov_epot = 0;
		doBondDerivatives(cov_epot );
		doAngleDerivatives( cov_epot);
		doDihedralDerivatives( cov_epot );

		float ratio = (float)Step/(float)Steps;
		TargetTemp = std::max(0.01,      ratio*endTemp + (1.0-ratio)*startTemp );
		applyForces_LangevinIntegration( TargetTemp , kin, temp );

		setPosePositionsFromCartesian( );

		if ( ( Step % 100 ) == 0 ) {
			pdbfile << "MODEL    " << I(5, mcount ) << std::endl;
			pose->dump_pdb( pdbfile );
			pdbfile << "ENDMDL" << std::endl;
			mcount++;
		}
		pdbfile.flush();
	}

	pdbfile.close();
}


void MolecularDynamics::testCartesianDerivatives( core::scoring::ScoreFunction const & scorefxn )
{
	using namespace core;

	// score it !
	Real start_score = scorefxn( *pose );
	scorefxn.show(std::cout, *pose);
	std::cout << "STARTSCORE: " << start_score << std::endl;
	/*
	Residue rsd_mod =  pose->residue(4);

	std::cout << rsd_mod.xyz(4).z() << std::endl;
	// move an atom
	Vector test = rsd_mod.xyz(4).xyz();
	Vector dd(0,0,0.1); yy
	Vector newvec;

	newvec = rsd_mod.xyz(4) + dd;
	pose->residue(4).atom(3).set_xyz(4, newvec);

	pose->set_residue
	pose->set_xyz( id, pose->xyz( id ) + dx );

	scorefxn( *pose );
	scorefxn.show(std::cout, *pose);
	std::cout << pose->residue(4).atom(3).xyz().z() << std::endl;
	*/
	const float dd = 0.00001;
	Vector dx(dd,0,0);
	Vector dy(0,dd,0);
	Vector dz(0,0,dd);

	utility::vector1< CartesianAtom > cartom;
	utility::vector1< Vector > numeriv;

	std::cout << "Calculating derivatives \n";

	createCartesianArray( );
	getCartesianDerivatives( scorefxn );

	//int cend = clock();
	pose->energies().reset_nblist();
	std::cout << "STARTSCORE: --------- " << std::endl;
	start_score = scorefxn( *pose );
	scorefxn.show(std::cout, *pose);

	//std::cout << "setup_for_scoring" << std::endl;
	Size const nres( pose->size() );

	for ( Size ir = 1; ir <= nres; ++ir ) {
		// get the appropriate residue from the pose->
		conformation::Residue const & rsd( pose->residue( ir ) );
		std::cout << "IR: " << ir << std::endl;
		// iterate across neighbors within 12 angstroms
		for ( Size i = 1; i <= rsd.natoms(); i++ ) {
			id::AtomID atom_id( i, ir );

			pose->energies().clear();
			start_score = scorefxn( *pose );

			Real new_score;
			Real deriv;
			Vector deriv_vector;
			Vector safepos =  pose->xyz( atom_id );
			pose->set_xyz( atom_id, safepos + dx );
			pose->energies().clear();
			new_score = scorefxn( *pose );
			deriv = (new_score - start_score)/dd;
			deriv_vector.x(deriv);

			pose->set_xyz( atom_id, safepos + dy );
			pose->energies().clear();
			new_score = scorefxn( *pose );
			deriv = (new_score - start_score)/dd;
			deriv_vector.y(deriv);

			pose->set_xyz( atom_id, safepos + dz );
			pose->energies().clear();
			new_score = scorefxn( *pose );
			deriv = (new_score - start_score)/dd;
			deriv_vector.z(deriv);
			pose->set_xyz( atom_id, safepos );

			numeriv.push_back( deriv_vector );
		}
	}


	for ( Size i = 1; i < cartom.size(); i++ ) {

		if (  ( fabs( cartom[i].force.x()  -  numeriv[i].x()     ) > 0.1 ) ||
				( fabs( cartom[i].force.y()  -  numeriv[i].y()     ) > 0.1 ) ||
				( fabs( cartom[i].force.z()  -  numeriv[i].z()     ) > 0.1 )  ) {
			std::cout << "DRV: "
				<< cartom[i].res   << "  "
				<< cartom[i].index  << "   "
				<< cartom[i].force.x() << "  "
				<< cartom[i].force.y() << "  "
				<< cartom[i].force.z() << "  "
				<< numeriv[i].x() << "  "
				<< numeriv[i].y() << "  "
				<< numeriv[i].z() << "  "
				<< std::endl;
		} else {
			std::cout << "DRVT: "
				<< cartom[i].res   << "  "
				<< cartom[i].index  << "   "
				<< cartom[i].force.x() << "  "
				<< cartom[i].force.y() << "  "
				<< cartom[i].force.z() << "  "
				<< numeriv[i].x() << "  "
				<< numeriv[i].y() << "  "
				<< numeriv[i].z() << "  "
				<< std::endl;
		}
	}
}


} // namespace optimization
} // namespace core
