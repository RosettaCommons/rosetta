// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/power_diagram/PowerDiagram.cc
/// @brief  Class for constructing and querying power diagrams
/// @author Jim Havranek (havranek@biochem.wustl.edu)

// Unit Headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/power_diagram/PowerDiagram.hh>

#include <boost/bind.hpp>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

// Project Headers
// Numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/constants.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <ObjexxFCL/format.hh>

#include <stack>
#include <iterator>
#include <functional>

static basic::Tracer TR( "core.scoring.power_diagram.PowerDiagram" );

using namespace ObjexxFCL::format;
using namespace core;

namespace core {
namespace scoring {
namespace power_diagram {

bool
matches_ptr( PDvertexOP elem, PDvertex * test ){
	return( elem.get() == test );
}


Real power_distance( Vector const & pt, PDsphere const * sph )
{
	return ( pt.distance_squared( sph->xyz() ) - sph->rad2() );
}


PowerDiagram::PowerDiagram(
	core::pose::Pose & pose
)
{
	construct_from_pose( pose );
}

void
PowerDiagram::construct_from_pose(
	core::pose::Pose & pose
)
{

	// Make sure all old info is wiped clean
	clear();

	Size total_res( pose.size() );

	Real max_radius( 0.0 );
	Vector min( pose.residue(1).xyz(1) );
	Vector max( pose.residue(1).xyz(1) );

	for ( Size ires = 1 ; ires <= total_res ; ++ires ) {
		Size const total_atm( pose.residue( ires ).natoms() );
		for ( Size iatm = 1 ; iatm <= total_atm ; ++iatm ) {
			Size const RADIUS_INDEX( pose.residue(1).atom_type_set().extra_parameter_index( "GK_RADIUS" ) );
			Real const Rwater( 1.4 );
			Real this_rad = pose.residue(ires).atom_type(iatm).extra_parameter( RADIUS_INDEX ) + Rwater;
			if ( this_rad > max_radius ) max_radius = this_rad;
			if ( pose.residue(ires).xyz(iatm).x() < min.x() ) {
				min.x() = pose.residue(ires).xyz(iatm).x();
			}
			if ( pose.residue(ires).xyz(iatm).y() < min.y() ) {
				min.y() = pose.residue(ires).xyz(iatm).y();
			}
			if ( pose.residue(ires).xyz(iatm).z() < min.z() ) {
				min.z() = pose.residue(ires).xyz(iatm).z();
			}
			if ( pose.residue(ires).xyz(iatm).x() > max.x() ) {
				max.x() = pose.residue(ires).xyz(iatm).x();
			}
			if ( pose.residue(ires).xyz(iatm).y() > max.y() ) {
				max.y() = pose.residue(ires).xyz(iatm).y();
			}
			if ( pose.residue(ires).xyz(iatm).z() > max.z() ) {
				max.z() = pose.residue(ires).xyz(iatm).z();
			}
		}
	}

	min -= 3.0*max_radius;
	max += 3.0*max_radius;
	// TR << "Max radius " << max_radius << std::endl;
	// TR << "box min " << min << std::endl;
	// TR << "box max " << max << std::endl;

	bool const use_virtual_box( true );

	if ( use_virtual_box ) {
		Vector average = 0.5*( max + min );

		// Make the six virtual spheres to confine the power diagram
		Vector virtual1 = average;
		virtual1.x() = max.x();
		Vector virtual2 = average;
		virtual2.y() = max.y();
		Vector virtual3 = average;
		virtual3.z() = max.z();
		Vector virtual4 = average;
		virtual4.x() = min.x();
		Vector virtual5 = average;
		virtual5.y() = min.y();
		Vector virtual6 = average;
		virtual6.z() = min.z();
		Real dummy_rad = 0.5*max_radius;
		PDsphereOP vsph1( make_new_sphere( virtual1, dummy_rad ) );
		PDsphereOP vsph2( make_new_sphere( virtual2, dummy_rad ) );
		PDsphereOP vsph3( make_new_sphere( virtual3, dummy_rad ) );
		PDsphereOP vsph4( make_new_sphere( virtual4, dummy_rad ) );
		PDsphereOP vsph5( make_new_sphere( virtual5, dummy_rad ) );
		PDsphereOP vsph6( make_new_sphere( virtual6, dummy_rad ) );

		spheres_.push_back( vsph1 );
		spheres_.push_back( vsph2 );
		spheres_.push_back( vsph3 );
		spheres_.push_back( vsph4 );
		make_initial_power_diagram();
		add_single_atom_to_power_diagram( vsph5 );
		spheres_.push_back( vsph5 );
		add_single_atom_to_power_diagram( vsph6 );
		spheres_.push_back( vsph6 );

	}

	// Resize lookup info based on pose info
	sphere_lookup_.resize( total_res );
	for ( Size ires = 1 ; ires <= total_res ; ++ires ) {
		sphere_lookup_[ ires ].resize( pose.residue( ires ).natoms() );
	}

	// Get index into the radii
	// Size const RADIUS_INDEX( pose.residue(1).atom_type_set().extra_parameter_index( "GK_RADIUS" ) );

	for ( Size ires = 1 ; ires <= total_res ; ++ires ) {
		Size const total_atm( pose.residue( ires ).natoms() );
		for ( Size iatm = 1 ; iatm <= total_atm ; ++iatm ) {
			PDsphereOP new_sph( make_new_sphere( pose, ires, iatm ) );
			if ( spheres_.size() < 4 && !use_virtual_box ) {
				store_new_sphere( new_sph );
				// If we have hit four spheres, we can bootstrap the starter power diagram
				if ( spheres_.size() == 4 ) {
					make_initial_power_diagram();
				}
			} else {
				//TR << "Adding atom for residue " << ires << " atomno " << iatm << std::endl;

				//    TR << "ATOMINFO " << new_sph->xyz()(1) << " " << new_sph->xyz()(2) << " " << new_sph->xyz()(3) << " " << new_sph->rad() << std::endl;

				add_single_atom_to_power_diagram( new_sph );
				store_new_sphere( new_sph );

				//print_vertices( finite_vertices_, infinite_vertices_ );

			}
		}
	}
	// Now get the cycles and associate them with the spheres

	for ( auto patom : spheres_ ) {

		std::list< PDinterOP > intersections( get_intersections_for_atom( patom.get() ) );

		for ( auto inter : intersections ) {
			find_common_intersection_atoms( inter );
		}
		patom->cycles() = get_cycles_from_intersections( intersections, patom.get() );

	}
}



void
PowerDiagram::clear()
{
	// reset vertex count
	vertex_count_ = 0;

	for ( auto this_sphere : spheres_ ) {
		this_sphere->vertices().clear();
		for ( auto cyc : this_sphere->cycles() ) {
			cyc.clear();
		}
	}

	for ( auto sl : sphere_lookup_ ) {
		sl.clear();
	}

	for ( auto fv : finite_vertices_ ) {
		fv->nonconst_partners().clear();
		fv->nonconst_generators().clear();
	}

	for ( auto iv : infinite_vertices_ ) {
		iv->nonconst_partners().clear();
		iv->nonconst_generators().clear();
	}

	// Clear the stored vertices
	spheres_.clear();
	finite_vertices_.clear();
	infinite_vertices_.clear();
	sphere_lookup_.clear();

}


Real
PowerDiagram::extract_sasa_for_atom( Size ires, Size iatm )
{

	Real surf_area( 0.0 );
	PDsphere* patom( sphere_lookup( ires, iatm ) );
	surf_area += get_sasa_from_cycles( patom->cycles(), patom );

	return surf_area;

}




PDsphereOP
PowerDiagram::make_new_sphere(
	core::pose::Pose & p,
	Size ires,
	Size iatm
)
{
	Size const RADIUS_INDEX( p.residue(1).atom_type_set().extra_parameter_index( "GK_RADIUS" ) );
	Real const Rwater( 1.4 );
	PDsphereOP new_sph( new PDsphere() );
	new_sph->nonconst_res() = ires;
	new_sph->nonconst_atom() = iatm;
	new_sph->nonconst_xyz() = p.residue(ires).xyz(iatm);
	new_sph->nonconst_rad() = p.residue(ires).atom_type(iatm).extra_parameter( RADIUS_INDEX ) + Rwater;
	new_sph->nonconst_rad2() = new_sph->rad()*new_sph->rad();

	return new_sph;
}

PDsphereOP
PowerDiagram::make_new_sphere(
	Vector & pos,
	Real radius
)
{
	PDsphereOP new_sph( new PDsphere() );
	new_sph->nonconst_res() = 9999;
	new_sph->nonconst_atom() = 9999;
	new_sph->nonconst_xyz() = pos;
	new_sph->nonconst_rad() = radius;
	new_sph->nonconst_rad2() = new_sph->rad()*new_sph->rad();

	return new_sph;
}

void
PowerDiagram::make_initial_power_diagram()
{
	// Make sure we have the required four spheres
	if ( spheres_.size() != 4 ) {
		TR << "Ill-formed call to make_initial_power_diagram():  exactly four spheres required!" << std::endl;
		std::exit(1);
	}

	// Each combo of three atoms defines an infinite vertex
	Vector dir1( (spheres_[1]->xyz() - spheres_[2]->xyz()).cross( spheres_[3]->xyz() - spheres_[2]->xyz() ).normalize() );
	Vector dir2( (spheres_[2]->xyz() - spheres_[3]->xyz()).cross( spheres_[4]->xyz() - spheres_[3]->xyz() ).normalize() );
	Vector dir3( (spheres_[3]->xyz() - spheres_[4]->xyz()).cross( spheres_[1]->xyz() - spheres_[4]->xyz() ).normalize() );
	Vector dir4( (spheres_[4]->xyz() - spheres_[1]->xyz()).cross( spheres_[2]->xyz() - spheres_[1]->xyz() ).normalize() );
	if ( dir1.dot( spheres_[4]->xyz() - spheres_[2]->xyz() ) > 0.0 ) dir1 *= -1.0;
	if ( dir2.dot( spheres_[1]->xyz() - spheres_[3]->xyz() ) > 0.0 ) dir2 *= -1.0;
	if ( dir3.dot( spheres_[2]->xyz() - spheres_[4]->xyz() ) > 0.0 ) dir3 *= -1.0;
	if ( dir4.dot( spheres_[3]->xyz() - spheres_[1]->xyz() ) > 0.0 ) dir4 *= -1.0;

	PDvertexOP vrt1( new PDvertex() );
	vrt1->set_finite( false );
	vrt1->set_live( true );
	vrt1->nonconst_id() = vertex_count_++;
	vrt1->nonconst_generators().push_back( spheres_[1].get() );
	vrt1->nonconst_generators().push_back( spheres_[2].get() );
	vrt1->nonconst_generators().push_back( spheres_[3].get() );
	std::sort( vrt1->nonconst_generators().begin(), vrt1->nonconst_generators().end() );
	vrt1->nonconst_direction() = dir1;
	link_vertex_to_generators( vrt1 );
	// TR << "Made initial infinite vertex with id " << vrt1->id() << std::endl;

	PDvertexOP vrt2( new PDvertex() );
	vrt2->set_finite( false );
	vrt2->set_live( true );
	vrt2->nonconst_id() = vertex_count_++;
	vrt2->nonconst_generators().push_back( spheres_[2].get() );
	vrt2->nonconst_generators().push_back( spheres_[3].get() );
	vrt2->nonconst_generators().push_back( spheres_[4].get() );
	std::sort( vrt2->nonconst_generators().begin(), vrt2->nonconst_generators().end() );
	vrt2->nonconst_direction() = dir2;
	link_vertex_to_generators( vrt2 );
	// TR << "Made initial infinite vertex with id " << vrt2->id() << std::endl;

	PDvertexOP vrt3( new PDvertex() );
	vrt3->set_finite( false );
	vrt3->set_live( true );
	vrt3->nonconst_id() = vertex_count_++;
	vrt3->nonconst_generators().push_back( spheres_[3].get() );
	vrt3->nonconst_generators().push_back( spheres_[4].get() );
	vrt3->nonconst_generators().push_back( spheres_[1].get() );
	std::sort( vrt3->nonconst_generators().begin(), vrt3->nonconst_generators().end() );
	vrt3->nonconst_direction() = dir3;
	link_vertex_to_generators( vrt3 );
	// TR << "Made initial infinite vertex with id " << vrt3->id() << std::endl;

	PDvertexOP vrt4( new PDvertex() );
	vrt4->set_finite( false );
	vrt4->set_live( true );
	vrt4->nonconst_id() = vertex_count_++;
	vrt4->nonconst_generators().push_back( spheres_[4].get() );
	vrt4->nonconst_generators().push_back( spheres_[1].get() );
	vrt4->nonconst_generators().push_back( spheres_[2].get() );
	std::sort( vrt4->nonconst_generators().begin(), vrt4->nonconst_generators().end() );
	vrt4->nonconst_direction() = dir4;
	link_vertex_to_generators( vrt4 );
	// TR << "Made initial infinite vertex with id " << vrt4->id() << std::endl;

	vrt1->nonconst_partners().push_back( vrt2.get() );
	vrt1->nonconst_partners().push_back( vrt3.get() );
	vrt1->nonconst_partners().push_back( vrt4.get() );
	vrt2->nonconst_partners().push_back( vrt3.get() );
	vrt2->nonconst_partners().push_back( vrt4.get() );
	vrt2->nonconst_partners().push_back( vrt1.get() );
	vrt3->nonconst_partners().push_back( vrt4.get() );
	vrt3->nonconst_partners().push_back( vrt1.get() );
	vrt3->nonconst_partners().push_back( vrt2.get() );
	vrt4->nonconst_partners().push_back( vrt1.get() );
	vrt4->nonconst_partners().push_back( vrt2.get() );
	vrt4->nonconst_partners().push_back( vrt3.get() );

	// Now get the initial finite vertex

	Vector const p_intersect( vertex_xyz_from_generators( spheres_[1].get(), spheres_[2].get(),
		spheres_[3].get(), spheres_[4].get() ) );

	PDvertexOP vrt5( new PDvertex() );
	vrt5->set_finite( true );
	vrt5->set_live( true );
	vrt5->nonconst_id() = vertex_count_++;
	vrt5->nonconst_power() = power_distance( p_intersect, spheres_[1].get() );
	vrt5->nonconst_generators().push_back( spheres_[1].get() );
	vrt5->nonconst_generators().push_back( spheres_[2].get() );
	vrt5->nonconst_generators().push_back( spheres_[3].get() );
	vrt5->nonconst_generators().push_back( spheres_[4].get() );
	std::sort( vrt5->nonconst_generators().begin(), vrt5->nonconst_generators().end() );
	vrt5->nonconst_xyz() = p_intersect;
	vrt5->nonconst_partners().push_back( vrt1.get() );
	vrt5->nonconst_partners().push_back( vrt2.get() );
	vrt5->nonconst_partners().push_back( vrt3.get() );
	vrt5->nonconst_partners().push_back( vrt4.get() );
	link_vertex_to_generators( vrt5 );
	// TR << "Made initial finite vertex with id " << vrt5->id() << std::endl;

	// Add this vertex to the infinite vertices
	vrt1->nonconst_partners().push_back( vrt5.get() );
	vrt2->nonconst_partners().push_back( vrt5.get() );
	vrt3->nonconst_partners().push_back( vrt5.get() );
	vrt4->nonconst_partners().push_back( vrt5.get() );

	// Set up the vertex lists
	infinite_vertices_.push_front( vrt1 );
	vrt1->my_itr() = infinite_vertices_.begin();
	infinite_vertices_.push_front( vrt2 );
	vrt2->my_itr() = infinite_vertices_.begin();
	infinite_vertices_.push_front( vrt3 );
	vrt3->my_itr() = infinite_vertices_.begin();
	infinite_vertices_.push_front( vrt4 );
	vrt4->my_itr() = infinite_vertices_.begin();

	finite_vertices_.push_front( vrt5 );
	vrt5->my_itr() = finite_vertices_.begin();
}




void
PowerDiagram::add_single_atom_to_power_diagram(
	PDsphereOP & new_sph
)
{
	////////////////////////////////////////////
	// Step 1 from Klenin paper: delete vertices
	////////////////////////////////////////////

	// TR << "*** STEP 1 ***" << std::endl;

	// Move finite vertices to be removed to the trash
	std::list< PDvertexOP > trash;

	bool const use_new_method( true );

	if ( use_new_method ) {

		/////////////////////////////////////////////////////////////
		/////// START NEW METHOD FOR DELETING VERTICES //////////////
		/////////////////////////////////////////////////////////////

		// Find a single vertex to be removed, then go from there
		// Optimization - start at end and follow lowest power

		// Get any finite vertex.
		PDvertex * srch_vrt( (*finite_vertices_.rbegin()).get() );
		PDvertex * next_vrt = nullptr;
		PDvertex * delete_vrt = nullptr;
		Real pd( power_distance( srch_vrt->xyz(), new_sph.get() ) );
		Real this_dist( pd - srch_vrt->power() );

		while ( this_dist > 0.0 && next_vrt != srch_vrt ) {
			if ( next_vrt != nullptr ) srch_vrt = next_vrt;
			next_vrt = find_next_vertex_with_smallest_dist( srch_vrt, new_sph, this_dist );
			//TR << "this_dist " << this_dist << " same_check " << (next_vrt == srch_vrt) << " null check " << (next_vrt == nullptr) << std::endl;
		}

		if ( this_dist < 0.0 ) {
			//TR << "Found finite vertex to delete!" << std::endl;
			if ( next_vrt == nullptr ) {
				delete_vrt = srch_vrt;
			} else {
				delete_vrt = next_vrt;
			}
		} else {
			//TR << "Stalled out at power_distance diff of " << this_dist << std::endl;
		}

		// If no finite vertices get deleted, try the infinite vertex partners
		if ( this_dist > 0.0 ) {

			auto itr = next_vrt->partners().begin();
			while ( itr != next_vrt->partners().end() && (*itr)->finite() ) ++itr;
			if ( itr != next_vrt->partners().end() ) {
				Real check_val( ( new_sph->xyz() - (*itr)->generators()[1]->xyz() ).dot( (*itr)->direction() ) );
				while ( ( check_val <= 0.0 ) && ( itr != next_vrt->partners().end() ) ) {
					++itr;
					if ( itr == next_vrt->partners().end() || (*itr)->finite() ) continue;
					check_val = ( new_sph->xyz() - (*itr)->generators()[1]->xyz() ).dot( (*itr)->direction() );
				}
				if ( itr != next_vrt->partners().end() ) {
					delete_vrt = (*itr);
					//TR << "Found infinite vertex to delete!" << std::endl;
				}
			}
		}

		//  if ( delete_vrt == nullptr ) {
		//   TR << "This sphere doesn't contribute!" << std::endl;
		//  }

		if ( delete_vrt != nullptr ) {
			// Ok, got one.  Now process this one and all its neighbors
			std::stack< PDvertex * > delete_me;
			delete_vrt->set_live( false );
			delete_me.push( delete_vrt );
			while ( !delete_me.empty() ) {
				// Get the next one to deldete
				PDvertex * this_vrt( delete_me.top() );
				delete_me.pop();

				// Check each partner vertex and push to deletion stack if necessary
				for ( auto partner_vrt : this_vrt->partners() ) {
					//    for ( Size i = 1 ; i <= this_vrt->partners().size() ; ++i ) {
					//     PDvertex * partner_vrt( this_vrt->nonconst_partners()[i] );
					//TR << "Checking partner at " << partner_vrt->xyz() << std::endl;
					// Different test for finite vs. infinite
					if ( partner_vrt->finite() ) {
						Real partner_pd( power_distance( partner_vrt->xyz(), new_sph.get() ) );
						if ( ( partner_pd < partner_vrt->power() ) && partner_vrt->live() ) {
							partner_vrt->set_live( false );
							delete_me.push( partner_vrt );
						}
					} else {
						Real const check_val( ( new_sph->xyz() - partner_vrt->generators()[1]->xyz() ).dot( partner_vrt->direction() ) );
						//TR << "Comparing to infinite vertex plane, check_val is " << check_val << std::endl;
						if ( ( check_val > 0.0 ) && partner_vrt->live() ) {
							//TR << "Removed an infinite vertex with id " << partner_vrt->id() << " at " << partner_vrt->xyz() << std::endl;
							partner_vrt->set_live( false );
							delete_me.push( partner_vrt );
						} else if ( check_val == 0.0 ) {
							// The new atom is on the generator plane.  Only remove
							// the infinite vertex if the lone finite partner vertex
							// will be deleted.
							//       TR << "New atom is on generator plane of infinite vertex" << std::endl;
							utility::vector1< PDvertex * >::const_iterator fv_itr( partner_vrt->partners().begin() );
							while ( !(*fv_itr)->finite() ) ++fv_itr;
							Real finite_partner_pd( power_distance( (*fv_itr)->xyz(), new_sph.get() ) );

							// Note we're just checking the finite partner here - we don't push it yet
							if ( finite_partner_pd < (*fv_itr)->power() && partner_vrt->live() ) {
								partner_vrt->set_live( false );
								delete_me.push( partner_vrt );
								//        TR << "Removed an infinite vertex after plane check because finite vertex partner is in trash" << std::endl;
							}
						}
					}
				}

				// Remove this vertex from the list and push onto the trash
				//    TR << "Deleting vertex at " << this_vrt->xyz() << std::endl;
				this_vrt->set_live( false );
				// This needs owning pointers, else the erase call below destroys the object
				trash.push_front( *(this_vrt->my_itr()) );

				if ( this_vrt->finite() ) {
					finite_vertices_.erase( this_vrt->my_itr() );
				} else {
					infinite_vertices_.erase( this_vrt->my_itr() );
				}
			}
		}

		//  TR << "Done with this part" << std::endl;

		/////////////////////////////////////////////////////////////
		/////// DONE NEW METHOD FOR DELETING VERTICES ///////////////
		/////////////////////////////////////////////////////////////

	}

	//TR << "Printing finite list" << std::endl;

	//for( std::list< PDvertexOP >::iterator fin_itr = finite_vertices_.begin() ; fin_itr != finite_vertices_.end() ; ++fin_itr ) {
	//  TR << "Finite vertex at " << (*fin_itr)->xyz() << std::endl;
	// }

	// TR << "Done finite list" << std::endl;
	//
	// TR << "Printing trash list" << std::endl;
	//
	// for( std::list< PDvertexOP >::iterator trash_itr = trash.begin() ; trash_itr != trash.end() ; ++trash_itr ) {
	//   TR << "Trashed vertex at " << (*trash_itr)->xyz() << std::endl;
	// }
	//
	// TR << "Done trash list" << std::endl;


	//////////////////////////////////////////////////////
	// Step 2 from Klenin paper:  Make new finite vertices
	//////////////////////////////////////////////////////

	// TR << "*** STEP 2 ***" << std::endl;

	std::list< PDvertexOP > new_vertices;

	for ( std::list< PDvertexOP >::iterator trash_itr( trash.begin() ) ;
			trash_itr != trash.end() ; ++trash_itr ) {

		// Remove this vertex from all of its generator atoms' lists
		for ( auto pvert : (*trash_itr)->nonconst_generators() ) {
			Size check_val( pvert->vertices().size() );

			pvert->vertices().remove( (*trash_itr).get() );

			debug_assert( pvert->vertices().size() == (check_val - 1)   );
			//   TR << "Atom res " << (*gen_itr)->res() << " atom " << (*gen_itr)->atom() << " lost a vertex" << std::endl;
		}

		PDvertex * trash_vrt( (*trash_itr).get() );

		for ( auto part_vrt : trash_vrt->nonconst_partners() ) {
			//   TR << "Examining trashed vertex " << trash_vrt->id() << " partner id " << (*part_itr)->id() << std::endl;

			// Skip any partners that are also marked for deletion
			if ( !(part_vrt->live()) ) continue;

			// Two infinite vertices don't form an edge, so no new finite vertex
			// in this case.  Just remove the trashed vertice from the partners
			// of the surviving infinite vertex.
			if ( !trash_vrt->finite() && !(part_vrt->finite()) ) {
				// For infinite vertex partners, remove the trashed vertex
				//    TR << "Infinite vertex id " << (*part_itr)->id() << " losing connection to trashed vertex " << trash_vrt->id() << std::endl;
				part_vrt->nonconst_partners().erase(
					std::remove(part_vrt->nonconst_partners().begin(),
					part_vrt->nonconst_partners().end(), trash_vrt ),
					part_vrt->nonconst_partners().end());
				continue;
			}

			// If we're here we need to make a new finite vertex
			PDvertexOP new_vrt( new PDvertex() );
			new_vrt->set_finite( true );
			new_vrt->set_live( true );
			new_vrt->nonconst_id() = vertex_count_++;
			//   TR << "Making a new finite vertex with id " << new_vrt->id() << std::endl;

			// Find the three generator atoms these two vertices have in common,
			// add in the new atom.
			// Sorting here isn't necessary.  Generators are sorted at vertex construction,
			// and they never change.
			std::set_intersection( part_vrt->nonconst_generators().begin(), part_vrt->nonconst_generators().end(),
				trash_vrt->nonconst_generators().begin(), trash_vrt->nonconst_generators().end(),
				std::back_inserter( new_vrt->nonconst_generators() ) );

			debug_assert( new_vrt->nonconst_generators().size() == 3 );
			//   TR << "Set intersection found " << new_vrt->nonconst_generators().size() << " in common " << std::endl;

			// Add in new atom
			new_vrt->nonconst_generators().push_back( new_sph.get() );
			std::sort( new_vrt->nonconst_generators().begin(), new_vrt->nonconst_generators().end() );
			link_vertex_to_generators( new_vrt );
			// Now we have enough info to get the xyz
			new_vrt->nonconst_xyz() = vertex_xyz_from_generators( new_vrt->nonconst_generators() );
			new_vrt->nonconst_power() = power_distance( new_vrt->xyz(), new_vrt->generators()[1] );

			//   TR << "Connecting new vertex id " << new_vrt->id() << " with old vertex id " << (*part_itr)->id() << std::endl;

			// Stitch up the connections on this edge
			new_vrt->nonconst_partners().push_back( part_vrt );
			// The trashed vertex is replaced by the new one in the surviving vertex
			debug_assert( std::find( part_vrt->nonconst_partners().begin(), part_vrt->nonconst_partners().end(),
				trash_vrt ) != part_vrt->nonconst_partners().end() );
			debug_assert( std::find( part_vrt->nonconst_partners().begin(), part_vrt->nonconst_partners().end(),
				new_vrt.get() ) == part_vrt->nonconst_partners().end() );
			std::replace( part_vrt->nonconst_partners().begin(), part_vrt->nonconst_partners().end(), trash_vrt, new_vrt.get() );
			debug_assert( std::find( part_vrt->nonconst_partners().begin(), part_vrt->nonconst_partners().end(),
				new_vrt.get() ) != part_vrt->nonconst_partners().end() );
			debug_assert( std::find( part_vrt->nonconst_partners().begin(), part_vrt->nonconst_partners().end(),
				trash_vrt ) == part_vrt->nonconst_partners().end() );

			// Move the new vertex to a temporary list - it still needs the
			// rest of its partners to be specified.
			new_vertices.push_front( new_vrt );
			new_vrt->my_itr() = new_vertices.begin();
		}
	}

	// Done with the trashed vertices.

	for ( auto tr_vrt : trash ) {
		tr_vrt->nonconst_partners().clear();
		tr_vrt->nonconst_generators().clear();
	}

	trash.clear();

	//////////////////////////////////////////////////////////////////////////////
	// Step 3 from Klenin paper:  Connect new vertices that share three generators
	//////////////////////////////////////////////////////////////////////////////

	// TR << "*** STEP 3 ***" << std::endl;

	// Note that we sorted the vectors of generator for all the new vertices,
	// so we can use set_intersection from the get-go.

	for ( std::list< PDvertexOP >::iterator itr1( new_vertices.begin() ) ;
			itr1 != new_vertices.end() ; ++itr1 ) {
		for ( std::list< PDvertexOP >::iterator itr2( std::next( itr1 ) ) ;
				itr2 != new_vertices.end() ; ++itr2 ) {
			utility::vector1< PDsphere* > shared_atoms;
			std::set_intersection( (*itr1)->nonconst_generators().begin(), (*itr1)->nonconst_generators().end(),
				(*itr2)->nonconst_generators().begin(), (*itr2)->nonconst_generators().end(),
				std::back_inserter( shared_atoms ) );
			if ( shared_atoms.size() == 3 ) {
				//    TR << "Connecting new finite vertex id " << (*itr1)->id() << " with other new finite vertex " << (*itr2)->id() <<std::endl;
				(*itr1)->nonconst_partners().push_back( (*itr2).get() );
				(*itr2)->nonconst_partners().push_back( (*itr1).get() );
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// Step 4 from Klenin paper:  Construct new infinite vertices
	/////////////////////////////////////////////////////////////

	// TR << "*** STEP 4 ***" << std::endl;

	// Any new finite vertices that have less than four partners
	// indicate that new infinite vertices must be made from these three.

	std::list< PDvertexOP > new_infinite_vertices;

	for ( auto nvp : new_vertices ) {
		if ( nvp->partners().size() < 4 ) {
			//   TR << "Checking generators of new finite vertex " << (*itr)->id()  << std::endl;
			// Check for each set of three (out of four) generators whether
			// any of the current partners shares all three.  If not, those
			// three are the generators of a new infinite vertex.

			// Make the four possible (sorted) sets of three
			for ( Size omit( 1 ) ; omit <= 4 ; ++omit ) {
				//    TR << "Working omit = " << omit << std::endl;
				// If prior passes through this loop have found all
				// of the required infinite vertices to make, we will
				// know because the partner count will have hit four,
				// and we can bag out.
				if ( nvp->partners().size() == 4 ) continue;
				utility::vector1< PDsphere * > only_three;
				for ( Size igen( 1 ) ; igen <= 4 ; igen++ ) {
					if ( igen != omit ) only_three.push_back( nvp->nonconst_generators()[ igen ] );
				}
				debug_assert( only_three.size() == 3 );
				// Check this set of three against each of the current partners
				bool build_infinite_vertex( true );
				for ( auto pvp : nvp->nonconst_partners() ) {
					if ( !build_infinite_vertex ) continue; // short circuit
					utility::vector1< PDsphere* > shared_gens;
					std::set_intersection( pvp->nonconst_generators().begin(), pvp->nonconst_generators().end(),
						only_three.begin(), only_three.end(),
						std::back_inserter( shared_gens ) );
					//     TR << "Infinite vertex build check shared gen atoms is " << shared_gens.size() << std::endl;
					if ( shared_gens.size() == 3 ) build_infinite_vertex = false;
				}
				if ( build_infinite_vertex ) {
					// We have the three generators to build an infinite vertex,
					// and the omitted atom tells us the opposite direction for the
					// plane normal.

					PDvertexOP new_inf_vrt( new PDvertex() );
					new_inf_vrt->set_finite( false );
					new_inf_vrt->set_live( true );
					new_inf_vrt->nonconst_id() = vertex_count_++;
					// Note these should already be sorted
					new_inf_vrt->nonconst_generators().push_back( only_three[1] );
					new_inf_vrt->nonconst_generators().push_back( only_three[2] );
					new_inf_vrt->nonconst_generators().push_back( only_three[3] );
					link_vertex_to_generators( new_inf_vrt );
					Vector dir( (only_three[1]->xyz() - only_three[2]->xyz()).cross( only_three[3]->xyz() - only_three[2]->xyz() ).normalize() );
					if ( dir.dot( nvp->generators()[omit]->xyz() - only_three[1]->xyz() ) > 0.0 ) dir *= -1.0;
					new_inf_vrt->nonconst_direction() = dir;
					new_infinite_vertices.push_front( new_inf_vrt );
					new_inf_vrt->my_itr() = new_infinite_vertices.begin();

					//     TR << "Making new infinite vertex with id " << new_inf_vrt->id() << std::endl;

					//     TR << "Generators are res " << only_three[1]->res() << " atom " << only_three[1]->atom() <<
					//           " res " << only_three[2]->res() << " atom " << only_three[2]->atom() <<
					//           " res " << only_three[3]->res() << " atom " << only_three[3]->atom() << std::endl;

					//       TR << "dir check atom is res " << (*itr)->generators()[omit]->res() << " atom " << (*itr)->generators()[omit]->atom() << std::endl;

					// Link up the connection
					//     TR << "Connecting new infinite vertex id " << new_inf_vrt->id() << " with old finite vertex id " << (*itr)->id() << std::endl;
					nvp->nonconst_partners().push_back( new_inf_vrt.get() );
					new_inf_vrt->nonconst_partners().push_back( nvp.get() );
				}
			}
		}
		debug_assert( nvp->partners().size() == 4 );
	}

	// At this point the new finite vertices are complete, and can be added
	// to the full set of live finite vertices
	finite_vertices_.splice( finite_vertices_.begin(), new_vertices );

	////////////////////////////////////////////////////////////////
	// Step 5 from Klenin paper:  Set partners for infinite vertices
	////////////////////////////////////////////////////////////////

	// TR << "*** STEP 5 ***" << std::endl;

	// First walk through the pre-existing infinite vertices.  If they
	// have four live partner vertices, they are fine.  Otherwise, they
	// need to find new partners.  They only need to look among the new
	// infinite vertices for matches, however.

	// TR << "There are currently " << infinite_vertices_.size() << " old infinite vertices" << std::endl;

	for ( auto ivp : infinite_vertices_ ) {
		// Check for dead partners
		while ( ivp->partners().size() != 4 ) {
			debug_assert( ivp->partners().size() <= 4 );
			for ( auto nivp : new_infinite_vertices ) {
				utility::vector1< PDsphere* > shared_atoms;
				std::set_intersection( ivp->nonconst_generators().begin(), ivp->nonconst_generators().end(),
					nivp->nonconst_generators().begin(), nivp->nonconst_generators().end(),
					std::back_inserter( shared_atoms ) );
				if ( shared_atoms.size() == 2 ) {
					// Link them up
					ivp->nonconst_partners().push_back( nivp.get() );
					nivp->nonconst_partners().push_back( ivp.get() );
				}
			}
		}
	}

	// Now the infinite vertices need to look amongst themselves to
	// find the last bunch of partners.

	for ( std::list< PDvertexOP >::iterator itr( new_infinite_vertices.begin() ) ;
			itr != new_infinite_vertices.end() ; ++itr ) {
		while ( (*itr)->partners().size() != 4 ) {

			debug_assert( (*itr)->partners().size() <= 4 );
			for ( std::list< PDvertexOP >::iterator new_itr( std::next( itr ) ) ;
					new_itr != new_infinite_vertices.end() ; ++new_itr ) {

				utility::vector1< PDsphere * > shared_atoms;
				std::set_intersection( (*itr)->nonconst_generators().begin(), (*itr)->nonconst_generators().end(),
					(*new_itr)->nonconst_generators().begin(), (*new_itr)->nonconst_generators().end(),
					std::back_inserter( shared_atoms ) );
				if ( shared_atoms.size() == 2 ) {
					// Link them up
					(*itr)->nonconst_partners().push_back( (*new_itr).get() );
					(*new_itr)->nonconst_partners().push_back( (*itr).get() );
				}
			}
		}
	}

	// Now merge the new infinite vertices in with the old
	infinite_vertices_.splice( infinite_vertices_.begin(), new_infinite_vertices );

	///////////////////////////////////////////////////////////////////
	// Step 6 (not from Klenin paper):  Perform a lot of error checking
	///////////////////////////////////////////////////////////////////

	// Check each finite vertex and make sure only the generators evaluate
	// to the stored power value for the vertex

#ifdef NOTDEF
			for( std::list< PDvertexOP >::iterator ec_itr( finite_vertices_.begin() ) ;
							ec_itr != finite_vertices_.end() ; ++ec_itr ) {
				Real const ec_power( (*ec_itr)->power() );
				for( utility::vector1< PDsphereOP >::iterator atm_itr =  atoms.begin() ;
								atm_itr != atoms.end() ; ++ atm_itr ) {
					if( std::find( (*ec_itr)->generators().begin(), (*ec_itr)->generators().end(), (*atm_itr) ) != (*ec_itr)->generators().end() ) continue;
					Real const gen_power( power_distance( (*ec_itr)->xyz(), (*atm_itr) ) );
					TR << "Comparing other power " << gen_power << " with vertex power " << ec_power << std::endl;
					if( gen_power < ec_power ) {
						TR.Error << "Calculating power between finite vertex id " << (*ec_itr)->id() << " and other atom less than generators!" << std::endl;
						TR.Error << "Calculated power " << gen_power << " less than generator power " << ec_power << std::endl;
					}
					debug_assert( gen_power > ec_power );
				}
			}

			// Check each infinite vertex and make sure all non-generators
			// are on the correct side of the generator plane

			for( std::list< PDvertexOP >::iterator ec_itr( infinite_vertices_.begin() ) ;
							ec_itr != infinite_vertices_.end() ; ++ec_itr ) {
				for( utility::vector1< PDsphereOP >::iterator atm_itr =  atoms.begin() ;
								atm_itr != atoms.end() ; ++ atm_itr ) {
					if( std::find( (*ec_itr)->generators().begin(), (*ec_itr)->generators().end(), (*atm_itr) ) != (*ec_itr)->generators().end() ) continue;
					Real const dir_dot( ( (*atm_itr)->xyz() - (*ec_itr)->generators()[1]->xyz() ).dot( (*ec_itr)->direction() ) );
					TR << "Comparing other dir dot product " << dir_dot << std::endl;
					if( dir_dot > 0.0 ) {
						TR.Error << "Calculating plane normal dot of infinite vertex id " << (*ec_itr)->id() << " with atom res " << (*atm_itr)->res() << " and atom " << (*atm_itr)->atom() << " on forbidden side!" << std::endl;
					}
					debug_assert( dir_dot < 0.0 );
				}
			}

	// Check all vertices to make sure they have the correct number of generators
	// and partners

	for( std::list< PDvertexOP >::iterator ec_itr( finite_vertices_.begin() ) ;
					ec_itr != finite_vertices_.end() ; ++ec_itr ) {
			if( (*ec_itr)->generators().size() != 4 ) {
				TR.Error << "finite vertex id " << (*ec_itr)->id() << " has wrong number of generators: " << (*ec_itr)->generators().size() << " - should be 4 " << std::endl;
			}
			debug_assert( (*ec_itr)->generators().size() == 4 );
	}

	for( std::list< PDvertexOP >::iterator ec_itr( infinite_vertices_.begin() ) ;
					ec_itr != infinite_vertices_.end() ; ++ec_itr ) {
			if( (*ec_itr)->generators().size() != 3 ) {
				TR.Error << "infinite vertex id " << (*ec_itr)->id() << " has wrong number of generators: " << (*ec_itr)->generators().size() << " - should be 3 " << std::endl;
			}
			debug_assert( (*ec_itr)->generators().size() == 3 );
	}


	for( std::list< PDvertexOP >::iterator ec_itr( finite_vertices_.begin() ) ;
					ec_itr != finite_vertices_.end() ; ++ec_itr ) {
			if( (*ec_itr)->partners().size() != 4 ) {
				TR.Error << "finite vertex id " << (*ec_itr)->id() << " has wrong number of partners: " << (*ec_itr)->partners().size() << " - should be 4 " << std::endl;
			}
			debug_assert( (*ec_itr)->partners().size() == 4 );
	}

	for( std::list< PDvertexOP >::iterator ec_itr( infinite_vertices_.begin() ) ;
					ec_itr != infinite_vertices_.end() ; ++ec_itr ) {
			if( (*ec_itr)->partners().size() != 4 ) {
				TR.Error << "infinite vertex id " << (*ec_itr)->id() << " has wrong number of partners: " << (*ec_itr)->partners().size() << " - should be 4 " << std::endl;
			}
			debug_assert( (*ec_itr)->partners().size() == 4 );
	}
#endif
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



std::list< PDinterOP >
PowerDiagram::get_intersections_for_atom( Size ires, Size iatm )
{

	PDsphere* psph( sphere_lookup_[ ires ][ iatm ] );
	return get_intersections_for_atom( psph );

}





std::list< PDinterOP >
PowerDiagram::get_intersections_for_atom( PDsphere* psph )
{

	std::list< PDinterOP > intersections;

	// Get circles of intersection.
	// Check through all sets of at least three finite
	// vertices that share a generator atom. Check whether
	// the convex polygon they form encloses the circle of intersection
	// where the plane slices the sphere.

	std::multimap< PDsphere*, PDvertex* > gen_dist2_map;
	utility::vector1< PDsphere* > valid_neighbors;
	// Load up a multimap with squared distances
	for ( auto vrt : psph->vertices() ) {
		if ( !vrt->finite() ) continue;
		// Check the dist between psph and each of the generators of vrt
		for ( auto osph : vrt->nonconst_generators() ) {
			if ( osph == psph ) continue;
			gen_dist2_map.insert( std::make_pair( osph, vrt ) );
			if ( std::find( valid_neighbors.begin(), valid_neighbors.end(), osph ) == valid_neighbors.end() ) {
				valid_neighbors.push_back( osph );
				//TR << "Pushing back " << valid_neighbors.size() << "-th atom with number " << osph->atom() << std::endl;
			}
		}
	}

	// Now iterate through the multimap to find enclosing neighbors
	for ( Size iat = 1 ; iat <= valid_neighbors.size() ; ++iat ) {
		PDsphere * osph( valid_neighbors[iat] );
		Size const vrt_count( gen_dist2_map.count( osph ) );
		//TR << "Testing other atom " << osph->atom() << " which has " << vrt_count << " shared finite vertices " <<std::endl;
		if ( vrt_count >= 3 ) {
			// Get iterators for these vertices
			std::pair< std::multimap< PDsphere*, PDvertex* >::iterator, std::multimap< PDsphere*, PDvertex* >::iterator > key_range;
			key_range = gen_dist2_map.equal_range( osph );

			// Quick check here to see if psph and osph have a common neighbor nsph
			// such that nsph intersects both.  If so, there is no circle.  Candidates
			// for nsph can be found by searching the vertices psph and osph have in
			// common and testing the generators for overlaps.

			if ( ( psph->xyz().distance( osph->xyz() ) > ( psph->rad() + osph->rad() ) ) ) {
				continue;
			}


			bool found_a_shared_overlap( false );
			for ( std::multimap< PDsphere*, PDvertex* >::iterator it = key_range.first ;
					it != key_range.second; ++it ) {
				if ( found_a_shared_overlap ) continue;
				PDvertex* nvrt( (*it).second );
				// Loop through the generators
				for ( auto nsph : nvrt->nonconst_generators() ) {
					if ( found_a_shared_overlap ) continue;
					if ( nsph == psph || nsph == osph ) continue;

					if ( ( psph->xyz().distance( nsph->xyz() ) < ( psph->rad() + nsph->rad() ) ) &&
							( osph->xyz().distance( nsph->xyz() ) < ( osph->rad() + nsph->rad() ) ) ) {
						found_a_shared_overlap = true;
					}
				}
			}

			if ( found_a_shared_overlap ) continue;

			//   TR << "Passed shared overlap check" << std::endl;

			// Need to determine if the polygon completely encloses the intersection
			// of the sphere with the plane of the polygon.

			Vector const v1( (*key_range.first).second->xyz() );
			Vector const v2( (*(std::next(key_range.first))).second->xyz() );
			Vector const v3( (*(std::next(std::next(key_range.first)))).second->xyz() );

			// Get the plane normal
			Vector const normal( (v2 - v1).cross( v3 - v1 ).normalized() );

			// Get the distance of the center of the sphere from the plane
			Vector const sphere_vec( psph->xyz() - v1 );
			//   Real const dist2_to_plane( sphere_vec.dot(normal)*sphere_vec.dot(normal) );
			Vector const center_in_plane( psph->xyz() - sphere_vec.dot(normal)*normal );

			//TR << "Comparing sphere distance^2 from plane of " << dist2_to_plane << " with rad^2 of " << psph->rad2() << std::endl;

			// If we have made it all the way here, we have a legitimate full
			// circle intersection.

			//   TR << "Found full circle intersection!" << std::endl;
			Vector new_inter_pt( center_in_plane );
			PDvertexCOP vrt1( (*key_range.first).second );
			PDinterOP new_inter( new PDinter( new_inter_pt, vrt1.get(), vrt1.get() ) );
			new_inter->add_atom( osph );
			new_inter->nonconst_circle() = true;
			intersections.push_back( new_inter );

		}
	}


	// Get points of intersections that correspond to triple atom overlaps
	for ( auto vrt1 : psph->vertices() ) {
		for ( auto vrt2 : vrt1->partners() ) {
			if ( vrt1->id() > vrt2->id() ) continue; // otherwise we do it twice
			// Make sure this other vertex is also in work_atom's list of vertices
			if ( find( psph->vertices().begin(), psph->vertices().end(), vrt2 ) !=
					psph->vertices().end() ) {

				// Handle by case - in any case you need a start point, a direction, and
				// a maximum extent if both vertices are finite.
				//TR << "Working with vertices " << vrt1->id() << " and " << vrt2->id() << std::endl;
				if ( vrt1->finite() && vrt2->finite() ) {
					Vector const start_pt( vrt1->xyz() );
					Vector const dir( ( vrt2->xyz() - vrt1->xyz() ).normalized() );
					Real const max_extent( vrt1->xyz().distance( vrt2->xyz() ) );
					//     Size const saved_size( intersections.size() );
					find_intersections( psph, vrt1, vrt2, start_pt, dir, max_extent, intersections );
					//TR << "Finite-finite pair found " << intersections.size() - saved_size << " intersections" << std::endl;

				} else if ( vrt1->finite() && !vrt2->finite() ) {
					Vector const start_pt( vrt1->xyz() );
					Vector const dir( vrt2->direction() );
					Real const max_extent( -1.0 );
					//     Size const saved_size( intersections.size() );
					find_intersections( psph, vrt1, vrt2, start_pt, dir, max_extent, intersections );
					//TR << "Finite-infinite pair found " << intersections.size() - saved_size << " intersections" << std::endl;

				} else if ( !vrt1->finite() && vrt2->finite() ) {
					Vector const start_pt( vrt2->xyz() );
					Vector const dir( vrt1->direction() );
					Real const max_extent( -1.0 );
					//     Size const saved_size( intersections.size() );
					find_intersections( psph, vrt1, vrt2, start_pt, dir, max_extent, intersections );
					//TR << "Infinite-finite pair found " << intersections.size() - saved_size << " intersections" << std::endl;

				} else {
					//TR << "Infinite-infinite pair - is this a thing?" << std::endl;
				} // Do nothing if both are infinite

			}
		}
	}

	return intersections;
}

void
find_intersections(
	PDsphere const * psph,
	PDvertex const * vrt1,
	PDvertex const * vrt2,
	Vector const & start_pt,
	Vector const & dir,
	Real const max_extent,
	std::list< PDinterOP > & intersections
)
{

	Vector const s( start_pt - psph->xyz() );
	Real const s_dot_dir( s.dot( dir ) );
	Real const s_mag_sqr( s.magnitude_squared() );
	Real const rad_sqr( psph->rad2() );
	Real const discrim( s_dot_dir*s_dot_dir - s_mag_sqr + rad_sqr );

	if ( discrim < 0.0 ) {
		//TR << "Discriminant is less than zero (value " << discrim << " ) - skipping" << std::endl;
		return; // No intersections
	}

	Real const discrim_sqrt( std::sqrt( discrim ) );

	Real const alpha1( -s_dot_dir + discrim_sqrt );
	Real const alpha2( -s_dot_dir - discrim_sqrt );

	//TR << "For comparison, alpha1 id " << alpha1 << " , alpha2 is " << alpha2 << " and max_extent is " << max_extent << std::endl;

	if ( alpha1 >= 0.0 && ( max_extent < 0.0 || alpha1 <= max_extent ) ) {
		Vector new_inter_pt( start_pt + alpha1*dir );
		//  PDsphereCOP other_atom( find_other_atom( psph, vrt1, vrt2 ) );
		PDinterOP new_inter( new PDinter( new_inter_pt, vrt1, vrt2 ) );
		intersections.push_back( new_inter );
		//TR << "Checky check atom radius " << psph->rad() << " compares with alpha1 point " << new_inter_pt.distance( psph->xyz() ) << std::endl;
	}

	if ( alpha2 >= 0.0 && ( max_extent < 0.0 || alpha2 <= max_extent ) ) {
		Vector new_inter_pt( start_pt + alpha2*dir );
		//  PDsphereCOP other_atom( find_other_atom( psph, vrt1, vrt2 ) );
		PDinterOP new_inter( new PDinter( new_inter_pt, vrt1, vrt2 ) );
		intersections.push_back( new_inter );
		//TR << "Checky check atom radius " << psph->rad() << " compares with alpha2 point " << new_inter_pt.distance( psph->xyz() ) << std::endl;
	}


	return;

}

bool
find_intersections(
	PDsphereCOP psph,
	Vector const & start_pt,
	Vector const & dir,
	Real const max_extent
)
{

	Vector const s( start_pt - psph->xyz() );
	Real const s_dot_dir( s.dot( dir ) );
	Real const s_mag_sqr( s.magnitude_squared() );
	Real const rad_sqr( psph->rad2() );
	Real const discrim( s_dot_dir*s_dot_dir - s_mag_sqr + rad_sqr );

	if ( discrim < 0.0 ) {
		//TR << "Discriminant is less than zero (value " << discrim << " ) - skipping" << std::endl;
		return false; // No intersections
	}

	Real const discrim_sqrt( std::sqrt( discrim ) );
	Real const alpha1( -s_dot_dir + discrim_sqrt );
	Real const alpha2( -s_dot_dir - discrim_sqrt );

	if ( alpha1 >= 0.0 && ( max_extent < 0.0 || alpha1 <= max_extent ) ) {
		Vector new_inter_pt( start_pt + alpha1*dir );
		//TR << "Checky check atom radius " << psph->rad() << " compares with alpha1 point " << new_inter_pt.distance( psph->xyz() ) << std::endl;
		return true;
	}

	if ( alpha2 >= 0.0 && ( max_extent < 0.0 || alpha2 <= max_extent ) ) {
		Vector new_inter_pt( start_pt + alpha2*dir );
		//TR << "Checky check atom radius " << psph->rad() << " compares with alpha2 point " << new_inter_pt.distance( psph->xyz() ) << std::endl;
		return true;
	}

	return false;

}

void
find_common_intersection_atoms(
	PDinterOP inter
)
{
	// Full circles already know their other atom
	if ( inter->circle() ) return;

	// All of the information is in the vertices' generator info
	PDvertex const * vrt1( inter->vrt1() );
	PDvertex const * vrt2( inter->vrt2() );

	if ( vrt1->generators().size() == 3 ) {
		//  TR << "Common atoms: ";
		//  for( utility::vector1< PDsphereOP const >::iterator itr( vrt1->generators().begin() ) ;
		for ( auto this_atom : vrt1->generators() ) {
			//   TR << " res " << atom->res() << " atom " << atom->atom() << "   ";
			inter->add_atom( this_atom );
		}
		//  TR << std::endl;
	} else if ( vrt2->generators().size() == 3 ) {
		//  TR << "Common atoms: ";
		for ( auto this_atom : vrt2->generators() ) {
			//   TR << " res " << atom->res() << " atom " << atom->atom() << "   ";
			inter->add_atom( this_atom );
		}
		//  TR << std::endl;
	} else {
		//  TR << "Common atoms: ";
		for ( auto atom1 : vrt1->generators() ) {
			if (  find( vrt2->generators().begin(), vrt2->generators().end(), atom1 ) != vrt2->generators().end() ) {
				//    TR << " res " << atom1->res() << " atom " << atom1->atom() << "   ";
				inter->add_atom( atom1 );
			}
		}
		//  TR << std::endl;
	}

}


utility::vector1< utility::vector1< SAnode > >
get_cycles_from_intersections(
	std::list< PDinterOP > & inters,
	PDsphere const * this_sph )
{
	Real const two_pi( 2.0*numeric::constants::d::pi );
	// Real const four_pi( 4.0*numeric::constants::d::pi );
	// Try out with just one cycle for now
	utility::vector1< utility::vector1< SAnode > > cycles;

	Size incorporated_inters( 0 );
	Size target_inters( inters.size() );

	// Process the intersections until none remain
	while ( incorporated_inters < target_inters ) {

		utility::vector1< SAnode > cycle;
		Size current_node_index(0);

		// Pick any intersection that isn't in a previous cycle and start from there

		PDinterOP i_start;

		for ( auto interp : inters ) {
			bool skip( false );
			// Check all cycles
			//   TR << "Checking potential new start intersection " << (*itr) << std::endl;
			for ( Size icyc = 1 ; icyc <= cycles.size() ; ++icyc ) {
				for ( Size ivrt = 1 ; ivrt <= cycles[icyc].size() ; ++ivrt ) {
					//     TR << "Against stored intersection " << ((cycles[icyc])[ivrt].inter()) << std::endl;
					if ( interp == ((cycles[icyc])[ivrt].inter()) ) skip = true;
				}
			}

			// Must be ok
			if ( !skip ) {
				i_start = interp;
				break;
			}
		}

		PDinterOP i_current( i_start );

		PDsphere const * other_sph( ( i_current->atoms()[1] == this_sph ?
			i_current->atoms()[2] :
			i_current->atoms()[1] ) );

		bool cycle_done( false );


		// Full circles are complete cycles on their own
		if ( i_current->circle() ) {
			// Circles only know the other atom
			PDsphere const * other_sph( i_current->atoms()[1] );
			SAnode this_node( i_current, other_sph );
			cycle.push_back( this_node );
			//   TR << "Pushing completed cycle with " << cycle.size() << " intersections" << std::endl;
			cycles.push_back( cycle ); // Maybe 'cycle' should hold pointers rather than SAnodes.  This copy could be costly.
			incorporated_inters += cycle.size();
			continue;
		}

		while ( !cycle_done ) {
			// Figure out the direction to search
			PDsphere const * limit_atom( ( i_current->atoms()[1] != this_sph && i_current->atoms()[1] != other_sph ) ?  i_current->atoms()[1] :
				( ( i_current->atoms()[2] != this_sph && i_current->atoms()[2] != other_sph ) ? i_current->atoms()[2] : i_current->atoms()[3] ) );

			//   TR << "Got the atoms to work with" << std::endl;
			//   TR << "This atom " << this_sph->atom() << std::endl;
			//   TR << "Other atom " << other_sph->atom() << std::endl;
			//   TR << "Limit atom " << limit_atom->atom() << std::endl;

			// Find the next intersection along this direction on this circle
			Vector axis( (other_sph->xyz() - this_sph->xyz()).normalized() );
			Vector temp( i_current->xyz() - this_sph->xyz() );
			Vector local_inter( temp - ((temp.dot( axis ))*axis) );
			Vector xaxis( local_inter.normalized() );
			Vector yaxis( axis.cross( xaxis ) );
			Vector limit_vec( i_current->xyz() - limit_atom->xyz() );
			//   if( yaxis.dot( limit_vec ) < 0.0 ) {
			//    TR << "Positive rotation brings me closer to the limit atom!" << std::endl;
			//   } else {
			//    TR << "Positive rotation takes me farther from the limit atom!" << std::endl;
			//   }

			// Ouch this is painful.  Another good place to stop and find a more
			// efficient way.  Basically, we want search positive to be false.
			// If it is currently true, switch the 'other' and 'limit' atoms.

			if ( yaxis.dot( limit_vec ) < 0.0 ) {
				//    TR << "Switching direction of path!" << std::endl;
				PDsphere const * tmp_switch( other_sph );
				other_sph = limit_atom;
				limit_atom = tmp_switch;
				axis = ( (other_sph->xyz() - this_sph->xyz()).normalized() );
				temp = ( i_current->xyz() - this_sph->xyz() );
				local_inter = ( temp - ((temp.dot( axis ))*axis) );
				xaxis = ( local_inter.normalized() );
				yaxis = ( axis.cross( xaxis ) );
				limit_vec = ( i_current->xyz() - limit_atom->xyz() );
				//    if( yaxis.dot( limit_vec ) < 0.0 ) {
				//     TR << "Positive rotation brings me closer to the limit atom!" << std::endl;
				//    } else {
				//     TR << "Positive rotation takes me farther from the limit atom!" << std::endl;
				//    }
			}

			SAnode this_node( i_current, other_sph );
			cycle.push_back( this_node );
			current_node_index++;

			//   TR << "Got the atoms to work with" << std::endl;
			//   TR << "This atom " << this_sph->atom() << std::endl;
			//   TR << "Other atom " << other_sph->atom() << std::endl;
			//   TR << "Limit atom " << limit_atom->atom() << std::endl;

			bool search_positive( yaxis.dot( limit_vec ) < 0.0 );
			Real best_val( -1.0 );
			PDinterOP next_inter;

			//   TR << "Search positive is " << search_positive << std::endl;

			// Now that we have it, check to see if it is the same as the first in this cycle.
			// If so, we have completed to cycle and we are done.

			for ( auto interp : inters ) {
				if ( interp == i_current ) {
					//     TR << "Skipping current intersection" << std::endl;
					continue;
				}

				// Full circles aren't part of other cycles.
				if ( interp->circle() ) continue;

				if ( share_axis_atoms( interp, this_sph, other_sph ) ) {
					// Calculate angle for this intersection point
					//     TR << "Checking angle for intersections" << std::endl;
					//     TR << "Current Inter atoms " << i_current->atoms()[1]->atom() << "   "
					//                << i_current->atoms()[2]->atom() << "   "
					//                << i_current->atoms()[3]->atom() << std::endl;
					//     TR << "Check Inter atoms " << (*itr)->atoms()[1]->atom() << "   "
					//                << (*itr)->atoms()[2]->atom() << "   "
					//                << (*itr)->atoms()[3]->atom() << std::endl;

					Vector temp2( interp->xyz() - this_sph->xyz() );
					Vector check( (temp2 - (temp2.dot(axis)*axis) ).normalized() );
					Real xcomp( check.dot( xaxis ) );
					Real ycomp( check.dot( yaxis ) );
					//     TR << "Xcomp, ycomp are:  " << xcomp << " and " << ycomp << std::endl;
					Real angle( std::atan2( ycomp, xcomp ) );
					//     TR << "Unadjusted angle " << angle << std::endl;
					if ( angle < 0.0 ) angle += two_pi;

					//     TR << "Found angle " << angle << std::endl;

					// See if this is better than the current best stored point
					if ( best_val < 0.0 ) {
						//      TR << "Found new best intersection" << std::endl;
						best_val = angle;
						next_inter = interp;
					} else if ( search_positive && ( angle > best_val ) ) {
						//      TR << "Found new best intersection" << std::endl;
						best_val = angle;
						next_inter = interp;
					} else if ( !search_positive && ( angle < best_val ) ) {
						//      TR << "Found new best intersection" << std::endl;
						best_val = angle;
						next_inter = interp;
					}
				} else {
					//     TR << "Skipping vertex - not on this great circle" << std::endl;
				}
			}

			if ( best_val < 0.0 ) {
				//    TR << "Didn't find any other intersection on this circle!" << std::endl;
				cycle_done = true;
			} else {
				i_current = next_inter;

				cycle[current_node_index].set_phi( best_val );
				Real b_i( this_sph->xyz().distance( other_sph->xyz() ) );
				cycle[current_node_index].set_cos_theta( 0.5*(b_i + ( this_sph->rad2() - other_sph->rad2() )/b_i ) / this_sph->rad() );

				if ( i_current->atoms()[1] != this_sph && i_current->atoms()[1] != other_sph ) {
					other_sph = i_current->atoms()[1];
				} else if ( i_current->atoms()[2] != this_sph && i_current->atoms()[2] != other_sph ) {
					other_sph = i_current->atoms()[2];
				} else {
					other_sph = i_current->atoms()[3];
				}
			}

			if ( next_inter == i_start ) {
				//    TR << "Looped back around to first intersection in cycle!" << std::endl;
				cycle_done = true;
			}
		} // end this cycle
		//  TR << "Pushing completed cycle with " << cycle.size() << " intersections" << std::endl;
		cycles.push_back( cycle ); // Maybe 'cycle' should hold pointers rather than SAnodes.  This copy could be costly.
		incorporated_inters += cycle.size();
	} // end while intersection list not empty

	// TR << "Found " << cycles.size() << " cycles" << std::endl;

	return cycles;

}

Real
get_sasa_from_cycles(
	utility::vector1< utility::vector1< SAnode > > & cycles,
	PDsphere* this_sph )
{

	// Real const two_pi( 2.0*numeric::constants::d::pi );
	Real const four_pi( 4.0*numeric::constants::d::pi );
	// Real contour_surface_area( four_pi * this_sph->rad2() ); // this is the exposed part
	Real const full_surf_area = four_pi * this_sph->rad2();
	Real contour_surface_area( 0.0 ); // this is the exposed part
	Real circle_surface_area( 0.0 ); // This part is all buried

	for ( Size i = 1 ; i <= cycles.size() ; ++i ) {
		//  for( utility::vector1< SAnode >::iterator itr = cycles[i].begin() ; itr != cycles[i].end() ; ++itr ) {
		//   TR << "Intersection with atoms " <<
		//    " atom " << (*itr).inter()->atoms()[1]->atom() <<
		//    " atom " << (*itr).inter()->atoms()[2]->atom() <<
		//    " atom " << (*itr).inter()->atoms()[3]->atom() << std::endl;
		//  }

		//  TR << "Found " << cycles[i].size() << " arcs in cycle" << std::endl;
		Real surf_area( get_area_from_cycle( this_sph, cycles[i] ) );
		if ( (*(cycles[i]).begin()).inter()->circle() ) {
			//   TR << "Got circle surf " << surf_area << std::endl;
			circle_surface_area += (full_surf_area - surf_area);
		} else {
			//   TR << "Got contour surf " << surf_area << std::endl;
			contour_surface_area += surf_area;
		}
	}
	// while( contour_surface_area > full_surf_area ) {
	//  contour_surface_area -= full_surf_area;
	// }
	Real surf_area( contour_surface_area - circle_surface_area );
	while ( surf_area < 0.0 ) {
		surf_area += full_surf_area;
	}
	// TR << "Surface area is " << surf_area << " out of full area " << full_surf_area << " with complement " << full_surf_area - surf_area << std::endl;

	return surf_area;

}



void
get_derivs_from_cycles(
	utility::vector1< utility::vector1< SAnode > > & cycles,
	PDsphereOP & this_atom,
	PDsphereOP & check_atom,
	Vector & f1,
	Vector & f2 )
{

	for ( Size i = 1 ; i <= cycles.size() ; ++i ) {
		get_derivs_from_cycle( cycles[i], this_atom, check_atom, f1, f2 );
	}

	return;
}

void
get_derivs_from_cycle(
	utility::vector1< SAnode > & cycle,
	PDsphereOP & this_atom,
	PDsphereOP & check_atom,
	Vector & f1,
	Vector & f2 )
{

	// Real const two_pi( 2.0*numeric::constants::d::pi );
	// Real const four_pi( 4.0*numeric::constants::d::pi );
	// Real const this_rad( this_atom->rad() );

	// Each cycle can be a complete circle or a series of arcs.

	// Handle full circles here
	if ( (*cycle.begin()).inter()->circle() ) {
		PDsphereCOP other_atom( (*cycle.begin()).inter()->atoms()[1] );
		if ( other_atom != check_atom ) { return; }

		//  TR << "Handling full circle derivative case" << std::endl;

		Real const this_rad( this_atom->rad() );
		Real const other_rad( other_atom->rad() );
		Real const atom_dist( this_atom->xyz().distance( other_atom->xyz() ) );

		Vector const p1( this_atom->xyz() );
		Vector const p2( other_atom->xyz() );
		Vector const r12( p2 - p1 );
		Vector const p2_cross_r12( -1.0*p2.cross_product( p1 ) );

		if ( atom_dist != atom_dist ) {
			TR << "Bad zero atom-atom distance" << std::endl;
			TR << "This atom xyz is " << p1 << std::endl;
			TR << "Other atom xyz is " << p2 << std::endl;
		}

		//  Real const para_fac_full( numeric::constants::d::pi*(this_rad+other_rad)*( 1.0 - ( this_rad - other_rad )*( this_rad - other_rad )/(atom_dist*atom_dist) )/atom_dist );
		Real const para_fac_full( numeric::constants::d::pi*this_rad*( 1.0 - ( this_rad - other_rad )*( this_rad - other_rad )/(atom_dist*atom_dist) )/atom_dist );

		f1 += para_fac_full*p2_cross_r12;
		f2 += -1.0*para_fac_full*r12;

		return;
	}

	// Loop through arcs, getting derivatives for each
	Size arc_num( 0 );
	for ( utility::vector1< SAnode >::iterator itr = cycle.begin() ; itr != cycle.end() ; ++itr ) {
		// What kind of information is useful here?
		// Certainly the 'other atom' that forms the arc with the atom of interest.
		// Maybe should store this.

		arc_num++;

		//Size arc_num_m1( arc_num - 1 );
		//if ( arc_num_m1 == 0 ) arc_num_m1 = cycle.size();
		//Size arc_num_p1( arc_num + 1 );
		//if ( arc_num_p1 > cycle.size() ) arc_num_p1 = 1;
		//  TR << "Working indices are " << arc_num << " and " << arc_num_m1 << std::endl;

		PDsphereCOP other_atom( cycle[ arc_num ].other_atom() );

		if ( other_atom != check_atom ) { continue; }

		//  TR << "Handling arc derivative case" << std::endl;

		Real const this_rad( this_atom->rad() );
		//  Real const other_rad( cycle[arc_num].inter()->atoms()[1]->rad() );
		Real const other_rad( other_atom->rad() );
		Real const atom_dist( this_atom->xyz().distance( other_atom->xyz() ) );
		Real const phi_i( cycle[arc_num].phi() );

		Vector const p1( this_atom->xyz() );
		Vector const p2( other_atom->xyz() );
		Vector const r12( p2 - p1 );

		Real const perp_dist_denom( atom_dist * atom_dist );

		Real const para_fac( 0.5*phi_i*this_rad*( 1.0 - ( this_rad*this_rad - other_rad*other_rad )/(atom_dist*atom_dist) )/atom_dist );
		Real const perp_fac( -1.0*this_rad/perp_dist_denom );

		PDinterCOP inter_i( (*itr).inter() );
		PDinterCOP inter_ip1( std::next( itr ) == cycle.end() ? (*cycle.begin()).inter() : (*std::next( itr )).inter() );

		Vector const arc_chord( inter_ip1->xyz() - inter_i->xyz() );
		Vector const chord_cross_r12( arc_chord.cross_product( r12 ) );
		Vector const p2_cross_r12( -1.0*p2.cross_product( p1 ) );


		//  Real const para_fac2( 0.5*phi_i*other_rad*( 1.0 - ( other_rad*other_rad - this_rad*this_rad )/(atom_dist*atom_dist) )/atom_dist );
		//  Real const perp_fac2( -1.0*other_rad/perp_dist_denom );
		Vector const r21( p1 - p2 );
		Vector const arc_chord2( -arc_chord );
		Vector const chord2_cross_r21( arc_chord2.cross_product( r21 ) );
		Vector const p1_cross_r21( -1.0*p1.cross_product( p2 ) );

		//  Real const para_fac_full( 0.5*phi_i*(this_rad+other_rad)*( 1.0 - ( this_rad - other_rad )*( this_rad - other_rad )/(atom_dist*atom_dist) )/atom_dist );
		//  Real const para_fac_full( 0.5*phi_i*this_rad*( 1.0 - ( this_rad - other_rad )*( this_rad - other_rad )/(atom_dist*atom_dist) )/atom_dist );


		// Get the effect on the surface area of second atom.

		f1 += para_fac*p2_cross_r12 + perp_fac*p2.cross_product(chord_cross_r12);
		f2 += -1.0*para_fac*r12 - perp_fac*chord_cross_r12;

		//  f1 += para_fac_full*p2_cross_r12 + (perp_fac-perp_fac2)*p2.cross_product(chord_cross_r12);
		//  f2 += -1.0*para_fac_full*r12 - (perp_fac-perp_fac2)*chord_cross_r12;


		continue;

	}




	return;
}

bool
share_axis_atoms(
	PDinterCOP v1,
	PDsphere const * a1,
	PDsphere const * a2
)
{
	Size same_count( 0 );
	// Brute force ok here
	if ( v1->atoms()[1] == a1 ) same_count++;
	if ( v1->atoms()[1] == a2 ) same_count++;
	if ( v1->atoms()[2] == a1 ) same_count++;
	if ( v1->atoms()[2] == a2 ) same_count++;
	if ( v1->atoms()[3] == a1 ) same_count++;
	if ( v1->atoms()[3] == a2 ) same_count++;

	// TR << "Intersection with atoms " <<
	//  " res " << v1->atoms()[1]->res() << " atom " << v1->atoms()[1]->atom() <<
	//  " res " << v1->atoms()[2]->res() << " atom " << v1->atoms()[2]->atom() <<
	//  " res " << v1->atoms()[3]->res() << " atom " << v1->atoms()[3]->atom() << std::endl;
	// TR << "Axis atom 1 res " << a1->res() << " atom " << a1->atom() << std::endl;
	// TR << "Axis atom 2 res " << a2->res() << " atom " << a2->atom() << std::endl;


	if ( same_count == 2 ) {
		//  TR << "Returning true" << std::endl;
		return true;
	} else {
		//  TR << "Returning false" << std::endl;
		return false;
	}

}


Real
get_area_from_cycle(
	PDsphere* this_atom,
	utility::vector1< SAnode > & cycle
)
{
	Real const two_pi( 2.0*numeric::constants::d::pi );
	Real const four_pi( 4.0*numeric::constants::d::pi );

	// Handle full circles here
	if ( (*cycle.begin()).inter()->circle() ) {
		Real const this_rad( this_atom->rad() );
		Real const other_rad( (*cycle.begin()).inter()->atoms()[1]->rad() );
		Real const atom_dist( this_atom->xyz().distance( (*cycle.begin()).inter()->atoms()[1]->xyz() ) );
		Real const cos_theta( (0.5/this_rad)*( atom_dist + ( this_rad*this_rad - other_rad*other_rad )/atom_dist ) );

		Real surf_area = four_pi*this_atom->rad2() - two_pi*( 1.0 - cos_theta )*this_atom->rad2();

		return surf_area;
	}


	// Try getting the surface area for a single atom.
	Real surf_area( two_pi );
	Size arc_num( 0 );
	for ( utility::vector1< SAnode >::iterator itr = cycle.begin() ; itr != cycle.end() ; ++itr ) {
		// What kind of information is useful here?
		// Certainly the 'other atom' that forms the arc with the atom of interest.
		// Maybe should store this.

		arc_num++;
		Size arc_num_m1( arc_num - 1 );
		if ( arc_num_m1 == 0 ) arc_num_m1 = cycle.size();

		//  TR << "Working indices are " << arc_num << " and " << arc_num_m1 << std::endl;

		PDsphere const * other_atom( cycle[ arc_num ].other_atom() );
		PDsphere const * other_atom_m1( cycle[ arc_num_m1 ].other_atom() );

		//  TR << "Other atom is " << other_atom->atom() << " and other minus one is " << other_atom_m1->atom() << std::endl;

		// Package all this into an inlined function later
		PDinterCOP inter_i( (*itr).inter() );
		PDinterCOP inter_ip1( std::next( itr ) == cycle.end() ? (*cycle.begin()).inter() : (*std::next( itr )).inter() );

		Vector e_i( ( other_atom->xyz() - this_atom->xyz() ).normalized() );
		Vector e_im1( ( other_atom_m1->xyz() - this_atom->xyz() ).normalized() );

		Vector v_i( ( inter_i->xyz() - this_atom->xyz() ).normalized() );
		Vector v_ip1( ( inter_ip1->xyz() - this_atom->xyz() ).normalized() );

		Vector outer_i( e_i.cross( v_i ) );
		Vector outer_im1( e_im1.cross( v_i ) );

		Real const cos_theta_i( cycle[arc_num].cos_theta() );
		Real const cos_theta_im1( cycle[arc_num_m1].cos_theta() );
		Real const sin_theta_i( std::sqrt( 1.0 - cos_theta_i*cos_theta_i ) );
		Real const sin_theta_im1( std::sqrt( 1.0 - cos_theta_im1*cos_theta_im1 ) );

		Real const phi_i( cycle[arc_num].phi() );

		//  Real const cos_phi_i( ( v_i.dot( v_ip1 ) - cos_theta_i*cos_theta_i ) / ( sin_theta_i*sin_theta_i ) );


		//  TR << "Cached phi is " << phi_i << " and paper version is " << calc_phi << std::endl;


		Real const pi_minus_psi( std::acos( ( e_i.dot( e_im1 ) - cos_theta_i*cos_theta_im1 )/(sin_theta_i*sin_theta_im1) ) );
		//  TR << "pi_minus_psi parts are " << e_i.dot( e_im1 ) << " cti " << cos_theta_i << " ctim1 " << cos_theta_im1 << " sti " << sin_theta_i << " stim1 " << sin_theta_im1 << " giving an argument to acos of: " <<  ( e_i.dot( e_im1 ) - cos_theta_i*cos_theta_im1 )/(sin_theta_i*sin_theta_im1) << std::endl;

		//  TR << "This contribution to the surf area is: " << ( phi_i*cos_theta_i - pi_minus_psi )  << " based on parts " << phi_i*cos_theta_i << " and " << pi_minus_psi << std::endl;

		surf_area += ( phi_i*cos_theta_i - pi_minus_psi );


		while ( surf_area > four_pi ) {
			surf_area -= four_pi;
		}
	}

	surf_area *= this_atom->rad2();

	return surf_area;

}


Vector
vertex_xyz_from_generators( PDsphere const * a1, PDsphere const * a2, PDsphere const * a3, PDsphere const * a4 )
{
	// Need three planes to interesect
	Vector const n1( (a2->xyz() - a1->xyz()).normalize() );
	Vector const n2( (a3->xyz() - a1->xyz()).normalize() );
	Vector const n3( (a4->xyz() - a1->xyz()).normalize() );

	// Calculate position along atom-atom vectors where plane is located
	Real const d_sqr1( (a2->xyz().distance_squared( a1->xyz() ) ) );
	Vector const p1( a1->xyz() + ( (a2->xyz() - a1->xyz()) *
		0.5*( 1.0 - (a2->rad2() - a1->rad2() )/d_sqr1 ) ) );

	Real const d_sqr2( (a3->xyz().distance_squared( a1->xyz() ) ) );
	Vector const p2( a1->xyz() + ( (a3->xyz() - a1->xyz()) *
		0.5*( 1.0 - (a3->rad2() - a1->rad2() )/d_sqr2 ) ) );

	Real const d_sqr3( (a4->xyz().distance_squared( a1->xyz() ) ) );
	Vector const p3( a1->xyz() + ( (a4->xyz() - a1->xyz()) *
		0.5*( 1.0 - (a4->rad2() - a1->rad2() )/d_sqr3 ) ) );

	// Check that the planes do indeed intersect
	Matrix const det_check( Matrix::cols( n1, n2, n3 ) );
	Real const det_val( det_check.det() );

	//TR << "Value of determinant is " << det_val << std::endl;

	// See, e.g. Plane-plane intersection on Wolfram Mathworld
	Vector p_intersect( ( p1.dot( n1 )*n2.cross( n3 ) +
		p2.dot( n2 )*n3.cross( n1 ) +
		p3.dot( n3 )*n1.cross( n2 ) )/det_val );

	return p_intersect;
}

Vector
vertex_xyz_from_generators( utility::vector1< PDsphere* > const & gen )
{
	return vertex_xyz_from_generators( gen[1], gen[2], gen[3], gen[4] );
}



// DIAGNOSTIC OUTPUT FOR POWER DIAGRAM STUFF

void
print_points( std::list< PDinterOP > & inters )
{

	Size counter( 0 );
	for ( std::list< PDinterOP >::iterator itr = inters.begin() ; itr != inters.end() ; ++itr ) {
		TR << "ATOM  " << I(5,++counter) << "  FE  VRT X" <<
			I(4,counter ) << "    " <<
			F(8,3,(*itr)->xyz()(1)) <<
			F(8,3,(*itr)->xyz()(2)) <<
			F(8,3,(*itr)->xyz()(3)) <<
			F(6,2,1.0) << F(6,2,1.0) << std::endl;
	}
	TR << "END" << std::endl;

	return;
}

void
print_partners( PDvertexCOP vrt )
{

	TR << "Partners: ";
	for ( std::vector< PDvertex * >::const_iterator itr( vrt->partners().begin() ) ;
			itr != vrt->partners().end() ; ++itr ) {
		TR << " id " << (*itr)->id() << "   ";
	}
	TR << std::endl;

	return;
}

void
print_generators( PDvertexCOP vrt )
{

	TR << "Generators: ";
	for ( auto gen : vrt->generators() ) {
		TR << " res " << gen->res() << " atom " << gen->atom() << "   ";
	}
	TR << std::endl;

	return;
}




void
print_vertices( std::list< PDvertexOP > & fv, std::list< PDvertexOP > & iv )
{

	for ( std::list< PDvertexOP >::iterator itr = fv.begin() ; itr != fv.end() ; ++itr ) {
		if ( !(*itr)->finite() ) continue;
		TR << "ATOM  " << I(5,(*itr)->id()) << "  FE  VRT X" <<
			I(4,(*itr)->id()) << "    " <<
			F(8,3,(*itr)->xyz()(1)) <<
			F(8,3,(*itr)->xyz()(2)) <<
			F(8,3,(*itr)->xyz()(3)) <<
			F(6,2,1.0) << F(6,2,1.0) << std::endl;
	}
	for ( std::list< PDvertexOP >::iterator itr = iv.begin() ; itr != iv.end() ; ++itr ) {
		PDvertex * finite_partner = nullptr;

		for ( std::vector< PDvertex * >::iterator part_itr( (*itr)->nonconst_partners().begin() ) ;
				part_itr != (*itr)->nonconst_partners().end(); ++part_itr ) {
			if ( (*part_itr)->finite() ) {
				finite_partner = (*part_itr);
			}
		}

		Vector vrt( finite_partner->xyz() + 5.0*(*itr)->direction() );
		TR << "ATOM  " << I(5,(*itr)->id()) << "  FE  VRT X" <<
			I(4,(*itr)->id() ) << "    " <<
			F(8,3,vrt(1)) <<
			F(8,3,vrt(2)) <<
			F(8,3,vrt(3)) <<
			F(6,2,1.0) << F(6,2,1.0) << std::endl;
	}
	TR << "END" << std::endl;

	return;
}

void
print_vertices( std::list< PDvertexOP > & v )
{

	Size finite_count( 0 );
	for ( std::list< PDvertexOP >::iterator itr = v.begin() ; itr != v.end() ; ++itr ) {
		if ( (*itr)->finite() ) {
			finite_count++;
			TR << "ATOM  " << I(5,(*itr)->id()) << "  FE  FIN X" <<
				I(4,(*itr)->id()) << "    " <<
				F(8,3,(*itr)->xyz()(1)) <<
				F(8,3,(*itr)->xyz()(2)) <<
				F(8,3,(*itr)->xyz()(3)) <<
				F(6,2,1.0) << F(6,2,1.0) << std::endl;
		} else {
			PDvertex * finite_partner = nullptr;

			for ( std::vector< PDvertex * >::iterator part_itr( (*itr)->nonconst_partners().begin() ) ;
					part_itr != (*itr)->nonconst_partners().end(); ++part_itr ) {
				if ( (*part_itr)->finite() ) {
					finite_partner = (*part_itr);
				}
			}

			//  TR << "There are " << finite_count << " finite vertices" << std::endl;

			Vector vrt( finite_partner->xyz() + 5.0*(*itr)->direction() );
			TR << "ATOM  " << I(5,(*itr)->id()) << "  FE  INF X" <<
				I(4,(*itr)->id() ) << "    " <<
				F(8,3,vrt(1)) <<
				F(8,3,vrt(2)) <<
				F(8,3,vrt(3)) <<
				F(6,2,1.0) << F(6,2,1.0) << std::endl;

		}
		print_generators( *itr );
		print_partners( *itr );
	}
	TR << "END" << std::endl;

	return;
}

PDvertex *
find_next_vertex_with_smallest_dist( PDvertex * srch_vrt, PDsphereOP & new_sph, Real & this_dist )
{
	PDvertex * return_vrt = nullptr;
	for ( auto itr = srch_vrt->partners().begin() ; itr != srch_vrt->partners().end() ; ++itr ) {
		// Only checking finite vertices
		if ( !(*itr)->finite() ) continue;
		Real check_dist( power_distance( (*itr)->xyz(), new_sph.get() ) - (*itr)->power() );
		if ( check_dist < this_dist ) {
			this_dist = check_dist;
			return_vrt = (*itr);
		}
	}
	if ( return_vrt == nullptr ) return_vrt = srch_vrt;
	return return_vrt;
}


} // power_diagram
} // scoring
} // core

