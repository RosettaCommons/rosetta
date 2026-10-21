// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ligand_docking/PlaceWaters.cc
/// @brief Places waters in the interface
/// @author Laura Shub (laurashub@gmail.com)

// Unit headers
#include <protocols/ligand_docking/PlaceWaters.hh>
#include <protocols/ligand_docking/PlaceWatersCreator.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueDatabaseIO.hh>
#include <core/chemical/residue_io.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <core/kinematics/Jump.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/grid/CartGrid.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/conversions.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <fstream>


// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/datacache/DataMap.hh>


static basic::Tracer TR( "protocols.ligand_docking.PlaceWaters" );

namespace protocols {
namespace ligand_docking {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PlaceWaters::PlaceWaters():
	Mover( PlaceWaters::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PlaceWaters::PlaceWaters( PlaceWaters const & ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PlaceWaters::~PlaceWaters(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
PlaceWaters::apply( core::pose::Pose & adding_water ){ //change "adding_water" to "pose" to make it look better.

	//core::pose::Pose adding_water(pose);
	debug_assert(ligand_chain_.size() == 1);
	core::Size chain_id = core::pose::get_chain_id_from_chain(ligand_chain_, adding_water);
	utility::vector1<core::Size> resnums( core::pose::get_resnums_for_chain (adding_water, chain_id));

	if ( resnums.size() > 1 ) {
		utility_exit_with_message( "Specified ligand chain has more than 1 residue.");
	}

	core::Size jump_id = core::pose::get_jump_id_from_chain_id(chain_id, adding_water);
	core::Size const begin(adding_water.conformation().chain_begin(chain_id));
	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(adding_water, jump_id));

	core::conformation::UltraLightResidueOP ligand(new core::conformation::UltraLightResidue ( adding_water.residue(begin).get_self_ptr()));

	if ( grid_size_ == 0 ) {
		grid_size_ = grid_size(ligand);
	}

	TR << "Using grid size " << grid_size_ << std::endl;
	initialize_grids( center, grid_size_, resolution_);

	fill_occupied(adding_water);

	utility::vector1<core::Vector> waters = empty_spaces(cutoff_);
	if ( output_water_coords_ ) {
		output_water_coords_file(waters, "1.emptyspaces_coords.pdb");
	}

	cluster_groups(waters, group_cutoff_, merge_cutoff_);
	TR << "Number of clusters created: " << waters.size() << std::endl;
	if ( output_water_coords_ ) {
		output_water_coords_file(waters, "2.clusters_coords.pdb");
	}

	//try to detect hydrogen bonding sites
	detect_hydrogens(adding_water);
	utility::vector1<water_info> final_waters = check_possible_hbonds(waters, hbonds_inner_ , hbonds_outer_);
	TR << "Number of final waters after hbond checks: " << final_waters.size() << std::endl;

	core::chemical::ChemicalManager * chem_mang = core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCAP atom_types = chem_mang->atom_type_set("fa_standard");
	core::chemical::ElementSetCAP elements = chem_mang->element_set("default");
	core::chemical::MMAtomTypeSetCAP mm_atom_types = chem_mang->mm_atom_type_set("fa_standard");
	core::chemical::orbitals::OrbitalTypeSetCAP orbital_types = chem_mang->orbital_type_set("fa_standard");
	// new water ResidueCOP where water_type option is given
	core::conformation::ResidueCOP new_water_resi( core::conformation::ResidueFactory::create_residue(
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( water_type_ ) ) );

	TR << "Adding " << final_waters.size() << " waters" << std::endl;
	// If outputting waters separately to PDB file outside of pose (with "output_water_coords" option)
	if ( output_water_coords_ ) {
		utility::vector1<core::Vector> ws ;
		if ( final_waters.size() > 0 ) {
			for ( core::Size i = 0; i <= final_waters.size() ; ++i ) {
				ws.push_back(final_waters[i].centroid);
			}
			output_water_coords_file(ws, "3.final_coords.pdb");
		}
	}

	core::Size water_num = 0;
	for ( water_info wi : final_waters ) {
		++water_num;

		//core::conformation::Residue new_water_resi (water_restype, true);
		adding_water.append_residue_by_jump(*new_water_resi, begin , "", "", true);

		//change to UltraLightResidue to transform
		core::Size water_begin(adding_water.conformation().chain_begin(adding_water.num_chains()));
		core::conformation::UltraLightResidue new_wat_ulr( adding_water.residue(water_begin).get_self_ptr() );
		/*
		//randomize rotation
		core::Real angle = 360;
		numeric::xyzMatrix<core::Real> rotation(
		numeric::z_rotation_matrix_degrees( angle*numeric::random::uniform() ) * (
		numeric::y_rotation_matrix_degrees( angle*numeric::random::uniform() ) *
		numeric::x_rotation_matrix_degrees( angle*numeric::random::uniform() ) ));
		*/

		numeric::xyzMatrix<core::Real> rotation( numeric::random::random_rotation() );

		//add to pose at calculated location
		TR << "Adding water at " << wi.centroid.to_string() << std::endl;

		new_wat_ulr.transform(rotation, wi.centroid, new_wat_ulr.center());
		new_wat_ulr.update_conformation(adding_water.conformation());
		adding_water.pdb_info()->chain(water_begin, water_chain_);
		adding_water.pdb_info()->number(water_begin, water_num);
	}
	adding_water.pdb_info()->obsolete(false); // ligand chain gets changed if this isn't set.
	//pose = adding_water;
} //apply()


utility::vector1<water_info>
PlaceWaters::check_possible_hbonds(utility::vector1<core::Vector> const & vs, core::Real inner_radius, core::Real outer_radius){

	utility::vector1<water_info> bridging;
	for ( core::Vector v : vs ) {
		water_info wi;
		wi.centroid = v;
		find_nearby_hbonds(h_grid_,v, wi, inner_radius, outer_radius);
		check_other_waters(v, vs, wi);
		if ( prot_lig_hbonds(wi) ) {
			bridging.push_back(wi);
			print_contacts(wi);
		}
	}

	return bridging;
}

void
PlaceWaters::output_water_coords_file(utility::vector1<core::Vector> const & waters, std::string outfile){
	std::ofstream fout;
	fout.open(outfile, std::ios::out);
	if ( !fout.is_open() ) {
		TR << "Can't open " << outfile << std::endl;
		utility_exit();
	}
	int resnum(0);
	for ( core::Vector water : waters ) {
		++resnum;
		core::Real x_val (water[0]); core::Real y_val (water[1]) ; core::Real z_val (water[2]);
		char line[100];
		std::sprintf(line,"%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n",
			"HETATM",  // record name
			resnum, // atom serial number
			"O",  // atom name
			"", // alt location
			"HOH", //residue name
			"W",  // chain id
			resnum,  // residue sequence number
			"",  // insertion code
			x_val, //x
			y_val, // y
			z_val,  // z
			1.00, // occupancy
			0.00,  // temp factor
			"O", // element symbol
			"" // charge
		);
		fout << line;
	} // waters loop
	fout.close();
}

void
PlaceWaters::print_contacts(water_info const & wi){
	TR << "Ligand Acceptors: ";
	for ( core::Vector v : wi.ligand_acceptors ) {
		TR << v.to_string() << " ";
	}
	TR << "" << std::endl;

	TR << "Ligand Donors: ";
	for ( core::Vector v : wi.ligand_donors ) {
		TR << v.to_string() << " ";
	}
	TR << "" << std::endl;

	TR << "Protein Acceptors: ";
	for ( core::Vector v : wi.protein_acceptors ) {
		TR << v.to_string() << " ";
	}
	TR << "" << std::endl;

	TR << "Protein Donors: ";
	for ( core::Vector v : wi.protein_donors ) {
		TR << v.to_string() << " ";
	}
	TR << "" << std::endl;

}


void
PlaceWaters::check_other_waters( core::Vector coords, const utility::vector1<core::Vector> & waters, water_info & wi){

	core::Real inner_radius = 2.25;
	core::Real outer_radius = 12.25;
	for ( core::Vector wat : waters ) {
		core::Real dist = coords.distance_squared(wat);
		if ( dist > inner_radius && dist < outer_radius ) {
			wi.other_waters.push_back(wat);
		}
	}
}

bool
PlaceWaters::prot_lig_hbonds(water_info const & wi){

	bool ligand_contacts = (wi.ligand_donors.size() || wi.ligand_acceptors.size());
	bool protein_contacts = (wi.protein_donors.size() || wi.protein_acceptors.size());

	return (ligand_contacts && protein_contacts);
}

void
PlaceWaters::find_nearby_hbonds(
	core::grid::CartGrid<core::Real> const & grid,
	core::Vector const coords,
	water_info & wi,
	core::Real inner_radius, //do we need this?
	core::Real outer_radius

){

	core::Real inner_radius2 = inner_radius*inner_radius;
	core::Real outer_radius2 = outer_radius*outer_radius;
	core::Size x_count(0);
	core::Size y_count(0);
	core::Size z_count(0);
	grid.getNumberOfPoints(x_count,y_count,z_count);
	debug_assert( x_count >0 && y_count > 0 && z_count > 0 );

	core::Vector vector_radius (outer_radius);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid.gridpt(coords - outer_radius);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid.gridpt(coords + outer_radius);

	for ( int x_index = std::max(0,grid_min.x()); x_index <= std::min(int(x_count)-1,grid_max.x()); ++x_index ) {
		for ( int y_index = std::max(0,grid_min.y()); y_index <= std::min(int(y_count)-1,grid_max.y()); ++y_index ) {
			for ( int z_index = std::max(0,grid_min.z()); z_index <= std::min(int(z_count)-1,grid_max.z()); ++z_index ) {
				core::grid::CartGrid<core::Real>::GridPt point(x_index, y_index, z_index);
				core::Vector box_center = grid.coords(point);
				if ( box_center.distance_squared(coords) >= inner_radius2 &&
						box_center.distance_squared(coords) <= outer_radius2 &&
						grid.getValue(box_center) != 0 ) {
					update_info(box_center, grid.getValue(box_center), wi);
				}
			}
		}
	}
}

void
PlaceWaters::update_info( core::Vector box_center, core::Real val, water_info & wi){
	if ( val == -1.0 ) wi.ligand_acceptors.push_back(box_center) ;
	else if ( val == -2.0 ) wi.protein_acceptors.push_back(box_center) ;
	else if  ( val == -3.0 ) wi.ligand_donors.push_back(box_center) ;
	else if ( val == -4.0 ) wi.protein_donors.push_back(box_center) ;
}

//returns true if all values in specified section  = value, false otherwise
bool
PlaceWaters::check_sphere(
	core::grid::CartGrid<core::Real> const & grid,
	core::Vector const & coords,
	core::Real radius,
	core::Real value,
	core::Vector & hit
){
	return check_ring( grid, coords, 0, radius, value, hit);
}

bool
PlaceWaters::check_ring(
	core::grid::CartGrid<core::Real> const & grid,
	core::Vector const & coords,
	core::Real inner_radius,
	core::Real outer_radius,
	core::Real value,
	core::Vector & hit
) {
	core::Real inner_radius2 = inner_radius*inner_radius;
	core::Real outer_radius2 = outer_radius*outer_radius;
	core::Size x_count(0);
	core::Size y_count(0);
	core::Size z_count(0);
	grid.getNumberOfPoints(x_count,y_count,z_count);
	debug_assert( x_count >0 && y_count > 0 && z_count > 0 );

	core::Vector vector_radius (outer_radius);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid.gridpt(coords - outer_radius);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid.gridpt(coords + outer_radius);

	for ( int x_index = std::max(0,grid_min.x()); x_index <= std::min(int(x_count)-1,grid_max.x()); ++x_index ) {
		for ( int y_index = std::max(0,grid_min.y()); y_index <= std::min(int(y_count)-1,grid_max.y()); ++y_index ) {
			for ( int z_index = std::max(0,grid_min.z()); z_index <= std::min(int(z_count)-1,grid_max.z()); ++z_index ) {
				core::grid::CartGrid<core::Real>::GridPt point(x_index, y_index, z_index);
				core::Vector box_center = grid.coords(point);
				if ( box_center.distance_squared(coords) >= inner_radius2 &&
						box_center.distance_squared(coords) <= outer_radius2 &&
						grid.getValue(box_center) != value ) {
					hit = box_center;
					return false;
				}
			}
		}
	}

	return true;

}

void
PlaceWaters::cluster_groups(utility::vector1<core::Vector> & vs, core::Real group_cutoff, core::Real merge_cutoff){

	utility::vector1<core::Vector> means;
	utility::vector1<utility::vector1 <core::Vector>> collections;
	core::Size num_groups(0);

	for ( core::Vector v : vs ) {
		core::Size closest = get_closest(means, v, group_cutoff);

		//didn't find group within cutoff distance
		if ( closest == 0 ) {
			num_groups++;
			closest = num_groups;
			means.push_back(v);
			utility::vector1<core::Vector> new_group;
			collections.push_back(new_group);
		}

		collections[closest].push_back(v);
		means[closest] = recalculate_mean(collections[closest]);
	}

	merge_close(means, merge_cutoff);
	vs = means;

	TR << "Number of void centers after merging: " << num_groups << std::endl;

}

void
PlaceWaters::merge_close(utility::vector1<core::Vector> & centers, core::Real merge_cutoff){

	//while there's another center within 2 angstroms, merge
	core::Vector a(0);
	core::Vector b(0);
	while ( too_close(centers, a, b, merge_cutoff) > 0 ) {
		TR << "Merging " << a.to_string() << " and " << b.to_string() << std::endl;
		utility::vector1<core::Vector> temp;
		core::Vector avg = (a + b) / 2.0;
		for ( core::Vector c : centers ) {
			if ( c != a && c != b ) {
				temp.push_back(c);
			}
		}
		temp.push_back(avg);
		centers = temp;
	}
}

core::Size
PlaceWaters::too_close(utility::vector1<core::Vector> & centers,
	core::Vector & a,
	core::Vector & b,
	core::Real cutoff){

	core::Real cutoff2 = cutoff*cutoff;
	core::Real min_dist = 10000;

	core::Size tc = 0;
	for ( core::Vector c : centers ) {
		for ( core::Vector o : centers ) {
			if ( c == o ) {
				continue;
			}

			core::Size cur_dist = c.distance_squared(o);
			if ( cur_dist < min_dist && cur_dist < cutoff2 ) {
				a = c;
				b = o;
				min_dist = cur_dist;
				tc++;
			}
		}
	}
	return tc;
}

core::Size
PlaceWaters::get_closest(utility::vector1<core::Vector> const & means, core::Vector const & v, core::Real group_cutoff){
	core::Size closest = 0;
	core::Size i = 1;
	core::Real dist = 100000;

	for ( core::Vector mean : means ) {
		core::Real mean_dist = std::sqrt(v.distance_squared(mean));
		if ( mean_dist < group_cutoff && mean_dist < dist ) {
			dist = mean_dist;
			closest = i;
		}
		i++;
	}

	return closest;

}

core::Vector
PlaceWaters::recalculate_mean(utility::vector1<core::Vector> col){

	core::Real x_sum = 0;
	core::Real y_sum = 0;
	core::Real z_sum = 0;

	for ( core::Vector v : col ) {
		x_sum += v.x();
		y_sum += v.y();
		z_sum += v.z();
	}

	return  core::Vector(x_sum/col.size(), y_sum/col.size(), z_sum/col.size() );
}

utility::vector1<core::Vector>
PlaceWaters::empty_spaces(core::Size cutoff){
	core::Size x_count(0);
	core::Size y_count(0);
	core::Size z_count(0);
	void_grid_.getNumberOfPoints(x_count,y_count,z_count);
	debug_assert( x_count >0 && y_count > 0 && z_count > 0 );

	utility::vector1<core::Vector> void_centers;
	core::Vector v;
	core::Size empty = 0;
	core::Size empty2 = 0;

	for ( int x_index = 0; x_index <= int(x_count)-1; ++x_index ) {
		for ( int y_index = 0; y_index <= int(y_count)-1; ++y_index ) {
			for ( int z_index = 0; z_index <= int(z_count)-1; ++z_index ) {
				core::grid::CartGrid<core::Real>::GridPt point(x_index, y_index, z_index);
				core::Vector cur_point = void_grid_.coords(point);
				if ( void_grid_.getValue(cur_point) == 0.0 ) {
					empty++;
				}
				if ( void_grid_.getValue(cur_point) == 0.0 && check_sphere(void_grid_, cur_point, cutoff, 0.0, v) ) {
					empty2++;
					void_centers.push_back(cur_point);
				}
			}
		}
	}

	TR << "Empty spaces in grid: " << empty << std::endl;
	TR << "Remaining coordinates: " << empty2 << std::endl;
	return void_centers;

}

void
PlaceWaters::initialize_grids( core::Vector const & center, core::Real width,core::Real resolution){

	auto num_pts = static_cast<core::Size>(width/resolution);

	core::Real const grid_halfwidth = width / 2.0;

	//void_grid_
	void_grid_.setBase(
		center.x() - grid_halfwidth,
		center.y() - grid_halfwidth,
		center.z() - grid_halfwidth
	);
	void_grid_.setDimensions(num_pts,num_pts,num_pts,resolution,resolution,resolution);
	void_grid_.setupZones();
	void_grid_.zero(); //set all to 0


	//h grid
	h_grid_.setBase(
		center.x() - grid_halfwidth,
		center.y() - grid_halfwidth,
		center.z() - grid_halfwidth
	);

	h_grid_.setDimensions(num_pts,num_pts,num_pts,resolution,resolution,resolution);
	h_grid_.setupZones();
	h_grid_.setFullOccupied(1.0); //set all to 1
}


//below functions shamelessly modified (stolen) from qsar/scoring_grid/SingleGrid
void
PlaceWaters::fill_occupied(core::pose::Pose const & pose){
	for ( core::Size residue_index = 1; residue_index <=pose.size(); ++residue_index ) {
		core::conformation::Residue const & residue(pose.residue(residue_index));
		for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
			core::chemical::AtomType const & cur_atom_type(residue.atom_type(atom_index));
			set_sphere(void_grid_, residue.xyz(atom_index),cur_atom_type.lj_radius(), 1.0);
		}
	}
}


//Go through heavy atoms and find regions where waters could be
void
PlaceWaters::detect_hydrogens(core::pose::Pose const & pose){
	core::Real sphere_radius = resolution_*0.5; //should ensure that only one coordinate is set
	for ( core::Size residue_index = 1; residue_index <=pose.size(); ++residue_index ) {
		bool is_ligand = false;
		core::conformation::Residue const & residue(pose.residue(residue_index));
		if ( !residue.is_protein() ) is_ligand = true;
		for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
			core::chemical::AtomType const & cur_atom_type(residue.atom_type(atom_index));
			//if acceptor found, value is -1
			if ( cur_atom_type.is_acceptor() ) {
				if ( is_ligand ) {
					set_sphere(h_grid_, residue.xyz(atom_index), sphere_radius,-1.0);
				} else {
					set_sphere(h_grid_, residue.xyz(atom_index), sphere_radius, -2.0 );
				}
			} else if ( cur_atom_type.is_polar_hydrogen() ) {
				if ( is_ligand ) {
					set_sphere(h_grid_, residue.xyz(atom_index), sphere_radius, -3.0);
				} else {
					set_sphere(h_grid_, residue.xyz(atom_index), sphere_radius, -4.0);
				}
			}
		}
	}
}

void
PlaceWaters::set_sphere(
	core::grid::CartGrid<core::Real> & grid,
	core::Vector const & coords,
	core::Real radius,
	core::Real value
){
	set_ring(grid, coords, 0, radius, value);
}

void
PlaceWaters::set_ring(
	core::grid::CartGrid<core::Real> & grid,
	core::Vector const & coords,
	core::Real inner_radius,
	core::Real outer_radius,
	core::Real value
){

	core::Real inner_radius2 = inner_radius*inner_radius;
	core::Real outer_radius2 = outer_radius*outer_radius;
	core::Size x_count(0);
	core::Size y_count(0);
	core::Size z_count(0);
	grid.getNumberOfPoints(x_count,y_count,z_count);
	debug_assert( x_count >0 && y_count > 0 && z_count > 0 );

	core::Vector vector_radius (outer_radius);
	core::grid::CartGrid<core::Real>::GridPt grid_min = grid.gridpt(coords - outer_radius);
	core::grid::CartGrid<core::Real>::GridPt grid_max = grid.gridpt(coords + outer_radius);

	for ( int x_index = std::max(0,grid_min.x()); x_index <= std::min(int(x_count)-1,grid_max.x()); ++x_index ) {
		for ( int y_index = std::max(0,grid_min.y()); y_index <= std::min(int(y_count)-1,grid_max.y()); ++y_index ) {
			for ( int z_index = std::max(0,grid_min.z()); z_index <= std::min(int(z_count)-1,grid_max.z()); ++z_index ) {
				core::grid::CartGrid<core::Real>::GridPt point(x_index, y_index, z_index);
				core::Vector box_center = grid.coords(point);
				if ( box_center.distance_squared(coords) >= inner_radius2
						&& box_center.distance_squared(coords) <= outer_radius2
						&& grid.getValue(box_center) != value
						) {
					grid.setValue(point, value);
				}
			}
		}
	}
}

core::Real
PlaceWaters::grid_size(core::conformation::UltraLightResidueOP lig){
	return 2.0*core::Real(lig->max_dist_to_center());
}

core::Real
PlaceWaters::get_grid_size(){
	return grid_size_;
}

core::Real
PlaceWaters::get_resolution(){
	return resolution_;
}

core::Real
PlaceWaters::get_cutoff(){
	return cutoff_;
}

core::Real
PlaceWaters::get_group_cutoff(){
	return group_cutoff_;
}

core::Real
PlaceWaters::get_merge_cutoff(){
	return merge_cutoff_;
}

core::Real
PlaceWaters::get_hbonds_inner(){
	return hbonds_inner_;
}

core::Real
PlaceWaters::get_hbonds_outer(){
	return hbonds_outer_;
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
PlaceWaters::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& /*datamap*/
) {

	if ( tag->getName() != "PlaceWaters" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");
	}
	if ( ! tag->hasOption("ligand_chain") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'PlaceWaters' mover requires 'ligand_chain' tag");
	if ( ! tag->hasOption("water_type") ) { water_type_ = "TP3"; } else { water_type_ = tag->getOption<std::string>("water_type"); }
	if ( tag->hasOption("output_water_coords") ) { output_water_coords_ = tag->getOption<bool>("output_water_coords", false); }

	ligand_chain_ = tag->getOption<std::string>("ligand_chain");
	water_chain_ = tag->getOption<char>("water_chain", 'W');
	resolution_ = tag->getOption<core::Real>("resolution", 0.15);
	cutoff_ = tag->getOption<core::Real>("cutoff", 1.5);
	grid_size_ = tag->getOption<core::Real>("grid_size", 0.0);
	cycles_ = tag->getOption<core::Size>("cycles", 100);
	group_cutoff_ = tag->getOption<core::Real>("group_cutoff", 3.75);
	merge_cutoff_ = tag->getOption<core::Real>("merge_cutoff", 3.0);
	hbonds_inner_ = tag->getOption<core::Real>("hbonds_inner", 1.5);
	hbonds_outer_ = tag->getOption<core::Real>("hbonds_outer", 5.0);
}

void PlaceWaters::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	attlist
		+ XMLSchemaAttribute::required_attribute("ligand_chain", xs_string, "The ligand chain, specified as the PDB chain ID")
		+ XMLSchemaAttribute::attribute_w_default("water_chain", xsct_char, "The chain ID the water will be added to--should be different than ligand_chain. W is recommended.", "W")
		+ XMLSchemaAttribute::attribute_w_default("water_type", xs_string, "Water type--this can either be HOH, TP3, TP2 from Rosetta database. TP3 is recommended.", "TP3")
		+ XMLSchemaAttribute::attribute_w_default("grid_size", xsct_real, "Grid size", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("resolution", xsct_real, "Grid resolution", "0.15")
		+ XMLSchemaAttribute::attribute_w_default("cutoff", xsct_real, "Remove all coordinates within this distance of the ligand or protein before clustering.", "1.0")
		+ XMLSchemaAttribute::attribute_w_default("cycles", xsct_non_negative_integer, "Cycles for water transforms, if optimizing waters.", "100")
		+ XMLSchemaAttribute::attribute_w_default("group_cutoff", xsct_real, "Cutoff to use for grouping void areas","3.5")
		+ XMLSchemaAttribute::attribute_w_default("merge_cutoff", xsct_real, "Cutoff to use for merging nearby void clusters", "2.0")
		+ XMLSchemaAttribute::attribute_w_default("hbonds_inner", xsct_real, "Inner cutoff for nearby hbonds", "1.5")
		+ XMLSchemaAttribute::attribute_w_default("hbonds_outer", xsct_real, "Outer cutoff for nearby hbonds", "3.5")
		+ XMLSchemaAttribute( "output_water_coords", xsct_rosetta_bool, "Output water coordinates at each step to PDB files for visual observation. "
		"Default is false. Files will be labeled 1.emptyspaces_coords.pdb, 2.clusters_coords.pdb, 3.final_coords.pdb .");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Places waters throughout the interface", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PlaceWaters::fresh_instance() const
{
	return protocols::moves::MoverOP( new PlaceWaters );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PlaceWaters::clone() const
{
	return protocols::moves::MoverOP( new PlaceWaters( *this ) );
}

std::string PlaceWaters::get_name() const {
	return mover_name();
}

std::string PlaceWaters::mover_name() {
	return "PlaceWaters";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
PlaceWatersCreator::create_mover() const
{
	return protocols::moves::MoverOP( new PlaceWaters );
}

std::string
PlaceWatersCreator::keyname() const
{
	return PlaceWaters::mover_name();
}

void PlaceWatersCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PlaceWaters::provide_xml_schema( xsd );
}


} //protocols
} //ligand_docking
