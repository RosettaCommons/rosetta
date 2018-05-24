// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/guffysl/zinc_statistic_generator.cc
/// @brief Builds zinc ions off coordinating residues in different positions to build a library of allowed zinc positions relative to the parent residue
/// @details The input pose will be a single residue which will be mutated to zinc-coordinating residues and tested in different rotamers
/// @author Sharon Guffy

//Devel
#include <devel/init.hh>
#include <protocols/sewing/hashing/LigandBindingResPlacer.hh>
//Protocols

//Core
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/Stub.hh>
#include <core/types.hh>

//Basic
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//Utility
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/ReferenceCount.hh>
//Numeric
#include <numeric/xyzVector.hh>



//tracers
using basic::Error;
using basic::Warning;
static basic::Tracer TR("apps.pilot.guffysl.zinc_statistic_generator");


namespace local {

//Local options
basic::options::BooleanOptionKey const allow_his( "allow_his" ); //default true
basic::options::BooleanOptionKey const allow_de( "allow_de" ); //default false
basic::options::StringOptionKey const secstruct( "secstruct" ); //default "H"
basic::options::RealOptionKey const sampling_range( "sampling_range" ); //default 1
basic::options::RealOptionKey const sampling_increment( "sampling_increment" ); //default 0.1
basic::options::StringOptionKey const output_name( "output_name" ); //default "stats.txt"

}

//Returns a stub for this residue
//core::pose::motif::get_backbone_reference_frame( N, CA, C)
//core::kinematics::Stub ligand_stub( CA, N, CA, C );
//Then call ligand_stub.global2local( zn_xyz ) to get the zinc's local coordinates for this rotamer


/*
struct ZnCoordInfo
{
numeric::xyzVector< core::Real > local_coords;
std::string coord_res_name; //This should distinguish HIS from HISD
std::string coord_atom_index;
utility::fixedsizearray1< core::Real, 4 > chis;
};
*/
namespace zinc_stats {

class ZincStatisticGenerator: public utility::pointer::ReferenceCount {

public:
	ZincStatisticGenerator(){
	}

	~ZincStatisticGenerator(){
	}

	void
	initialize_from_options(){

		sampling_range_ = basic::options::option[ local::sampling_range ].value();
		sampling_increment_ = basic::options::option[ local::sampling_increment ].value();
		//If these are negative, it will all crash! Make sure that doesn't happen
		if ( sampling_range_ < 0 ) sampling_range_ = sampling_range_ * -1;
		if ( sampling_increment_ < 0 ) sampling_increment_ = sampling_increment_ * -1;

		allow_his_ = basic::options::option[ local::allow_his ].value();

		allow_de_ = basic::options::option[ local::allow_de ].value();

		allowed_secstruct_ = basic::options::option[ local::secstruct ].value();

		base_filename_ = basic::options::option[ local::output_name ].value();

	}

	numeric::xyzVector< core::Real >
	get_local_zn_coords( core::conformation::Residue const & res, numeric::xyzVector< core::Real > zn_coords ){
		//First make the stub
		core::kinematics::Stub ligand_stub( res.xyz( 2 ), res.xyz( 1 ), res.xyz( 2 ), res.xyz( 3 ) );
		return ligand_stub.global2local( zn_coords );
	}

	utility::vector1< protocols::sewing::hashing::LigandCoordInfo >
	get_local_zn_positions_for_residue( core::pose::Pose & pose, core::Size ligand_resnum, char secstruct ){

		//Set up variables
		core::pack::rotamers::SingleResidueRotamerLibraryCOP rotlib;
		utility::vector1< protocols::sewing::hashing::LigandCoordInfo > local_zn_positions;

		//Determine backbone torsion angles to sample
		utility::vector1< utility::fixedsizearray1< core::Real, 5 > > secstruct_phi_psi;
		//Ideal phi/psi:
		//Helix: -60, -60
		//Strand: -135, +135
		//Loop? alpha+beta plus:
		//* +60, +60
		//* -60, -30
		//* -90, 0
		//* -60, 120

		if ( secstruct == 'E' ) {
			utility::fixedsizearray1< core::Real, 5 > bb_angles;
			bb_angles[ 1 ] = -135;
			bb_angles[ 2 ] = 135;
			secstruct_phi_psi.push_back( bb_angles );
		} else if ( secstruct == 'H' ) {
			utility::fixedsizearray1< core::Real, 5 > bb_angles;
			bb_angles[ 1 ] = -60;
			bb_angles[ 2 ] = -60;
			secstruct_phi_psi.push_back( bb_angles );
		} else { //Treat as a loop, check all possible combinations
			utility::fixedsizearray1< core::Real, 5 > bb_angles;
			bb_angles[ 1 ] = -135;
			bb_angles[ 2 ] = 135;
			secstruct_phi_psi.push_back( bb_angles );
			bb_angles[ 1 ] = -60;
			bb_angles[ 2 ] = -60;
			secstruct_phi_psi.push_back( bb_angles );
			bb_angles[ 1 ] = 60;
			bb_angles[ 2 ] = 60;
			secstruct_phi_psi.push_back( bb_angles );
			bb_angles[ 1 ] = -60;
			bb_angles[ 2 ] = -30;
			secstruct_phi_psi.push_back( bb_angles );
			bb_angles[ 1 ] = -90;
			bb_angles[ 2 ] = 0;
			secstruct_phi_psi.push_back( bb_angles );
			bb_angles[ 1 ] = -60;
			bb_angles[ 2 ] = 120;
			secstruct_phi_psi.push_back( bb_angles );
		}

		if ( allow_his_ ) {


			//Try everything for histidine (ND1 coordinating)
			mutate_residue( pose, ligand_resnum, "HIS" );
			//Get rotamer library
			rotlib = core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue( ligand_resnum ).type() );
			for ( utility::fixedsizearray1< core::Real, 5 > bb_angles: secstruct_phi_psi ) {
				utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > rotamers = std::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const >( rotlib )->get_all_rotamer_samples(
					bb_angles );
				//Loop through most likely rotamers of the specified secstruct
				for ( core::pack::dunbrack::DunbrackRotamerSampleData rotamer: rotamers ) {
					//If the rotamer probability is too low, ignore it
					if ( rotamer.probability() < 0.05 ) {
						continue;
					}
					//Within these rotamers, loop through +- specified # of standard deviations
					//The number of chi angles will determine the number of nested loops--recursion will be the cleanest way to deal with this
					utility::fixedsizearray1< core::Real, 4 > all_chi_angle_values;
					recursively_model_rotamer_chis( pose, ligand_resnum, 1, all_chi_angle_values, rotamer, local_zn_positions );
				} //End iterate over rotamers
			} //End iterate over phi/psi


			//Try everything for HISD
			mutate_residue( pose, ligand_resnum, "HIS_D" );
			//Get rotamer library
			rotlib = core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue( ligand_resnum ).type() );
			for ( utility::fixedsizearray1< core::Real, 5 > bb_angles: secstruct_phi_psi ) {
				utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > rotamers = std::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const >( rotlib )->get_all_rotamer_samples(
					bb_angles );
				//Loop through most likely rotamers of the specified secstruct
				for ( core::pack::dunbrack::DunbrackRotamerSampleData rotamer: rotamers ) {
					//If the rotamer probability is too low, ignore it
					if ( rotamer.probability() < 0.05 ) {
						continue;
					}
					//Within these rotamers, loop through +- specified # of standard deviations
					//The number of chi angles will determine the number of nested loops--recursion will be the cleanest way to deal with this
					utility::fixedsizearray1< core::Real, 4 > all_chi_angle_values;
					recursively_model_rotamer_chis( pose, ligand_resnum, 1, all_chi_angle_values, rotamer, local_zn_positions );
				} //End iterate over rotamers
			} //End iterate over phi/psi

		}

		//If asp/glu not allowed, we're done
		if ( !allow_de_ ) {
			return local_zn_positions;
		}

		//Otherwise we'll continue with both asp and glu
		mutate_residue( pose, ligand_resnum, "ASP" );
		//Get rotamer library
		rotlib = core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue( ligand_resnum ).type() );
		for ( utility::fixedsizearray1< core::Real, 5 > bb_angles: secstruct_phi_psi ) {
			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > rotamers = std::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const >( rotlib )->get_all_rotamer_samples(
				bb_angles );
			//Loop through most likely rotamers of the specified secstruct
			for ( core::pack::dunbrack::DunbrackRotamerSampleData rotamer: rotamers ) {
				//If the rotamer probability is too low, ignore it
				if ( rotamer.probability() < 0.05 ) {
					continue;
				}
				//Within these rotamers, loop through +- specified # of standard deviations
				//The number of chi angles will determine the number of nested loops--recursion will be the cleanest way to deal with this
				utility::fixedsizearray1< core::Real, 4 > all_chi_angle_values;
				recursively_model_rotamer_chis(pose, ligand_resnum, 1, all_chi_angle_values, rotamer, local_zn_positions );
			} //End iterate over rotamers
		} //End iterate over phi/psi
		//Build zinc off outer positions of both oxygens (not bidentate position)
		mutate_residue( pose, ligand_resnum, "GLU" );
		//Get rotamer library
		rotlib = core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue( ligand_resnum ).type() );
		for ( utility::fixedsizearray1< core::Real, 5 > bb_angles: secstruct_phi_psi ) {
			utility::vector1< core::pack::dunbrack::DunbrackRotamerSampleData > rotamers = std::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const >( rotlib )->get_all_rotamer_samples(
				bb_angles );
			//Loop through most likely rotamers of the specified secstruct
			for ( core::pack::dunbrack::DunbrackRotamerSampleData rotamer: rotamers ) {
				//If the rotamer probability is too low, ignore it
				if ( rotamer.probability() < 0.05 ) {
					continue;
				}
				//Within these rotamers, loop through +- specified # of standard deviations
				//The number of chi angles will determine the number of nested loops--recursion will be the cleanest way to deal with this
				utility::fixedsizearray1< core::Real, 4 > all_chi_angle_values;
				recursively_model_rotamer_chis( pose, ligand_resnum, 1, all_chi_angle_values, rotamer,  local_zn_positions );
			} //End iterate over rotamers
		} //End iterate over phi/psi
		//Build zinc off outer positions of both oxygens (not bidentate position)
		return local_zn_positions;
	}

	void
	recursively_model_rotamer_chis(
		core::pose::Pose & pose,
		core::Size resnum,
		core::Size current_chi,
		utility::fixedsizearray1< core::Real, 4 > & all_chi_angle_values, //Reference prevents an explosion of copies (just reassign the value)
		core::pack::dunbrack::DunbrackRotamerSampleData current_rotamer,
		utility::vector1< protocols::sewing::hashing::LigandCoordInfo > & output )
	{
		if ( current_chi > current_rotamer.nchi() ) { //All the chi angles have been computed
			//Set all chi angles that this residue doesn't have to 0 to avoid weird output.
			for ( core::Size i = current_chi; i <= 4; ++i ) {
				all_chi_angle_values[ i ] = 0;
			}
			//First, set chi angles
			for ( core::Size i = 1; i <= pose.residue( resnum ).nchi(); ++ i ) {
				pose.set_chi( i, resnum, all_chi_angle_values[ i ] );
			}
			//For each, build ideal position of zinc off each side chain nitrogen (actually we don't have to build it, we just need to know its XYZ coords)
			//If we're looking at ASP/GLU, there will be 2 possible zinc positions (we'll only look at the 2 outer positions and ignore the bidentate position)
			//If we're looking at HIS or HIS_D, there will also be 2 possible zinc positions (one from each nitrogen)
			if ( pose.residue( resnum ).type().base_name() == "HIS" ) {
				protocols::sewing::hashing::LigandCoordInfo nd1_zinc_coord_info;

				nd1_zinc_coord_info.coord_res_name = "HIS";
				nd1_zinc_coord_info.coord_atom_index = pose.residue( resnum ).atom_index( "ND1" );
				nd1_zinc_coord_info.chis = all_chi_angle_values;
				nd1_zinc_coord_info.ligand_atom_index = 1;

				numeric::xyzVector< core::Real > nd1_coords = pose.residue( resnum ).atom( "ND1" ).xyz();
				numeric::xyzVector< core::Real > cg_coords = pose.residue( resnum ).atom( "CG" ).xyz();
				numeric::xyzVector< core::Real > ce1_coords = pose.residue( resnum ).atom( "CE1" ).xyz();
				numeric::xyzVector< core::Real > cd2_coords = pose.residue( resnum ).atom( "CD2" ).xyz();

				//In this case, we can determine the plane of the ring from ND1, CG, and CE1
				numeric::xyzVector< core::Real > plane_normal = ( cg_coords - nd1_coords ).cross( ce1_coords - nd1_coords );
				//I think the plane normal was computed correctly


				//THis part is probably wrong
				//NE2: Subtract midpoint of CD2 and CE1
				//ND1: Subtract midpoint of CG and CE1
				numeric::xyzVector< core::Real > nd1_zinc = nd1_coords - numeric::midpoint( cg_coords, ce1_coords );

				nd1_zinc = nd1_zinc - ( nd1_zinc.projected_parallel( plane_normal ) );
				nd1_zinc = nd1_zinc.normalized( 2.2 );
				//So I think we actually need to add this to the original coords to get the right answer
				nd1_zinc = nd1_zinc + nd1_coords;
				//OK, I think this may be the problem in this vein

				nd1_zinc_coord_info.local_coords = get_local_zn_coords( pose.residue( resnum ), nd1_zinc );
				output.push_back( nd1_zinc_coord_info );

			} else if ( pose.residue( resnum ).type().base_name() == "HIS_D" ) {
				protocols::sewing::hashing::LigandCoordInfo ne2_zinc_coord_info;
				ne2_zinc_coord_info.coord_res_name = "HIS_D";
				ne2_zinc_coord_info.coord_atom_index = pose.residue( resnum ).atom_index( "NE2" );
				ne2_zinc_coord_info.chis = all_chi_angle_values;
				ne2_zinc_coord_info.ligand_atom_index = 1;


				numeric::xyzVector< core::Real > nd1_coords = pose.residue( resnum ).atom( "ND1" ).xyz();
				numeric::xyzVector< core::Real > ne2_coords = pose.residue( resnum ).atom( "NE2" ).xyz();
				numeric::xyzVector< core::Real > cg_coords = pose.residue( resnum ).atom( "CG" ).xyz();
				numeric::xyzVector< core::Real > ce1_coords = pose.residue( resnum ).atom( "CE1" ).xyz();
				numeric::xyzVector< core::Real > cd2_coords = pose.residue( resnum ).atom( "CD2" ).xyz();


				numeric::xyzVector< core::Real > plane_normal = ( cg_coords - nd1_coords ).cross( ce1_coords - nd1_coords );



				numeric::xyzVector< core::Real > ne2_zinc = ne2_coords - numeric::midpoint( cd2_coords, ce1_coords );
				ne2_zinc = ne2_zinc - ( ne2_zinc.projected_parallel( plane_normal ) );
				ne2_zinc = ne2_zinc.normalized( 2.2 );
				ne2_zinc = ne2_zinc + ne2_coords;
				ne2_zinc_coord_info.local_coords = get_local_zn_coords( pose.residue( resnum ), ne2_zinc );
				output.push_back( ne2_zinc_coord_info );

			} else if ( pose.residue( resnum ).type().base_name() == "ASP" ) {
				protocols::sewing::hashing::LigandCoordInfo od1_zinc_coord_info;
				protocols::sewing::hashing::LigandCoordInfo od2_zinc_coord_info;

				od1_zinc_coord_info.coord_res_name = "ASP";
				od2_zinc_coord_info.coord_res_name = "ASP";
				od1_zinc_coord_info.coord_atom_index = pose.residue( resnum ).atom_index( "OD1" );
				od2_zinc_coord_info.coord_atom_index = pose.residue( resnum ).atom_index( "OD2" );
				od1_zinc_coord_info.chis = all_chi_angle_values;
				od2_zinc_coord_info.chis = all_chi_angle_values;


				od1_zinc_coord_info.ligand_atom_index = 1;
				od2_zinc_coord_info.ligand_atom_index = 1;

				numeric::xyzVector< core::Real > od1_coords = pose.residue( resnum ).atom( "OD1" ).xyz();
				numeric::xyzVector< core::Real > od2_coords = pose.residue( resnum ).atom( "OD2" ).xyz();

				//We'll need to calculate the plane normal vector from CG, OD1, and OD2
				numeric::xyzVector< core::Real > cg_coords = pose.residue( resnum ).atom( "CG" ).xyz();
				numeric::xyzVector< core::Real > plane_normal = ( od1_coords - cg_coords ).cross( od2_coords - cg_coords );

				//and then subtract the projection onto that vector from the zinc vectors
				numeric::xyzVector< core::Real > od1_zinc = ( od1_coords - od2_coords );
				numeric::xyzVector< core::Real > od2_zinc = ( od2_coords - od1_coords );

				od1_zinc = od1_zinc - ( od1_zinc.projected_parallel( plane_normal ) );
				od2_zinc = od2_zinc - ( od2_zinc.projected_parallel( plane_normal ) );
				od1_zinc = od1_zinc.normalized( 2.2 );
				od2_zinc = od2_zinc.normalized( 2.2 );
				od1_zinc = od1_zinc + od1_coords;
				od2_zinc = od2_zinc + od2_coords;

				od1_zinc_coord_info.local_coords = get_local_zn_coords( pose.residue( resnum ), od1_zinc );
				od2_zinc_coord_info.local_coords = get_local_zn_coords( pose.residue( resnum ), od2_zinc );
				output.push_back( od1_zinc_coord_info );
				output.push_back( od2_zinc_coord_info );
			} else if ( pose.residue( resnum ).type().base_name() == "GLU" ) {
				protocols::sewing::hashing::LigandCoordInfo oe1_zinc_coord_info;
				protocols::sewing::hashing::LigandCoordInfo oe2_zinc_coord_info;

				oe1_zinc_coord_info.coord_res_name = "GLU";
				oe2_zinc_coord_info.coord_res_name = "GLU";
				oe1_zinc_coord_info.coord_atom_index = pose.residue( resnum ).atom_index( "OE1" );
				oe2_zinc_coord_info.coord_atom_index = pose.residue( resnum ).atom_index( "OE2" );
				oe1_zinc_coord_info.chis = all_chi_angle_values;
				oe2_zinc_coord_info.chis = all_chi_angle_values;


				oe1_zinc_coord_info.ligand_atom_index = 1;
				oe2_zinc_coord_info.ligand_atom_index = 1;

				numeric::xyzVector< core::Real > oe1_coords = pose.residue( resnum ).atom( "OE1" ).xyz();
				numeric::xyzVector< core::Real > oe2_coords = pose.residue( resnum ).atom( "OE2" ).xyz();

				//We'll need to calculate the plane normal vector from CD, OE1, and OE2
				numeric::xyzVector< core::Real > cd_coords = pose.residue( resnum ).atom( "CD" ).xyz();
				numeric::xyzVector< core::Real > plane_normal = ( oe1_coords - cd_coords ).cross( oe2_coords - cd_coords );

				//and then subtract the projection onto that vector from the zinc vectors
				numeric::xyzVector< core::Real > oe1_zinc = ( oe1_coords - oe2_coords );
				numeric::xyzVector< core::Real > oe2_zinc = ( oe2_coords - oe1_coords );



				//Don't know about this
				oe1_zinc = oe1_zinc - ( oe1_zinc.projected_parallel( plane_normal ) );
				oe2_zinc = oe2_zinc - ( oe2_zinc.projected_parallel( plane_normal ) );
				oe1_zinc = oe1_zinc.normalized( 2.2 );
				oe2_zinc = oe2_zinc.normalized( 2.2 );
				oe1_zinc = oe1_zinc + oe2_coords;
				oe2_zinc = oe2_zinc + oe2_coords;

				oe1_zinc_coord_info.local_coords = get_local_zn_coords( pose.residue( resnum ), oe1_zinc );
				oe2_zinc_coord_info.local_coords = get_local_zn_coords( pose.residue( resnum ), oe2_zinc );
				output.push_back( oe1_zinc_coord_info );
				output.push_back( oe2_zinc_coord_info );

			}
			//For each of those atoms, build a vector of length 2.2 that is 1) in the plane of the residue (-180/0/180) and b) has a 120 degree (+-?) angle about the atom
			//Check old match constraints for exact values
		} else {
			for ( core::Real zscore = -1 * sampling_range_; zscore <= sampling_range_; zscore = zscore + sampling_increment_ ) {
				//Set this chi angle to all the allowed values
				core::Real new_chi = current_rotamer.chi_mean()[ current_chi ] + ( zscore * current_rotamer.chi_sd()[ current_chi ] );
				all_chi_angle_values[ current_chi ] = new_chi;
				recursively_model_rotamer_chis( pose, resnum, current_chi + 1, all_chi_angle_values, current_rotamer, output );
			}
		}
	}

	void
	mutate_residue( core::pose::Pose & pose, core::Size res_position, std::string target_res_name ){
		core::chemical::ResidueTypeSetCOP restype_set( pose.conformation().residue_type_set_for_conf() );
		core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue(
			restype_set->name_map( target_res_name ), pose.residue( res_position ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( res_position ), *new_res, pose.conformation(), true );
		pose.replace_residue( res_position, *new_res, false );
		pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only( res_position );
	}

	std::map< char, utility::vector1< protocols::sewing::hashing::LigandCoordInfo > >
	iterate_over_residues( core::pose::Pose & pose ){
		std::map< char, utility::vector1< protocols::sewing::hashing::LigandCoordInfo > > secstruct_coords;
		for ( core::Size resnum = 2; resnum < pose.total_residue(); ++resnum ) { //Note that in the current implementation the pose only has one internal residue anyway
			//These are strings and therefore use 0 indexing
			for ( core::Size i = 0; i < allowed_secstruct_.size(); ++i ) {
				char secstruct = allowed_secstruct_.at( i );
				//utility::vector1< LigandCoordInfo > res_coord = get_local_zn_positions_for_residue( pose, resnum, secstruct );
				//secstruct_coords[ secstruct ].insert( secstruct_coords[ secstruct ].end(), res_coord.begin(), res_coord.end() );
				secstruct_coords[ secstruct ] = get_local_zn_positions_for_residue( pose, resnum, secstruct );
			}
		}
		return secstruct_coords;
	}

	void
	output_file( std::string prefix, utility::vector1< protocols::sewing::hashing::LigandCoordInfo > coords ){
		std::string output_filename = prefix + base_filename_;
		//If the output file exists, print a warning message and append to the file

		utility::io::ozstream output;
		//Print headers (tab-separated)
		//std::stringstream header;
		//header << "Coord_res_name\tCoord_atom_number\tLigand_atom_number\tLocal_x\tLocal_y\tLocal_z\tChi1\tChi2\tChi3\tChi4" << std::endl;

		if ( utility::file::file_exists( output_filename ) ) {
			TR << "WARNING: A file with the specified name exists! Overwriting file " << output_filename << std::endl;
		} else {
			//If the file is new, it will be created and the header will be printed
			TR << "Creating new file " << output_filename << std::endl;
		}
		//output.open_append_if_existed( output_filename, header ); //Header will only be printed if the file didn't already exist
		output.open( output_filename );
		output << prefix.substr( 0, 1 ) << std::endl;
		output << "ZN" << std::endl;
		output << "Coord_res_name\tCoord_atom_number\tLigand_atom_number\tLocal_x\tLocal_y\tLocal_z\tChi1\tChi2\tChi3\tChi4" << std::endl;
		//Iterate through coords
		for ( protocols::sewing::hashing::LigandCoordInfo coord: coords ) {
			//Output data in order: coord_res_name coord_atom_number ligand_atom_number local_coords chis (because the # of chis may vary)
			output << coord.coord_res_name << "\t" << coord.coord_atom_index << "\t" << "1" << "\t";
			//Print coordinates
			numeric::xyzVector< core::Real > lc = coord.local_coords;
			output << utility::to_string( lc.x() ) << "\t" << utility::to_string( lc.y() ) << "\t" << utility::to_string( lc.z() ) << "\t";

			//Print chi angles
			utility::fixedsizearray1< core::Real, 4 > chi = coord.chis;
			output << utility::to_string( chi[ 1 ] ) << "\t" << utility::to_string( chi[ 2 ] ) << "\t" << utility::to_string( chi[ 3 ] ) << "\t" << utility::to_string( chi[ 4 ] ) << std::endl;
		}

		//Close file
		output.flush();
		output.close();

	}

private:
	core::Real sampling_range_ = 1;
	core::Real sampling_increment_ = 0.1;
	bool allow_de_ = false;
	bool allow_his_ = true;
	std::string allowed_secstruct_ = "H";
	std::string base_filename_ = "stats.txt";

};


}//namespace zinc_stats
int main( int argc, char* argv[] ){
	try{

		//Set up options
		basic::options::option.add( local::allow_his, "Calculate local zinc coordinates for HIS and HIS_D rotamers?" ).def( false );
		basic::options::option.add( local::allow_de, "Calculate local zinc coordinates for aspartate and glutamate rotamers?" ).def( false );
		basic::options::option.add( local::secstruct, "String containing DSSP codes for which to compute local zinc coordinates (e.g. \"HL\" for loops and helices" ).def( "H" );
		basic::options::option.add( local::sampling_range, "How many standard deviations from the mean to sample chi angles?" ).def( 1.0 );
		basic::options::option.add( local::sampling_increment, "How fine ( in standard deviations) should the sampling of chi angles be?" ).def( 0.1 );
		basic::options::option.add( local::output_name, "Base name of output file (DSSP codes will be prepended) " ).def( "stats.txt" );

		devel::init( argc, argv );

		//BODY OF APP
		//Initialize variables
		//std::map< char, utility::vector1< LigandCoordInfo> > secstruct_coords;

		//Construct the pose ( will be a single residue )
		//in core/pose/annotated_sequence.hh
		core::pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "AAA", core::chemical::FA_STANDARD );

		//Run iterate_over_residues
		zinc_stats::ZincStatisticGenerator zincstat;
		zincstat.initialize_from_options();
		std::map< char, utility::vector1< protocols::sewing::hashing::LigandCoordInfo > > secstruct_coords = zincstat.iterate_over_residues( pose );
		//Output data
		for ( std::pair< char, utility::vector1< protocols::sewing::hashing::LigandCoordInfo > > coords: secstruct_coords ) {
			TR << "Writing file for " << utility::to_string( coords.first ) << std::endl;
			std::string prefix = utility::to_string( coords.first ) + "_";
			zincstat.output_file( prefix, coords.second );
		}
		//END BODY OF APP

		return 0;
	}
catch ( utility::excn::Exception const &e ){
	std::cout << "Caught exception " << e.msg() << std::endl;
	return -1;
}
}
