// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk

// Project Headers
#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/constants.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>

#include <numeric/constants.hh>
#include <numeric/conversions.hh>

//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>

//Auto Headers
#include <core/pose/annotated_sequence.hh>


using namespace core;
using namespace core::scoring;
using namespace core::scoring::hbonds;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( String, atom_type )
OPT_KEY( Real, occluding_radius )
OPT_KEY( Real, water_entropy_scaling_factor )

static basic::Tracer TR( "apps.pilot.johnk_parametrize_geometric_solvation.main" );

struct GridInfo
{
	core::Size xnum_points, ynum_points, znum_points;
	core::Real xstep, ystep, zstep;
	core::Real xorigin, yorigin, zorigin;
};


std::string setup_pdb_line( core::Vector const & atom_pos, core::Size const & atomnum, core::Real const & Bfac ) {
	std::stringstream pdbline;
	pdbline << "HETATM";
	if ( atomnum<10 ) pdbline << "    ";
	else if ( atomnum<100 ) pdbline << "   ";
	else if ( atomnum<1000 ) pdbline << "  ";
	else pdbline << " ";
	pdbline << atomnum << "  GR  GRD";
	pdbline<<" X   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<atom_pos.x()<<std::setw(8)<<atom_pos.y()<<std::setw(8)<<atom_pos.z();
	pdbline<<"  1.00"<<std::setw(6)<<std::fixed<<std::setprecision(2)<<Bfac;
	return pdbline.str();
}

// Dump the grid to a PDB file, with desired grid info as the B-factor
void dump_grid_pdb( std::string const & pdb_fname, GridInfo const & grid_info,
	std::vector < std::vector < std::vector <core::Real> > > const & grid_values ) {
	utility::io::ozstream outpdb_stream;
	outpdb_stream.open( pdb_fname, std::ios::out);
	core::Size atomnum = 1;

	core::Vector atom_pos(grid_info.xorigin,grid_info.yorigin,grid_info.zorigin);

	// for (core::Size tx=0;tx<grid_info.xnum_points;tx++){
	// jk just print a single x slice (at the middle)
	core::Size tx = ( grid_info.xnum_points + 1 ) / 2;
	atom_pos.x() = 0;

	//  atom_pos.x() += grid_info.xstep;
	atom_pos.y() = grid_info.yorigin;
	for ( core::Size ty=0; ty<grid_info.ynum_points; ty++ ) {
		atom_pos.y() += grid_info.ystep;
		atom_pos.z() = grid_info.zorigin;
		for ( core::Size tz=0; tz<grid_info.znum_points; tz++ ) {
			atom_pos.z() += grid_info.zstep;
			//    if ( ( ( tx % 40 ) == 1 ) && ( ( ty % 40 ) == 1 ) && ( ( tz % 32 ) == 1 ) ) {
			core::Real const curr_grid_val = grid_values[tx][ty][tz];
			if ( curr_grid_val > 0.0001 ) {
				outpdb_stream << setup_pdb_line( atom_pos, atomnum, curr_grid_val ) << '\n';
				++atomnum;
			}
			//    }
		}
	}
	// }

	outpdb_stream.close();
	outpdb_stream.clear();
}


// Dump the grid to a datafile, to be read later
void dump_water_grid_file( std::string const & fname, GridInfo const & grid_info,
	std::vector < std::vector < std::vector <core::Real> > > const & grid_values ) {
	utility::io::ozstream outstream;
	outstream.open( fname, std::ios::out);

	core::Vector atom_pos(grid_info.xorigin,grid_info.yorigin,grid_info.zorigin);

	for ( core::Size tx=0; tx<grid_info.xnum_points; tx++ ) {
		atom_pos.x() += grid_info.xstep;
		atom_pos.y() = grid_info.yorigin;
		for ( core::Size ty=0; ty<grid_info.ynum_points; ty++ ) {
			atom_pos.y() += grid_info.ystep;
			atom_pos.z() = grid_info.zorigin;
			for ( core::Size tz=0; tz<grid_info.znum_points; tz++ ) {
				atom_pos.z() += grid_info.zstep;
				core::Real const curr_grid_val = grid_values[tx][ty][tz];
				if ( curr_grid_val > 0.0001 ) {
					outstream << atom_pos.x() << ' ' << atom_pos.y() << ' ' << atom_pos.z() << ' ' << curr_grid_val << "\n";
				}
			}
		}
	}
	outstream.close();
	outstream.clear();

}


/// General testing code
int
main( int argc, char * argv [] )
{
	try {
		NEW_OPT( atom_type, "The donor / acceptor type to parametrize", "" );
		NEW_OPT( occluding_radius, "The radius of the occluding atom (default is 2.0, eg. for a carbon)", 2. );
		NEW_OPT( water_entropy_scaling_factor, "water_entropy_scaling_factor", 1.0 );

		// note: radii are: 2.0 for C, 1.75 for N, 1.55 for O, 1.9 for S, 2.15 for P, 1.0 for polar H, and 1.2 for non-polar H

		devel::init(argc, argv);

		std::string const atom_type_to_parametrize = option[ atom_type ];
		core::Real const radius_of_occluding_atom = option[ occluding_radius ];
		core::Real const entropy_scaling = option[ water_entropy_scaling_factor ];

		TR << "jk geometric solvation parametrization" << std::endl;

		core::Real const kT = 0.593;
		core::Real const water_O_H_distance = 0.958;

		// Build a peptide with multiple aa's, to use as a template for measuring bond lengths and angles...
		pose::Pose pose;
		core::pose::make_pose_from_sequence( pose, "AACDHKNRSW", *( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD )));
		// for ( Size i = 1; i <= pose.size(); i++ ) {
		//  if ( ! pose.residue(i).is_protein() ) continue;
		//  pose.set_phi( i, -180 ); // note: -150 for typical extended peptide, eg. protocols/topology_broker/SequenceClaimer.cc
		//  pose.set_psi( i, 180); // note: 150 for typical extended peptide, eg. protocols/topology_broker/SequenceClaimer.cc
		//  pose.set_omega( i, 180 );
		// }
		core::Size const ala_resnum = 2;
		// core::Size const cys_resnum = 3;
		core::Size const asp_resnum = 4;
		core::Size const his_resnum = 5;
		core::Size const lys_resnum = 6;
		core::Size const asn_resnum = 7;
		core::Size const arg_resnum = 8;
		core::Size const ser_resnum = 9;
		core::Size const trp_resnum = 10;

		// Setup Etable
		core::scoring::etable::EtableOP etable_ptr( new core::scoring::etable::Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ), core::scoring::etable::EtableOptions() ) );
		// core::Real const methyl_radius = etable_ptr->lj_radius( pose.residue( ala_resnum ).atom_type_index( id::AtomID( pose.residue( ala_resnum ).atom_index("CB") , ala_resnum ).atomno() ) );
		core::Real const water_radius = 1.4;

		// Things that depend on which donor/acceptor we're trying to parametrize
		bool WATER_IS_DONOR;
		HBEvalTuple hbond_eval_tuple;
		core::Real max_possible_LK;
		id::AtomID outer_atom_id, base_atom_id;
		core::Size template_resnum;

		// note: we can't use switch/case for strings, so we'll do if/else...

		// from atom_properties.txt, we find the following cases:
		// 1 -> Nbb donor, eg. backbone NH
		// 2 -> Ntrp donor, eg. Trp
		// 3 -> NH2O donor, eg. Asn/Gln
		// 4 -> Nlys donor, eg. Lys
		// 5 -> Narg donor, eg. Arg
		// 6 -> Npro donor, eg. Pro
		// 7 -> OH donor, eg. Ser/Thr/Tyr
		// 8 -> explicit water donor
		// 9 -> OCbb acceptor, eg. backbone CO
		// 10 -> Nhis acceptor, eg. His
		// 11 -> OH acceptor, eg. Ser/Thr/Tyr
		// 12 -> ONH2 acceptor, eg. Asn/Gln
		// 13 -> OOC acceptor, eg. Asp/Glu
		// 14 -> Oaro acceptor, eg. furan
		// 15 -> explicit water acceptor

		// Note: 6, 8, 14, and 15 are not yet set up

		if ( atom_type_to_parametrize == "" ) {
			TR << "Need to specify atom_type on the command-line" << std::endl;
			exit(1);

		} else if ( atom_type_to_parametrize == "donor_Nbb" ) {
			TR << "Parametrizing for Nbb donor, eg. backbone NH" << std::endl;
			WATER_IS_DONOR = false;
			// Water is an acceptor (hbacc_SP3), set the donor and base types here
			template_resnum = ala_resnum;
			outer_atom_id = id::AtomID( pose.residue( ala_resnum ).atom_index("H") , ala_resnum );
			base_atom_id = id::AtomID( pose.residue( ala_resnum ).atom_index("N") , ala_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( ala_resnum ).atom_type_index( base_atom_id.atomno() ) );

		} else if ( atom_type_to_parametrize == "donor_Ntrp" ) {
			TR << "Parametrizing for Ntrp donor, eg. His/Trp" << std::endl;
			WATER_IS_DONOR = false;
			// Water is an acceptor (hbacc_SP3), set the donor and base types here
			template_resnum = trp_resnum;
			outer_atom_id = id::AtomID( pose.residue( trp_resnum ).atom_index("HE1") , trp_resnum );
			base_atom_id = id::AtomID( pose.residue( trp_resnum ).atom_index("NE1") , trp_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( trp_resnum ).atom_type_index( base_atom_id.atomno() ) );

		} else if ( atom_type_to_parametrize == "donor_NH2O" ) {
			TR << "Parametrizing for NH2O donor, eg. Asn/Gln" << std::endl;
			WATER_IS_DONOR = false;
			// Water is an acceptor (hbacc_SP3), set the donor and base types here
			template_resnum = asn_resnum;
			outer_atom_id = id::AtomID( pose.residue( asn_resnum ).atom_index("1HD2") , asn_resnum );
			base_atom_id = id::AtomID( pose.residue( asn_resnum ).atom_index("ND2") , asn_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( asn_resnum ).atom_type_index( base_atom_id.atomno() ) ) / 2;

		} else if ( atom_type_to_parametrize == "donor_Nlys" ) {
			TR << "Parametrizing for Nlys donor, eg. Lys" << std::endl;
			WATER_IS_DONOR = false;
			// Water is an acceptor (hbacc_SP3), set the donor and base types here
			template_resnum = lys_resnum;
			outer_atom_id = id::AtomID( pose.residue( lys_resnum ).atom_index("1HZ") , lys_resnum );
			base_atom_id = id::AtomID( pose.residue( lys_resnum ).atom_index("NZ") , lys_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( lys_resnum ).atom_type_index( base_atom_id.atomno() ) ) / 3;

		} else if ( atom_type_to_parametrize == "donor_Narg" ) {
			TR << "Parametrizing for Narg donor, eg. Arg" << std::endl;
			WATER_IS_DONOR = false;
			// Water is an acceptor (hbacc_SP3), set the donor and base types here
			template_resnum = arg_resnum;
			outer_atom_id = id::AtomID( pose.residue( arg_resnum ).atom_index("1HH1") , arg_resnum );
			base_atom_id = id::AtomID( pose.residue( arg_resnum ).atom_index("NH1") , arg_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( arg_resnum ).atom_type_index( base_atom_id.atomno() ) ) / 2;

		} else if ( atom_type_to_parametrize == "donor_Npro" ) {
			TR << "Parametrizing for Npro donor, eg. Pro" << std::endl;
			TR << "Skipping this donor, since usage is unclear..." << std::endl;
			exit(1);

		} else if ( atom_type_to_parametrize == "donor_OH" ) {
			TR << "Parametrizing for OH donor, eg. Ser/Thr/Tyr" << std::endl;
			WATER_IS_DONOR = false;
			// Water is an acceptor (hbacc_SP3), set the donor and base types here
			template_resnum = ser_resnum;
			outer_atom_id = id::AtomID( pose.residue( ser_resnum ).atom_index("HG") , ser_resnum );
			base_atom_id = id::AtomID( pose.residue( ser_resnum ).atom_index("OG") , ser_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( ser_resnum ).atom_type_index( base_atom_id.atomno() ) );

		} else if ( atom_type_to_parametrize == "acceptor_OCbb" ) {
			TR << "Parametrizing for OCbb acceptor, eg. backbone CO" << std::endl;
			WATER_IS_DONOR = true;
			// set the acceptor and base types here
			template_resnum = ala_resnum;
			outer_atom_id = id::AtomID( pose.residue( ala_resnum ).atom_index("O") , ala_resnum );
			base_atom_id = id::AtomID( pose.residue( ala_resnum ).atom_index("C") , ala_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( ala_resnum ).atom_type_index( outer_atom_id.atomno() ) );

		} else if ( atom_type_to_parametrize == "acceptor_Nhis" ) {
			TR << "Parametrizing for Nhis acceptor, eg. His" << std::endl;
			WATER_IS_DONOR = true;
			// set the acceptor and base types here
			template_resnum = his_resnum;
			outer_atom_id = id::AtomID( pose.residue( his_resnum ).atom_index("ND1") , his_resnum );
			base_atom_id = id::AtomID( pose.residue( his_resnum ).atom_index("CG") , his_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( his_resnum ).atom_type_index( outer_atom_id.atomno() ) );

		} else if ( atom_type_to_parametrize == "acceptor_OH" ) {
			TR << "Parametrizing for OH acceptor, eg. Ser/Thr/Tyr" << std::endl;
			WATER_IS_DONOR = true;
			// set the acceptor and base types here
			template_resnum = ser_resnum;
			outer_atom_id = id::AtomID( pose.residue( ser_resnum ).atom_index("OG") , ser_resnum );
			base_atom_id = id::AtomID( pose.residue( ser_resnum ).atom_index("CB") , ser_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( ser_resnum ).atom_type_index( outer_atom_id.atomno() ) );

		} else if ( atom_type_to_parametrize == "acceptor_ONH2" ) {
			TR << "Parametrizing for ONH2 acceptor, eg. Asn/Gln" << std::endl;
			WATER_IS_DONOR = true;
			// set the acceptor and base types here
			template_resnum = asn_resnum;
			outer_atom_id = id::AtomID( pose.residue( asn_resnum ).atom_index("OD1") , asn_resnum );
			base_atom_id = id::AtomID( pose.residue( asn_resnum ).atom_index("CG") , asn_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( asn_resnum ).atom_type_index( outer_atom_id.atomno() ) );

		} else if ( atom_type_to_parametrize == "acceptor_OOC" ) {
			TR << "Parametrizing for OOC acceptor, eg. Asp/Glu" << std::endl;
			WATER_IS_DONOR = true;
			// set the acceptor and base types here
			template_resnum = asp_resnum;
			outer_atom_id = id::AtomID( pose.residue( asp_resnum ).atom_index("OD1") , asp_resnum );
			base_atom_id = id::AtomID( pose.residue( asp_resnum ).atom_index("CG") , asp_resnum );
			max_possible_LK = etable_ptr->lk_dgfree( pose.residue( asp_resnum ).atom_type_index( outer_atom_id.atomno() ) );

		} else if ( atom_type_to_parametrize == "acceptor_Oaro" ) {
			TR << "Parametrizing Oaro acceptor, eg. furan" << std::endl;
			TR << "Skipping this acceptor for now..." << std::endl;
			exit(1);

		} else {
			TR << "Unknown atom type for parametrizing: " << atom_type_to_parametrize << std::endl;
			exit(1);
		}

		TR << "Max possible LK solvation energy is " << max_possible_LK << std::endl;

		max_possible_LK = -5.;
		TR << "Changing max_possible_LK to -5." << std::endl;

		if ( WATER_IS_DONOR ) {
			hbond_eval_tuple = HBEval_lookup( hbdon_H2O, get_hb_acc_chem_type( outer_atom_id.atomno(), pose.residue( template_resnum ) ), seq_sep_other);
		} else {
			hbond_eval_tuple = HBEval_lookup( get_hb_don_chem_type( outer_atom_id.atomno(), pose.residue( template_resnum ) ), hbacc_H2O, seq_sep_other);
		}

		// params independent of Hbond donor/acceptor type
		core::Vector outer_atom_position(0,0,0);  // note: this is the proton for a donor, or the oxygen for an acceptor
		core::Vector base_atom_position(0,0,0);
		base_atom_position.z() = - pose.conformation().bond_length( base_atom_id, outer_atom_id );
		core::Real occ_dist_step, occ_angle_step, max_occ_dist, max_occ_angle, min_occ_dist, min_occ_angle;
		min_occ_dist = radius_of_occluding_atom + etable_ptr->lj_radius( pose.residue( template_resnum ).atom_type_index( outer_atom_id.atomno() ) );
		min_occ_angle = 0.;
		max_occ_dist = 7.;
		occ_dist_step = 0.1;
		max_occ_angle = 150.;
		occ_angle_step = 2.;

		// Build a grid for water weights
		// Grid params - grid should be ~4 A deep and 8 A across...?
		GridInfo water_weight_grid_info;

		// JK FOR CHUNKY PICTURE TO ILLUSTRATE GRID
		core::Real water_grid_width = 8.;
		core::Real water_grid_depth = 8.;
		water_weight_grid_info.xnum_points = 17;
		water_weight_grid_info.ynum_points = 17;
		water_weight_grid_info.znum_points = 17;


		// core::Real water_grid_width = 10.;
		// core::Real water_grid_depth = 8.;
		// water_weight_grid_info.xnum_points = 401;
		// water_weight_grid_info.ynum_points = 401;
		// water_weight_grid_info.znum_points = 321;
		water_weight_grid_info.xstep = water_grid_width / ( water_weight_grid_info.xnum_points - 1 );
		water_weight_grid_info.ystep = water_grid_width / ( water_weight_grid_info.ynum_points - 1 );
		water_weight_grid_info.zstep = water_grid_depth / ( water_weight_grid_info.znum_points - 1 );
		// Note: the point at the origin will NOT be considered in calculations - the grid starts AFTER the origin!!
		water_weight_grid_info.xorigin = -water_weight_grid_info.xstep * ( water_weight_grid_info.xnum_points + 1) / 2;
		water_weight_grid_info.yorigin = -water_weight_grid_info.ystep * ( water_weight_grid_info.ynum_points + 1) / 2;
		water_weight_grid_info.zorigin = 0;
		// Setup grid, initialize to zero
		std::vector < std::vector < std::vector <core::Real> > > grid_water_weights;
		grid_water_weights.resize(water_weight_grid_info.xnum_points);
		for ( core::Size tx=0; tx<water_weight_grid_info.xnum_points; tx++ ) {
			grid_water_weights[tx].resize(water_weight_grid_info.ynum_points);
			for ( core::Size ty=0; ty<water_weight_grid_info.ynum_points; ty++ ) {
				grid_water_weights[tx][ty].resize(water_weight_grid_info.znum_points, 0.);
			}
		}


		// Fill the water weights in the grid
		core::Real sum_grid_water_weight = 0.;
		core::Vector water_position(water_weight_grid_info.xorigin,water_weight_grid_info.yorigin,water_weight_grid_info.zorigin);
		hbonds::HBondDatabaseCOP hb_database = HBondDatabase::get_database();

		for ( core::Size tx=0; tx<water_weight_grid_info.xnum_points; tx++ ) {
			water_position.x() += water_weight_grid_info.xstep;
			water_position.y() = water_weight_grid_info.yorigin;
			for ( core::Size ty=0; ty<water_weight_grid_info.ynum_points; ty++ ) {
				water_position.y() += water_weight_grid_info.ystep;
				water_position.z() = water_weight_grid_info.zorigin;
				for ( core::Size tz=0; tz<water_weight_grid_info.znum_points; tz++ ) {
					water_position.z() += water_weight_grid_info.zstep;

					// Compute the current geometry
					core::Real AHdis, xD, xH;
					if ( WATER_IS_DONOR ) {

						// water is the donor, give it perfect geometry
						xD = 0.9999;

						// compute the distance to the accepting water proton
						// subtract the water's OH distance to get the AHdis,
						// since the distance computed was from the acceptor to the water oxygen
						// note: water proton lies on the line between the acceptor and the water oxygen
						AHdis = distance ( outer_atom_position, water_position );
						AHdis -= water_O_H_distance; // water O-H distance

						// find cosine of the base-acceptor-water_proton angle (xH)
						// note: this is the same as the base-acceptor-water_oxygen angle
						xH = dot( (outer_atom_position - base_atom_position).normalize(),  (water_position - outer_atom_position).normalize() );

					} else {

						// water is the acceptor, give it perfect geometry
						xH = 1./3.;  // perfect geometry is cos( 180 - 109.5 degrees), which is 1/3

						// compute the distance to the accepting water
						AHdis = distance ( outer_atom_position, water_position );

						// find the cosine of the base-proton-water angle (xD)
						xD = dot( (outer_atom_position - base_atom_position).normalize(),  (water_position - outer_atom_position).normalize() );

					}

					if ( xH < MIN_xH ) continue;
					if ( xH > MAX_xH ) continue;
					if ( xD < MIN_xD ) continue;
					if ( xD > MAX_xD ) continue;
					if ( AHdis < MIN_R ) continue;
					if ( AHdis > MAX_R ) continue;

					// Get the Hbond energy
					core::Real curr_water_hbond;
					core::scoring::hbonds::hbond_compute_energy(hb_database, hbond_eval_type, AHdis, xD, xH, curr_water_hbond);

					// Save the Hbond energy
					curr_water_hbond *= entropy_scaling;
					if ( curr_water_hbond < 0 ) {
						core::Real curr_water_weight = exp( - curr_water_hbond / kT );
						grid_water_weights[tx][ty][tz] = curr_water_weight;
						sum_grid_water_weight += curr_water_weight;
					}
				}
			}
		}

		// Dump the grid to a PDB file, with the water weight as the B-factor
		dump_grid_pdb("water_weight_grid.pdb", water_weight_grid_info, grid_water_weights);

		// JK FOR CHUNKY PICTURE TO ILLUSTRATE GRID
		exit(0);

		// Dump the grid to a datafile, for use later
		// dump_water_grid_file("grid_"+atom_type_to_parametrize+".txt", water_weight_grid_info, grid_water_weights);

		core::Real Emax_weight = exp( max_possible_LK / kT );
		core::Real Ebulk_weight = ( sum_grid_water_weight * Emax_weight ) / ( 1. - Emax_weight);

		// core::Real Ebulk = -kT * log( Ebulk_weight );
		// TR << "jk for the sake of curiosity, Ebulk is " << Ebulk << std::endl;

		// This grid constant is the denominator in computing solvation energies,
		// it depends on the grid dimensions, and sets the max possible solvation energy (in this case to match LK)
		core::Real grid_constant = sum_grid_water_weight + Ebulk_weight;

		TR << "jk finished computing water grid" << std::endl;

		utility::io::ozstream out_energy_pdb_stream;
		out_energy_pdb_stream.open( "energy_of_occlusion_grid.pdb", std::ios::out);
		core::Size grid_atomnum = 1;

		// Print the axes, for display and fitting with octave
		utility::io::ozstream out_dist_octave;
		out_dist_octave.open("octave_dist.out", std::ios::out);
		core::Size occ_dist_nsteps(0);
		for ( core::Real occ_dist = min_occ_dist; occ_dist <= max_occ_dist ; occ_dist += occ_dist_step ) {
			out_dist_octave << std::setprecision(3) << occ_dist << '\n';
			++occ_dist_nsteps;
		}
		out_dist_octave.close();
		out_dist_octave.clear();

		utility::io::ozstream out_cos_angle_octave;
		out_cos_angle_octave.open("octave_cos_angle.out", std::ios::out);
		core::Size occ_angle_nsteps(0);
		for ( core::Real occ_angle = min_occ_angle; occ_angle <= max_occ_angle ; occ_angle += occ_angle_step ) {
			out_cos_angle_octave << std::setprecision(3) << std::cos( numeric::conversions::radians(occ_angle) ) << '\n';
			++occ_angle_nsteps;
		}
		out_cos_angle_octave.close();
		out_cos_angle_octave.clear();

		// Store occlusion energies here, for fitting
		std::vector < std::vector <core::Real> > geo_sol_desired;
		geo_sol_desired.resize(occ_dist_nsteps);
		for ( core::Size td=0; td<occ_dist_nsteps; td++ ) {
			geo_sol_desired[td].resize(occ_angle_nsteps);
		}

		// Compute occlusion energies, writing to file as we go
		// Anything solvent at dist_cut angstroms away is "occluded" by an atom of the specified radius here
		core::Real const sq_dist_cut = ( radius_of_occluding_atom + water_radius ) * ( radius_of_occluding_atom + water_radius );
		utility::io::ozstream out_energy_list;
		out_energy_list.open("solvation_energies.out", std::ios::out);

		utility::io::ozstream out_energy_list_octave;
		out_energy_list_octave.open("octave_solvation_energies.out", std::ios::out);

		core::Vector occ_atom_position(0,0,0);

		core::Size occ_dist_curr_step_num = 0;
		for ( core::Real occ_dist = min_occ_dist; occ_dist <= max_occ_dist ; occ_dist += occ_dist_step ) {
			core::Size occ_angle_curr_step_num = 0;
			for ( core::Real occ_angle = min_occ_angle; occ_angle <= max_occ_angle ; occ_angle += occ_angle_step ) {

				occ_atom_position.x() = 0.;
				occ_atom_position.y() = occ_dist * std::sin( numeric::conversions::radians(occ_angle) );
				occ_atom_position.z() = occ_dist * std::cos( numeric::conversions::radians(occ_angle) );

				// Loop over all water positions, sum contributions from those occluded
				core::Real sum_occluded_water_weights = 0.;
				core::Vector water_position(water_weight_grid_info.xorigin,water_weight_grid_info.yorigin,water_weight_grid_info.zorigin);
				for ( core::Size wx=0; wx<water_weight_grid_info.xnum_points; wx++ ) {
					water_position.x() += water_weight_grid_info.xstep;
					core::Real sq_xdist = ( water_position.x() - occ_atom_position.x() ) * ( water_position.x() - occ_atom_position.x() );
					if ( sq_xdist > sq_dist_cut ) continue;
					water_position.y() = water_weight_grid_info.yorigin;
					for ( core::Size wy=0; wy<water_weight_grid_info.ynum_points; wy++ ) {
						water_position.y() += water_weight_grid_info.ystep;
						core::Real sq_ydist = ( water_position.y() - occ_atom_position.y() ) * ( water_position.y() - occ_atom_position.y() );
						if ( sq_ydist > sq_dist_cut ) continue;
						water_position.z() = water_weight_grid_info.zorigin;
						for ( core::Size wz=0; wz<water_weight_grid_info.znum_points; wz++ ) {
							water_position.z() += water_weight_grid_info.zstep;
							core::Real sq_zdist = ( water_position.z() - occ_atom_position.z() ) * ( water_position.z() - occ_atom_position.z() );
							if ( sq_zdist > sq_dist_cut ) continue;
							core::Real sq_curr_dist = sq_xdist + sq_ydist + sq_zdist;
							if ( sq_curr_dist < sq_dist_cut ) sum_occluded_water_weights += grid_water_weights[wx][wy][wz];
						}
					}
				}

				// compute solvation energy for occluding atom at this point, store data both as line of pdb file and separate data file
				core::Real geometric_solvation_energy = - kT * log( 1 - ( sum_occluded_water_weights / grid_constant ) );
				geo_sol_desired[occ_dist_curr_step_num][occ_angle_curr_step_num] = geometric_solvation_energy;
				out_energy_list_octave << std::setprecision(6) << geometric_solvation_energy << ' ';
				out_energy_list << std::setprecision(3) << occ_dist << ' ' << std::setprecision(3) << occ_angle << ' '
					<< std::setprecision(6) << geometric_solvation_energy << '\n';

				out_energy_pdb_stream << setup_pdb_line( occ_atom_position, grid_atomnum, geometric_solvation_energy ) << '\n';
				++grid_atomnum;
				++occ_angle_curr_step_num;

			}

			out_energy_list_octave << '\n';
			++occ_dist_curr_step_num;

		}
		out_energy_list.close();
		out_energy_list.clear();
		out_energy_list_octave.close();
		out_energy_list_octave.clear();
		out_energy_pdb_stream.close();
		out_energy_pdb_stream.clear();

		TR << "jk finished computing desired geometric solvation potential" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

