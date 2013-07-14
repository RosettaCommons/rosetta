// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_Mg_KnowledgeBasedPotential.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Jump.hh>

#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/viewer/viewers.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef
// option key includes

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <numeric/xyz.functions.hh>

using namespace core;
using namespace ObjexxFCL::fmt;
using ObjexxFCL::FArray3D;
using namespace basic::options::OptionKeys;

using utility::vector1;
using io::pdb::dump_pdb;


typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( Real, xyz_step )
//OPT_KEY( Real, elec_min_dis )
//OPT_KEY( Real, elec_weight )
OPT_KEY( Real, score_cut )
OPT_KEY( Real, score_cut_PDB )
OPT_KEY( Boolean, brute_force )

////////////////////////////////////////////////////
void	strip_out_magnesiums( pose::Pose & pose ){

	utility::vector1< Size > slice_res;

	for (Size n = 1; n <= pose.total_residue(); n++ ){
		if ( pose.residue(n).name3() == " MG" ) continue;
		slice_res.push_back( n );
	}

	protocols::swa::pdbslice( pose, slice_res );

}

///////////////////////////////////////////////////////////////////////////////
void
figure_out_box_bounds( pose::Pose const & pose,
											 Real & xmin, Real & xmax,
											 Real & ymin, Real & ymax,
											 Real & zmin, Real & zmax ){

	// determine bounds of scan.
	bool init( false );
	for ( Size i = 1; i <= pose.total_residue(); i++ ){

		if ( pose.residue( i ).is_virtual_residue() ) continue;

		if ( pose.residue( i ).name3() == " MG" ) continue;

		for ( Size j = 1; j <= pose.residue( i ).natoms(); j++ ){
			Vector pos = pose.residue( i ).xyz( j );

			if ( !init ) {
				xmin = pos.x();
				xmax = pos.x();
				ymin = pos.y();
				ymax = pos.y();
				zmin = pos.z();
				zmax = pos.z();
				init = true;
			}

			if ( pos.x() < xmin ) xmin = pos.x();
			if ( pos.x() > xmax ) xmax = pos.x();
			if ( pos.y() < ymin ) ymin = pos.y();
			if ( pos.y() > ymax ) ymax = pos.y();
			if ( pos.z() < zmin ) zmin = pos.z();
			if ( pos.z() > zmax ) zmax = pos.z();

		}
	}

	// a little padding to be safe, and round to integer.
	xmin = int( xmin ) - 1.0;
	ymin = int( ymin ) - 1.0;
	zmin = int( zmin ) - 1.0;

	xmax = int( xmax ) + 1.0;
	ymax = int( ymax ) + 1.0;
	zmax = int( zmax ) + 1.0;

}

///////////////////////////////////////////////////////////////////////////////
void
create_grid( Real const xmin,
						 Real const xmax,
						 Real const ymin,
						 Real const ymax,
						 Real const zmin,
						 Real const zmax,
						 Real const xyz_increment,
						 FArray3D< Real > & energy_grid ){

	Size xgridsize = int( ( xmax - xmin ) / xyz_increment + 0.5 );
	Size ygridsize = int( ( ymax - ymin ) / xyz_increment + 0.5 );
	Size zgridsize = int( ( zmax - zmin ) / xyz_increment + 0.5 );
	energy_grid.dimension( xgridsize, ygridsize, zgridsize );
	energy_grid = 0.0;
}

///////////////////////////////////////////////////////////////////////////////
void
define_bins( Real const x,
						 Real const subgrid_radius,
						 Real const xmin,
						 Size const xgridsize,
						 Real const xyz_increment,
						 Size & xbinmin,
						 Size & xbinmax ){

	xbinmin = std::max( 1,              int( ( x - subgrid_radius - xmin )/xyz_increment ) + 1 ); // The +1 is because we are indexing by 1.
	xbinmax = std::min( int(xgridsize), int( ( x + subgrid_radius - xmin )/xyz_increment ) + 1 );
}

///////////////////////////////////////////////////////////////////////////////
Real
get_position( Size const xbin, Real const xmin, Real const xyz_increment ){
  return (  xmin + ( xbin-0.5 )*xyz_increment ); // The -0.5 is to get to the center of the bin.
}


///////////////////////////////////////////////////////////////////////////////
Real
get_score( pose::Pose & pose,
					 Vector const mg_position,
					 scoring::ScoreFunctionOP scorefxn ){

	id::AtomID mg_id( 1, pose.total_residue() );
	pose.set_xyz( mg_id, mg_position );

	Real const score = (*scorefxn)( pose );
	return score;
}

///////////////////////////////////////////////////////////////////////////////
void
scan_magnesium( pose::Pose & pose,
								utility::vector1< Vector > & mg_positions,
								utility::vector1< Real > & mg_energies,
								core::scoring::ScoreFunctionOP scorefxn ){

	using namespace core::scoring;
	using namespace basic::options;

	rna::RNA_Mg_KnowledgeBasedPotential rna_mg_low_res_potential;
	rna_mg_low_res_potential.setup_info_for_mg_calculation( pose );

	rna::RNA_ScoringInfo const & rna_scoring_info( rna::rna_scoring_info_from_pose( pose ) );
	utility::vector1< utility::vector1< Size > > const
			atom_numbers_for_mg_calculation( rna_scoring_info.atom_numbers_for_mg_calculation() );

	Real xmin,xmax,ymin,ymax,zmin,zmax; // initialized below
	figure_out_box_bounds( pose, xmin, xmax, ymin, ymax, zmin,zmax );

	// there's probably a nice object somewhere in Rosetta (e.g., map) for doing this.
	// Anyway, FArray will work for now -- ask Fang later if he has a useful class.
	FArray3D< Real > energy_grid;
	Real const xyz_increment = option[ xyz_step ]();
	create_grid( xmin, xmax, ymin, ymax, zmin, zmax, xyz_increment, energy_grid );

	Size xgridsize = energy_grid.size1();
	Size ygridsize = energy_grid.size2();
	Size zgridsize = energy_grid.size3();
	FArray3D< bool > energy_assigned( xgridsize, ygridsize, zgridsize );

	// vdw calculation takes a while, actually -- do a faster calculation...
	ScoreFunctionOP scorefxn_fast = new ScoreFunction;
	scorefxn_fast->set_weight( rna_mg,  1.0 );
	scorefxn_fast->set_weight( rna_mg_indirect,  1.0 );

	Real init_score_fast = get_score( pose,  Vector(xmax + (xmax-xmin), ymax, zmax), scorefxn_fast );
	Real init_score      = get_score( pose,  Vector(xmax + (xmax-xmin), ymax, zmax), scorefxn );

	Real const score_cutoff_fast( -6.0 );
	Real const score_cutoff( option[ score_cut ]() );

	//	std::cout << "BOUNDS: " << xmin << ' ' << xgridsize << ' ' << ymin << ' ' << ygridsize << ' ' << zmin << ' ' << zgridsize << std::endl;

	pose::Pose pose_fast = pose; // 'scratch' pose for fast checks on mg/ligand binding.

	utility::vector1< Size > scan_res;
	if ( option[ in::file::input_res ].user() ) {
		scan_res = option[ in::file::input_res ]();
	} else {
		for ( Size n = 1; n <= pose.total_residue(); n++ ) scan_res.push_back( n );
	}

	for ( Size q = 1; q <= scan_res.size(); q++ ){
		Size const n = scan_res[ q ];

		if ( n > pose.total_residue() ) continue;

		utility::vector1< Size > const & atom_numbers1   ( atom_numbers_for_mg_calculation[ n ]  );
		core::conformation::Residue const & rsd1 = pose.residue( n );

		for ( Size m = 1; m <= atom_numbers1.size(); ++m ) {

			Size const i = atom_numbers1[ m ];
			if (rsd1.is_virtual(i) ) continue;
			if (!rsd1.heavyatom_is_an_acceptor(i) ) continue;

			Vector const & i_xyz( rsd1.xyz(i) );

			Real const subgrid_radius = 3.0;
			Real const min_radius = 1.5;

			Size xbinmin, xbinmax, ybinmin, ybinmax, zbinmin, zbinmax;
			define_bins( i_xyz.x(), subgrid_radius, xmin, xgridsize, xyz_increment, xbinmin, xbinmax );
			define_bins( i_xyz.y(), subgrid_radius, ymin, ygridsize, xyz_increment, ybinmin, ybinmax );
			define_bins( i_xyz.z(), subgrid_radius, zmin, zgridsize, xyz_increment, zbinmin, zbinmax );

			//			std::cout << "SUB-BOUNDS: " << xbinmin << ' ' << xbinmax << ' ' << ybinmin << ' ' << ybinmax << ' ' << zbinmin << ' ' << zbinmax << std::endl;

			for ( Size xbin = xbinmin; xbin <= xbinmax; xbin++ ){
				for ( Size ybin = ybinmin; ybin <= ybinmax; ybin++ ){
					for ( Size zbin = zbinmin; zbin <= zbinmax; zbin++ ){

						if ( energy_assigned( xbin,ybin,zbin ) ) continue;

						Real const x = get_position( xbin, xmin, xyz_increment );
						Real const y = get_position( ybin, ymin, xyz_increment );
						Real const z = get_position( zbin, zmin, xyz_increment );

						Vector const mg_position = Vector(x,y,z);
						Real const r = ( mg_position - i_xyz ).length() ;
						if ( r > subgrid_radius ) continue;
						if ( r < min_radius ) continue;

						energy_assigned( xbin, ybin, zbin ) = true; // prevents recalculation

						// this is a quick filter. another way to accelerate things would be to do a check for
						// steric problems *first* and immediately discard anything with a bad clash.
						Real const score_fast = get_score( pose_fast, mg_position, scorefxn_fast ) - init_score_fast;
						if ( score_fast > score_cutoff_fast ) continue;

						//						std::cout << x << " " << y << " " << z << ":  " << score_fast << std::endl;

						// Now really score...
						Real const score = get_score( pose, mg_position, scorefxn ) - init_score;

						energy_grid( xbin, ybin, zbin ) = score; // this isn't actually used?

						//scorefxn->show( std::cout, pose );

						if ( score <= score_cutoff ){
							std::cout << x << " " << y << " " << z << ":  " << score << std::endl;
							mg_positions.push_back( mg_position );
							mg_energies.push_back( score );
						}


					}
				}
			}
		}

	}

}

///////////////////////////////////////////////////////////////////////////////
void
add_single_magnesium( pose::Pose & pose )
{
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();
	core::chemical::ResidueTypeCOPs const & rsd_type_list ( residue_set.name3_map ( " MG" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	pose.append_residue_by_jump ( *new_res , pose.total_residue() );
}


///////////////////////////////////////////////////////////////////////////////
// There's a much smarter/faster way to do this, but I want this to be a simple, transparent brute-force search.
void
scan_magnesium_SLOW( pose::Pose & pose,
										 utility::vector1< Vector > & mg_positions,
										 utility::vector1< Real > & mg_energies,
										 core::scoring::ScoreFunctionOP scorefxn ){

	using namespace core::scoring;
	using namespace core::id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Real xmin,xmax,ymin,ymax,zmin,zmax; // initialized below
	figure_out_box_bounds( pose, xmin, xmax, ymin, ymax, zmin,zmax );

	std::cout << "Carrying out scan in following box: " << std::endl;
	std::cout << "  X: " << xmin << " to " << xmax << std::endl;
	std::cout << "  Y: " << ymin << " to " << ymax << std::endl;
	std::cout << "  Z: " << zmin << " to " << zmax << std::endl;

	Distance const xyz_increment = option[ xyz_step ]();

	ScoreFunctionOP scorefxn_fast = new ScoreFunction;
	scorefxn_fast->set_weight( rna_mg,  1.0 );
	//	scorefxn->set_weight( rna_vdw, 1.0 );

	// initially move mg(2+) really far away.
	Real init_score = get_score( pose,  Vector(xmax + (xmax-xmin), ymax, zmax), scorefxn_fast );
	//	runtime_assert( init_score == 0.0 );

	Real const score_cutoff( -8.0 );

	scorefxn->show( std::cout, pose );
	// first scan without VDW (which takes a long time -- not sure why)
	for ( Real x = xmin; x <= xmax; x += xyz_increment ){
		for ( Real y = ymin; y <= ymax; y += xyz_increment ){
			for ( Real z = zmin; z <= zmax; z += xyz_increment ){

				Vector const xyz = Vector(x,y,z);
				Real const score = 	get_score( pose, xyz, scorefxn_fast ) - init_score;

				if ( score < score_cutoff ) {
					std::cout << x << " " << y << " " << z << ":  " << score << std::endl;
					mg_positions.push_back( xyz );
					mg_energies.push_back( score );
				}
			}
		}
	}


	//scorefxn->set_weight( rna_vdw, 1.0 );
	Real const score_cutoff_with_vdw( 999.0 );

	init_score = get_score( pose,  Vector(xmax + (xmax-xmin), ymax, zmax), scorefxn );

	utility::vector1< Vector > mg_positions_filter;
	utility::vector1< Real > mg_energies_filter;


	for ( Size n = 1; n <= mg_positions.size(); n++ ){

		Vector const xyz = mg_positions[ n ];
		Real const score = 	get_score( pose, xyz, scorefxn ) - init_score;

		if ( score < score_cutoff_with_vdw ) {
			std::cout << "filter ==> " << xyz.x() << " " << xyz.y() << " " << xyz.z() << ":  " << score << std::endl;
			mg_positions_filter.push_back( xyz );
			mg_energies_filter.push_back( score );
		}

	}

	mg_positions = mg_positions_filter;
	mg_energies = mg_energies_filter;

}

///////////////////////////////////////////////////////////////////////////////
void
cluster_mg( utility::vector1< Vector > & mg_positions,
						utility::vector1< Real > & mg_energies ){

	// this is such a pain. Need a list to do the sorting...
	std::list< std::pair< Real, Vector > > mg_energy_position_list;
	for ( Size n = 1; n <= mg_positions.size(); n++ ){
		mg_energy_position_list.push_back( std::make_pair( mg_energies[n], mg_positions[n] ) );
	}
	mg_energy_position_list.sort();

	utility::vector1< Vector > mg_positions_cluster;
	utility::vector1< Real > mg_energies_cluster;

	//	Real const CLUSTER_DISTANCE_CUTOFF( 4.0 );
	Real const CLUSTER_DISTANCE_CUTOFF( 2.0 );

	// iterate...
	for ( std::list< std::pair< Real, Vector > >::const_iterator iter = mg_energy_position_list.begin();
				iter != mg_energy_position_list.end(); iter++ ){

		Vector mg_position = iter->second;

		bool too_close( false );
		for ( Size n = 1; n <= mg_positions_cluster.size(); n++ ){
			if ( ( mg_position - mg_positions_cluster[n] ).length() < CLUSTER_DISTANCE_CUTOFF ){
				too_close = true; break;
			}
		}
		if ( !too_close ){
			mg_positions_cluster.push_back( mg_position );
			Real const mg_energy = iter->first;
			mg_energies_cluster.push_back( mg_energy );
			std::cout << mg_positions_cluster.size() << ": " << mg_position.x() << " " << mg_position.y() << " " << mg_position.z() << "  " << mg_energy << std::endl;
		}

	}

	mg_positions = mg_positions_cluster;
	mg_energies = mg_energies_cluster;

}


///////////////////////////////////////////////////////////////////////////////
void
output_mg_to_silent_file( utility::vector1< Vector > & mg_positions,
													utility::vector1< Real > & mg_energies,
													pose::Pose & pose,
													pose::Pose const & reference_pose,
													std::string const & silent_file,
													core::scoring::ScoreFunctionOP scorefxn ){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::io::silent;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );
	core::chemical::ResidueTypeCOPs const & rsd_type_list ( rsd_set->name3_map ( " MG" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );

	pose::Pose pose_mg;
	pose_mg.append_residue_by_bond( *new_res );

	SilentFileData silent_file_data;

	for ( Size n = 1; n <= mg_positions.size(); n++ ){

		Vector mg_position = mg_positions[ n ];
		pose_mg.set_xyz( AtomID( 1, 1 ), mg_position );

		std::string const out_file_tag = "S_" + ObjexxFCL::string_of( n );
		BinaryRNASilentStruct s( pose_mg, out_file_tag );

		// this is that working pose with RNA. move the mg(2+) and rescore it
		get_score( pose, mg_position, scorefxn );
		s.energies_from_pose( pose );

		s.add_energy( "mg_score", mg_energies[ n ] );

		Distance min_dist( 0.0 );
		for ( Size m = 1; m <= reference_pose.total_residue(); m++ ){
			if ( reference_pose.residue(m).name3() != " MG"  ) continue;
			Distance dist = ( reference_pose.residue(m).xyz(1) - mg_position ).length();
			if ( dist < min_dist || min_dist == 0.0 ) min_dist = dist;
		}
		s.add_energy( "rms", min_dist );

		silent_file_data.write_silent_struct( s, silent_file, false /*score_only*/ );

	}
}

/////////////////////////////////////////////////////////////////////////////
void
output_mg_to_PDB( utility::vector1< Vector > & mg_positions,
									utility::vector1< Real > & mg_energies,
									pose::Pose const & pose,
									std::string const pdbfilename,
									Real const score_cut )
{
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::io::silent;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );
	core::chemical::ResidueTypeCOPs const & rsd_type_list ( rsd_set->name3_map ( " MG" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );

	pose::Pose pose_mg = pose;

	for ( Size n = 1; n <= mg_positions.size(); n++ ){

		if ( mg_energies[n] > score_cut ) continue;

		pose_mg.append_residue_by_jump( *(new_res->clone()), 1 );

		Vector mg_position = mg_positions[ n ];
		pose_mg.set_xyz( AtomID( 1, pose_mg.total_residue() ), mg_position );

	}

	pose_mg.dump_pdb( pdbfilename );


}


///////////////////////////////////////////////////////////////////////////////
// gridsearch for potential mg(2+) locations.
///////////////////////////////////////////////////////////////////////////////
void
mg_scan_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::methods;



	std::string const file_path( option[ in::path::pdb ]( 1 ) );
	std::string pdb_file = option[ in::file::s ]()[1];
	std::string silent_file = option[ out::file::silent ]();

	std::string const new_prefix = ".no_cluster.out";
	std::string::size_type pos = silent_file.find( ".out", 0 );
	std::string silent_file_no_cluster = silent_file;
	silent_file_no_cluster.replace( pos, new_prefix.length(), new_prefix );

	Size count( 0 );

	pose::Pose pose;
	//Size total_residues( 0 );  // unused ~Labonte

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	import_pose::pose_from_pdb( pose, *rsd_set,  file_path + '/' + pdb_file );

	if ( count == 0 ) protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );

	/////////////////////////////////////////////
	ScoreFunctionOP scorefxn = new ScoreFunction;
	if ( option[ score::weights ].user() ){
		scorefxn = getScoreFunction();
	} else {
		EnergyMethodOptions energy_options;
		//		energy_options.elec_min_dis(  option[ elec_min_dis ]() );
		scorefxn->set_energy_method_options( energy_options );
		scorefxn->set_weight( rna_mg,  1.0 );
		scorefxn->set_weight( fa_rep,  1.0 );
		//		scorefxn->set_weight( fa_elec,  option[ elec_weight ]() );
	}


	pose::Pose reference_pose = pose; //keep this around, using Mg(2+) as 'references'.

	strip_out_magnesiums( pose );

	add_single_magnesium( pose );

	utility::vector1< Vector > mg_positions;
	utility::vector1< Real > mg_energies;
	if ( option[ brute_force ]() ){
		scan_magnesium_SLOW( pose, mg_positions, mg_energies, scorefxn ); // brute-force & slow.
	} else {
		scan_magnesium( pose, mg_positions, mg_energies, scorefxn ); // brute-force & slow.
	}


	output_mg_to_silent_file( mg_positions, mg_energies, pose, reference_pose, silent_file_no_cluster, scorefxn ); // just save Mg(2+). Can extract later.

	cluster_mg( mg_positions, mg_energies );

	output_mg_to_silent_file( mg_positions, mg_energies, pose, reference_pose, silent_file, scorefxn ); // just save Mg(2+). Can extract later.

	if ( option[ out::file::o ].user() ){
		output_mg_to_PDB( mg_positions, mg_energies, pose, option[out::file::o](), option[ score_cut_PDB ]() );
	}

}


///////////////////////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	mg_scan_test();
  protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


	using namespace basic::options;

	NEW_OPT( xyz_step, "increment in Angstroms for xyz scan", 0.50 );
	NEW_OPT( brute_force, "slow xyz scan", false );
	NEW_OPT( score_cut, "score cut for silent output", -8.0 );
	NEW_OPT( score_cut_PDB, "score cut for PDB output", 0.0 );
	//	NEW_OPT( elec_weight, "elec_weight", 0.5 );
	option.add_relevant( in::file::input_res );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
