// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/kkappel/calculate_distances.cc
/// @brief Calculate distances in RNA/protein structures for score function development
/// @author Kalli Kappel kappel@stanford.edu


//Unit Headers
//Package Headers
//Project Headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/util/SwitchResidueTypeSet.hh>
//Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
//Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/xyzMatrix.hh>
//C++ Headers
#include <iostream>
#include <fstream>

#include <devel/init.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// New options for this application
using namespace basic::options::OptionKeys;

//OPT_KEY( String, mutfile )
//OPT_KEY( Boolean, move_backbone )
//OPT_KEY( Boolean, move_protein_backbone )
//OPT_KEY( Boolean, min_jumps )
//OPT_KEY( String, out_prefix )
//OPT_KEY( String, unbound_protein )
//OPT_KEY( String, unbound_RNA )
//OPT_KEY( Filecore::Vector, unbound_wtRNA_ens )
//OPT_KEY( Boolean, mutate_unbound_RNA )
//OPT_KEY( Integer, prot_offset_num )
//OPT_KEY( Integer, RNA_offset_num )
//OPT_KEY( Integer, protein_start_res )
//OPT_KEY( Integer, protein_end_res )
//OPT_KEY( Real, relax_cutoff_dist )
//OPT_KEY( Boolean, min_only )
//OPT_KEY( Boolean, dump_bound_rna )
//OPT_KEY( String, bound_rna_dump_tag )
//OPT_KEY( Integer, protein_pack_reps )
//OPT_KEY( Integer, Nreps )
OPT_KEY( Boolean, get_vdw )
OPT_KEY( Boolean, get_dist )
OPT_KEY( Boolean, get_restypes )
OPT_KEY( Boolean, use_CEN )
OPT_KEY( Boolean, dump_pdbs )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.kkappel.calculate_distances" );

void get_distances( core::pose::Pose const & pose, std::string struct_name, std::ofstream & out_file ) {

	// Figure out distances
	// Loop through the residues, figure out if they are RNA or protein
	for ( core::Size rsd1 = 1; rsd1 <= pose.total_residue(); ++rsd1 ) {
		// Is it RNA or protein?
		bool is_rsd1_protein( pose.residue( rsd1 ).is_protein() );
		// Loop through the remaining residues
		for ( core::Size rsd2 = (rsd1 + 1); rsd2 <= pose.total_residue(); ++rsd2 ) {
			bool is_rsd2_protein( pose.residue( rsd2 ).is_protein() );
			// Only look at distances between RNA and protein
			// (don't care about protein/protein or RNA/RNA distances)
			if ( (is_rsd1_protein && is_rsd2_protein) || (!is_rsd1_protein && !is_rsd2_protein) ) {
				continue;
			}
			// Calculate all atom/atom distances
			for ( core::Size atom1 = 1; atom1 <= pose.residue( rsd1 ).natoms(); ++atom1 ) {
				// only want to look at heavy atoms
				if ( !pose.residue( rsd1 ).atom_type(atom1).is_heavyatom() ||
						pose.residue( rsd1 ).atom_name(atom1)==" H5'" ||
						pose.residue( rsd1 ).atom_name(atom1)=="H5''" ) continue;
				core::Vector xyz1( pose.residue( rsd1 ).xyz( atom1 ));
				for ( core::Size atom2 = 1; atom2 <= pose.residue( rsd2 ).natoms(); ++atom2 ) {
					// only want to look at heavy atoms
					if ( !pose.residue( rsd2 ).atom_type(atom2).is_heavyatom() ||
							pose.residue( rsd2 ).atom_name(atom2)==" H5'" ||
							pose.residue( rsd2 ).atom_name(atom2)=="H5''" ) continue;
					core::Vector xyz2( pose.residue( rsd2 ).xyz( atom2 ));
					core::Vector delta_xyz = xyz1 - xyz2;
					core::Real distance = ( delta_xyz ).length();
					// Print out resnum1, residue1, atom1, resnum2, residue2, atom2, resnum2, distance
					// For testing
					//if ( distance < 10.0 ) {
					out_file << struct_name << " " << rsd1 << " " << pose.residue(rsd1).name3() << " "
						<< pose.residue(rsd1).atom_name(atom1)
						<< " " << rsd2 << " " << pose.residue(rsd2).name3() << " "
						<< pose.residue(rsd2).atom_name(atom2)
						<< " distance: " << distance << "\n";
					//}
				}
			}
		}
	}
}

void get_atom_vdw( core::pose::Pose const & poseFA,
	//std::string const & struct_name,
	std::ofstream & out_file,
	utility::vector1< std::string > const & RNA_atoms_common,
	utility::vector1< std::string > const & protein_atoms,
	utility::vector1< std::string > const & RNA_atoms_A,
	utility::vector1< std::string > const & RNA_atoms_C,
	utility::vector1< std::string > const & RNA_atoms_G,
	utility::vector1< std::string > const & RNA_atoms_U )
{
	// Convert to centroid for protein part
	core::pose::Pose pose = poseFA;
	if ( basic::options::option[ use_CEN ]() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, true, true );
	}


	char seq_a = 'a';
	char seq_g = 'g';
	char seq_GLY = 'G';
	char seq_c = 'c';
	char seq_u = 'u';

	for ( core::Size rsd1 = 1; rsd1 <= pose.total_residue(); ++rsd1 ) {
		if ( !pose.residue( rsd1 ).is_RNA() ) continue;
		//core::Vector P, C5p, C1p, C3p, N1, O2p;
		//core::Vector N6, N7, N3, N4, C6, O2, O6, N7, N2, O4, C6, O2;

		for ( core::Size rsd2 = 1; rsd2 <= pose.total_residue(); ++rsd2 ) {
			if ( !pose.residue( rsd2 ).is_protein() ) continue;
			// Check how close they are
			if ( basic::options::option[ use_CEN ]() ) {
				if ( (pose.residue( rsd1 ).xyz( " P  " ) - pose.residue( rsd2 ).xyz( "CEN" )).length() > 20.0 ) continue;
			} else {
				if ( (pose.residue( rsd1 ).xyz( " P  " ) - pose.residue( rsd2 ).actcoord()).length() > 20.0 ) continue;
			}
			//std::cout << "Residue 1 " << pose.residue( rsd1 ).name1() << " Residue 2 " << pose.residue( rsd2 ).name1() << std::endl;
			for ( core::Size prot_atom = 1; prot_atom <= protein_atoms.size(); ++prot_atom ) {
				// CB is the same as CEN in glycine residues
				if ( basic::options::option[ use_CEN ]() ) {
					if ( pose.residue( rsd2 ).name1() == seq_GLY && protein_atoms[prot_atom] == " CB " ) continue;
				}
				for ( core::Size RNA_atom = 1; RNA_atom <= RNA_atoms_common.size(); ++RNA_atom ) {
					//std::cout << "RNA atom " << RNA_atoms_common[ RNA_atom ] << std::endl;
					//std::cout << "protein atom " << protein_atoms[ prot_atom ] << std::endl;

					core::Vector rna_to_prot;
					if ( !basic::options::option[ use_CEN ]() && protein_atoms[prot_atom] == "CEN" ) {
						rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_common[ RNA_atom ] ) -
							pose.residue( rsd2 ).actcoord();
					} else {
						rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_common[ RNA_atom ] ) -
							pose.residue( rsd2 ).xyz( protein_atoms[ prot_atom ] );
					}
					core::Real distance = rna_to_prot.length();
					out_file << pose.residue( rsd1 ).name1() << " "
						<< pose.residue( rsd2 ).name1() << " "
						<< RNA_atoms_common[ RNA_atom ] << " "
						<< protein_atoms[ prot_atom ] << " "
						<< distance << "\n";
				}

				if ( pose.residue( rsd1 ).name1() == seq_a ) {
					for ( core::Size RNA_atom = 1; RNA_atom <= RNA_atoms_A.size(); ++RNA_atom ) {
						core::Vector rna_to_prot;
						if ( !basic::options::option[ use_CEN ]() && protein_atoms[prot_atom] == "CEN" ) {
							rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_A[ RNA_atom ] ) -
								pose.residue( rsd2 ).actcoord();
						} else {
							rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_A[ RNA_atom ] ) -
								pose.residue( rsd2 ).xyz( protein_atoms[ prot_atom ] );
						}
						core::Real distance = rna_to_prot.length();
						out_file << pose.residue( rsd1 ).name1() << " "
							<< pose.residue( rsd2 ).name1() << " "
							<< RNA_atoms_A[ RNA_atom ] << " "
							<< protein_atoms[ prot_atom ] << " "
							<< distance << "\n";
					}
				}
				if ( pose.residue( rsd1 ).name1() == seq_g ) {
					for ( core::Size RNA_atom = 1; RNA_atom <= RNA_atoms_G.size(); ++RNA_atom ) {
						core::Vector rna_to_prot;
						if ( !basic::options::option[ use_CEN ]() && protein_atoms[prot_atom] == "CEN" ) {
							rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_G[ RNA_atom ] ) -
								pose.residue( rsd2 ).actcoord();
						} else {
							rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_G[ RNA_atom ] ) -
								pose.residue( rsd2 ).xyz( protein_atoms[ prot_atom ] );
						}
						core::Real distance = rna_to_prot.length();
						out_file << pose.residue( rsd1 ).name1() << " "
							<< pose.residue( rsd2 ).name1() << " "
							<< RNA_atoms_G[ RNA_atom ] << " "
							<< protein_atoms[ prot_atom ] << " "
							<< distance << "\n";
					}
				}
				if ( pose.residue( rsd1 ).name1() == seq_c ) {
					for ( core::Size RNA_atom = 1; RNA_atom <= RNA_atoms_C.size(); ++RNA_atom ) {
						core::Vector rna_to_prot;
						if ( !basic::options::option[ use_CEN ]() && protein_atoms[prot_atom] == "CEN" ) {
							rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_C[ RNA_atom ] ) -
								pose.residue( rsd2 ).actcoord();
						} else {
							rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_C[ RNA_atom ] ) -
								pose.residue( rsd2 ).xyz( protein_atoms[ prot_atom ] );
						}
						core::Real distance = rna_to_prot.length();
						out_file << pose.residue( rsd1 ).name1() << " "
							<< pose.residue( rsd2 ).name1() << " "
							<< RNA_atoms_C[ RNA_atom ] << " "
							<< protein_atoms[ prot_atom ] << " "
							<< distance << "\n";
					}
				}
				if ( pose.residue( rsd1 ).name1() == seq_u ) {
					for ( core::Size RNA_atom = 1; RNA_atom <= RNA_atoms_U.size(); ++RNA_atom ) {
						core::Vector rna_to_prot;
						if ( !basic::options::option[ use_CEN ]() && protein_atoms[prot_atom] == "CEN" ) {
							rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_U[ RNA_atom ] ) -
								pose.residue( rsd2 ).actcoord();
						} else {
							rna_to_prot = pose.residue( rsd1 ).xyz( RNA_atoms_U[ RNA_atom ] ) -
								pose.residue( rsd2 ).xyz( protein_atoms[ prot_atom ] );
						}
						core::Real distance = rna_to_prot.length();
						out_file << pose.residue( rsd1 ).name1() << " "
							<< pose.residue( rsd2 ).name1() << " "
							<< RNA_atoms_U[ RNA_atom ] << " "
							<< protein_atoms[ prot_atom ] << " "
							<< distance << "\n";
					}
				}


			}

		}
	}
}

void get_distance_around_rna_base( core::pose::Pose const & poseFA,
	std::string const & struct_name,
	std::ofstream & out_file ) {
	// Get the distance from the rna base centroid to the centroid of protein residues
	// Want delta x, y, and z (not just distance)
	// Use the RNA base coordinate system to get these values
	// Loop through all residue pairs
	// Copy the pose and switch it to CENTROID for protein residues
	core::pose::Pose pose = poseFA;
	if ( basic::options::option[ use_CEN ]() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, true, true );
	}
	//core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, false, true );
	for ( core::Size rsd1 = 1; rsd1 <= pose.total_residue(); ++rsd1 ) {
		if ( !pose.residue(rsd1).is_RNA() ) continue;
		numeric::xyzMatrix< core::Real > rna_base_coord_sys;
		core::Vector rna_centroid;
		rna_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd1), false /*verbose*/);
		rna_base_coord_sys = core::chemical::rna::get_rna_base_coordinate_system( pose.residue(rsd1), rna_centroid );
		core::Vector x_1 = rna_base_coord_sys.col_x();
		core::Vector y_1 = rna_base_coord_sys.col_y();
		core::Vector z_1 = rna_base_coord_sys.col_z();
		for ( core::Size rsd2 = 1; rsd2 <= pose.total_residue(); ++rsd2 ) {
			if ( !pose.residue(rsd2).is_protein() ) continue;
			core::Vector protein_centroid;
			// For using the C-beta atom
			//if (pose.residue(rsd2).name3() != "GLY"){
			// protein_centroid = pose.residue( rsd2 ).xyz( " CB " );
			//} else {
			// protein_centroid = pose.residue( rsd2 ).xyz( " CA " );
			//}
			// Use the actual centroid
			//protein_centroid = pose.residue( rsd2 ).actcoord();
			if ( basic::options::option[ use_CEN ]() ) {
				protein_centroid = pose.residue( rsd2 ).xyz( "CEN" );
			} else {
				protein_centroid = pose.residue( rsd2 ).actcoord();
			}
			// Figure out the delta_x, delta_y, and delta_z in the RNA coordinate system
			core::Vector dist_1_2 = protein_centroid - rna_centroid;
			core::Real distance( dist_1_2.length() );
			if ( distance > 20.0 ) continue;
			core::Real const dist_x = dot_product(dist_1_2, x_1);
			core::Real const dist_y = dot_product(dist_1_2, y_1);
			core::Real const dist_z = dot_product(dist_1_2, z_1);
			// Write it to the file
			out_file << struct_name << " " << rsd1 << " "
				<< pose.residue(rsd1).name3() << " "
				<< rsd2 << " " << pose.residue(rsd2).name3() << " "
				<< dist_x << " " << dist_y << " " << dist_z << " "
				<< distance << "\n";
		}
	}

}

void get_distances_around_rna_backbone( core::pose::Pose const & poseFA,
	std::string const & struct_name,
	std::ofstream & out_file ) {

	// Get the distance from the backbone phosphate to the centroid of protein residues
	// For now just look at the distance (not x, y, z, b/c unclear what the coordinate
	// system would be)
	// Loop through all residue pairs
	core::pose::Pose pose = poseFA;
	if ( basic::options::option[ use_CEN ]() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, true, true );
	}

	for ( core::Size rsd1 = 1; rsd1 <= pose.total_residue(); ++rsd1 ) {
		if ( !pose.residue(rsd1).is_RNA() ) continue;
		core::Vector rna_phosphate;
		rna_phosphate = pose.residue(rsd1).xyz(" P  ");
		for ( core::Size rsd2 = 1; rsd2 <= pose.total_residue(); ++rsd2 ) {
			if ( !pose.residue(rsd2).is_protein() ) continue;
			core::Vector protein_centroid;
			// For using the C-beta atom
			//if (pose.residue(rsd2).name3() != "GLY"){
			// protein_centroid = pose.residue( rsd2 ).xyz( " CB " );
			//} else {
			// protein_centroid = pose.residue( rsd2 ).xyz( " CA " );
			//}
			// Use the actual centroid
			if ( basic::options::option[ use_CEN ]() ) {
				protein_centroid = pose.residue( rsd2 ).xyz( "CEN" );
			} else {
				protein_centroid = pose.residue( rsd2 ).actcoord();
			}
			// Figure out the delta_x, delta_y, and delta_z in the RNA coordinate system
			core::Vector dist_1_2 = protein_centroid - rna_phosphate;
			core::Real distance( dist_1_2.length() );
			if ( distance > 20.0 ) continue;
			// Write it to the file
			out_file << struct_name << " " << rsd1 << " "
				<< pose.residue(rsd1).name3() << " "
				<< rsd2 << " " << pose.residue(rsd2).name3() << " "
				<< distance << "\n";
		}
	}

}

utility::vector1< std::pair< core::Size, bool > > get_interface_type( core::pose::Pose const & poseFA ) {
	utility::vector1< std::pair< core::Size, bool > > is_interface;
	core::Real const INTERFACE_CUTOFF( 10.0 );

	core::pose::Pose pose = poseFA;
	if ( basic::options::option[ use_CEN ]() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, true, true );
	}

	for ( core::Size rsd = 1; rsd <= pose.total_residue(); ++rsd ) {
		core::Size nbrs = 0;
		// Get the rsd_centroid
		// use the rna base centroid for RNA
		core::Vector rsd_centroid;
		bool const is_rsd_protein( pose.residue(rsd).is_protein() );
		if ( pose.residue(rsd).is_RNA() ) {
			rsd_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd), false /*verbos*/ );
		} else {
			if ( !pose.residue(rsd).is_protein() ) continue; // only look at RNA and protein
			//rsd_centroid = pose.residue(rsd).actcoord();
			if ( basic::options::option[ use_CEN ]() ) {
				rsd_centroid = pose.residue(rsd).xyz("CEN");
			} else {
				rsd_centroid = pose.residue(rsd).actcoord();
			}
		}
		for ( core::Size rsd2 = 1; rsd2 <= pose.total_residue(); ++rsd2 ) {
			if ( rsd == rsd2 ) continue;
			core::Vector rsd2_centroid;
			bool const is_rsd2_protein( pose.residue(rsd2).is_protein() );
			// Only look at RNA/protein and protein/RNA distances here
			if ( (is_rsd2_protein && is_rsd_protein) || (!is_rsd2_protein && !is_rsd_protein) ) continue;

			if ( pose.residue(rsd2).is_RNA() ) {
				rsd2_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd2), false /*verbos*/ );
			} else {
				if ( !pose.residue(rsd2).is_protein() ) continue;
				if ( basic::options::option[ use_CEN ]() ) {
					rsd2_centroid = pose.residue(rsd2).xyz("CEN");
				} else {
					rsd2_centroid = pose.residue(rsd2).actcoord();
				}
			}
			core::Real distance = ( rsd_centroid - rsd2_centroid ).length();
			if ( distance < INTERFACE_CUTOFF ) {
				++nbrs;
			}
		}
		// Write out the number of neighbors for this residue
		if ( nbrs > 0 ) {
			//std::cout << rsd << " INTERFACE " << std::endl;
			is_interface.push_back( std::make_pair( rsd, true ));
		} else {
			//std::cout << rsd << " NON-INTERFACE " << std::endl;
			is_interface.push_back( std::make_pair( rsd, false ));
		}
	}
	return is_interface;
}

utility::vector1< std::pair< core::Size, core::Size > > get_num_nbrs( core::pose::Pose const & poseFA ) {
	utility::vector1< std::pair< core::Size, core::Size > > neighbors;
	core::Real const CUTOFF( 10.0 );

	core::pose::Pose pose = poseFA;
	if ( basic::options::option[ use_CEN ]() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, true, true );
	}

	for ( core::Size rsd = 1; rsd <= pose.total_residue(); ++rsd ) {
		core::Size nbrs = 0;
		// Get the rsd_centroid
		// use the rna base centroid for RNA
		core::Vector rsd_centroid;
		bool const is_rsd_protein( pose.residue(rsd).is_protein() );
		if ( pose.residue(rsd).is_RNA() ) {
			rsd_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd), false /*verbos*/ );
		} else {
			if ( !pose.residue(rsd).is_protein() ) continue;
			if ( basic::options::option[ use_CEN ]() ) {
				rsd_centroid = pose.residue(rsd).xyz("CEN");
			} else {
				rsd_centroid = pose.residue(rsd).actcoord();
			}
		}
		for ( core::Size rsd2 = 1; rsd2 <= pose.total_residue(); ++rsd2 ) {
			if ( rsd == rsd2 ) continue;
			core::Vector rsd2_centroid;
			bool const is_rsd2_protein( pose.residue(rsd2).is_protein() );
			// Only look at RNA/RNA and protein/protein distances here
			if ( (is_rsd2_protein && !is_rsd_protein) || (!is_rsd2_protein && is_rsd_protein) ) continue;

			if ( pose.residue(rsd2).is_RNA() ) {
				rsd2_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd2), false /*verbos*/ );
			} else {
				if ( !pose.residue(rsd2).is_protein() ) continue;
				if ( basic::options::option[ use_CEN ]() ) {
					rsd2_centroid = pose.residue(rsd2).xyz("CEN");
				} else {
					rsd2_centroid = pose.residue(rsd2).actcoord();
				}
			}
			core::Real distance = ( rsd_centroid - rsd2_centroid ).length();
			if ( distance < CUTOFF ) {
				++nbrs;
			}
		}
		// Write out the number of neighbors for this residue
		//std::cout << rsd << " has " << nbrs << " neighbors" << std::endl;
		neighbors.push_back( std::make_pair( rsd, nbrs ) );
	}
	return neighbors;
}

void write_to_restype_file( std::ofstream & restype_file,
	utility::vector1< std::pair< core::Size, bool > > is_interface ,
	utility::vector1< std::pair< core::Size, core::Size > > num_nbrs,
	std::string const & pdb_name ) {
	for ( core::Size i = 1; i<= is_interface.size(); ++i ) {
		restype_file << pdb_name;
		restype_file << " ";
		restype_file << i;
		if ( is_interface[i].second ) {
			restype_file << " INTERFACE ";
		} else {
			restype_file << " NON-INTERFACE ";
		}
		restype_file << num_nbrs[i].second;
		restype_file << "\n";
	}

}

void test_actcoord( core::pose::Pose const & pose ) {
	// Just do this for 3 residues
	core::Size i = 0;
	for ( core::Size rsd = 1; rsd <= pose.total_residue(); ++rsd ) {
		if ( i > 3 ) break;
		if ( !pose.residue( rsd ).is_protein() ) continue;
		std::cout << pose.residue(rsd).name3() << std::endl;
		core::Vector centroid( 0.0, 0.0, 0.0 );
		core::Size num_atoms = 0;
		for ( core::Size atom = 1; atom <= pose.residue(rsd).natoms(); ++atom ) {
			if ( !pose.residue( rsd ).atom_type( atom ).is_heavyatom() ) continue;
			centroid += pose.residue(rsd).xyz(atom);
			std::cout << "centroid " << atom << std::endl;
			num_atoms += 1;
		}
		for ( core::Size s = 1; s<=pose.residue(rsd).actcoord_atoms().size(); ++s ) {
			core::Size actcoord = pose.residue(rsd).actcoord_atoms()[s];
			std::cout << "actcoord " << pose.residue(rsd).actcoord_atoms()[s] << std::endl;
			std::cout << pose.residue(rsd).atom_name( actcoord ) << std::endl;

		}
		std::cout << "Centroid: " << (centroid/num_atoms).length() << std::endl;
		std::cout << "Actcoord: " << (pose.residue(rsd).actcoord()).length() << std::endl;
		//std::cout << pose.residue(rsd).actcoord_atoms().size() << std::endl;
		++i;
	}
}

void calculate_distances( ) {

	using namespace basic::options;

	// For the atom VDW stuff
	utility::vector1< std::string > RNA_atoms_common;
	utility::vector1< std::string > protein_atoms;
	utility::vector1< std::string > RNA_atoms_A;
	utility::vector1< std::string > RNA_atoms_C;
	utility::vector1< std::string > RNA_atoms_G;
	utility::vector1< std::string > RNA_atoms_U;
	RNA_atoms_common.push_back( " P  ");
	RNA_atoms_common.push_back( " C5'");
	RNA_atoms_common.push_back( " C1'");
	RNA_atoms_common.push_back( " C3'");
	RNA_atoms_common.push_back( " N1 ");
	RNA_atoms_A.push_back(" N6 ");
	RNA_atoms_A.push_back(" N7 ");
	RNA_atoms_A.push_back(" N3 ");
	RNA_atoms_A.push_back(" O2'");
	RNA_atoms_C.push_back(" N4 ");
	RNA_atoms_C.push_back(" C6 ");
	RNA_atoms_C.push_back(" O2 ");
	RNA_atoms_C.push_back(" C2'");
	RNA_atoms_G.push_back(" O6 ");
	RNA_atoms_G.push_back(" N7 ");
	RNA_atoms_G.push_back(" N2 ");
	RNA_atoms_G.push_back(" O2'");
	RNA_atoms_U.push_back(" O4 ");
	RNA_atoms_U.push_back(" C6 ");
	RNA_atoms_U.push_back(" O2 ");
	RNA_atoms_U.push_back(" C2'");
	protein_atoms.push_back("CEN");
	protein_atoms.push_back(" CA ");
	protein_atoms.push_back(" CB ");
	protein_atoms.push_back(" C  ");
	protein_atoms.push_back(" O  ");
	protein_atoms.push_back(" N  ");

	core::pose::Pose pose;

	std::ofstream distance_file;
	std::ofstream restype_file;
	std::ofstream distance_from_base_file;
	std::ofstream distance_from_bb_file;

	if ( basic::options::option[ get_restypes ]() ) {
		if ( basic::options::option[ use_CEN ]() ) {
			restype_file.open("RNP_restypes.txt");
		} else {
			restype_file.open("RNP_restypes_actualCEN.txt");
		}
	}

	if ( basic::options::option[ get_dist ]() ) {
		if ( basic::options::option[ use_CEN ]() ) {
			distance_file.open("RNP_distances.txt");
			distance_from_base_file.open("RNP_distances_from_base.txt");
			distance_from_bb_file.open("RNP_distances_from_backbone.txt");
		} else {
			distance_file.open("RNP_distances_actualCEN.txt");
			distance_from_base_file.open("RNP_distances_from_base_actualCEN.txt");
			distance_from_bb_file.open("RNP_distances_from_backbone_actualCEN.txt");
		}
	}

	std::ofstream atom_vdw_file;
	if ( basic::options::option[ get_vdw ]() ) {
		if ( basic::options::option[ use_CEN ]() ) {
			atom_vdw_file.open("RNP_atom_vdw_distances.txt");
		} else {
			atom_vdw_file.open("RNP_atom_vdw_distances_actualCEN.txt");
		}
	}
	//core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	//pose = *(protocols::stepwise::setup::initialize_pose_and_other_poses_from_command_line( rsd_set ));
	// Can make handling PDB vs silent files nicer later, just want it to work for now
	// Read in the PDB files
	using namespace core::import_pose::pose_stream;
	if ( option[ in::file::s ].user() ) {
		utility::vector1< std::string > input_structs = basic::options::option[ basic::options::OptionKeys::in::file::s ]();
		for ( core::Size s = 1; s <= input_structs.size(); ++s ) {
			std::cout << "Working on " << input_structs[ s ] << " " << s
				<< " out of " << input_structs.size() << std::endl;
			core::import_pose::pose_from_file( pose, input_structs[ s ] );
			//std::cout << pose.sequence() << std::endl;
			//test_actcoord( pose );
			if ( basic::options::option[ get_dist ]() ) {
				//get_distances( pose, input_structs[ s ], distance_file );
				get_distance_around_rna_base( pose, input_structs[ s ], distance_from_base_file );
				get_distances_around_rna_backbone( pose, input_structs[ s ], distance_from_bb_file );
			}
			if ( basic::options::option[ get_restypes ]() ) {
				utility::vector1< std::pair< core::Size, bool > > is_interface; // resnum, interface?
				utility::vector1< std::pair< core::Size, core::Size > > num_nbrs; // resnum, nbrs
				is_interface = get_interface_type( pose );
				num_nbrs = get_num_nbrs( pose );
				write_to_restype_file( restype_file, is_interface, num_nbrs, input_structs[ s ] );
			}
			if ( basic::options::option[ get_vdw ]() ) {
				get_atom_vdw( pose, atom_vdw_file, RNA_atoms_common, protein_atoms, RNA_atoms_A, RNA_atoms_C, RNA_atoms_G, RNA_atoms_U );
				//get_atom_vdw( pose, input_structs[ s ], atom_vdw_file, RNA_atoms_common, protein_atoms, RNA_atoms_A, RNA_atoms_C, RNA_atoms_G, RNA_atoms_U );
			}
			if ( basic::options::option[ dump_pdbs ]() ) {
				std::string tag = "rosetta_after_calc_dist_pdbs/";
				tag += input_structs[ s ];
				pose.dump_pdb( tag );
			}

		}
	} else if ( option[ in::file::silent ].user() ) {
		// Read in silent file

		// setup input stream
		PoseInputStreamOP input;
		if ( option[ in::file::tags ].user() ) {
			input = PoseInputStreamOP( new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
				));
		} else {
			input = PoseInputStreamOP( new SilentFilePoseInputStream(
				option[ in::file::silent ]()
				));
		}

		core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		while ( input->has_another_pose() ) {
			// fill/score start_pose
			input->fill_pose( pose, *rsd_set );
			std::string tag = tag_from_pose( pose );

			std::cout << "Working on " << tag << std::endl;
			if ( option[ get_dist ]() ) {
				//get_distances( pose, input_structs[ s ], distance_file );
				get_distance_around_rna_base( pose, tag, distance_from_base_file );
				get_distances_around_rna_backbone( pose, tag, distance_from_bb_file );
			}
			if ( basic::options::option[ get_restypes ]() ) {
				utility::vector1< std::pair< core::Size, bool > > is_interface; // resnum, interface?
				utility::vector1< std::pair< core::Size, core::Size > > num_nbrs; // resnum, nbrs
				is_interface = get_interface_type( pose );
				num_nbrs = get_num_nbrs( pose );
				write_to_restype_file( restype_file, is_interface, num_nbrs, tag );
			}
			if ( basic::options::option[ get_vdw ]() ) {
				get_atom_vdw( pose, atom_vdw_file, RNA_atoms_common, protein_atoms, RNA_atoms_A, RNA_atoms_C, RNA_atoms_G, RNA_atoms_U );
				//get_atom_vdw( pose, tag, atom_vdw_file, RNA_atoms_common, protein_atoms, RNA_atoms_A, RNA_atoms_C, RNA_atoms_G, RNA_atoms_U );
			}
			if ( basic::options::option[ dump_pdbs ]() ) {
				std::string output_tag = "rosetta_after_calc_dist_pdbs/";
				output_tag += tag;
				output_tag += ".pdb";
				pose.dump_pdb( output_tag );
			}


		}

	}


	// // what's up with my heavy atom calculation?
	// for (core::Size rsd1 = 1; rsd1 <= pose.total_residue(); ++rsd1) {
	//  for (core::Size atom =1; atom<= pose.residue(rsd1).natoms(); ++atom) {
	//   if ( !pose.residue( rsd1).atom_type(atom).is_heavyatom() ){
	//    std::cout << pose.residue( rsd1).atom_name(atom) << std::endl;
	//   }
	//  }
	// }
	// // H5' and H5'' are considered heavy atoms?! according to above!

	// Also want to figure out the environment: interface or non-interface
	// and for proteins: buried or exposed
	// and for RNA: base-paired or not base-paired

	// Interface or non-interface I can figure out from the distances
	// Buried or exposed I should calculate here
	// Base-paired or not base paired I should also calculate here
	// Loop through all the residues, for proteins, count all the protein residues around it
	// print out number of residues within e.g. 5A?
	// For now actually let's just do this for both RNA and protein, look at the centroid of the residue
	// figure out how many residue centroid are within a certain distance of this, write this out
	// might not be natural enough for RNA, but I think it's ok to start
	// core::Real const CUTOFF( 10.0 );
	//
	// for ( core::Size rsd = 1; rsd <= pose.total_residue(); ++rsd ) {
	//  core::Size nbrs = 0;
	//  // Get the rsd_centroid
	//  // use the rna base centroid for RNA
	//  core::Vector rsd_centroid;
	//  bool const is_rsd_protein( pose.residue(rsd).is_protein() );
	//  if (pose.residue(rsd).is_RNA()) {
	//   rsd_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd), false /*verbos*/ );
	//  } else {
	//   rsd_centroid = pose.residue(rsd).actcoord();
	//  }
	//  for ( core::Size rsd2 = 1; rsd2 <= pose.total_residue(); ++rsd2 ) {
	//   if ( rsd == rsd2 ) continue;
	//   core::Vector rsd2_centroid;
	//   bool const is_rsd2_protein( pose.residue(rsd2).is_protein() );
	//   // Only look at RNA/RNA and protein/protein distances here
	//   if ((is_rsd2_protein && !is_rsd_protein) || (!is_rsd2_protein && is_rsd_protein)) continue;
	//
	//   if (pose.residue(rsd2).is_RNA()) {
	//    rsd2_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd2), false /*verbos*/ );
	//   } else {
	//    rsd2_centroid = pose.residue(rsd2).actcoord();
	//   }
	//   core::Real distance = ( rsd_centroid - rsd2_centroid ).length();
	//   if ( distance < CUTOFF ) {
	//    ++nbrs;
	//   }
	//  }
	//  // Write out the number of neighbors for this residue
	//  std::cout << rsd << " has " << nbrs << " neighbors" << std::endl;
	// }
	//
	// // Figure out if it's an interface or non-interface residue
	//
	// core::Real const INTERFACE_CUTOFF( 10.0 );
	//
	// for ( core::Size rsd = 1; rsd <= pose.total_residue(); ++rsd ) {
	//  core::Size nbrs = 0;
	//  // Get the rsd_centroid
	//  // use the rna base centroid for RNA
	//  core::Vector rsd_centroid;
	//  bool const is_rsd_protein( pose.residue(rsd).is_protein() );
	//  if (pose.residue(rsd).is_RNA()) {
	//   rsd_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd), false /*verbos*/ );
	//  } else {
	//   rsd_centroid = pose.residue(rsd).actcoord();
	//  }
	//  for ( core::Size rsd2 = 1; rsd2 <= pose.total_residue(); ++rsd2 ) {
	//   if ( rsd == rsd2 ) continue;
	//   core::Vector rsd2_centroid;
	//   bool const is_rsd2_protein( pose.residue(rsd2).is_protein() );
	//   // Only look at RNA/protein and protein/RNA distances here
	//   if ((is_rsd2_protein && is_rsd_protein) || (!is_rsd2_protein && !is_rsd_protein)) continue;
	//
	//   if (pose.residue(rsd2).is_RNA()) {
	//    rsd2_centroid = core::chemical::rna::get_rna_base_centroid( pose.residue(rsd2), false /*verbos*/ );
	//   } else {
	//    rsd2_centroid = pose.residue(rsd2).actcoord();
	//   }
	//   core::Real distance = ( rsd_centroid - rsd2_centroid ).length();
	//   if ( distance < INTERFACE_CUTOFF ) {
	//    ++nbrs;
	//   }
	//  }
	//  // Write out the number of neighbors for this residue
	//  if ( nbrs > 0 ) {
	//   std::cout << rsd << " INTERFACE " << std::endl;
	//  } else {
	//   std::cout << rsd << " NON-INTERFACE " << std::endl;
	//  }
	// }

	// Write out 2 files, one that lists every residue in every structure and
	// lists the number of neighbors and whether it's interface or non-interface
	// the second file lists every heavy atom in every structure and the distance
	// to every other atom

	if ( basic::options::option[ get_vdw ]() ) {
		atom_vdw_file.close();
	}

	if ( basic::options::option[ get_dist ]() ) {
		distance_file.close();
		distance_from_base_file.close();
		distance_from_bb_file.close();
	}
	if ( basic::options::option[ get_restypes ]() ) {
		restype_file.close();
	}
}

int main( int argc, char ** argv ) {

	try {
		using namespace basic::options;
		//NEW_OPT( mutfile, "File listing mutations", "test_mutfile" );
		//NEW_OPT( move_backbone, "Move the backbone?", false );
		//NEW_OPT( move_protein_backbone, "Move the protein backbone?", false );
		//NEW_OPT( min_jumps, "Minimize the jumps?", false ); // Added 1-12-16: before this jumps were minimized
		//NEW_OPT( out_prefix, "Prefix for the output files", "" );
		//NEW_OPT( unbound_protein, "Structure to use for the unbound protein", "" );
		//NEW_OPT( unbound_RNA, "Structure to use for the unbound RNA", "" );
		//NEW_OPT( unbound_wtRNA_ens, "An ensemble of unbound RNA structures for the wt sequence", "" );
		//NEW_OPT( mutate_unbound_RNA, "Should the unbound RNA get mutated to the new sequence?", false );
		//NEW_OPT( prot_offset_num, "Offset to apply for mutfile numbering, how many protein residues come before RNA to mutate?", 0 );
		//NEW_OPT( RNA_offset_num, "Offset to apply for mutfile numbering for RNA, does RNA numbering start at 1?", 0 );
		//NEW_OPT( protein_start_res, "protein start residue", 0 );
		//NEW_OPT( protein_end_res, "protein end residue", 0 );
		//NEW_OPT( relax_cutoff_dist, "Residues within this cutoff distance of mutations will be relaxed", 20);
		//NEW_OPT( min_only, "Minimize only? (Don't re-pack)", false );
		//NEW_OPT( dump_bound_rna, "Dump a pdb containing only the bound RNA (not minimized alone)", false);
		//NEW_OPT( bound_rna_dump_tag, "Tag for the structure from dump_bound_rna", "");
		//NEW_OPT( protein_pack_reps, "Number of times to repack the unbound protein structure then avg", 1);
		//NEW_OPT( Nreps, "Number of times to calculate complex score (top 20% will be averaged)", 10);
		NEW_OPT( get_vdw, "Do the vdw atom calculation", false );
		NEW_OPT( get_dist, "Calculate distances around the base", false );
		NEW_OPT( get_restypes, "Figure out the residue environments", false );
		NEW_OPT( use_CEN, "Use CEN atom for the protein centroid (otherwise use the actual centroid)", true );
		NEW_OPT( dump_pdbs, "dump a pdb of the pose", false );

		devel::init( argc, argv );
		calculate_distances();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
