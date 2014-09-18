// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/motifs/MotifSearch.cc
/// @brief Motif searching protocol
/// @author sthyme (sthyme@gmail.com)

// Unit Headers
#include <protocols/motifs/MotifSearch.hh>

// Package Headers
#include <protocols/motifs/BuildPosition.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifHit.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/motif_utils.hh>

// Project Headers (protocols)
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/DnaInterfaceFinder.hh>
// AUTO-REMOVED #include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/dna/util.hh>
#include <protocols/simple_moves/MinMover.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/dna/setup.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/string_util.hh>

#include <numeric/xyzVector.hh>

// C++ Headers
#include <iostream>

// Option Key Includes
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace motifs {

static thread_local basic::Tracer ms_tr( "protocols.motifs.MotifSearch", basic::t_info );

MotifSearch::~MotifSearch()
{}

MotifSearch::MotifSearch()
	: motif_library_(),
		target_positions_(),
		build_positionOPs_(0),
		target_conformers_map_(),
		// Option flags/parameters: default to command line options
		ztest_cutoff_1_(	basic::options::option[ basic::options::OptionKeys::motifs::z1 ]() ),
		ztest_cutoff_2_(	basic::options::option[ basic::options::OptionKeys::motifs::z2 ]() ),
		rmsd_cutoff_1_(	basic::options::option[ basic::options::OptionKeys::motifs::r1 ]() ),
		rmsd_cutoff_2_(	basic::options::option[ basic::options::OptionKeys::motifs::r2 ]() ),
		dtest_cutoff_(	basic::options::option[ basic::options::OptionKeys::motifs::dtest ]() ),
		rot_level_(	basic::options::option[ basic::options::OptionKeys::motifs::rotlevel ]() ),
		minimize_(	basic::options::option[ basic::options::OptionKeys::motifs::minimize ]() )
{
	init_options();
}

MotifSearch::MotifSearch( MotifSearch const & src ) :
utility::pointer::ReferenceCount( src )
{
	(*this) = src;
}

MotifSearch const &
MotifSearch::operator = ( MotifSearch const & src )
{
	if( this != &src ) {
		motif_library_ = src.motif_library_ ;
		target_positions_ = src.target_positions_ ;
		build_positionOPs_ = src.build_positionOPs_ ;
		target_conformers_map_ = src.target_conformers_map_ ;
		ztest_cutoff_1_ = src.ztest_cutoff_1_ ;
		ztest_cutoff_2_ = src.ztest_cutoff_2_ ;
		rmsd_cutoff_1_ = src.rmsd_cutoff_1_ ;
		rmsd_cutoff_2_ = src.rmsd_cutoff_2_ ;
		dtest_cutoff_ = src.dtest_cutoff_ ;
		rot_level_ = src.rot_level_ ;
		minimize_ = src.minimize_;
		bpdata_ = src.bpdata_;
		bpdata_filename_ = src.bpdata_filename_;
		output_ = src.output_;
		output_filename_ = src.output_filename_;
		data_ = src.data_;
		data_filename_ = src.data_filename_;
		quick_and_dirty_ = src.quick_and_dirty_;
		dump_motifs_ = src.dump_motifs_;
		clear_bprots_ = src.clear_bprots_;
		rots2add_ = src.rots2add_;
		restrict_to_wt_ = src.restrict_to_wt_;
		rerun_motifsearch_ = src.rerun_motifsearch_;
		bump_check_ = src.bump_check_;
	}
	return *this;
}

void
MotifSearch::run(
	Pose & pose,
	utility::vector1< Size > & input_BPs
)
{
	initialize( pose, input_BPs );
	incorporate_motifs( pose );
}

void
MotifSearch::initialize(
	Pose & pose,
	utility::vector1< Size > & input_BPs
)
{

	// Obtain all necessary user input
	if( motif_library_.empty() ) {
		MotifLibrary motifs( get_MotifLibrary_user() );
		MotifCOPs motifcops = motifs.library();
		motif_library_ = motifcops;
	} // if it's not empty that means that the app must have filled the motif_library_
	core::conformation::ResidueOPs conformers( get_targetconformers_user() );
	target_conformers_map_ = setup_conformer_map( conformers );
	if( input_BPs.empty() ) {
		target_positions_ = get_target_position_map_make_dna_mutations( pose );
	}

		// This section is only relevant if you are using defs to get amino acid build positions (unlikely)
		// Temporary container for build positions until I can associate with targets
		DnaDesignDefOPs build_position_defs;
		build_position_defs = get_motif_build_position_defs_user();
		position_vector_setup( pose );

	// If there are no input target positions wild-type DNA is assumed
	// I am not sure if the interface finder is the best way to get positions,
	// I am worried that it is restrictive? change parameter?
	if( build_position_defs.empty() && target_positions_.empty() ) {
		ms_tr << "No input build or target positions, will be doing a motif search using every protein position in the interface as a build position." << std::endl;
		for( utility::vector1< Size >::const_iterator it( dna_positions_.begin() ), end( dna_positions_.end() );
				it != end; ++it ) {
			std::set< std::string > names;
			names.insert( protocols::dna::dna_full_name3( (pose.residue( (*it) ).name3() ) ) );
			target_positions_[*it] = names;
		}
		if( ! input_BPs.empty() ) {
			for( Size i(1); i <= input_BPs.size(); ++i ) {
				BuildPosition_from_Size( pose, input_BPs[i] );
			}
		}	else {
			identify_motif_BuildPositions( pose );
		}
	} else if( build_position_defs.empty() ) {
		ms_tr << "Identifying build positions based on input target positions (or they're from DnaInterfacePacker)." << std::endl;
		if ( ! input_BPs.empty() ) {
			for( Size i(1); i <= input_BPs.size(); ++i ) {
				BuildPosition_from_Size( pose, input_BPs[i] );
			}
		}	else {
			identify_motif_BuildPositions( pose );
		}
		// If the BuildPosition vector remains empty then you're allowed to use any positions
		// If statement below (in search fxn itself) would be if( contains name3 of current motif || .empty );
	} else if ( target_positions_.empty() ) {
		ms_tr << "No target positions given, will be using closest position to placed motif." << std::endl;
		for ( utility::vector1< Size >::const_iterator it( dna_positions_.begin() ), end( dna_positions_.end() );
				it != end; ++it ) {
			std::set< std::string > names;
			names.insert( protocols::dna::dna_full_name3( (pose.residue( (*it) ).name3() ) ) );
			target_positions_[*it] = names;
		}
		defs2BuildPositions_findts( pose, build_position_defs );
	} else {
		defs2BuildPositions( pose, build_position_defs );
		// When there is both an input for build positions and target positions,
		// the input is forced (each build position uses all input target sites)
	}
	if( bpdata_ ) {
		for( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
				ir != end_ir; ++ir ) {
			load_build_position_data( **ir, bpdata_filename_, pose );
			if( clear_bprots_ ) {
				// I could avoid adding them in the first place, but that loading function is more complicated than just clearing and I might want them later
				(*ir)->clear_rots();
			}
		}
		if( ! motif_library_.empty() ) {
			ms_tr << "WARNING!! Loading BPData, but also loaded a MotifLibrary of all motifs" << std::endl;
		}
	}
}

void
MotifSearch::incorporate_motifs(
	Pose const & pose
)
{
	core::pose::Pose posecopy( pose );
	core::pose::Pose posecopy2( pose );

	// Setup output files for motifs and rotamers
	utility::io::ozstream motif_output_file;
	utility::io::ozstream data_output_file;
	utility::io::ozstream qd_output_file;
	if( output_ ) {
		motif_output_file.open_append( output_filename_ );
	}
	if( data_ ) {
		data_output_file.open_append( data_filename_ );
	}
	if( quick_and_dirty_ ) {
		std::string motif_filename( pose.pdb_info()->name() );
		motif_filename.erase( motif_filename.end()-4, motif_filename.end() );
		std::string qd_output_filename = motif_filename + ".qd_motifs";
		qd_output_file.open_append( qd_output_filename );
	}

	// for every protein backbone position (motif build position)
	for( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
			ir != end_ir; ++ir ) {
		// Map of all of the very best residues for each amino acid type, to make sure I don't add 2,000 Args and only 2 Tyrs
		std::map< std::string, std::map< Real, MotifHitOP > > best_mhits_all;

		// If we have rotamers coming in from files and they weren't cleared in initialization, then the search won't happen on this BuildPosition
		if ( ( ! ((*ir)->best_rotamers()).empty() ) && ( ! rerun_motifsearch_ ) ) continue;
		ms_tr << "WORKING ON PROTEIN POSITION " << (*ir)->seqpos() << std::endl;
		MotifCOPs bp_best_motifs( (*ir)->best_motifs() );
		core::pack::rotamer_set::Rotamers bp_best_rotamers( (*ir)->best_rotamers() );
		// Need to clear the best_rotamers and best_motifs before I collect new ones
		(*ir)->clear_data();

		Size seqpos( (*ir)->seqpos() );
		std::stringstream firstline;
		firstline << "POSITION " << seqpos;
		if( output_ ) {
			motif_output_file << firstline.str() << "\n";
		}
		if( data_ ) {
			data_output_file << firstline.str() << "\n";
		}
		if( quick_and_dirty_ ) {
			qd_output_file << firstline.str() << "\n";
		}

		MotifCOPs motif_library;
		// If the BP has motifs in it at this point then they came from an input file from the cmd line
		if( ! bp_best_motifs.empty() ) {
				motif_library = bp_best_motifs;
		} else {
			motif_library = motif_library_;
		}

		std::map< core::Size, core::pack::rotamer_set::RotamerSetOP > rotamer_sets;
		if( bp_best_rotamers.empty() ) {
			for( Size i(1); i <= core::chemical::num_canonical_aas; ++i ) {
				utility::vector1< bool > aa_info( core::chemical::num_canonical_aas, false );
				aa_info[i] = true;
				core::pack::rotamer_set::RotamerSetOP rotset = build_rotamers_lite( posecopy, seqpos, aa_info, rot_level_, bump_check_ );
				rotamer_sets[i] = rotset;
			}
		} else {
				Size bp_rots( bp_best_rotamers.size() );
				for( Size i(1); i <= core::chemical::num_canonical_aas; ++i ) {
					core::pack::rotamer_set::RotamerSetFactory rsf;
					core::pack::rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( posecopy.residue((*ir)->seqpos()) );
					for( Size r(1); r <= bp_rots; ++r ) {
						if( bp_best_rotamers[r]->name3() == core::chemical::name_from_aa(core::chemical::AA(i)) ) {
							rotset->add_rotamer( *((bp_best_rotamers)[r]) );
						}
					}
					rotamer_sets[i] = rotset;
				}
		}

		// For every motif in the motif_library_
		Real dtest_cutoff_sq = ( dtest_cutoff_ * dtest_cutoff_ );
		for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = motif_library.begin(), end_itr = motif_library.end();
				motifcop_itr != end_itr; ++motifcop_itr ) {
			// WARNING: everything in this code assumes that residue 1 in the motif is the protein position and residue 2 is the dna
			bool passed_quick_and_dirty(false);
			protocols::motifs::MotifCOP motifcop( *motifcop_itr );
			//ms_tr << "WORKING ON MOTIF: " << motifcop->remark() << std::endl;

			// The BuildPosition may at some point have the ability to restrict ahead of time
			// currently I am not sure how to get restrictions from resfiles and so on . . .
			// it may not be an issue because there are ways other than resfiles, such as the defs for protein positions
			// and the main reason I would want this functionality is for homology models and theoretically the input pose
			// should already have all of the position types set as wild-type and be relaxed and minimized
			std::set< std::string > allowedtypes( (*ir)->allowed_types() );

			// Tighter cutoffs for residues with 3 and 4 chi angles
			Real ztest_cutoff_1 = ztest_cutoff_1_;
			Real ztest_cutoff_2 = ztest_cutoff_2_;
			Real rmsd_cutoff_1 = rmsd_cutoff_1_;
			Real rmsd_cutoff_2 = rmsd_cutoff_2_;
			Real dtest_cutoff = dtest_cutoff_sq;
			if( (motifcop->restype_name1() == "ARG") && ( ! quick_and_dirty_ ) ) {
				rmsd_cutoff_1 = (rmsd_cutoff_1 / 2.0 );
				rmsd_cutoff_2 = (rmsd_cutoff_2 / 2.0 );
				dtest_cutoff = (dtest_cutoff / 2.0 );
			}
			if( ( (motifcop->restype_name1() == "MET") || (motifcop->restype_name1() == "LYS") || (motifcop->restype_name1() == "GLU") ||  (motifcop->restype_name1() == "GLN") ) && ( ! quick_and_dirty_ ) ) {
				rmsd_cutoff_1 = (rmsd_cutoff_1 / 1.33 );
				rmsd_cutoff_2 = (rmsd_cutoff_2 / 1.33 );
				dtest_cutoff = (dtest_cutoff / 1.33 );
			}

			// Atoms required for the parallel base test
			// maybe these vectors should be stored as a part of this class?
			// for non-DNA motif searches I won't be using the z-test at all
			utility::vector1< std::string > atoms;
			if( protocols::dna::dna_full_name3( motifcop->restype_name2() ) == "GUA" || protocols::dna::dna_full_name3( motifcop->restype_name2() ) == "ADE" ) {
				atoms.push_back("C5");
				atoms.push_back("C6");
				atoms.push_back("N3");
				atoms.push_back("C2");
				atoms.push_back("N1");
				atoms.push_back("C4");
			} else if( protocols::dna::dna_full_name3( motifcop->restype_name2() ) == "CYT" || protocols::dna::dna_full_name3( motifcop->restype_name2() ) == "THY" ) {
				atoms.push_back("C5");
				atoms.push_back("C4");
				atoms.push_back("N1");
				atoms.push_back("C2");
				atoms.push_back("N3");
				atoms.push_back("C6");
			} else {
				ms_tr << "Residue you are planning to do parallel base test with is not a DNA base!" << std::endl;
			}

			bool automorphism(false);
			bool passed_automorphism(false);
			if(
				motifcop->restype_name1() == "ASP" ||
				motifcop->restype_name1() == "GLU" ||
				motifcop->restype_name1() == "PHE" ||
				motifcop->restype_name1() == "LEU" ||
				motifcop->restype_name1() == "ARG" ||
				motifcop->restype_name1() == "TYR" ||
				motifcop->restype_name1() == "VAL"
			)
			{
				automorphism = true;
			}

			// Related to only using motifs that include an allowed amino acid at the BuildPosition
			// At this point the allowedtypes is always going to be empty
			bool allowed(false);
			if( allowedtypes.empty() ) {
				allowed = true;
			}
			for( std::set< std::string >::const_iterator ir2(allowedtypes.begin() ), end_ir2( allowedtypes.end() );
					ir2 != end_ir2; ++ir2 ) {
				if( (*ir2) == motifcop->restype_name1() ) {
					allowed = true;
				}
			}
			if( ! allowed ) continue;

			// Create a standard dna base of type in motif
			// For ligand motifs we will need to create a residue of the type we are trying to target
			// Since there won't be a second residue in the motif
			// restype_name2() will need to be set to something in the motif?
			std::string basetype( motifcop->restype_name2() );
			// If there are conformers
			core::conformation::ResidueOPs DNAResidueOPs( target_conformers_map_[basetype] );
			bool noconformers( false );
			if( DNAResidueOPs.empty() ) {
				noconformers = true;
			}
			core::conformation::ResidueOP baseres2 = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( basetype ) );
			// Get the rotamer set with type matching motif at the motif build position
			core::pack::rotamer_set::RotamerSetOP rotset = rotamer_sets[ core::chemical::aa_from_name(motifcop->restype_name1()) ];

			Real final(100);
			std::pair< core::conformation::ResidueOP, core::conformation::ResidueOP > bestpair;
			bool b_bestpair( false );
			Size rs1( rotset->num_rotamers() );

			// For every rotamer in the rotamer set
			for( Size ir2(1); ir2 <= rs1; ++ir2 ) {
				Real ztest_ir2(0.0);
				Real rmsdtest_ir2(100.0);

				// This atom type is only C1' because that atom is common
				// to all DNA bases and it is basically in the plane of the base
				core::conformation::Atom atm( baseres2->atom("C1'") );
				core::conformation::Atom auto_atm( baseres2->atom("C1'") );
				motifcop->place_atom( *(rotset->nonconst_rotamer(ir2)), *baseres2, atm );
				if( automorphism ) {
					motifcop->place_atom( *(rotset->nonconst_rotamer(ir2)), *baseres2, auto_atm, false );
				}

				// Important in case of the very strange situation where more than one
				// target_position is close enough to pass the tests
				bool tftest(false);
				Size bestpos(0);
				Real test(1000);

				Sizes target_positions( (*ir)->target_positions() );
				for( Sizes::const_iterator bpos( target_positions.begin() ), end( target_positions.end() );
						bpos != end; ++bpos ) {
					std::set< std::string > allowed_types( target_positions_[*bpos] );
					for(	std::set< std::string >::const_iterator type(	allowed_types.begin()), atend( allowed_types.end() );
							type != atend; ++type ) {
						if( basetype == *type ) {
							// Should have a test to ensure that it has the atom types I'll be using?
							Real dtest1( atm.xyz().distance_squared( posecopy.residue( *bpos ).xyz( "C1'" ) ) );
							Real dtest1_auto(100);
							if( automorphism ) {
								dtest1_auto = ( auto_atm.xyz().distance_squared( posecopy.residue( *bpos ).xyz( "C1'" ) ) );
							}

							if( ! automorphism ) {
								if( dtest1 > dtest_cutoff ) continue;
							} else {
								if( dtest1 > dtest_cutoff && dtest1_auto > dtest_cutoff ) continue;
								if( dtest1 < dtest_cutoff && dtest1_auto < dtest_cutoff ) {
									if( dtest1_auto < dtest1 ) {
										passed_automorphism = true;
									} else {
										passed_automorphism = false;
									}
								} else if( dtest1_auto < dtest_cutoff ) {
									passed_automorphism = true;
								} else {
									passed_automorphism = false;
								}
							}
							core::conformation::ResidueOP posebase = new core::conformation::Residue( posecopy.residue( *bpos ) );
							if( passed_automorphism ) {
								motifcop->place_atoms( *(rotset->nonconst_rotamer(ir2)), *posebase, atoms, false );
							} else {
								motifcop->place_atoms( *(rotset->nonconst_rotamer(ir2)), *posebase, atoms );
							}
							Real ztest = protocols::motifs::parallel_base_test( *posebase, posecopy.residue( *bpos ) );

							if( ztest < ztest_cutoff_1 ) continue;
							ztest_ir2 = ztest;
							if( quick_and_dirty_ ) {
								qd_output_file << *motifcop;
								passed_quick_and_dirty = true;
								break;
							}
							if( protocols::dna::dna_full_name3( posecopy.residue(*bpos).name3() ) != basetype ) {
								make_base_pair_mutation( posecopy, *bpos , core::chemical::aa_from_name( basetype ) );
							}
							Real rmsdtest = atom_specific_rms( *posebase, posecopy.residue(*bpos), atoms );

							if( rmsdtest > rmsd_cutoff_1 ) continue;
							rmsdtest_ir2 = rmsdtest;
							if( rmsdtest < test ) {
								test = rmsdtest;
								bestpos = *bpos;
								posecopy2 = posecopy;
								tftest = true;
							} // if( rmsdtest < test )
						}
					} // loop over allowed_types for the allowed DNA types
					if( passed_quick_and_dirty ) break;
				} // loop over potential target DNA positions

				if( passed_quick_and_dirty ) break;
				if( ! tftest ) continue;

				MotifHitOP motifhit = new MotifHit( *motifcop, bestpos, passed_automorphism );

				// If there are no conformers will use the base from the pose, so rmsd can be 0 theoretically
				if( noconformers ) {
					core::conformation::ResidueOP posebase2 = new core::conformation::Residue( posecopy2.residue( bestpos) );
					if( passed_automorphism ) {
						motifcop->place_residue(*(rotset->nonconst_rotamer(ir2)), *posebase2, false );
					} else {
						motifcop->place_residue(*(rotset->nonconst_rotamer(ir2)), *posebase2 );
					}
					// The RMSD is def different if you do the entire base instead of the few atoms as above
					//Real ztest2 = parallel_base_test( *posebase2, posecopy2.residue(bestpos) );
					//Real rmsdtest2 = core::scoring::automorphic_rmsd( *posebase2, posecopy2.residue(bestpos), false );
					//Real finaltest2 = ( ( rmsdtest2 * 100 ) / ( ztest2 * 100 ) );
					Real finaltest = ( ( rmsdtest_ir2 * 100 ) / ( ztest_ir2 * 100 ) );
					//NOTE: Do I want to keep any statistics about percentages of passing certain cutoffs??
					if( (rmsdtest_ir2 < rmsd_cutoff_2) && (ztest_ir2 > ztest_cutoff_2) ) {
						ms_tr << "Passed! RMSD between DNA resi (rosetta #), no conformers " << bestpos << " and motif DNA = " << rmsdtest_ir2 << " and Z-test = " << ztest_ir2 << " and combined score = " << finaltest << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
						//ms_tr << "Passed! RMSD between DNA resi (rosetta #) " << bestpos << " and motif DNA = " << rmsdtest2 << " and Z-test = " << ztest2 << " and combined score = " << finaltest2 << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
						if( data_ ) {
							data_output_file << "Passed! RMSD between DNA resi (rosetta #), no conformers " << bestpos << " and motif DNA = " << rmsdtest_ir2 << " and Z-test = " << ztest_ir2 << " and combined score = " << finaltest << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
						}
						// Fix to make it so that rotamers with different hydrogen positions are added, even though they are not part of the motif and thus have the same finaltest value
						if ( ( motifcop->restype_name1() == "SER" ) || ( motifcop->restype_name1() == "THR" ) || ( motifcop->restype_name1() == "TYR" ) || ( motifcop->restype_name1() == "ILE" ) || ( motifcop->restype_name1() == "CYS" ) ) {
							finaltest = finaltest + ( 0.00001 * ir2 );
						}
						motifhit->final_test( finaltest );
						motifhit->build_rotamer( *(rotset->nonconst_rotamer(ir2) ) );
						motifhit->target_conformer( *posebase2 );
						best_mhits_all[motifcop->restype_name1()][finaltest] = motifhit->clone();
						if ( finaltest < final ) {
							b_bestpair = true;
							bestpair = std::make_pair( (rotset->nonconst_rotamer(ir2))->clone(), posebase2->clone() );
							final = finaltest;
						}
					} // if passed the second round of tests
				} else {
						for( core::conformation::ResidueOPs::const_iterator resop( DNAResidueOPs.begin() ), end( DNAResidueOPs.end() );
								resop != end; ++resop ) {
							if( passed_automorphism ) {
								motifcop->place_residue(*(rotset->nonconst_rotamer(ir2)), **resop, false );
							} else {
								motifcop->place_residue(*(rotset->nonconst_rotamer(ir2)), **resop );
							}
							Real ztest2 = parallel_base_test( **resop, posecopy2.residue(bestpos) );
							if( ztest2 < ztest_cutoff_2 ) continue;
							Real rmsdtest2 = core::scoring::automorphic_rmsd( **resop, posecopy2.residue(bestpos), false );
							if( rmsdtest2 > rmsd_cutoff_2 ) continue;
							Real finaltest = ( ( rmsdtest2 * 100 ) / ( ztest2 * 100 ) );
							Real finaltestc = ( ( rmsdtest_ir2 * 100 ) / ( ztest_ir2 * 100 ) );
							//NOTE: Do I want to keep any statistics about percentages of passing certain cutoffs??
							if( (rmsdtest2 < rmsd_cutoff_2) && (ztest2 > ztest_cutoff_2) ) {
								ms_tr << "Passed! RMSD between DNA resi (rosetta #) " << bestpos << " and motif DNA = " << rmsdtest2 << " and Z-test = " << ztest2 << " and combined score = " << finaltest << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
								if( data_ ) {
									data_output_file << "Passed! RMSD between DNA resi (rosetta #) " << bestpos << " and motif DNA = " << rmsdtest2 << " and Z-test = " << ztest2 << " and combined score = " << finaltest << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
								}
								// Fix to make it so that rotamers with different hydrogen positions are added, even though they are not part of the motif and thus have the same finaltest value
								if ( ( motifcop->restype_name1() == "SER" ) || ( motifcop->restype_name1() == "THR" ) || ( motifcop->restype_name1() == "TYR" ) || ( motifcop->restype_name1() == "ILE" ) || ( motifcop->restype_name1() == "CYS" ) ) {
									finaltest = finaltest + ( 0.00001 * ir2 );
								}
								motifhit->final_test( finaltest );
								motifhit->build_rotamer( *(rotset->nonconst_rotamer(ir2) ) );
								motifhit->target_conformer( **resop );
								best_mhits_all[motifcop->restype_name1()][finaltestc] = motifhit->clone();
								if ( finaltest < final ) {
									b_bestpair = true;
									bestpair = std::make_pair( (rotset->nonconst_rotamer(ir2))->clone(), (*resop)->clone() );
									final = finaltest;
								}
							} // if conformer passed the second round of tests
						} // loop over conformers
					} // if not noconfomers
			}  // loop over first rotamer set

			// Dump the best rotamer and conformer pair for each motif, mainly for debugging purposes
			if( dump_motifs_ ) {
				if( b_bestpair ) {
						core::pose::Pose pose_dump2;
						pose_dump2.append_residue_by_jump( *(bestpair.second), 1);
						pose_dump2.append_residue_by_jump( *(bestpair.first), 1);
						std::stringstream pose2_name_full;
						if( passed_automorphism ) {
							pose2_name_full << "Test_auto_" << motifcop->restype_name2()[0] << "_" << (*ir)->seqpos() << motifcop->restype_name1() << "_" << motifcop->remark() << ".pdb";
						} else {
							pose2_name_full << "Test_" << motifcop->restype_name2()[0] << "_" << (*ir)->seqpos() << motifcop->restype_name1()  << "_" << motifcop->remark() << ".pdb";
						}
						core::io::pdb::dump_pdb( pose_dump2, pose2_name_full.str() );
				}
			}
		}
		if( ! best_mhits_all.empty() ) {
			core::pose::Pose pose_dump( pose );
			for( std::map< std::string, std::map< Real, MotifHitOP > >::const_iterator bh( best_mhits_all.begin() ),
					end( best_mhits_all.end() ); bh != end; ++bh ) {
				Size hits = 0;
				for( std::map< Real, MotifHitOP >::const_iterator bh2( (bh->second).begin() ),
						end2( (bh->second).end() ); bh2 != end2; ++bh2 ) {
					MotifHitOP motifhitop( bh2->second );
					if( ! minimize_ ) {
						(*ir)->keep_rotamer( *(motifhitop->build_rotamer()) );
						(*ir)->keep_motif( *(motifhitop->motifcop()) );
						(*ir)->keep_motifhit( *(motifhitop) );
						++hits;
						if( output_ ) {
							motif_output_file << *(motifhitop->motifcop());
							motif_output_file << "RESIDUE " << (*(motifhitop->build_rotamer()));
						}
						if( hits > rots2add_ ) break;
					} else {
						using namespace core::scoring;
						// Need a copy of the pose to generate the constraints
						// Maybe there's a way to rewrite the constraint making code to avoid this?
						// No, because the minimize itself needs the residues to be a part of the pose?
						// Or maybe instead of copying the pose you could save the residue you are replacing and replace it back?
						//core::pose::Pose pose_dump( pose );
						core::conformation::ResidueOP build_rotamer = new core::conformation::Residue( *(motifhitop->build_rotamer()) );
						pose_dump.replace_residue( (*ir)->seqpos(), *build_rotamer, false );
						// The residue is probably already placed . . .
					//	if( passed_automorphism ) {
					//		motifcop->place_residue( motifhitop->build_rotamer(), motifhitop->target_conformer(), false );
					//	} else {
					//		motifcop->place_residue( motifhitop->build_rotamer(), motifhitop->target_conformer() );
				//		}
						if ( protocols::dna::dna_full_name3( pose_dump.residue(motifhitop->vbpos()).name3() ) != protocols::dna::dna_full_name3( (motifhitop->target_conformer())->name3() ) ) {
							make_base_pair_mutation( pose_dump, motifhitop->vbpos(), core::chemical::aa_from_name( protocols::dna::dna_full_name3( motifhitop->target_conformer()->name3() ) ) );
						}
						if( motifhitop->passed_automorphism() ) {
							(motifhitop->motifcop())->place_residue( pose_dump.residue( motifhitop->vbpos() ), *build_rotamer, false );
						} else {
							(motifhitop->motifcop())->place_residue( pose_dump.residue( motifhitop->vbpos() ), *build_rotamer );
						}
						/*core::pose::Pose pose_dump2( pose );
						pose_dump2.replace_residue( (*ir)->seqpos(), *build_rotamer, false );
						std::stringstream pose2_name_full;
						std::stringstream pose_name_full;
						std::stringstream pose3_name_full;
						pose2_name_full << "AfterReplace_" << bh2->first << ".pdb";
						pose_name_full << "BeforeReplace_" << bh2->first << ".pdb";
						pose3_name_full << "AfterMinBeforeReplace_" << bh2->first << ".pdb";
						core::io::pdb::dump_pdb( pose_dump2, pose2_name_full.str() );
						core::io::pdb::dump_pdb( pose_dump, pose_name_full.str() );*/

						constraints::ConstraintSetOP sc_cst_set( new constraints::ConstraintSet() );
						add_motif_sc_constraints( sc_cst_set, pose_dump, (*ir)->seqpos(), *build_rotamer, motifhitop->motifcop(), false );
						//add_motif_sc_constraints( sc_cst_set, pose_dump2, (*ir)->seqpos(), *build_rotamer, motifhitop->motifcop(), false );
						ScoreFunctionOP score_fxn( get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );
						methods::EnergyMethodOptions options( score_fxn->energy_method_options() );
						score_fxn->set_energy_method_options( options );
						score_fxn->set_weight( coordinate_constraint, 10.0 );
						pose_dump.constraint_set( sc_cst_set );
						//pose_dump2.constraint_set( sc_cst_set );
						//core::Real pre_sc_constraint_check( pose_dump.energies().total_energies()[ coordinate_constraint ] );
						//core::Real pre_sc_constraint_check2( pose_dump2.energies().total_energies()[ coordinate_constraint ] );
						//ms_tr << "Before sidechain refinement constraints score is " << pre_sc_constraint_check << std::endl;
						//ms_tr << "2Before sidechain refinement constraints score is " << pre_sc_constraint_check2 << std::endl;
						/*if( data_ ) {
							data_output_file << "Before sidechain refinement constraints score is " << pre_sc_constraint_check << std::endl;
						}*/
						core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
						movemap->set_chi( (*ir)->seqpos(), true );
						protocols::simple_moves::MinMoverOP minmover = new protocols::simple_moves::MinMover( movemap, score_fxn, "dfpmin_armijo_nonmonotone_atol", 0.000001, true );
						minmover->apply( pose_dump );
						//core::io::pdb::dump_pdb( pose_dump, pose3_name_full.str() );
						core::Real sc_constraint_check( pose_dump.energies().total_energies()[ coordinate_constraint ] );
						ms_tr << "After sidechain refinement constraints score is " << sc_constraint_check << std::endl;
						/*if( data_ ) {
							data_output_file << "After sidechain refinement constraints score is " << sc_constraint_check << std::endl;
						}*/
						/*if( sc_constraint_check < 10.0 && sc_constraint_check < pre_sc_constraint_check ) {
							(*ir)->keep_rotamer( (pose_dump.residue((*ir)->seqpos())) );
							++hits;
							if( output_ ) {
								motif_output_file << *(motifhitop->motifcop());
								motif_output_file << "RESIDUE " << (pose_dump.residue((*ir)->seqpos()));
							}
						} else if( pre_sc_constraint_check < 10.0 && sc_constraint_check > pre_sc_constraint_check ) {
							pose_dump.replace_residue( (*ir)->seqpos(), *(motifhitop->build_rotamer()), true );
							(*ir)->keep_rotamer( (pose_dump.residue((*ir)->seqpos())) );
							++hits;
							if( output_ ) {
								motif_output_file << *(motifhitop->motifcop());
								motif_output_file << "RESIDUE " << (pose_dump.residue((*ir)->seqpos()));
							}
						}*/
						if( sc_constraint_check < 10.0 ) {
							(*ir)->keep_rotamer( (pose_dump.residue((*ir)->seqpos())) );
							(*ir)->keep_motif( *(motifhitop->motifcop()) );
							(*ir)->keep_motifhit( *(motifhitop) );
							++hits;
							if( output_ ) {
								motif_output_file << *(motifhitop->motifcop());
								motif_output_file << "RESIDUE " << (pose_dump.residue((*ir)->seqpos()));
							}
						}
						pose_dump.replace_residue( (*ir)->seqpos(), *(motifhitop->build_rotamer()), true );
						(*ir)->keep_rotamer( (pose_dump.residue((*ir)->seqpos())) );
						(*ir)->keep_motif( *(motifhitop->motifcop()) );
						(*ir)->keep_motifhit( *(motifhitop) );
						if( output_ ) {
							motif_output_file << *(motifhitop->motifcop());
							motif_output_file << "RESIDUE " << (pose_dump.residue((*ir)->seqpos()));
						}
						if( hits > rots2add_ ) break;
					}
				}
			}
		}
	}
	if( output_ ) {
		motif_output_file.close();
	}
	if( data_ ) {
		data_output_file.close();
	}
	if( quick_and_dirty_ ) {
		qd_output_file.close();
	}
}

core::pack::rotamer_set::Rotamers
MotifSearch::bp_rotamers(
	Size const seqpos
)
{
	core::pack::rotamer_set::Rotamers best_rotamers;
	for ( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
			ir != end_ir; ++ir ) {
		if( (*ir)->seqpos() != seqpos ) continue;
			if ( ! ((*ir)->best_rotamers()).empty() ) {
				std::set< std::string > allowedtypes( (*ir)->allowed_types() );
				if( allowedtypes.empty() ) {
					best_rotamers = (*ir)->best_rotamers();
				}
				Size rs( ((*ir)->best_rotamers()).size() );
				for ( Size r(1) ; r <= rs; ++r ) {
					for( std::set< std::string >::const_iterator ir2(allowedtypes.begin() ), end_ir2( allowedtypes.end() );
							ir2 != end_ir2; ++ir2 ) {
						if( (*ir2) == ((*ir)->best_rotamers()[r])->name3() ) {
							best_rotamers.push_back( (*ir)->best_rotamers()[r] );
						}
					}
				}
				//best_rotamers = (*ir)->best_rotamers();
			} else {
				ms_tr << "There were no rotamers to be included for position " << seqpos << std::endl;
			}
		}
	return best_rotamers;
}

protocols::motifs::MotifHitCOPs
MotifSearch::bp_motifhits(
	Size const seqpos
)
{
	protocols::motifs::MotifHitCOPs motifhitcops;
	for ( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
			ir != end_ir; ++ir ) {
		if( (*ir)->seqpos() != seqpos ) continue;
		if ( ! ((*ir)->best_motifhits()).empty() ) {
			Size rs( ((*ir)->best_motifhits()).size() );
			for ( Size r(1) ; r <= rs; ++r ) {
				motifhitcops.push_back( (*ir)->best_motifhits()[r] );
			}
		} else {
			ms_tr << "There were no motif hits for " << seqpos << ". Check to be sure that MotifSearch protocol actually ran. Use flag rerun_motifsearch to ensure that it runs even with input rotamers from BPData flag." << std::endl;
		}
	}
	return motifhitcops;
}

// Maybe this belongs with Motif.cc or MotifLibrary.cc, or both, but not here
bool
MotifSearch::protein_dna_motif()
{
	using namespace core::chemical;
	// This function only works for a motif that is made up of only two residues
	bool protein_dna( false );
	// motif_library_ has to be filled with the input MotifLibrary before you can call this fxn
	if ( !	motif_library_.empty() ) {
		// Check to see if motif has a protein component and a dna component
		//THIS IS RIDICULOUS, DON'T MAKE RESIDUES, MAKE RESIDUE TYPE FROM NAME3
			core::conformation::ResidueOP res1 = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map( protocols::dna::dna_full_name3( motif_library_[1]->restype_name1() ) ) );
			core::conformation::ResidueOP res2 = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map( protocols::dna::dna_full_name3( motif_library_[1]->restype_name2() ) ) );
		if	( ( res1->is_protein() && res2->is_DNA() ) || ( res2->is_protein() && res1->is_DNA() ) ) {
			protein_dna = true;
		}
	} else {
		ms_tr << "MotifLibrary has not been initialized yet, cannot yet identify the type of motifs being used, assuming protein-DNA with possible disastrous consequences." << std::endl;
		protein_dna = true;
	}
	return protein_dna;
}

// Will make this more general!  Should keep the chains separate
// Maybe make a map of vectors, keeping track of the different chains
void
MotifSearch::position_vector_setup(
	Pose const & pose
)
{
	for ( Size i(1), end( pose.total_residue() ); i <= end; ++i ) {
		if ( pose.residue_type(i).is_protein() ) {
			protein_positions_.push_back(i);
		}
		if ( pose.residue_type(i).is_DNA() ) {
			dna_positions_.push_back(i);
		}
	}
}

void
MotifSearch::identify_motif_build_positions(
	Pose const & pose,
	utility::vector1< Size > & build_positions
)
{
	if ( protein_dna_motif() ) {
		protein_DNA_motif_build_positions_JA( pose, build_positions, dna_positions_ ); //could easily change this line to take one of Phil's DNA interface fxns that also fills a vector
	} else {
		ms_tr << "ERROR! These motifs are not protein-DNA, need to add another else if statement that allows for a different type of interface finding function." << std::endl;
	}
}

void
MotifSearch::identify_motif_BuildPositions(
	Pose const & pose
)
{
	Sizes positions(0);
	if ( protein_dna_motif() ) {
		Sizes target_positions( map2keyvector( target_positions_ ) );
		protein_DNA_motif_build_positions_JA( pose, positions, target_positions );
		for ( Sizes::const_iterator pos( positions.begin() ), end( positions.end() );
				pos != end; ++pos ) {
			Size seqpos(*pos);
			Sizes short_target_positions( shorten_target_list( pose, seqpos, target_positions ) );
			std::set< std::string > allowed_types; // this vector will remain empty in this sitatuation since there is no input Def to limit the types of amino acids allowed
			if( restrict_to_wt_ ) {
				allowed_types.insert( pose.residue( seqpos ).name3() );
			}
			BuildPositionOP build_position = new BuildPosition( seqpos, short_target_positions, allowed_types );
			build_positionOPs_.push_back( build_position );
		}
	} else {
		ms_tr << "ERROR! These motifs are not protein-DNA, need to add another else if statement that allows for a different type of interface finding function." << std::endl;
	}
}

void
MotifSearch::BuildPosition_from_Size(
	Pose const & pose,
	Size const input_BP
)
{
	Sizes target_positions( map2keyvector( target_positions_ ) );
	Sizes short_target_positions( shorten_target_list( pose, input_BP, target_positions ) );
	std::set< std::string > allowed_types; // this set will remain empty in this sitatuation since there is no input Def to limit the types of amino acids allowed
	if( restrict_to_wt_ ) {
		allowed_types.insert( pose.residue( input_BP ).name3() );
	}
	BuildPositionOP build_position = new BuildPosition( input_BP, short_target_positions, allowed_types );
	build_positionOPs_.push_back( build_position );
}

void
MotifSearch::defs2BuildPositions(
	Pose const & pose,
	DnaDesignDefOPs const & defs
)
{
	if ( protein_dna_motif() ) {
		Sizes full_tl( map2keyvector( target_positions_ ) );
		std::map< Size, std::set< std::string > > mappositions( bpdefs2map( pose, defs ) );
		for ( std::map<Size, std::set< std::string > >::const_iterator it( mappositions.begin() ),
				end( mappositions.end() ); it != end; ++it ) {
			BuildPositionOP build_position = new BuildPosition( it->first, full_tl, it->second );
			build_positionOPs_.push_back( build_position );
		}
	} else {
		ms_tr << "ERROR! These motifs are not protein-DNA, need to add another else if statement that allows for a different type of interface finding function." << std::endl;
	}
}

void
MotifSearch::defs2BuildPositions_findts(
	Pose const & pose,
	DnaDesignDefOPs const & defs
)
{
	if ( protein_dna_motif() ) {
		Sizes full_tl( map2keyvector( target_positions_ ) );
		std::map< Size, std::set< std::string > > mappositions( bpdefs2map( pose, defs ) );
		for ( std::map<Size, std::set< std::string > >::const_iterator it( mappositions.begin() ),
				end( mappositions.end() ); it != end; ++it ) {
			Size test(it->first);
			Sizes short_tl( shorten_target_list( pose, test, full_tl ) );
			BuildPositionOP build_position = new BuildPosition( it->first, short_tl, it->second );
			build_positionOPs_.push_back( build_position );
		}
	} else {
		ms_tr << "ERROR! These motifs are not protein-DNA, need to add another else if statement that allows for a different type of interface finding function." << std::endl;
	}
}

utility::vector1< core::Size >
MotifSearch::map2keyvector(
	std::map< Size, std::set< std::string > > mappositions
)
{
	Sizes positions(0);
	for ( std::map<Size, std::set< std::string > >::const_iterator it( mappositions.begin() ),
				end( mappositions.end() ); it != end; ++it ) {
		positions.push_back( it->first );
	}
	return positions;
}

utility::vector1< core::Size >
MotifSearch::shorten_target_list(
	Pose const & pose,
	Size const bp,
	Sizes & full_tl
)
{
	Sizes short_tl(0);
	Sizes bps(0);
	bps.push_back( bp );
	protein_DNA_motif_target_positions_JA( pose, bps, full_tl, short_tl );
	return short_tl;
}

void
MotifSearch::protein_DNA_motif_build_positions_JA(
	Pose const & pose,
	Sizes & build_positions,
	Sizes & target_positions
)
{
	protocols::dna::DnaInterfaceFinderOP interface = new protocols::dna::DnaInterfaceFinder( 10*10, 3.9*3.9, 6., true );
	if ( ! target_positions.empty() ) { // again, this won't work well if there are multiple target positions in vector
		interface->determine_protein_interface( pose, protein_positions_, target_positions ); // unless target_positions_ is empty - will actually deal with that later, will fill with all DNA in initialize if no protein pos or dna pos are given
		protocols::dna::DnaNeighbors protein_neighbors = interface->protein_neighbors();
		for ( protocols::dna::DnaNeighbors::const_iterator itr( protein_neighbors.begin() ),
				end( protein_neighbors.end() ); itr != end; ++itr ) {
			if ( (*itr).second.contact() ) {
				build_positions.push_back( itr->first );
				ms_tr << "Positions being targeted for motif design " << itr->first << std::endl;
			}
		}
		ms_tr << "Attempting to identify build positions when there are no target positions." << std::endl;
	}
}

void
MotifSearch::protein_DNA_motif_target_positions_JA(
	Pose const & pose,
	Sizes & build_positions,
	Sizes & target_positions,
	Sizes & short_tl
)
{
	// JA used 3.7, I picked 3.9 to access a position I knew was important
	protocols::dna::DnaInterfaceFinderOP interface = new protocols::dna::DnaInterfaceFinder( 10*10, 3.9*3.9, 6., true );
	if ( ! build_positions.empty() ) {
		interface->determine_dna_interface( pose, build_positions, target_positions );
		protocols::dna::DnaNeighbors dna_neighbors = interface->dna_neighbors();
		for ( protocols::dna::DnaNeighbors::const_iterator itr( dna_neighbors.begin() ),
				end( dna_neighbors.end() ); itr != end; ++itr ) {
			if ( (*itr).second.contact() ) {
				short_tl.push_back( itr->first );
				ms_tr << "Positions (DNA) being targeted for motif design " << itr->first << std::endl;
			}
		}
		ms_tr << "Attempting to identify build positions when there are no target positions." << std::endl;
	}
}

void
MotifSearch::override_option_input(
	Real const & r1,
	Real const & z1,
	Real const & r2,
	Real const & z2,
	Real const & d1,
	Size const & rlevel,
	bool const bpdata,
	bool const bump_check
)
{
	ztest_cutoff_1_ = z1;
	ztest_cutoff_2_ = z2;
	rmsd_cutoff_1_ = r1;
	rmsd_cutoff_2_ = r2;
	dtest_cutoff_ = d1;
	rot_level_ = rlevel;
	bpdata_ = bpdata;
	bump_check_ = bump_check;
}

void
MotifSearch::reset_option_input()
{
	ztest_cutoff_1_ = basic::options::option[ basic::options::OptionKeys::motifs::z1 ]();
	ztest_cutoff_2_ =	basic::options::option[ basic::options::OptionKeys::motifs::z2 ]();
	rmsd_cutoff_1_ =	basic::options::option[ basic::options::OptionKeys::motifs::r1 ]();
	rmsd_cutoff_2_ =	basic::options::option[ basic::options::OptionKeys::motifs::r2 ]();
	dtest_cutoff_ =	basic::options::option[ basic::options::OptionKeys::motifs::dtest ]();
	rot_level_ =	basic::options::option[ basic::options::OptionKeys::motifs::rotlevel ]();
}

void
MotifSearch::set_motif_library(
	MotifLibrary & motiflibrary
)
{
	for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = motiflibrary.begin(), end_itr = motiflibrary.end();
			motifcop_itr != end_itr; ++motifcop_itr ) {
		protocols::motifs::MotifCOP motifcop( *motifcop_itr );
		motif_library_.push_back( motifcop );
	}
}

// should probably be a private fxn, so should other ones . . . need to organize code
void
MotifSearch::init_options()
{
	if( basic::options::option[ basic::options::OptionKeys::motifs::BPData ].user() ) {
		bpdata_ = true;
		bpdata_filename_ = basic::options::option[ basic::options::OptionKeys::motifs::BPData ]();
	} else {
		bpdata_ = false;
	}
	if( basic::options::option[ basic::options::OptionKeys::motifs::output_file ].user() ) {
		output_filename_ = basic::options::option[ basic::options::OptionKeys::motifs::output_file ]();
		output_ = true;
	} else {
		output_ = false;
	}
	if( basic::options::option[ basic::options::OptionKeys::motifs::data_file ].user() ) {
		data_filename_ = basic::options::option[ basic::options::OptionKeys::motifs::data_file ]();
		data_ = true;
	} else {
		data_ = false;
	}
	if( (basic::options::option[ basic::options::OptionKeys::motifs::quick_and_dirty ]).user() ) {
		quick_and_dirty_ = true;
	} else {
		quick_and_dirty_ = false;
	}
	if( (basic::options::option[ basic::options::OptionKeys::motifs::dump_motifs ]).user() ) {
		dump_motifs_ = true;
	} else {
		dump_motifs_ = false;
	}
	if( basic::options::option[ basic::options::OptionKeys::motifs::clear_bprots ].user() ) {
		clear_bprots_ = true;
	} else {
		clear_bprots_ = false;
	}
	if( basic::options::option[ basic::options::OptionKeys::motifs::rots2add ].user() ) {
		rots2add_ = basic::options::option[ basic::options::OptionKeys::motifs::rots2add ]();
	} else {
		rots2add_ = 100;
	}
	if( (basic::options::option[ basic::options::OptionKeys::motifs::restrict_to_wt ]).user() ) {
		restrict_to_wt_ = true;
	} else {
		restrict_to_wt_ = false;
	}
	if( (basic::options::option[ basic::options::OptionKeys::motifs::rerun_motifsearch ]).user() ) {
		rerun_motifsearch_ = true;
	} else {
		rerun_motifsearch_ = false;
	}
	if( (basic::options::option[ basic::options::OptionKeys::motifs::no_rotamer_bump ]).user() ) {
		bump_check_ = false;
	} else {
		bump_check_ = true;
	}


}

} // motifs
} // protocols
