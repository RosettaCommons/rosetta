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

// phil headers
#include <apps/pilot/phil/capri15_scoring.hh>
#include <apps/pilot/phil/loop_model.hh>


// libRosetta headers
#include <devel/dna/protocols.hh>
#include <devel/dna/util.hh>
#include <devel/dna/ProteinDNA_Relax.hh>
//#include <devel/dna/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/frags/TorsionFragment.hh>
#include <utility/excn/Exceptions.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/DME_FilterMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
//#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
//#include <core/scoring/elec/FA_ElecEnergy.hh>
//#include <core/scoring/etable/EtableEnergy.hh>
//#include <core/scoring/etable/count_pair/CountPairAll.hh>
//#include <core/scoring/etable/count_pair/CountPairFunction.hh>
//#include <core/scoring/etable/count_pair/CountPairFactory.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/VariantType.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>

#include <core/kinematics/FoldTree.hh>
#include <protocols/viewer/visualize.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/id/AtomID_Map.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
//#include <basic/options/after_opts.hh>

#include <basic/prof.hh> // profiling
#include <basic/basic.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/util.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// phil headers


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <sstream>
#include <math.h>

//silly using/typedef


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/OK.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>


using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;
using namespace protocols;

using utility::vector1;
using std::string;
using std::cout;
using std::endl;



static THREAD_LOCAL basic::Tracer TR( "apps.pilot.phil.loop_model" );

////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace id;
namespace OK = OptionKeys;
using utility::vector1;
using std::string;

// silly typedef
typedef std::map< std::string, std::string > SS_Map;


/// FWD DEC
void
capri15_relax(
							pose::Pose & pose,
							protocols::loops::Loops const & loops
							);

ScoreFunctionOP
get_relax_scorefxn();

void
setup_sam_constraints_for_t033(
															 pose::Pose & pose
															 );

SS_Map
get_t033_ss();

///////////////////////////////////////////////////////////////////////////////

void
get_sequence_and_secstruct_from_dssp(
																		 std::string const & filename,
																		 std::string & sequence,
																		 std::string & secstruct
																		 )
{
	sequence.clear();
	secstruct.clear();

	std::ifstream data( filename.c_str() );

	bool started( false );
	std::string line;
	while ( getline( data,line) ) {
		if ( !started ) {
			if ( line.substr(0,8) == "  #  RES" ) started = true;
		} else {
			sequence += line[13];
			char const sec = line[16];
			Size const pos = int_of( line.substr(0,5) );

			if ( sec == 'E' ) {
				secstruct.push_back( 'E' );
			} else if ( sec == 'H' || sec == 'I' || sec == 'G' ) {
				secstruct.push_back( 'H' );
			} else {
				secstruct.push_back( 'L' ); //all other B, S, and T
			}

			assert( pos == sequence.size() );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
std::string
filebase( std::string const & file )
{
	size_t found = file.find_last_of("/\\");
	if ( found == std::string::npos ) return file;
	else return file.substr(found+1);
}


/**

	 centroid score for docking / loop building:

	 "cyan" experimental data
	 reward RNA contacts to conserved residues in protein

	 for the latter, seems like we need an alignment to t033_.fasta ? In case we've done some trimming??

	 K/R     centroid contacts to phosphate backbone OP2/OP1
	 S/T/N/Q centroid contacts to phosphate backbone OP2/OP1
	 D/E/H   centroid contacts to O2'

	 vdw/hybrid vdw

	 backbone O to O2'
	 backbone N to OP2, OP1

	 distance between SAM CE and rGU N1

**/


/// @details  Sets up a trimmed t033 protein using the 1p91 structure with 1qao SAM in place and rna added
/// @note -- minimally trimmed... can be further trimmed stochastically as needed


/// @details  Given a pose with protein + SAM + rna, does random re-orientation and slide into contact to get
///  distribution of docked conformations


/// @details  Trims back a full-length t033 protein+sam pose to include only aligned residues
/// Used for trimming Rhiju's homology models back to only the aligned positions
///
/// @note  Uses trim_back_sequence_mapping hence results are stochastic


// returns a mapping from t033 fasta to the output pose sequence
id::SequenceMapping
trim_back_t033_pose( pose::Pose & pose )
{
	Size const min_loop_size( 5 );


	/// read a mapping from 1p91a sequence to t033 protein sequence
	id::SequenceMapping mapping;
	std::string source_seq, t033_seq;
	read_alignment_file( "input/alignment_for_trimming.txt", source_seq, t033_seq, mapping );
	assert( t033_seq == pose.sequence().substr(0,t033_seq.size() ) );


	// this is stochastic
	protocols::loops::trim_back_sequence_mapping( mapping, source_seq, t033_seq, min_loop_size );

	/// now go through and delete residues that are unaligned in t033_seq
	mapping.reverse();

	Size i_end(0);
	while( mapping[i_end+1 ] == 0 ) ++i_end;
	if ( i_end ) {
		TR.Warning << "NOT TRIMMING BACK AN N-TERMINAL LOOP!! " << i_end << std::endl;
	}

	Size i_begin( mapping.size1() );
	while( mapping[ i_begin ] == 0 ) --i_begin;
	if ( i_begin < mapping.size1() ) {
		TR.Warning << "NOT TRIMMING BACK A C-TERMINAL LOOP!! " << mapping.size1() - i_begin << std::endl;
	}

	id::SequenceMapping mapping_from_t033_to_pose( id::SequenceMapping::identity( t033_seq.size() ) );

	for ( Size i=i_begin; i>i_end; --i ) {
		if ( mapping[i] == 0 ) {
			mapping.delete_source_residue( i );
			assert( pose.residue(i).name1() == t033_seq[ i-1 ] );
			pose.conformation().delete_residue_slow( i );
			mapping_from_t033_to_pose.delete_target_residue( i );
		}
	}

	for ( Size i=1; i<= t033_seq.size(); ++i ) {
		if ( mapping_from_t033_to_pose[i] ) {
			assert( t033_seq[ i-1 ] == pose.residue( mapping_from_t033_to_pose[i] ).name1() );
		}
	}
	return mapping_from_t033_to_pose;
}


void
setup_bonded_protein_rna_pose(
															pose::Pose & pose,
															pose::Pose rna_pose // make a local copy
															)
{

	Size const sam_pos( pose.total_residue() ); // SAM comes last
	assert( rna_pose.total_residue() == 74 - 8 );
	Size const rna_root_pos( 54 - 6 );

	assert( rna_pose.residue( rna_root_pos ).aa() == na_rgu );

	core::pose::add_variant_type_to_pose_residue( rna_pose, "RGU_H1_DELETION", rna_root_pos );

	/// now append the rna by a bond to the methyl group of SAM
	///
	assert( pose.residue( sam_pos ).name3() == "SAM" );


	core::pose::add_variant_type_to_pose_residue( pose, "SAM_CE_CONNECT", sam_pos );


	std::string const anchor_atom_name( "CE" );
	std::string const rna_root_atom( "N1" );

	Residue const & root_rsd( rna_pose.residue( rna_root_pos ) );
	Residue const & anchor_rsd( pose.residue( sam_pos ) );
	assert(   root_rsd.n_non_polymeric_residue_connections() == 1 );
	assert( anchor_rsd.n_non_polymeric_residue_connections() == 1 );
	pose.append_residue_by_bond( root_rsd, false, root_rsd.n_possible_residue_connections(), sam_pos,
															 anchor_rsd.n_possible_residue_connections() );


	// now add the other rna residues
	for ( Size i=rna_root_pos-1; i>= 1; --i ) {
		pose.prepend_polymer_residue_before_seqpos( rna_pose.residue(i), pose.total_residue() - ( rna_root_pos-i )+1,
																								false );
	}
	for ( Size i=rna_root_pos+1; i<= rna_pose.total_residue(); ++i ) {
		pose.append_polymer_residue_after_seqpos( rna_pose.residue(i), pose.total_residue(), false );
	}
}

/// @details  This should be kept in sync with the next routine
utility::vector1< DOF_ID >
get_bonded_rna_dof_ids( pose::Pose const & pose )
{
	/// find the sequence positions of the SAM and the rna root residue
	Size sam_pos( 0 ), rna_root_pos( 0 );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).name3() == "SAM" ) sam_pos = i;
		else if ( pose.residue(i).is_RNA() ) {
			assert( sam_pos );
			if ( pose.residue(i).is_bonded( pose.residue(sam_pos) ) ) {
				rna_root_pos = i;
				break;
			}
		}
	}

	assert( sam_pos && rna_root_pos );

	Residue const & anchor_rsd( pose.residue( sam_pos ) );
	Residue const & root_rsd( pose.residue( rna_root_pos ) );


	AtomID const     SD( anchor_rsd.atom_index("SD"), sam_pos );
	AtomID const     CE( anchor_rsd.atom_index("CE"), sam_pos );
	AtomID const     N1(   root_rsd.atom_index("N1"), rna_root_pos );
	AtomID const  CHILD( pose.atom_tree().atom( N1 ).child(0)->id() );
	AtomID const GCHILD( pose.atom_tree().atom( CHILD ).child(0)->id() );

	utility::vector1< DOF_ID > dof_ids;
	Real tmpoffset;
	dof_ids.push_back( pose.atom_tree().bond_length_dof_id( CE, N1 ) );
	dof_ids.push_back( pose.atom_tree().bond_angle_dof_id( SD, CE, N1, tmpoffset ) );
	dof_ids.push_back( pose.atom_tree().bond_angle_dof_id( CE, N1, CHILD, tmpoffset ) );
	dof_ids.push_back( DOF_ID( N1    , id::PHI ) );
	dof_ids.push_back( DOF_ID( CHILD , id::PHI ) );
	dof_ids.push_back( DOF_ID( GCHILD, id::PHI ) );
	return dof_ids;
}

void
make_bonded_rna_move(
										 pose::Pose & pose,
										 Real const scale,
										 bool const reorient
										 )
{
	using numeric::conversions::radians;
	using numeric::random::gaussian;
	using numeric::random::uniform;

	/// find the sequence positions of the SAM and the rna root residue
	Size sam_pos( 0 ), rna_root_pos( 0 );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).name3() == "SAM" ) sam_pos = i;
		else if ( pose.residue(i).is_RNA() ) {
			assert( sam_pos );
			if ( pose.residue(i).is_bonded( pose.residue(sam_pos) ) ) {
				rna_root_pos = i;
				break;
			}
		}
	}

	assert( sam_pos && rna_root_pos );

	Residue const & anchor_rsd( pose.residue( sam_pos ) );
	Residue const & root_rsd( pose.residue( rna_root_pos ) );


	AtomID const     SD( anchor_rsd.atom_index("SD"), sam_pos );
	AtomID const     CE( anchor_rsd.atom_index("CE"), sam_pos );
	AtomID const     N1(   root_rsd.atom_index("N1"), rna_root_pos );
	AtomID const  CHILD( pose.atom_tree().atom( N1 ).child(0)->id() );
	AtomID const GCHILD( pose.atom_tree().atom( CHILD ).child(0)->id() );

	Real bond1, angle1, angle2, phi1, phi2, phi3;
	if ( reorient ) {
		// randomize
		bond1 = 3.5 + std::abs( gaussian() * 2.0 );
		angle1 = radians( 180.0  - std::abs( gaussian() * 30.0 ) );
		angle2 = radians( 120.0 + gaussian() * 45.0 );
		phi1 = radians( uniform() * 360.0 );
		phi2 = radians( uniform() * 360.0 );
		phi3 = radians( 0.0 /* inspection */+ gaussian() * 60.0 );
	} else {
		// perturb existing values
		bond1  = pose.conformation().bond_length( CE, N1 )        + scale * gaussian() * 0.5;
		angle1 = pose.conformation().bond_angle ( SD, CE,    N1 ) + scale * radians( gaussian() * 3.0 );
		angle2 = pose.conformation().bond_angle ( CE, N1, CHILD ) + scale * radians( gaussian() * 5.0 );
		phi1 = pose.dof( DOF_ID( N1    , id::PHI ) )              + scale * radians( gaussian() * 5.0 );
		phi2 = pose.dof( DOF_ID( CHILD , id::PHI ) )              + scale * radians( gaussian() * 5.0 );
		phi3 = pose.dof( DOF_ID( GCHILD, id::PHI ) )              + scale * radians( gaussian() * 2.5 );
	}

	/// set the new values
	pose.conformation().set_bond_length( CE, N1, bond1 );
	pose.conformation().set_bond_angle( SD, CE,    N1, angle1 );
	pose.conformation().set_bond_angle( CE, N1, CHILD, angle2 );
	pose.set_dof( DOF_ID( N1    , id::PHI ), phi1 );
	pose.set_dof( DOF_ID( CHILD , id::PHI ), phi2 );
	pose.set_dof( DOF_ID( GCHILD, id::PHI ), phi3 );
}


class BondedRNA_Mover : public protocols::moves::Mover {

public:

	BondedRNA_Mover( Real const scale ):
		protocols::moves::Mover( "BondedRNA_Mover" ),
		scale_( scale )
	{}


	virtual
	void
	apply( pose::Pose & pose )
	{
		make_bonded_rna_move( pose, scale_, false );
	}

private:
	Real scale_;
};


ScoreFunctionOP
get_loop_scorefxn()
{
	// loop building scorefxn
	//std::string function_tag("cen_std"), patch_tag("score4L");
	//ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(function_tag, patch_tag) );
	ScoreFunctionOP scorefxn( new ScoreFunction() );

	//	scorefxn->energy_method_options().atom_vdw_atom_type_set_name( chemical::HYBRID_FA_STANDARD_CENTROID );
	// Safer:
	methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.atom_vdw_atom_type_set_name( chemical::HYBRID_FA_STANDARD_CENTROID );
	scorefxn->set_energy_method_options( options );


	// censtd
	scorefxn->set_weight( env, 1.0 );
	scorefxn->set_weight( pair, 1.0 );
	scorefxn->set_weight( cbeta, 1.0 );
	scorefxn->set_weight( vdw, 1.0 );
	// score4l
	scorefxn->set_weight( rama, 0.1 );
	scorefxn->set_weight( rg, 2.0 );
	scorefxn->set_weight( chainbreak, 1.0 );
	scorefxn->set_weight( hbond_lr_bb, 1.0 );
	scorefxn->set_weight( hbond_sr_bb, 1.0 );

	scorefxn->set_weight( hybrid_vdw, 1.0 );

	std::map< ScoreType, Real > wts;
	wts[ capri_cen  ] = 0.4;
	wts[ capri_bb   ] = 0.4;
	scorefxn->add_extra_method( wts, CapriTwoBodyEnergy() );
	wts.clear();
	wts[ capri_blue ] = 2.0;
	wts[ capri_cyan ] = 2.0;
	wts[ capri_pink ] = 2.0;
	wts[ capri_red  ] = 2.0;
	wts[ capri_cons ] = 1.0;
	wts[ capri_dist ] = 1.0;
	scorefxn->add_extra_method( wts, CapriTotalEnergy() );
	return scorefxn;
}

ScoreFunctionOP
get_loop_scorefxn_old()
{
	// loop building scorefxn
	//std::string function_tag("cen_std"), patch_tag("score4L");
	//ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(function_tag, patch_tag) );
	ScoreFunctionOP scorefxn( new ScoreFunction() );

	//	scorefxn->energy_method_options().atom_vdw_atom_type_set_name( chemical::HYBRID_FA_STANDARD_CENTROID );
	//Safer:
	methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.atom_vdw_atom_type_set_name( chemical::HYBRID_FA_STANDARD_CENTROID );
	scorefxn->set_energy_method_options( options );

	// censtd
	scorefxn->set_weight( env, 1.0 );
	scorefxn->set_weight( pair, 1.0 );
	scorefxn->set_weight( cbeta, 1.0 );
	scorefxn->set_weight( vdw, 1.0 );
	// score4l
	scorefxn->set_weight( rama, 0.1 );
	scorefxn->set_weight( rg, 2.0 );
	scorefxn->set_weight( chainbreak, 1.0 );
	scorefxn->set_weight( hbond_lr_bb, 1.0 );
	scorefxn->set_weight( hbond_sr_bb, 1.0 );

	scorefxn->set_weight( hybrid_vdw, 1.0 );

	std::map< ScoreType, Real > wts;
	wts[ capri_cen  ] = 2.0;
	wts[ capri_bb   ] = 2.0;
	scorefxn->add_extra_method( wts, CapriTwoBodyEnergy() );
	wts.clear();
	wts[ capri_cyan ] = 2.0;
	wts[ capri_cons ] = 0.5;
	wts[ capri_dist ] = 2.0;
	scorefxn->add_extra_method( wts, CapriTotalEnergy() );
	return scorefxn;
}

ScoreFunctionOP
get_centroid_dock_scorefxn()
{
	ScoreFunctionOP centroid_scorefxn( new ScoreFunction() );

	centroid_scorefxn->set_weight( hybrid_vdw, 0.4 );

	std::map< ScoreType, Real > wts;
	wts[ capri_cen  ] = 0.4;
	wts[ capri_bb   ] = 0.4;
	centroid_scorefxn->add_extra_method( wts, CapriTwoBodyEnergy() );
	wts.clear();

	wts[ capri_blue ] = 2.0;
	wts[ capri_cyan ] = 2.0;
	wts[ capri_pink ] = 2.0;
	wts[ capri_red  ] = 2.0;
	wts[ capri_cons ] = 1.0;
	wts[ capri_dist ] = 1.0;
	centroid_scorefxn->add_extra_method( wts, CapriTotalEnergy() );

	return centroid_scorefxn;
}

ScoreFunctionOP
get_centroid_dock_scorefxn_old()
{
	ScoreFunctionOP centroid_scorefxn( new ScoreFunction() );

	centroid_scorefxn->set_weight( hybrid_vdw, 0.4 );

	std::map< ScoreType, Real > wts;
	wts[ capri_cen  ] = 2.0;
	wts[ capri_bb   ] = 2.0;
	centroid_scorefxn->add_extra_method( wts, CapriTwoBodyEnergy() );
	wts.clear();
	wts[ capri_cyan ] = 2.0;
	wts[ capri_cons ] = 0.5;
	wts[ capri_dist ] = 2.0;
	centroid_scorefxn->add_extra_method( wts, CapriTotalEnergy() );
	return centroid_scorefxn;
}

///////////////////////////////////////////////////////////////////////////////
/// @details  Take a trimmed, docked model and add the loops
protocols::loops::Loops
rebuild_trimmed_loops(
											pose::Pose & pose,
											std::string const & t033_seq,
											id::SequenceMapping const & mapping_from_t033_to_pose, // just to the protein rsds
											SS_Map const & ss_map
											)
{

	id::SequenceMapping mapping( mapping_from_t033_to_pose );
	mapping.reverse();

	std::string source_seq, target_seq( t033_seq );
	for ( Size i=1; i<= mapping.size1(); ++i ) {
		assert( mapping[i] );
		assert( pose.residue( i ).name1() == t033_seq[ mapping[i] - 1 ] );
		source_seq += t033_seq[ mapping[i] - 1 ];
	}


	protocols::loops::extend_sequence_mapping( pose /* const */, mapping, source_seq, target_seq );

	assert( source_seq == pose.sequence() && mapping.size1() == source_seq.size() &&
					mapping.size2()==target_seq.size() );

	protocols::loops::apply_sequence_mapping( pose, target_seq, mapping );

	assert( pose.sequence() == target_seq );
	assert( pose.sequence().substr(0,t033_seq.size()) == t033_seq );
	assert( t033_seq.size() == ss_map.begin()->second.size() );

	// now the protein sequence has been extended, have to update the conservation info stored in the pose
	id::SequenceMapping identity_mapping( id::SequenceMapping::identity( t033_seq.size() ) );
	setup_capri_data( pose, option[ OK::dna::specificity::tf ], identity_mapping );

	// now try loop building/refinement
	protocols::loops::Loops loops;

	// stochastic:
	setup_loops_from_mapping( pose, mapping, loops, true );


	std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
	std::map< Size, bool > frag_libs_init;

	// setup the ss quotas
	SS_Quotas ss_quotas;
	ss_quotas.push_back( std::make_pair( ss_map.find( "1p91A"      )->second, 40 ) );
	ss_quotas.push_back( std::make_pair( ss_map.find( "psipred"    )->second, 80 ) );
	ss_quotas.push_back( std::make_pair( ss_map.find( "sam_dssp"   )->second, 40 ) );
	ss_quotas.push_back( std::make_pair( ss_map.find( "sam_stride" )->second, 40 ) );

	setup_frags_from_vall( pose, loops, frag_libs, frag_libs_init, ss_quotas, true /*pick 6mers*/ );

	// dont skip long frags:
	protocols::loops::perturb_loops_with_ccd( pose, loops, frag_libs, get_loop_scorefxn() );

	return loops;
}
/// @details  Doesnt do any modifying of protein or rna internal coords
/// Just attaches rna to protein at sam and does random sampling of protein-RNA dofs


void
sample_rna_dofs_new(
										pose::Pose const & start_pose,
										std::string const & output_tag,
										Size const nstruct,
										std::string const & t033_seq,
										id::SequenceMapping const & mapping_from_t033_to_pose,
										SS_Map const & ss_map,
										bool const fa_relax
										)
{
// 	Size const n1( 5 );
// 	Size const n2_outer( 1 );
// 	Size const n2_inner( 10 );
 	Size const n1( 50 );
 	Size const n2_outer( 5 );
 	Size const n2_inner( 100 );
	Real const scale_factor( 0.2 );

	bool const fa_output( true );
	//bool const fa_relax( true );
	bool const rebuild_loops( true );

	retrieve_capri_data_from_pose( start_pose ); // assert that this has been setup

	for ( Size nn=1; nn<= nstruct; ++nn ) {
		std::string filename( output_tag + lead_zero_string_of( nn, 6 ) + ".pdb" );
		Pose pose( start_pose );

		ScoreFunctionOP centroid_scorefxn( get_centroid_dock_scorefxn() );


		protocols::moves::MonteCarlo mc( pose, *centroid_scorefxn, 2.0 );

		// random moves
		mc.set_temperature( 0.0 );
		for ( Size i=1; i<= n1; ++i ) {
			make_bonded_rna_move( pose, 0.0, true );
			mc.boltzmann( pose, "reorient" );
		}

		// perturb moves
		mc.set_temperature( 2.0 );
		for ( Size n=1; n<= n2_outer; ++n ) {
			Real const scale( scale_factor * ( 1.0 + Real( n2_outer-n ) / n2_outer ) );
			std::string const movetype( "pert_scale="+string_of( scale ) );
			for ( Size m=1; m<= n2_inner; ++m ) {
				make_bonded_rna_move( pose, scale, false );
				mc.boltzmann( pose, movetype );
			}
			mc.recover_low( pose );
			mc.show_counters();
		}

		Real finalscore( (*centroid_scorefxn)( pose ) );

		std::cout << "final_scores: " << finalscore << ' ' << filename << ' ' <<
			pose.energies().total_energies().weighted_string_of( centroid_scorefxn->weights() ) << std::endl;

		//pose.dump_pdb( filename );

		if ( rebuild_loops ) {

			if ( fa_output ) filename += "_loops_FArelax.pdb";
			else filename += "_loops.pdb";

			protocols::loops::Loops loops( rebuild_trimmed_loops( pose, t033_seq, mapping_from_t033_to_pose, ss_map ) );

			ScoreFunctionOP loop_scorefxn( get_loop_scorefxn() );
			finalscore = (*loop_scorefxn)(pose);
			std::cout << "final_loop_scores: " << finalscore << ' ' << filename << ' ' <<
				pose.energies().total_energies().weighted_string_of( loop_scorefxn->weights() ) << std::endl;

			if ( !fa_output ) pose.dump_pdb( filename );
			else {
				core::util::switch_to_residue_type_set( pose, FA_STANDARD );

				if ( fa_relax ) {
					setup_sam_constraints_for_t033( pose );
					capri15_relax( pose, loops );
					ScoreFunctionOP fa_scorefxn( get_relax_scorefxn() );
					finalscore = (*fa_scorefxn)( pose );
					std::cout << "final_FA_scores: " << finalscore << ' ' << filename << ' ' <<
						pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << std::endl;

					pose.dump_scored_pdb( filename, *fa_scorefxn );

				} else {

					setup_sam_constraints_for_t033( pose );
					ScoreFunctionOP fa_scorefxn( get_relax_scorefxn() );

					full_protein_repack( pose, *fa_scorefxn );

					finalscore = (*fa_scorefxn)(pose);

					std::cout << "final_FA_scores: " << finalscore << ' ' << filename << ' ' <<
						pose.energies().total_energies().weighted_string_of( fa_scorefxn->weights() ) << std::endl;

					pose.dump_scored_pdb( filename, *fa_scorefxn );
				}
			}
		}
	} // nstruct
}


///////////////////////////////////////////////////////////////////////////////
/// @details  Loops over starting structures, homology models for t033 including sam
void
capri_t033_centroid_trim_dock_test()
{
	/// cmdline args
	Size const nloop( option[ OK::dna::specificity::n_outer ] );
	Size const nstruct( option[ OK::out::nstruct ] ); // passed to rna-dof-sampling routine
	std::string const output_tag( option[ OK::out::output_tag ] );


	vector1< std::string > protein_files( start_files() );

	std::string const pssm_file( option[ OK::dna::specificity::tf ] );

	SS_Map const ss_map( get_t033_ss() );

	for ( Size k=1; k<= nloop; ++k ) {

		/// now the big loop
		for ( Size n=1; n<= protein_files.size(); ++n ) {

			/// read the homology model
			Pose protein_pose;
			core::import_pose::pose_from_file( protein_pose, protein_files[n] , core::import_pose::PDB_file);

			Pose rna_pose;
			core::import_pose::pose_from_file( rna_pose, "input/bound_rna.pdb" , core::import_pose::PDB_file);


			std::string const tag( output_tag +
														 "_protein" + string_of(n) +
														 "_loop" + string_of(k) );

			Pose pose( protein_pose );

			/// trims back unaligned positions, stochastic
			id::SequenceMapping mapping_from_t033_to_pose( trim_back_t033_pose( pose ) );

			std::string const t033_seq( protein_pose.sequence().substr(0,mapping_from_t033_to_pose.size1() ) );

			setup_bonded_protein_rna_pose( pose, rna_pose );

			setup_capri_data( pose, pssm_file, mapping_from_t033_to_pose );

			core::util::switch_to_residue_type_set( pose, chemical::HYBRID_FA_STANDARD_CENTROID );

			//pose.dump_pdb( "switched.pdb" );

			sample_rna_dofs_new( pose, tag, nstruct, t033_seq, mapping_from_t033_to_pose, ss_map, false );
		}
	}
}

void
trim_dock_rebuild_relax_test_rhiju()
{
	/// cmdline args
	//Size const nloop( option[ OK::dna::specificity::n_outer ] );
	Size const nstruct( option[ OK::out::nstruct ] ); // passed to rna-dof-sampling routine
	std::string const output_tag( option[ OK::out::output_tag ] );


	vector1< std::string > protein_files( start_files() );
	//vector1< std::string > rna_files( option[ OK::dna::specificity::frag_files ]() );

	std::string const pssm_file( option[ OK::dna::specificity::tf ] );

	SS_Map const ss_map( get_t033_ss() );

	for ( Size k=1; k<= nstruct; ++k ) {

		/// now the big loop
		for ( Size n=1; n<= protein_files.size(); ++n ) {

			/// read the homology model
			Pose protein_pose;
			core::import_pose::pose_from_file( protein_pose, protein_files[n] , core::import_pose::PDB_file);

			Pose rna_pose;
			core::import_pose::pose_from_file( rna_pose, "input/bound_rna.pdb" , core::import_pose::PDB_file);


			std::string const tag( output_tag +
														 "_protein" + string_of(n) +
														 "_loop" + string_of(k) );

			Pose pose( protein_pose );

			/// trims back unaligned positions, stochastic
			id::SequenceMapping mapping_from_t033_to_pose( trim_back_t033_pose( pose ) );

			std::string const t033_seq( protein_pose.sequence().substr(0,mapping_from_t033_to_pose.size1() ) );

			setup_bonded_protein_rna_pose( pose, rna_pose );

			setup_capri_data( pose, pssm_file, mapping_from_t033_to_pose );

			core::util::switch_to_residue_type_set( pose, chemical::HYBRID_FA_STANDARD_CENTROID );

			sample_rna_dofs_new( pose, tag, 10, t033_seq, mapping_from_t033_to_pose, ss_map, true );
		}
	}
}


/// @details  Loops over the loop-diversified models
/// Each one has t033 sequence plus SAM minus unaligned residues
///
/// MISNOMER AS WE ARE NOT DOING RELAX (NIGHT BEFORE TARGET IS DUE)
void
trim_dock_rebuild_relax_test()
{
	/// cmdline args
	Size const nstruct( option[ OK::out::nstruct ] );
	std::string const output_tag( option[ OK::out::output_tag ] );


	vector1< std::string > protein_files( start_files() );


	std::string const pssm_file( option[ OK::dna::specificity::tf ] );

	SS_Map const ss_map( get_t033_ss() );

	/// shuffle file lists
	utility::vector1< Size > protein_index;
	for ( Size i=1; i<= protein_files.size(); ++i ) protein_index.push_back(i);
	numeric::random::random_permutation( protein_index, numeric::random::rg() );

	for ( Size k=1; k<= nstruct; ++k ) {

		/// now the big loop
		for ( Size nn=1; nn<= protein_files.size(); ++nn ) {
			Size const n( protein_index[ nn ] );

			/// read the homology model
			Pose protein_pose;
			core::import_pose::pose_from_file( protein_pose, protein_files[ n ] , core::import_pose::PDB_file);

			Pose rna_pose;
			core::import_pose::pose_from_file( rna_pose, "input/bound_rna.pdb" , core::import_pose::PDB_file);


			std::string const tag( output_tag + "relax"
														 "_protein" + string_of(n) +
														 "_loop" + string_of(k) );

			Pose pose( protein_pose );


			//// setup a mapping from this starting model to t033
			id::SequenceMapping mapping_to_t033;
			std::string t033_seq;
			{
				id::SequenceMapping tmp_mapping;
				std::string source_seq, target_seq;
				read_alignment_file( "input/alignment_for_trimming.txt", source_seq, target_seq, tmp_mapping );
				t033_seq = target_seq;
				Size const nres_t033( t033_seq.size() );

				mapping_to_t033 = id::SequenceMapping::identity( nres_t033 );


				/// now go through and delete residues that are unaligned in target_seq
				tmp_mapping.reverse();
				// tmp_mapping is now a mapping from t033 to p91a
				// mapping_to_t033
				assert( tmp_mapping.size1() == nres_t033 );
				for ( Size i=nres_t033; i>=1; --i ) {
					if ( tmp_mapping[i] == 0 ) {
						tmp_mapping.delete_source_residue( i );
						mapping_to_t033.delete_source_residue( i );
					}
				}
			} // scope


			{ // trim back mapping to allow loop closure after docking
				std::string const protein_seq( protein_pose.sequence().substr(0,protein_pose.total_residue()-1 ) );
				assert( protein_seq.size() == mapping_to_t033.size1() );
				for ( Size i=1; i<= protein_seq.size(); ++i ) assert(pose.residue(i).name1() == t033_seq[mapping_to_t033[i]-1]);
				protocols::loops::trim_back_sequence_mapping( mapping_to_t033, protein_seq, t033_seq, 5 ); // seqs are const
			}

			{ // now delete residues from the pose that aren't aligned in this trimmed mapping
				assert( pose.total_residue() == protein_pose.total_residue() );
				for ( Size i=pose.total_residue()-1; i>0; --i ) {
					if ( mapping_to_t033[i] == 0 ) {
						pose.conformation().delete_residue_slow( i );
						mapping_to_t033.delete_source_residue( i );
					}
				}
				for ( Size i=1; i<= pose.total_residue()-1; ++i ) {
					assert( pose.residue(i).name1() == t033_seq[ mapping_to_t033[i]-1 ]);
				}
			}

			setup_bonded_protein_rna_pose( pose, rna_pose );

			// construct the reversed mapping for the next two routines
			id::SequenceMapping mapping_from_t033_to_pose( mapping_to_t033 );
			mapping_from_t033_to_pose.reverse();

			setup_capri_data( pose, pssm_file, mapping_from_t033_to_pose );

			core::util::switch_to_residue_type_set( pose, chemical::HYBRID_FA_STANDARD_CENTROID );

			sample_rna_dofs_new( pose, tag, 1, t033_seq, mapping_from_t033_to_pose, ss_map, false );

		} // protein_files
	} // nstruct
}

/// @details  Silly helper function
Size
pose_pos_from_pdb_pos( int const pdb_pos, char const pdb_chain, pose::Pose const & pose )
{
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.pdb_info()->number(i) == pdb_pos && pose.pdb_info()->chain(i) == pdb_chain ) return i;
	}
	TR.Fatal << "no such pdb pos+chain " << pdb_pos << ' ' << pdb_chain << std::endl;
	utility_exit();
	return 0;
}


void
setup_sam_constraints(
											pose::Pose & pose,
											id::SequenceMapping const & mapping_from_1p91A_to_pose,
											pose::Pose const & p91A_pose // for pdb numbering
											)
{
	/// get some important positions -- use the mapping to get sequence numbers
	char const P91_chain( 'A' );
	Size const pos93 ( mapping_from_1p91A_to_pose[ pose_pos_from_pdb_pos(  93, P91_chain, p91A_pose ) ]);
	Size const pos116( mapping_from_1p91A_to_pose[ pose_pos_from_pdb_pos( 116, P91_chain, p91A_pose ) ]);
	Size const pos117( mapping_from_1p91A_to_pose[ pose_pos_from_pdb_pos( 117, P91_chain, p91A_pose ) ]);
	Size const pos138( mapping_from_1p91A_to_pose[ pose_pos_from_pdb_pos( 138, P91_chain, p91A_pose ) ]);
	Size const pos139( mapping_from_1p91A_to_pose[ pose_pos_from_pdb_pos( 139, P91_chain, p91A_pose ) ]);
	Size const pos156( mapping_from_1p91A_to_pose[ pose_pos_from_pdb_pos( 156, P91_chain, p91A_pose ) ]);
	Size const pos158( mapping_from_1p91A_to_pose[ pose_pos_from_pdb_pos( 158, P91_chain, p91A_pose ) ]);
	Size sam_pos(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) if ( pose.residue(i).name3() == "SAM" ) sam_pos= i;
	assert( sam_pos );

	assert( pose.residue( pos93  ).aa() == aa_gly ); // hbond from O to SAM N
	assert( pose.residue( pos116 ).aa() == aa_asp ); // 2 hbonds from carboxyl Os to SAM ribose O2' O3'
	assert( pose.residue( pos117 ).aa() == aa_ile ); // pack against SAM base
	assert( pose.residue( pos138 ).aa() == aa_asp ); // hbond from carboxyl O1/2 to SAM N6
	assert( pose.residue( pos139 ).aa() == aa_leu ); // hbond from N to SAM N1
	assert( pose.residue( pos156 ).aa() == aa_phe ); // carbonyl O near SAM SD (3-4A?)
	assert( pose.residue( pos158 ).aa() == aa_pro ); // pack against SAM cleft between base and ribose
	assert( pose.residue( sam_pos ).name3() == "SAM" );

	using namespace scoring::constraints;
	ConstraintSetOP cst_set( new ConstraintSet() );
	typedef AtomPairConstraint APC;
	typedef HarmonicFunc HF;
	// pid just returns the AtomID of the given atom
	cst_set->add_constraint( new APC( pid(  "O", pos93 , pose ), pid(  "N" , sam_pos, pose ), new HF( 2.9, 0.25 ))); // G
	cst_set->add_constraint( new APC( pid( "CG", pos116, pose ), pid( "O2'", sam_pos, pose ), new HF( 3.3, 0.25 ))); // D
	cst_set->add_constraint( new APC( pid( "CG", pos116, pose ), pid( "O3'", sam_pos, pose ), new HF( 3.3, 0.25 ))); // "
	cst_set->add_constraint( new APC( pid( "CB", pos117, pose ), pid( "C2" , sam_pos, pose ), new HF( 3.7, 0.50 ))); // I
	cst_set->add_constraint( new APC( pid( "CB", pos117, pose ), pid( "C4" , sam_pos, pose ), new HF( 3.7, 0.50 ))); // "
	cst_set->add_constraint( new APC( pid( "CB", pos117, pose ), pid( "C6" , sam_pos, pose ), new HF( 3.7, 0.50 ))); // "
	cst_set->add_constraint( new APC( pid( "CG", pos138, pose ), pid( "N6" , sam_pos, pose ), new HF( 3.6, 0.50 ))); // D
	cst_set->add_constraint( new APC( pid(  "N", pos139, pose ), pid( "N1" , sam_pos, pose ), new HF( 2.9, 0.25 ))); // L
	cst_set->add_constraint( new APC( pid(  "O", pos156, pose ), pid( "SD" , sam_pos, pose ), new HF( 3.7, 0.50 ))); // F
	cst_set->add_constraint( new APC( pid( "CB", pos158, pose ), pid( "C8" , sam_pos, pose ), new HF( 4.5, 0.50 ))); // P
	cst_set->add_constraint( new APC( pid( "CB", pos158, pose ), pid( "C5'", sam_pos, pose ), new HF( 4.5, 0.50 ))); // "

	pose.constraint_set( cst_set );
}

void
setup_sam_constraints_for_t033(
															 pose::Pose & pose
															 )
{
	static bool init( false );
	static Pose p91A_pose;
	static id::SequenceMapping mapping_from_p91A_to_t033;

	if ( !init ) {
		init = true;

		core::import_pose::pose_from_file( p91A_pose, "input/1P91_chainA_w_1qao_SAM.pdb" , core::import_pose::PDB_file);

		std::string source_seq, target_seq;
		read_alignment_file( "input/alignment_for_trimming.txt", source_seq, target_seq, mapping_from_p91A_to_t033 );
		for ( Size i=1; i<= target_seq.size(); ++i ) assert( pose.residue(i).name1() == target_seq[ i-1 ] );
	}

	setup_sam_constraints( pose, mapping_from_p91A_to_t033, p91A_pose );

}

/// begins with full protein repack with soft soft rep, so can call immediately after core::util::switch_to_residue_type_set(FA)

void
juke_sam_pos(
						 pose::Pose & pose,
						 id::SequenceMapping const & mapping_from_1p91A_to_pose,
						 pose::Pose const & p91A_pose // for pdb numbering
						 )
{

	/// get some important positions -- use the mapping to get sequence numbers
	setup_sam_constraints( pose, mapping_from_1p91A_to_pose, p91A_pose );
	Size const sam_pos( pose.total_residue() );

	/// soft rep packer wts plus constraints
	ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( SOFT_REP_WTS ) );
	Real const save_fa_rep_weight( (*scorefxn)[ fa_rep ] );
	scorefxn->set_weight( atom_pair_constraint, 1.0 );


	/// rb-jump movemap + chimin!
	Size const sam_jump( pose.fold_tree().get_residue_edge( sam_pos ).label() );
	assert( Size( pose.fold_tree().downstream_jump_residue( sam_jump ) ) == sam_pos );
	kinematics::MoveMap mm;
	mm.set_jump( sam_jump, true );
	mm.set_chi( true );

	/// first pack and minimize with a lower rep weight //////////////////////////////////////////
	scorefxn->set_weight( fa_rep, 0.1 );

	// repack
	full_protein_repack( pose, *scorefxn );

	// minimize
	optimization::AtomTreeMinimizer().run( pose, mm, *scorefxn, optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true));

	/// now pack and minimize with a full strenght (soft) rep weight //////////////////////////////
	scorefxn->set_weight( fa_rep, save_fa_rep_weight );

	/// repack
	full_protein_repack( pose, *scorefxn );

	/// minimize
	optimization::AtomTreeMinimizer().run( pose, mm, *scorefxn, optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true));

	/// now try again with hard rep scorefxn ////////////////////////////////////////////////////
	scorefxn = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	scorefxn->set_weight( atom_pair_constraint, 1.0 );

	/// repack
	full_protein_repack( pose, *scorefxn );

	/// minimize
	optimization::AtomTreeMinimizer().run( pose, mm, *scorefxn, optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true));


}


ScoreFunctionOP
get_relax_scorefxn()
{


	ScoreFunctionOP scorefxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );

	// protein-NA mods
	scorefxn->set_weight( fa_pair, 0.0 );
	scorefxn->set_weight( fa_elec, 0.5 );

	// to hold the SAM in place
	scorefxn->set_weight( atom_pair_constraint, 1.0 );

	// to keep the loops closed
	scorefxn->set_weight( chainbreak, 1.0 );

	// these need to be optimized:
	std::map< ScoreType, Real > wts;
	wts[ capri_blue ] = 1.0;
	wts[ capri_cyan ] = 1.0;
	wts[ capri_pink ] = 1.0;
	wts[ capri_red  ] = 1.0;
	wts[ capri_cons ] = 0.5;
	wts[ capri_dist ] = 2.0;
	scorefxn->add_extra_method( wts, CapriTotalEnergy() );

	return scorefxn;
}

/**

	 start with a docked, loop-rebuilt pose containing protein+sam+rna

	 -- switch to fullatom residue set

	 -- ramp up the repulsive

	 -- scorefxn:

	    * standard wts
			* atompair constraints holding the sam
			* capri15 scores

	 -- dofs:

	    * protein-SAM jump
			* sam-rna dofs (6)
			* protein chi angles

	 -- moves: ( use nbr-atom dme filter as in zf-relax )

	    * sam-jump move
			* sam-rna dof move as in centroid docking, but with smaller scale
			* loop small moves

**/

void
capri15_relax(
							pose::Pose & pose,
							protocols::loops::Loops const & loops
							)
{

	// params
	Size const ramping_cycles( 2 );
	Real const ramping_initial_weight( 0.02 );
	Size const n_inner( option[ OK::dna::specificity::n_inner ] );
	Size const n_outer( option[ OK::dna::specificity::n_outer ] );
	Real const small_mover_temperature( 0.5 );
	Size const small_mover_nmoves( 3 );
	Real const dme_threshold( 0.12 );
	Size const max_tries( 20 );
	Real const bonded_rna_mover_scale( 0.35 );
	Real const trans_mag( 0.25 );
	Real const rot_mag( 3.0 );
	Real const energycut( 0.1 );
	Real const ramping_min_atol( 2.0 );
	Real const min_atol( 0.25 );
	Real const final_min_atol( 0.1 );

	ScoreFunctionOP scorefxn( get_relax_scorefxn() );

	//// setup task and move maps
	// pack_task:
	//  just repack protein residues, include current
	//
	// rotamer_trials_task: same as pack_task but maybe with more rotamers
	//
	pack::task::PackerTaskOP pack_task( pack::task::TaskFactory::create_packer_task( pose ));
	pack_task->initialize_from_command_line();
	pack_task->or_include_current( true );
	Size const nres( pose.total_residue() );	for ( Size i = 1; i <= nres; ++i ) {
		if ( pose.residue(i).is_protein() ) {
			pack_task->nonconst_residue_task( i ).restrict_to_repacking();
		} else {
			pack_task->nonconst_residue_task( i ).prevent_repacking();
		}
	}
	pack::task::PackerTaskOP rottrial_task( pack_task->clone() );
	for ( Size i = 1; i <= nres; ++i ) {
		if ( pose.residue(i).is_protein() ) {
			rottrial_task->nonconst_residue_task(i).or_ex1( true );
			rottrial_task->nonconst_residue_task(i).or_ex2( true );
			rottrial_task->nonconst_residue_task(i).or_ex1aro( true );
			rottrial_task->nonconst_residue_task(i).or_ex1aro_sample_level( pack::task::EX_THREE_THIRD_STEP_STDDEVS );
			rottrial_task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		}
	}

	// movemap:
	//  - protein chi
	//  - protein->SAM jump
	//  - SAM-RNA dofs (6)
	//  - loop backbone torsions
	//
	kinematics::MoveMapOP mm( new kinematics::MoveMap ); // default state is DOF FIXED
	Size sam_pos(0);
	for ( Size i = 1; i <= nres; ++i ) {
		if ( pose.residue(i).is_protein() ) {
			mm->set_chi( i, true );
			if ( loops.is_loop_residue( i ) ) mm->set_bb( i, true );
		} else if ( pose.residue(i).name3() == "SAM" ) sam_pos = i;
	}
	assert( sam_pos );
	Size const sam_jump( pose.fold_tree().get_residue_edge( sam_pos ).label() );
	mm->set_jump( sam_jump, true );
	// add the sam -> rna dofs
	{
		utility::vector1< DOF_ID > dof_ids( get_bonded_rna_dof_ids( pose ) );
		assert( dof_ids.size() == 6 );
		for ( Size i=1; i<= 6; ++i ) {
			mm->set( dof_ids[i], true );
		}
	}


	// 	{ // test mm dof min
	// 		mm->set_bb( false );
	// 		mm->set_chi( false );
	// 		mm->set_jump( false );

	// 		/// confirm that nangles == 6 here!
	// 		pose.dump_pdb( "start.pdb");
	// 		optimization::AtomTreeMinimizer().run( pose, *mm, *scorefxn, optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true));

	// 		pose.dump_pdb( "final.pdb");

	// 		exit(0);
	// 	}

	//// MonteCarlo object
	using namespace protocols::moves;
	using devel::dna::RB_Mover;
	using devel::dna::setup_MCM_trial;
	MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, 0.8 ) );


	//////////////////
	// setup movers:

	// packmover
	protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover( scorefxn, pack_task, 25 ) );

	// min mover
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( mm, scorefxn, "lbfgs_armijo_nonmonotone_atol", min_atol, true ) );

	// rb mover
	MoverOP sam_mover( new protocols::simple_moves::DME_FilterMover( new RB_Mover( mm, trans_mag, rot_mag ), dme_threshold, max_tries ) );

	// rna mover
	MoverOP rna_mover( new protocols::simple_moves::DME_FilterMover( new BondedRNA_Mover( bonded_rna_mover_scale ), dme_threshold, max_tries ) );

	// small mover
	MoverOP small_mover( new protocols::simple_moves::SmallMover( mm, small_mover_temperature, small_mover_nmoves ) );

	// rotamer trials w/ energycut
	protocols::simple_moves::EnergyCutRotamerTrialsMoverOP rottrial_mover
		( new protocols::simple_moves::EnergyCutRotamerTrialsMover( scorefxn, *rottrial_task, mc, energycut ) );

	// trials:
	TrialMoverOP sam_min_trial = setup_MCM_trial( sam_mover, rottrial_mover, min_mover, mc );

	TrialMoverOP rna_min_trial = setup_MCM_trial( rna_mover, rottrial_mover, min_mover, mc );

	TrialMoverOP small_min_trial = setup_MCM_trial( small_mover, rottrial_mover, min_mover, mc );

	TrialMoverOP pack_trial = new TrialMover( pack_mover, mc );

	TrialMoverOP min_trial = new TrialMover( min_mover, mc );

	// now the inner cycles mover
	MoverOP inner_cycle  = new SequenceMover( sam_min_trial, rna_min_trial, small_min_trial );
	MoverOP inner_cycles = new RepeatMover( inner_cycle, n_inner );

	// ramp up repulsive, dunbrack(?)
	//
	// this guy to the ramping_cycles power should be 1.0 / initial_weight
	//
	// the first ramping_cycles cycles are done at lower repulsive
	//
	core::Real const ramping_multiplier
		( ( ramping_cycles < 1 ) ? 1.0 : std::exp( std::log( 1.0 / ramping_initial_weight ) / ( ramping_cycles)));

	assert( std::abs( std::pow( ramping_multiplier, Real(ramping_cycles) ) - 1.0 / ramping_initial_weight ) < 1e-3 );

	core::Real ramping_weight_factor( ramping_initial_weight / ramping_multiplier );
	core::Real const final_fa_rep_weight( (*scorefxn)[ scoring::fa_rep ] );

	basic::prof_reset();
	for ( Size m1=1; m1<= n_outer; ++m1 ) {

		if ( ramping_cycles ) {
			if ( m1 <= ramping_cycles+1 ) {
				// ramp the weights
				ramping_weight_factor *= ramping_multiplier;
				scorefxn->set_weight( scoring::fa_rep, final_fa_rep_weight * ramping_weight_factor );
				mc->score_function( *scorefxn );
			}
			// set the min tolerance
			if ( m1 <= ramping_cycles ) min_mover->tolerance( ramping_min_atol );
			else min_mover->tolerance( min_atol );
			// sanity check:
			Real const current_rep_weight( scorefxn->get_weight( fa_rep ) );
			if ( m1 > ramping_cycles ) assert( std::abs( current_rep_weight - final_fa_rep_weight ) < 1e-3 );
			else assert( std::abs( current_rep_weight - final_fa_rep_weight ) > 1e-3 );
		}

		pack_trial->apply( pose );

		min_trial->apply( pose );

		inner_cycles->apply( pose );

		mc->recover_low( pose );

		mc->show_counters();

		basic::prof_show();
	}

	pack_trial->apply( pose );

	min_mover->tolerance( final_min_atol );
	min_trial->apply( pose );

	mc->recover_low( pose );
	basic::prof_show();


}


/**

	 read in a full-length model from recent simulations

	 make chemical bond modifications

	 setup foldtree:: attach sam by jump, attach rna by bond from sam

	 read alignment_for_trimming, setup an extended, trimmed mapping from 1p91a pdb to pose

	 setup loops from this mapping


**/

void
relax_test()
{
	ScoreFunctionOP scorefxn( get_relax_scorefxn() );

	Pose pose;
	core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file);


	// identify root and anchor positions
	Size const nres( pose.total_residue() );
	Size sam_pos(0);
	for ( Size i=1; i<= nres; ++i ) if ( pose.residue(i).name3() == "SAM" ) sam_pos = i;
	Size const rna_root_pos( sam_pos + 54 - 6 );
	Size const sam_anchor_pos( 102 ); // the GLY whose OC is hbonding to SAM N
	assert( sam_pos && ( nres == sam_pos + ( 74 - 8 ) ) );
	assert( pose.residue( rna_root_pos ).aa() == na_rgu );
	assert( pose.residue( sam_anchor_pos ).aa() == aa_gly );

	// chemically modify as necessary
	core::pose::add_variant_type_to_pose_residue( pose, "RGU_H1_DELETION", rna_root_pos );
	core::pose::add_variant_type_to_pose_residue( pose, "SAM_CE_CONNECT", sam_pos );
	pose.conformation().declare_chemical_bond( sam_pos, "CE", rna_root_pos, "N1" );
	// now want to set the jumps + bonds properly
	{
		// jump from pose to sam
		// chemical edge from sam to rna
		kinematics::FoldTree f( pose.total_residue() );
		//
		f.new_jump( sam_anchor_pos, sam_pos, sam_pos-1 );
		f.new_chemical_bond( sam_pos, rna_root_pos, "CE", "N1", sam_pos );
		pose.fold_tree( f );
		std::cout << "f2: " << f << std::endl;
		assert( pose.residue(      sam_pos ).is_bonded( pose.residue( rna_root_pos ) ) );
		assert( pose.residue( rna_root_pos ).is_bonded( pose.residue( sam_pos      ) ) );
	}

	// fill the capri15 data!
	id::SequenceMapping identity_mapping( id::SequenceMapping::identity( nres ) );
	setup_capri_data( pose, option[ OK::dna::specificity::tf ], identity_mapping );

	// score test
	(*scorefxn)(pose);

	//// setup loops object with guess as to loop definitions
	protocols::loops::Loops loops;
	{
		Pose p91A_pose;
		core::import_pose::pose_from_file( p91A_pose, "input/1P91_chainA_w_1qao_SAM.pdb" , core::import_pose::PDB_file);

		id::SequenceMapping mapping;
		std::string source_seq, target_seq;
		read_alignment_file( "input/alignment_for_trimming.txt", source_seq, target_seq, mapping );

		// stochastic:
		protocols::loops::trim_back_sequence_mapping( mapping, source_seq, target_seq, 5 ); // seqs are const

		protocols::loops::extend_sequence_mapping( p91A_pose /* const */, mapping, source_seq, target_seq );

		assert( source_seq == p91A_pose.sequence() && mapping.size1() == source_seq.size() &&
						mapping.size2()==target_seq.size() );

		// stochastic:
		// this will add cutpoints+jumps to the pose's foldtree
		// need to allow for not having rna sequence included in the mapping!
		setup_loops_from_mapping( pose, mapping, loops, false /*start_extended*/ );

		assert( !loops.is_loop_residue( sam_anchor_pos ) );
		std::cout << "prerelax foldtree: " << pose.fold_tree() << std::endl;
	} // scope

	// for graphics:
	protocols::viewer::add_conformation_viewer( pose.conformation(), "relax_pose" );

	capri15_relax( pose, loops );

}

/// @details  For calibrating the new centroid scores, rescore some models made with the old protocol
void
centroid_rescore_test()
{
	ScoreFunction scorefxn;
	std::map< ScoreType, Real > wts;
	wts[ capri_cen  ] = 1.0;
	wts[ capri_bb   ] = 1.0;
	scorefxn.add_extra_method( wts, CapriTwoBodyEnergy() );
	wts.clear();
	wts[ capri_blue ] = 1.0;
	wts[ capri_cyan ] = 1.0;
	wts[ capri_pink ] = 1.0;
	wts[ capri_red  ] = 1.0;
	wts[ capri_cons ] = 1.0;
	wts[ capri_dist ] = 1.0;
	scorefxn.add_extra_method( wts, CapriTotalEnergy() );

	ScoreFunctionOP loop_scorefxn( get_loop_scorefxn() );

	for ( Size n=1; n<= start_files().size(); ++n ) {

		Pose pose;
		core::import_pose::pose_from_file( pose, start_files()[n] , core::import_pose::PDB_file);

		// identify root and anchor positions
		Size const nres( pose.total_residue() );
		Size sam_pos(0);
		for ( Size i=1; i<= nres; ++i ) if ( pose.residue(i).name3() == "SAM" ) sam_pos = i;
		Size const rna_root_pos( sam_pos + 54 - 6 );
		Size const sam_anchor_pos( 102 ); // the GLY whose OC is hbonding to SAM N
		assert( sam_pos && ( nres == sam_pos + ( 74 - 8 ) ) );
		assert( pose.residue( rna_root_pos ).aa() == na_rgu );
		assert( pose.residue( sam_anchor_pos ).aa() == aa_gly );

		// chemically modify as necessary
		core::pose::add_variant_type_to_pose_residue( pose, "RGU_H1_DELETION", rna_root_pos );
		core::pose::add_variant_type_to_pose_residue( pose, "SAM_CE_CONNECT", sam_pos );
		pose.conformation().declare_chemical_bond( sam_pos, "CE", rna_root_pos, "N1" );
		// now want to set the jumps + bonds properly
		{
			// jump from pose to sam
			// chemical edge from sam to rna
			kinematics::FoldTree f( pose.total_residue() );
			//
			f.new_jump( sam_anchor_pos, sam_pos, sam_pos-1 );
			f.new_chemical_bond( sam_pos, rna_root_pos, "CE", "N1", sam_pos );
			pose.fold_tree( f );
			std::cout << "f2: " << f << std::endl;
			assert( pose.residue(      sam_pos ).is_bonded( pose.residue( rna_root_pos ) ) );
			assert( pose.residue( rna_root_pos ).is_bonded( pose.residue( sam_pos      ) ) );
		}

		// fill the capri15 data!
		id::SequenceMapping identity_mapping( id::SequenceMapping::identity( nres ) );
		setup_capri_data( pose, option[ OK::dna::specificity::tf ], identity_mapping );


		Real score( scorefxn( pose ) );

		// count protein-RNA nbrs:
		Size n_protein_rna_nbrs( 0 );
		{
			EnergyGraph const & energy_graph( pose.energies().energy_graph() );
			for ( Size seqpos=1; seqpos<= pose.total_residue(); ++seqpos ) {
				if ( !pose.residue(seqpos).is_protein() ) continue;
				for ( graph::Graph::EdgeListConstIter
								iru  = energy_graph.get_node( seqpos )->const_edge_list_begin(),
								irue = energy_graph.get_node( seqpos )->const_edge_list_end();
							iru != irue; ++iru ) {
					EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
					Size const other_pos( edge->get_other_ind( seqpos ) );
					if ( pose.residue( other_pos ).is_RNA() ) ++n_protein_rna_nbrs;
				}
			}
		}

		std::cout << "rescore_noloop: " << start_files()[n] << ' ' << score << ' ' << n_protein_rna_nbrs << ' ' <<
			pose.energies().total_energies().weighted_string_of( scorefxn.weights() ) << std::endl;


		// now rescore with loop scorefxn, 1st try
		core::util::switch_to_residue_type_set( pose, HYBRID_FA_STANDARD_CENTROID );
		score = (*loop_scorefxn)(pose);
		std::cout << "rescore_loop: " << start_files()[n] << ' ' << score << ' ' << n_protein_rna_nbrs << ' ' <<
			pose.energies().total_energies().weighted_string_of( loop_scorefxn->weights() ) << std::endl;


	}


}


/// @details  Diversify the conformation of the loop that passes beneath the SAM
///  protocol:
/// ** read 1P91A /w 1qao SAM
/// ** read alignment for trimming, delete unaligned target residues, extend, apply
/// ** setup loops file for the loop under SAM: "IFSPANY" in rlm2 sequence
/// ** delete SAM
/// ** convert to centroid
/// ** centroid perturb
/// ** convert to fullatom
/// ** copy in SAM
/// ** SAM-juke using Jim's constraints
///
/// ** full protein repack


void
diversify_sam_loop_test()
{
	Size const nstruct( option[ OK::out::nstruct ] );
	std::string const output_tag( option[ OK::out::output_tag ] );

	std::string const loopseq( "IFSPANY" );
	Pose p91A_pose, pose;
	core::import_pose::pose_from_file( p91A_pose, "input/1P91_chainA_w_1qao_SAM.pdb" , core::import_pose::PDB_file);

	pose = p91A_pose;

	id::SequenceMapping mapping;
	std::string source_seq, target_seq;
	read_alignment_file( "input/alignment_for_trimming.txt", source_seq, target_seq, mapping );

	/// now go through and delete residues that are unaligned in target_seq
	mapping.reverse();
	for ( Size i=mapping.size1(); i>=1; --i ) {
		if ( mapping[i] == 0 ) {
			mapping.delete_source_residue( i );
			target_seq.erase( i-1, 1 );
		}
	}
	mapping.reverse();

	protocols::loops::extend_sequence_mapping( pose /* pose is const */, mapping, source_seq, target_seq );

	assert( source_seq == pose.sequence() && mapping.size1() == source_seq.size() &&
					mapping.size2()==target_seq.size() );

	protocols::loops::apply_sequence_mapping( pose, target_seq, mapping );

	assert( pose.sequence() == target_seq );

	//pose.dump_pdb("trimmed.pdb");

	Pose const start_pose( pose );

	// setup fragment files
	assert( target_seq.find( loopseq ) != std::string::npos );
	Size const loop_size( loopseq.size() );
	Size const loop_begin( target_seq.find( loopseq )+ 1 );
	Size const loop_end( loop_begin + loop_size - 1 );
	std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
	std::map< Size, bool > frag_libs_init;

	// nstruct loop here -- stochastically choose cutpoint for rebuilding

	for ( Size n=1; n<= nstruct; ++n ) {

		pose = start_pose;

		//// convert to centroid mode
		Size const sam_pos( pose.total_residue() );
		assert( pose.residue( sam_pos ).name3() == "SAM" );
		// delete sam
		Pose sam_pose;
		sam_pose.append_residue_by_bond( pose.residue( sam_pos ) );
		pose.conformation().delete_residue_slow( sam_pos );
		for ( Size i=1; i<= pose.total_residue(); ++i ) assert( pose.residue(i).is_protein() );

		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

		//// setup for rebuilding the loop under the sam
		Size const cutpoint( loop_begin - 1 + static_cast< int >( numeric::random::uniform() * (loop_size+1)) );
		//Size const cutpoint( loop_begin+2 );
		//assert( pose.residue(cutpoint+1).aa() == aa_pro );
		assert( cutpoint >= loop_begin - 1 && cutpoint <= loop_end );

		// setup the loops object
		protocols::loops::Loops loops;
		bool const start_extended( numeric::random::uniform() < 0.25 );

		loops.add_loop( loop_begin, loop_end, cutpoint, 0.0, 0 /*offset*/, start_extended );
		if ( n == 1 ) setup_frags_from_vall( pose, loops, frag_libs, frag_libs_init );


		basic::prof_reset();

		protocols::loops::set_loop_cutpoint_in_pose_fold_tree( cutpoint, pose, loop_begin, loop_end );
		protocols::loops::perturb_loops_with_ccd( pose, loops, frag_libs );

		basic::prof_show();

		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );

		Size sam_anchor_pos( loops.begin()->stop() + 1 );
		assert(  loops.is_loop_residue( loops.begin()->start() ) );
		assert( !loops.is_loop_residue( sam_anchor_pos ) );

		pose.append_residue_by_jump( sam_pose.residue(1), sam_anchor_pos );

		juke_sam_pos( pose, mapping, p91A_pose );


		std::string filename( output_tag + "loop_juke" );
		if ( start_extended ) filename += "_extended"+ lead_zero_string_of( n, 4 ) + ".pdb";
		else filename += "_folded"+ lead_zero_string_of( n, 4 ) + ".pdb";


		ScoreFunctionOP scorefxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );

		Real const final_score( (*scorefxn)(pose) );
		std::cout << "final_scores: "<< final_score << ' ' << filename << ' ' <<
			pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;
		pose.dump_pdb( filename );
	}

}


void
capri_t033_trim_dock_test()
{
	// read the pose -- 1P91.pdb, chain A  with SAM group
	Pose pose;
	core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file);


	/// read in the original full crazy alignment
	id::SequenceMapping mapping;
	std::string source_seq, target_seq;
	read_alignment_file( option[ OK::phil::align_file ], source_seq, target_seq, mapping );

	/// now go through and delete residues that are unaligned in target_seq
	mapping.reverse();
	for ( Size i=mapping.size1(); i>=1; --i ) {
		if ( mapping[i] == 0 ) {
			mapping.delete_source_residue( i );
			target_seq.erase( i-1, 1 );
		}
	}
	mapping.reverse();

	protocols::loops::extend_sequence_mapping( pose /* pose is const */, mapping, source_seq, target_seq );

	assert( source_seq == pose.sequence() && mapping.size1() == source_seq.size() &&
					mapping.size2()==target_seq.size() );

	protocols::loops::apply_sequence_mapping( pose, target_seq, mapping );

	assert( pose.sequence() == target_seq );

	//pose.dump_pdb("trimmed.pdb");

	/// get some important positions -- use the mapping to get sequence numbers
	char const P91_chain( 'A' );
	Size const pos93 ( mapping[ pose_pos_from_pdb_pos(  93, P91_chain, pose ) ]);
	Size const pos116( mapping[ pose_pos_from_pdb_pos( 116, P91_chain, pose ) ]);
	Size const pos117( mapping[ pose_pos_from_pdb_pos( 117, P91_chain, pose ) ]);
	Size const pos138( mapping[ pose_pos_from_pdb_pos( 138, P91_chain, pose ) ]);
	Size const pos139( mapping[ pose_pos_from_pdb_pos( 139, P91_chain, pose ) ]);
	Size const pos156( mapping[ pose_pos_from_pdb_pos( 156, P91_chain, pose ) ]);
	Size const pos158( mapping[ pose_pos_from_pdb_pos( 158, P91_chain, pose ) ]);
	Size const sam_pos( pose.total_residue() );

	assert( pose.residue( pos93  ).aa() == aa_gly ); // hbond from O to SAM N
	assert( pose.residue( pos116 ).aa() == aa_asp ); // 2 hbonds from carboxyl Os to SAM ribose O2' O3'
	assert( pose.residue( pos117 ).aa() == aa_ile ); // pack against SAM base
	assert( pose.residue( pos138 ).aa() == aa_asp ); // hbond from carboxyl O1/2 to SAM N6
	assert( pose.residue( pos139 ).aa() == aa_leu ); // hbond from N to SAM N1
	assert( pose.residue( pos156 ).aa() == aa_phe ); // carbonyl O near SAM SD (3-4A?)
	assert( pose.residue( pos158 ).aa() == aa_pro ); // pack against SAM cleft between base and ribose


	/// try optimizing the sam position subject to some constraints
	{
		using namespace scoring::constraints;
		ConstraintSetOP cst_set( new ConstraintSet() );
		typedef AtomPairConstraint APC;
		typedef HarmonicFunc HF;
		// pid just returns the AtomID of the given atom
		cst_set->add_constraint( new APC( pid(  "O", pos93 , pose ), pid(  "N" , sam_pos, pose ), new HF( 2.9, 0.25 ))); // G
		cst_set->add_constraint( new APC( pid( "CG", pos116, pose ), pid( "O2'", sam_pos, pose ), new HF( 3.3, 0.25 ))); // D
		cst_set->add_constraint( new APC( pid( "CG", pos116, pose ), pid( "O3'", sam_pos, pose ), new HF( 3.3, 0.25 ))); // "
		cst_set->add_constraint( new APC( pid( "CB", pos117, pose ), pid( "C2" , sam_pos, pose ), new HF( 3.7, 0.50 ))); // I
		cst_set->add_constraint( new APC( pid( "CB", pos117, pose ), pid( "C4" , sam_pos, pose ), new HF( 3.7, 0.50 ))); // "
		cst_set->add_constraint( new APC( pid( "CB", pos117, pose ), pid( "C6" , sam_pos, pose ), new HF( 3.7, 0.50 ))); // "
		cst_set->add_constraint( new APC( pid( "CG", pos138, pose ), pid( "N6" , sam_pos, pose ), new HF( 3.6, 0.50 ))); // D
		cst_set->add_constraint( new APC( pid(  "N", pos139, pose ), pid( "N1" , sam_pos, pose ), new HF( 2.9, 0.25 ))); // L
		cst_set->add_constraint( new APC( pid(  "O", pos156, pose ), pid( "SD" , sam_pos, pose ), new HF( 3.7, 0.50 ))); // F
		cst_set->add_constraint( new APC( pid( "CB", pos158, pose ), pid( "C8" , sam_pos, pose ), new HF( 4.5, 0.50 ))); // P
		cst_set->add_constraint( new APC( pid( "CB", pos158, pose ), pid( "C5'", sam_pos, pose ), new HF( 4.5, 0.50 ))); // "


		pose.constraint_set( cst_set );

		/// soft rep packer wts plus constraints
		ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( SOFT_REP_WTS ) );
		Real const save_fa_rep_weight( (*scorefxn)[ fa_rep ] );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );


		/// rb-jump movemap + chimin!
		assert( pose.num_jump() == 1 );
		kinematics::MoveMap mm;
		mm.set_jump( true );
		mm.set_chi( true );

		/// first pack and minimize with a lower rep weight //////////////////////////////////////////
		scorefxn->set_weight( fa_rep, 0.1 );

		// repack
		pose.dump_scored_pdb("prepack1.pdb", *scorefxn );
		full_protein_repack( pose, *scorefxn );
		pose.dump_scored_pdb("afterpack1.pdb", *scorefxn );

		// minimize
		optimization::AtomTreeMinimizer().run( pose, mm, *scorefxn, optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true));
		pose.dump_scored_pdb("aftermin1.pdb", *scorefxn );

		/// now pack and minimize with a full strenght (soft) rep weight //////////////////////////////
		scorefxn->set_weight( fa_rep, save_fa_rep_weight );

		/// repack
		pose.dump_scored_pdb("prepack2.pdb", *scorefxn );
		full_protein_repack( pose, *scorefxn );
		pose.dump_scored_pdb("afterpack2.pdb", *scorefxn );

		/// minimize
		optimization::AtomTreeMinimizer().run( pose, mm, *scorefxn, optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true));
		pose.dump_scored_pdb("aftermin2.pdb", *scorefxn );

		/// now try again with hard rep scorefxn ////////////////////////////////////////////////////
		scorefxn = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );

		/// repack
		pose.dump_scored_pdb("prepack3.pdb", *scorefxn );
		full_protein_repack( pose, *scorefxn );
		pose.dump_scored_pdb("afterpack3.pdb", *scorefxn );

		/// minimize
		optimization::AtomTreeMinimizer().run( pose, mm, *scorefxn, optimization::MinimizerOptions("lbfgs_armijo_nonmonotone",0.001,true));
		pose.dump_scored_pdb("aftermin3.pdb", *scorefxn );


		///////////////////////////////////// all done
	}

	/// read rhiju's model for the hairpin
	Pose rna_pose;
	core::import_pose::pose_from_file( rna_pose, "rna_full_model.pdb", core::import_pose::PDB_file);
	//core::import_pose::pose_from_file( rna_pose, "rna_hairpin_model.pdb", core::import_pose::PDB_file);

	assert( false );
	Size const rna_root_pos( 54 );
// 	Size const rna_root_pos( 5 );

	core::pose::add_variant_type_to_pose_residue( rna_pose, "RGU_H1_DELETION", rna_root_pos );

	/// now append the rna by a bond to the methyl group of SAM
	///
	assert( pose.residue( sam_pos ).name3() == "SAM" );


	core::pose::add_variant_type_to_pose_residue( pose, "SAM_CE_CONNECT", sam_pos );


	std::string const anchor_atom_name( "CE" );
	assert( rna_pose.residue( rna_root_pos ).aa() == na_rgu );
	std::string const rna_root_atom( "N1" );

	{
		Residue const & root_rsd( rna_pose.residue( rna_root_pos ) );
		Residue const & anchor_rsd( pose.residue( sam_pos ) );
		assert(   root_rsd.n_non_polymeric_residue_connections() == 1 );
		assert( anchor_rsd.n_non_polymeric_residue_connections() == 1 );
		pose.append_residue_by_bond( root_rsd, false, root_rsd.n_possible_residue_connections(), sam_pos,
																 anchor_rsd.n_possible_residue_connections() );


		// now add the other rna residues
		for ( Size i=rna_root_pos-1; i>= 1; --i ) {
			pose.prepend_polymer_residue_before_seqpos( rna_pose.residue(i), pose.total_residue() - ( rna_root_pos-i )+1,
																									false );
		}
		for ( Size i=rna_root_pos+1; i<= rna_pose.total_residue(); ++i ) {
			pose.append_polymer_residue_after_seqpos( rna_pose.residue(i), pose.total_residue(), false );
		}


		// set things up a bit
		//
		// want to set the distance from CE to N1 to be ~3.5?
		// the SD-CE-N1 angle to ~180
		// the CE-N1-1stchild angle to be ~120??
		// searching over the SD-CE-N1-1st child bond torsion PHI
		//
		Size const new_rna_root_pos( pose.total_residue() - ( rna_pose.total_residue() - rna_root_pos ) );
		assert( pose.residue( new_rna_root_pos ).aa() == na_rgu );
		assert( pose.residue( new_rna_root_pos ).type().has_variant_type( "RGU_H1_DELETION" ) );
		using numeric::conversions::radians;
		AtomID const     SD( anchor_rsd.atom_index("SD"), sam_pos );
		AtomID const     CE( anchor_rsd.atom_index("CE"), sam_pos );
		AtomID const     N1(   root_rsd.atom_index("N1"), new_rna_root_pos );
		AtomID const  CHILD( pose.atom_tree().atom( N1 ).child(0)->id() );
		AtomID const GCHILD( pose.atom_tree().atom( CHILD ).child(0)->id() );

		std::cout << N1 << std::endl;
		std::cout << CHILD << std::endl;
		std::cout << GCHILD << std::endl;

		pose.conformation().set_bond_length( CE, N1, 3.5 );
		pose.conformation().set_bond_angle( SD, CE, N1, radians( 180.0 ) );
		pose.conformation().set_bond_angle( CE, N1, CHILD, radians( 120.0 ) );

		{
			DOF_ID tor_id( GCHILD, id::PHI );
			pose.set_dof( tor_id, 0.0 ); // by inspection
		}

		//conformation::show_atom_tree( pose.atom_tree().atom( CE ), pose.conformation(), std::cout );
		//pose.atom_tree().atom( CE ).show();


		/// now try exploring different values of the rna dofs
		///
		ScoreFunction soft_rep_scorefxn;
		soft_rep_scorefxn.set_etable( FA_STANDARD_SOFT );
		soft_rep_scorefxn.set_weight( fa_rep, 1.0 );

		Size const nstruct( 99999 );
		Real min_combo_score( 1e9 );
		for ( Size n=1; n<= nstruct; ++n ) {
			using namespace numeric::random;
			pose.conformation().set_bond_length( CE, N1, 3.5 + std::abs( gaussian() * 2.0 ) );
			pose.conformation().set_bond_angle( SD, CE, N1, radians( 180.0  - std::abs( gaussian() * 30.0 ) ) );
			pose.conformation().set_bond_angle( CE, N1, CHILD, radians( 120.0 + gaussian() * 30.0 ) );
			pose.set_dof( DOF_ID( GCHILD, id::PHI ), radians( 0.0 + gaussian() * 60.0 ) );
			pose.set_dof( DOF_ID( CHILD , id::PHI ), radians( uniform() * 360.0 ) );

			Real const rep_score( soft_rep_scorefxn( pose ) );
			assert( false );
			Real const contact_score( 0.0 ); //eval_protein_rna_contact_dev( pose ) );
			Real const combo_score( contact_score * 100 + rep_score );
			if ( combo_score < min_combo_score ) {
				min_combo_score = combo_score;
				pose.dump_pdb( "contact"+string_of(n)+".pdb" );
			}
			std::cout << "SCORES: " << n << " combo= " << combo_score << " rep= " << rep_score << " contact= " <<
				contact_score << std::endl;
		}
	}

	assert( false );
	exit(0);
}


///////////////////////////////////////////////////////////////////////////////

void
capri_t033_loop_test()
{

	// read the pose
	Pose pdb_pose;
	core::import_pose::pose_from_file( pdb_pose, basic::options::start_file() , core::import_pose::PDB_file);

	// starting point for each nstruct loop
	Pose const start_pose( pdb_pose );

	Size const nstruct( option[ OK::out::nstruct ] );
	std::string const output_tag( option[ OK::out::output_tag ] );


	for ( Size n=1; n<= nstruct; ++n ) {
		Pose pose( start_pose );

		// read the alignment file
		id::SequenceMapping mapping;
		std::string source_seq, target_seq;
		read_alignment_file( option[ OK::phil::align_file ], source_seq, target_seq, mapping );

		std::string const filename( output_tag+"loopmodel"+lead_zero_string_of( n, 4 )+".pdb" );

		// this is stochastic
		protocols::loops::trim_back_sequence_mapping( mapping, source_seq, target_seq, 5 ); //seqs are const refs

		protocols::loops::extend_sequence_mapping( pose /* const */, mapping, source_seq, target_seq );

		assert( source_seq == pose.sequence() && mapping.size1() == source_seq.size() &&
						mapping.size2()==target_seq.size() );

		protocols::loops::apply_sequence_mapping( pose, target_seq, mapping );

		assert( pose.sequence() == target_seq );


		// remove non-protein residues from pose, store them in sam_pose, preserve a mapping from sam_pose residues
		// to their original locations in pose
		//
		//id::SequenceMapping sam_pose_mapping;
		Pose sam_pose;
		{ // hacky short-term thing
			Size sam_pos(0);

			for ( Size i=pose.total_residue(); i>= 1; --i ) {
				if ( !pose.residue(i).is_protein() ) {
					assert( !sam_pos );
					sam_pos = i;
					sam_pose.append_residue_by_jump( pose.residue(i), 1 );
					//sam_pose.insert_residue_by_jump( pose.residue(i),
					pose.conformation().delete_residue_slow( i );
					mapping.delete_target_residue( i );
				}
			}
			assert( sam_pose.total_residue() == 1 );
		}


		// now try loop building/refinement
		protocols::loops::Loops loops;

		setup_loops_from_mapping( pose, mapping, loops, true );


		std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
		std::map< Size, bool > frag_libs_init;
		setup_frags_from_vall( pose, loops, frag_libs, frag_libs_init );

		{ // now do the centroid modeling

			// delete non-protein positions
			for ( Size i=pose.total_residue(); i>= 1; --i ) {
				if ( !pose.residue(i).is_protein() ) pose.conformation().delete_residue_slow( i );
			}

			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

			// for graphics:
			//protocols::viewer::add_conformation_viewer( pose.conformation(), "loops_pose" );

			basic::prof_reset();

			protocols::loops::perturb_loops_with_ccd( pose, loops, frag_libs );

			basic::prof_show();

			core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );

			// create a fullatom scorefxn
			ScoreFunctionOP scorefxn( get_score_function() );
			scorefxn->set_weight( chainbreak, 1.0 ); // unnecessary

			Size sam_anchor_pos( loops.begin()->stop() + 1 );
			assert( loops.is_loop_residue( loops.begin()->start() ) );
			assert( !loops.is_loop_residue( sam_anchor_pos ) );

			pose.append_residue_by_jump( sam_pose.residue(1), sam_anchor_pos );

			full_protein_repack( pose, *scorefxn );

			refine_loops_with_ccd( pose, loops );

			Real const final_score( (*scorefxn)(pose) );
			std::cout << "final_score: " << filename << ' ' << final_score << ' ' <<
				pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;

			pose.dump_scored_pdb( filename, *scorefxn );

			basic::prof_show();
		}
	} // nstruct loop

	exit(0);

}

char
map_sec( char const sec )
{
	if ( sec == 'E' ) {
		return 'E';
	} else if ( sec == 'H' || sec == 'I' || sec == 'G' ) {
		return 'H';
	} else {
		std::cout << "mapped to L: " << sec << std::endl;
		return 'L';
	}
}

std::string
map_secstruct( std::string const & ss )
{
	std::string tmp;
	for ( Size i=0; i< ss.size(); ++i ) {
		tmp += map_sec( ss[i] );
	}
	return tmp;
}


/// @details   Returns a map from type to ss, indexed by: "1p91A", "psipred", "sam_stride", "sam_dssp"


SS_Map
get_t033_ss()
{


	std::string t033_psipred    = map_secstruct( "CCCCCHHHHHHCCHHHCCCCCCCCHHEECCCEEECCCCCCCCCCCCCEEECCCCCCCCCCCCHHHHHHHHHHHHCCCCHHHHHHHHHHHHHCCCCCEEEEECCCCCHHHHHHHHHCCCCEEEEEECCHHHHHHHHHHCCCCCCEEEEECHHHCCCCCCCEEEEEECCCHHHHHHHHHHCCCCEEEEEEECCCHHHHHHHHHHHCCCCCCCCCHHHHHHHHHCCCCEEEEEEEEEEEECCHHHHHHHHHCCHHHCCCCHHHHHHCCCCCEEEEEEEEEEEEECC" );

	std::string t033_sam_dssp   = map_secstruct( "CCCCCCCHCHHHHHHHHHCTTTHHHHEHETTEEEECTTCEEEHHHHCCEECCCCCCCCTTCHHHHHHHHHHHHHHTCCHHHHHHHHHHHHTCCTTCEEEEECCCHHHHHHHHHHHCTTCEEEEEECCHHHHHHHHHHHTTTCCEEEEECCTGSCCCTSCCEEEEECSCCCCHHHHHHTCCTTEEEEEEECCTTCHHHHHHHHHHHHHCCCCCCHHHHHHHHHHHHHHHEEEEEEEEECCHHHHHHHHHHCHHHHTCCHHHHHHHHHTTCCEEEEEEEEEEEEC" );

	std::string t033_sam_stride = map_secstruct( "CCCHHHHHHHHHHHHHHTTTTTTTTTEETTTEEEETTTTEETTTTEEEECCCCTTCCTTTTCHHHHHHHHHHHHHHCHHHHHHHHHHHHHHTTTTCEEEECCTTHHHHHHHHHHHTTTTEEEEEECCHHHHHHHHHHHHTTCCEEEEEEGGGCCTTTTTEEEEEECTTTTCHHHHHHHHTTTEEEEEEECTTTTHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHHHEEEEEEEECCCHHHHHHHHHHCHHHHHCCHHHHHHHHHHCCTCEEEEEEEEEEEC" );


	// read the dssp ss for template pdb
	std::string t033_ss_from_1p91A;
	{
		std::string sequence, secstruct_1p91A;

		get_sequence_and_secstruct_from_dssp( "input/1p91A.dssp", sequence, secstruct_1p91A );

		id::SequenceMapping mapping;
		std::string source_seq, t033_seq;
		read_alignment_file( "input/alignment_for_trimming.txt", source_seq, t033_seq, mapping );
		assert( t033_seq.size() == t033_sam_stride.size() );

		// now copy over ss
		mapping.reverse(); // now goes from t033 to 1p91A
		t033_ss_from_1p91A = t033_seq;
		bool in_nterminal_loop( true );
		for ( Size i=1; i<= t033_seq.size(); ++i ) {
			char sec;
			if ( mapping[i] ) {
				sec = secstruct_1p91A[ mapping[i] - 1 ];
				in_nterminal_loop = false;
			} else {
				sec = 'L';
			}
			t033_ss_from_1p91A[ i-1 ] = sec;
		}
	} // scope

	SS_Map ss_map;
	ss_map[ "1p91A"      ] = t033_ss_from_1p91A;
	ss_map[ "psipred"    ] = t033_psipred;
	ss_map[ "sam_dssp"   ] = t033_sam_dssp;
	ss_map[ "sam_stride" ] = t033_sam_stride;
	return ss_map;
}


void
capri15_test()
{
	get_t033_ss();

	exit(0);


}

void*
my_main( void* )
{

 	std::string const mode( option[ OK::dna::specificity::mode ].value() );
	if ( mode == "loop" ) {
		capri_t033_loop_test();
		exit(0);
	} else if ( mode == "trim_dock" ) {
		capri_t033_trim_dock_test();
// 	} else if ( mode == "rhiju_trim_dock" ) {
// 		capri_t033_rhiju_trim_dock_test();
	} else if ( mode == "centroid_trim_dock" ) {
		capri_t033_centroid_trim_dock_test();
	} else if ( mode == "diversify_sam_loop" ) {
		diversify_sam_loop_test();
	} else if ( mode == "relax" ) {
		relax_test();
	} else if ( mode == "trim_dock_rebuild_relax" ) {
		trim_dock_rebuild_relax_test();
	} else if ( mode == "trim_dock_rebuild_relax_rhiju" ) {
		trim_dock_rebuild_relax_test_rhiju();
	} else if ( mode == "centroid_rescore" ) {
		centroid_rescore_test();
	} else if ( mode == "test" ) {
		capri15_test();
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
