// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilot/chrisk/rotamer_analysis.cc
/// @brief

//core library
#include <math.h>
#include <stdlib.h>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/chemical/util.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/types.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/types.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/sasa.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/constants.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/tools/make_map.hh>
#include <utility/string_util.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rtmin.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>

//protocols library (Movers)
#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

//calculator stuff
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>
#include <basic/MetricValue.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

//utilities

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rot_anl.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

//local options
namespace basic { namespace options { namespace OptionKeys {
}}}//basic::options::OptionKeys

////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace basic;
using namespace id;
using namespace pose;
using namespace pack;
using namespace conformation;
using namespace chemical;
using namespace kinematics;
using namespace scoring;
using namespace options;
using namespace basic::options::OptionKeys;
using namespace optimization;
namespace OK = OptionKeys;
using utility::vector1;
using utility::to_string;
using std::string;
using import_pose::pose_from_file;
using io::pdb::old_dump_pdb; // deprecated though
using namespace ObjexxFCL;
using basic::T;
using basic::Warning;
using basic::Error;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.chrisk.hbscan" );


//local options
namespace hbscan
{
basic::options::StringOptionKey tag( "hbscan:tag" );
basic::options::RealOptionKey repeat( "hbscan:repeat" );
basic::options::IntegerOptionKey nloop_hbscan( "hbscan:nloop_hbscan" );
basic::options::RealOptionKey hbond_sc_cutoff( "hbscan:hbond_sc_cutoff" );
basic::options::RealOptionKey fa_rep_cutoff( "hbscan:fa_rep_cutoff" );
basic::options::IntegerOptionKey lig_seqpos( "hbscan:lig_seqpos" );
basic::options::StringOptionKey restypes( "hbscan:restypes" );
basic::options::IntegerOptionKey lig_min_atom_sep( "hbscan:lig_min_atom_sep" );
basic::options::RealOptionKey hbond_wtd_energy_cutoff( "hbscan:hbond_wtd_energy_cutoff" );
basic::options::BooleanOptionKey write_cst_file( "hbscan:write_cst_file" );
}

//check for clashes based on atom distances
//only check atoms in AtomID vector
//should be way faster than calculating entire score
bool
fast_clash_check(
	Pose const & pose,
	vector1< id::AtomID > const check_atids,
	Real const clash_dist_cut
)
{
	Real const clash_dist2_cut( clash_dist_cut * clash_dist_cut );
	for ( Size iatid = 1; iatid <= check_atids.size(); ++iatid ) {
		Vector const at1_xyz( pose.xyz( check_atids[ iatid ] ) );
		for ( Size res2 = 1; res2 <= pose.total_residue(); ++res2 ) {
			for ( Size at2 = 1; at2 <= pose.residue( res2 ).natoms(); ++at2 ) {
				//skip virtual atoms!
				if ( pose.residue( res2 ).atom_type( at2 ).lj_wdepth() == 0.0 ) continue;
				id::AtomID atid2( at2, res2 );
				//skip if atid2 is in check_atids
				bool skip_at2( false );
				for ( Size jatid = 1; jatid <= check_atids.size(); ++jatid ) {
					if ( atid2 == check_atids[ jatid ] ) { skip_at2 = true; break; }
				}
				if ( skip_at2 ) continue;
				Real const dist2( at1_xyz.distance_squared( pose.xyz( atid2 ) ) );
				if ( dist2 < clash_dist2_cut ) {
					//TR << "CLASH!: " << check_atids[ iatid ] << " - " << atid2 <<
					//  " = " << dist2 << std::endl;
					return true;
				}
			}
		}
	}
	return false;
}

Size
max_atom_hbonds(
	chemical::ResidueType const & restype,
	Size const & atomno
) {
	chemical::AtomType const atomtype( restype.atom_type( atomno ) );
	if ( atomtype.is_polar_hydrogen() ) return Size( 1 );
	Size n_bonded_atoms( restype.number_bonded_hydrogens( atomno ) + restype.number_bonded_heavyatoms( atomno ) );
	switch( atomtype.hybridization() ){
	case chemical::SP2_HYBRID : return( Size( 3 - n_bonded_atoms ) );
	case chemical::SP3_HYBRID : return( Size( 4 - n_bonded_atoms ) );
	case chemical::RING_HYBRID : return( Size( 3 - n_bonded_atoms ) );
	case chemical::UNKNOWN_HYBRID : return( Size( 999 - n_bonded_atoms ) );
	default : return Size( 999 );
	}
}

//get # hbonds by type
// just looks bt 2 residues
void
get_n_hbonds_from_pair(
	Pose const pose,
	Size const rec_seqpos,
	Size const lig_seqpos,
	int const & lig_min_atom_sep,
	Real const & hbond_wtd_energy_cutoff,
	Size & n_sc_hbonds,
	Size & n_bb_hbonds,
	Size & n_hbonds,
	std::string & hbond_pairs_str
)
{
	using namespace core::scoring::hbonds;
	HBondSet hbset;
	hbset.setup_for_residue_pair_energies( pose, false, false );
	conformation::Residue rec_rsd( pose.residue( rec_seqpos ) );
	conformation::Residue lig_rsd( pose.residue( lig_seqpos ) );
	//count n hbonds for donors and acceptors
	vector1< Size > rec_atom_hbonds( rec_rsd.natoms(), 0 );
	vector1< Size > lig_atom_hbonds( lig_rsd.natoms(), 0 );
	for ( Size i = 1; i <= hbset.nhbonds(); ++i ) {
		HBond hb( hbset.hbond( i ) );
		//skip if hbond not good enough
		if ( hb.energy() * hb.weight() > hbond_wtd_energy_cutoff ) continue;
		if ( ( hb.don_res() == rec_seqpos && hb.acc_res() == lig_seqpos ) ||
				( hb.acc_res() == rec_seqpos && hb.don_res() == lig_seqpos ) ) {
			Size rec_atom( 0 );
			Size lig_atom( 0 );
			//if ligand is acceptor
			if ( hb.don_res() == rec_seqpos && hb.acc_res() == lig_seqpos ) {
				rec_atom = hb.don_hatm();
				lig_atom = hb.acc_atm();
			} else {
				rec_atom = hb.acc_atm();
				lig_atom = hb.don_hatm();
			}
			//atoms can only make so many simultaneous Hbonds
			if ( lig_atom_hbonds[ lig_atom ] >= max_atom_hbonds( lig_rsd.type(), lig_atom ) ) continue;
			if ( rec_atom_hbonds[ rec_atom ] >= max_atom_hbonds( rec_rsd.type(), rec_atom ) ) continue;
			//skip this hbond if this lig funxl group is already bonded
			bool skip_hbond( false );
			for ( Size iat = 1; iat <= lig_atom_hbonds.size(); ++iat ) {
				//    TR << "path_distance " + lig_rsd.atom_name( lig_atom ) + " -> " + lig_rsd.atom_name( iat ) + " = " << lig_rsd.type().path_distance( lig_atom, iat ) << "\n";
				if ( lig_atom_hbonds[ iat ] >= 1 &&
						lig_rsd.type().path_distance( lig_atom, iat ) < lig_min_atom_sep ) {
					//     TR << "skipping hbond " + rec_rsd.atom_name( rec_atom ) + "-" + lig_rsd.atom_name( lig_atom ) + "\n";
					skip_hbond = true;
					break;
				}
			}
			if ( skip_hbond ) continue;

			//then add to count if hbond not redundant
			lig_atom_hbonds[ lig_atom ] = ( lig_atom_hbonds[ lig_atom ] + 1 );
			rec_atom_hbonds[ rec_atom ] = ( rec_atom_hbonds[ rec_atom ] + 1 );
			hbond_pairs_str += ( " " + rec_rsd.atom_name( rec_atom ) + "-" + lig_rsd.atom_name( lig_atom ) );
			//if receptor res is backbone atom
			if ( hb.don_res() == rec_seqpos && hb.don_hatm_is_protein_backbone() ) ++n_bb_hbonds;
			else if ( hb.acc_res() == rec_seqpos && hb.acc_atm_is_protein_backbone() ) ++n_bb_hbonds;
			//else is sidechain
			else ++n_sc_hbonds;
		}
	}
	n_hbonds = n_sc_hbonds + n_bb_hbonds;
	/*
	//calc unsat hbonds a la BuriedUnsatisfied calculator
	for( Size atm = 1; atm <= rsd.nheavyatoms(); ++atm ){
	if( !( rsd.atom_type( atm ).is_acceptor() || rsd.atom_type( atm ).is_donor() ) ) continue;
	//we need at least one hbond
	if( n_atom_hbonds[ atm ] == 0 ){
	n_unsat += 1; continue;
	}
	//may need > 1 hbond, see BuriedUnsatCalc
	//behavior copied from Florian!
	std::string atom_type( rsd.type().atom_type( atm ).name() );
	Size satisfac_cut = 3;
	if( atom_type == "OH" || atom_type == "OCbb" || atom_type == "S" ){
	satisfac_cut = 2;
	}
	Size bonded_heavyatoms = rsd.n_bonded_neighbor_all_res( atm )
	- rsd.type().number_bonded_hydrogens( atm );
	if( bonded_heavyatoms + n_atom_hbonds[ atm ] < satisfac_cut ){
	n_unsat += 1;
	}
	}
	*/
}

/*
//get automorphic rmsd
Real
get_sc_automorphic_rmsd(
Pose pose,
Pose ref_pose,
Size const seqpos,
Size const ref_seqpos,
bool const super
)
{
//return bogus value if mismatch
if( ref_seqpos == 0 ||
pose.residue( seqpos ).name3().compare( ref_pose.residue( ref_seqpos ).name3() ) != 0 ||
pose.residue( seqpos ).natoms() != ref_pose.residue( ref_seqpos ).natoms()
) return 9999.;
// if( pose.residue( seqpos ).name3().compare( ref_pose.residue( ref_seqpos ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
//add variant type to virt bb atoms
core::pose::add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_BB, seqpos );
core::pose::add_variant_type_to_pose_residue( ref_pose, core::chemical::VIRTUAL_BB, ref_seqpos );

Residue rsd( pose.residue( seqpos ) );
Residue ref_rsd( ref_pose.residue( ref_seqpos ) );

Real sc_rmsd( scoring::automorphic_rmsd( pose.residue( seqpos ), ref_pose.residue( ref_seqpos ), super ) );
return sc_rmsd;
}
*/

void
print_cst_file(
	std::string const filename,
	std::string const & r1, char const & r2,
	std::string const & r1a1, std::string const & r1a2, std::string const & r1a3,
	std::string const & r2a1, std::string const & r2a2, std::string const & r2a3,
	Real const & distab,
	Real const & anga, Real const & angb,
	Real const & tora, Real const & torab, Real const & torb
) {
	utility::io::ozstream out( filename, std::ios::out );

	out << "CST::BEGIN\n\n";
	out << "\tTEMPLATE::   ATOM_MAP: 1 atom_name: " << r1a1 << " " << r1a2 << " " << r1a3 << "\n";
	out << "\tTEMPLATE::   ATOM_MAP: 1 residue3: " << r1 << "\n\n";
	out << "\tTEMPLATE::   ATOM_MAP: 2 atom_name: " << r2a1 << " " << r2a2 << " " << r2a3 << "\n";
	out << "\tTEMPLATE::   ATOM_MAP: 2 residue1: " << r2 << "\n";
	out << "\n\n";
	out << "\tCONSTRAINT:: distanceAB:\t" << distab << "\t0.20\t100.00\t0\t1\n";
	out << "\tCONSTRAINT:: angle_A:\t" << anga << "\t10.00\t100.00\t360.00\t2\n";
	out << "\tCONSTRAINT:: angle_B:\t" << angb << "\t10.00\t100.00\t360.00\t2\n";
	out << "\tCONSTRAINT:: torsion_A:\t" << tora << "\t10.00\t50.00\t360.00\t1\n";
	out << "\tCONSTRAINT:: torsion_AB:\t" << torab << "\t10.00\t50.00\t360.00\t1\n";
	out << "\tCONSTRAINT:: torsion_B:\t" << torb << "\t10.00\t50.00\t360.00\t1\n";
	out << "\n";
	out << "\tALGORITHM_INFO:: match\n";
	out << "\t\tCHI_STRATEGY:: CHI 1 EX_TWO_HALF_STEP_STDDEVS\n";
	out << "\t\tCHI_STRATEGY:: CHI 2 EX_TWO_HALF_STEP_STDDEVS\n";
	out << "\tALGORITHM_INFO::END\n";
	out << "\n";
	out << "CST::BEGIN\n\n";
	out << std::endl;
	out.close();
}


core::scoring::hbonds::HBondDatabaseCOP hb_database_;

void
scan_hbond_jumps(
	pose::Pose pose,
	scoring::ScoreFunctionOP scorefxn,
	Size seqpos,
	Size atomno,
	conformation::Residue new_rsd,
	Size new_atomno
)
{
	using namespace id;
	using namespace conformation;
	using namespace scoring;
	using namespace scoring::hbonds;
	using namespace kinematics;
	using numeric::constants::f::pi;
	using numeric::constants::f::pi_2;

	Residue rsd( pose.residue( seqpos ) );
	Pose ref_pose( pose );
	vector1< std::string > filename_in_split( utility::string_split( core::pose::tag_from_pose( pose ), '.' ) );
	filename_in_split.pop_back();
	std::string const pdbname_out_tag( utility::join( filename_in_split, "." ) );

	//append by jump from seqpos atomno to new_rsd atom 1
	//  pose.append_residue_by_jump( new_rsd, seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), true );
	pose.append_residue_by_jump( new_rsd, seqpos );
	Size new_seqpos( pose.total_residue() );
	Size jump_number( pose.fold_tree().num_jump() );
	Jump jump( pose.jump( jump_number ) );

	//give ideal conformation to new residue by repacking only with dunbrack term
	if ( pose.residue( new_seqpos ).is_protein() ) {
		pose.set_phi( new_seqpos, -57.0 );
		pose.set_psi( new_seqpos, -47.0 );
		pose.set_omega( new_seqpos, 180.0 );
		pack::task::TaskFactoryOP task_factory( new pack::task::TaskFactory );
		pack::task::operation::RestrictResidueToRepackingOP restrict_to_repack_taskop( new pack::task::operation::RestrictResidueToRepacking() );
		pack::task::operation::PreventRepackingOP prevent_repack_taskop( new pack::task::operation::PreventRepacking() );
		//prevent repacking all original seqpos and repack new rsd
		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( i == new_seqpos ) {
				restrict_to_repack_taskop->include_residue( i );
			} else {
				prevent_repack_taskop->include_residue( i );
			}
		}
		task_factory->push_back( restrict_to_repack_taskop );
		task_factory->push_back( prevent_repack_taskop );
		pack::task::PackerTaskOP task( task_factory->create_task_and_apply_taskoperations( pose ));
		core::scoring::ScoreFunctionOP fadun_scorefxn( new core::scoring::ScoreFunction() );
		fadun_scorefxn->set_weight( fa_dun, 1.0 );
		protocols::simple_moves::PackRotamersMoverOP pack( new protocols::simple_moves::PackRotamersMover( fadun_scorefxn, task, 1 ) );
		pack->apply( pose );
	}
	//\end give ideal conformation...

	//store  atom ids for clash check
	vector1< id::AtomID > clash_check_atids;
	for ( Size iat = 1; iat <= new_rsd.natoms(); ++iat ) {
		clash_check_atids.push_back( id::AtomID( iat, new_seqpos ) );
	}

	//which is acceptor, donor?
	Size aatm( 0 ), acc_pos( 0 ), hatm( 0 ), don_pos( 0 );
	bool new_rsd_is_acc( false );
	// is acceptor
	if ( rsd.atom_type( atomno ).is_polar_hydrogen() &&
			new_rsd.atom_type( new_atomno ).is_acceptor() ) {
		don_pos = seqpos;
		hatm = atomno;
		acc_pos = new_seqpos;
		aatm = new_atomno;
		new_rsd_is_acc = true;
	} else if ( new_rsd.atom_type( new_atomno ).is_polar_hydrogen() &&
			//or  is donor
			rsd.atom_type( atomno ).is_acceptor() ) {
		don_pos = new_seqpos;
		hatm = new_atomno;
		acc_pos = seqpos;
		aatm = atomno;
	} else { utility_exit_with_message( "ERROR: res " + to_string( seqpos ) + " atom " + to_string( atomno ) +
		" res " + to_string( new_seqpos ) + " atom " + to_string( new_atomno ) + " is not HB don/acc pair!!\n" );
	}

	//if less than 4 atoms in res, add vrt res so final torsion exists
	Size n_new_rsd( 1 );
	if ( new_rsd.natoms() < 4 ) {
		chemical::ResidueTypeSetCOP rsd_set( rsd.residue_type_set() );
		conformation::ResidueOP vrt_rsd( conformation::ResidueFactory::create_residue( rsd_set->name_map( "VRT" ) ) );
		pose.append_residue_by_jump( *vrt_rsd, pose.total_residue() );
		n_new_rsd += 1;
	}

	FoldTree f_jump( pose.fold_tree() );
	//just min the new jump
	//MoveMapOP mm = new MoveMap;
	//mm->set_jump( jump_number, true );
	//protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm, scorefxn, "dfpmin", 0.01, true );

	//new naive fold tree
	FoldTree f_rot( pose.total_residue() );
	//switch to chem bond so can use bond angle defs directly
	f_rot.new_chemical_bond( seqpos, new_seqpos, rsd.atom_name( atomno ), new_rsd.atom_name( new_atomno ), pose.total_residue() - n_new_rsd );
	pose.fold_tree( f_rot );

	//now get their base atoms to get datm and batm
	Size datm( pose.residue( don_pos ).atom_base( hatm ) );
	Size batm( pose.residue( acc_pos ).atom_base( aatm ) );
	//  Size b2atm( pose.residue( acc_pos ).abase2( aatm ) );
	Size bbatm( pose.residue( acc_pos ).atom_base( batm ) ); //base of base needed for defined torsion DOF, is sometimes diff than abase2
	Size dbatm( pose.residue( don_pos ).atom_base( datm ) ); //hpol base2

	Size hb_states_tot( 0 );
	Size hb_states_good( 0 );
	hbonds::HBEvalTuple hbe_type( datm, pose.residue( don_pos ), aatm, pose.residue( acc_pos ) );
	//granularity of enumeration
	//must import from hbonds/constants.hh
	// in case you're wondering, we're sampling in angle space instead of cos(angle)
	// space so we sample more equally in polar coord space (not more densely when cos(pi-angle)->1)
	// NOTE: these min/max angle vals have offsets to avoid 0,pi boundary issues that make torsions ambiguous
	Real AHdist_min(MIN_R), AHdist_max(MAX_R); Size AHdist_steps( 6 );
	Real BAHang_min( pi - std::acos( MIN_xH ) + 0.01 ), BAHang_max( pi - std::acos( MAX_xH ) - 0.01 ); Size BAHang_steps( 5 );
	Real AHDang_min( pi - std::acos( MIN_xD ) + 0.01 ), AHDang_max( pi - std::acos( MAX_xD ) - 0.01 ); Size AHDang_steps( 5 );
	Real other_tor_min( 0.01 ), other_tor_max( pi_2 - 0.01 ); Size other_tor_steps( 10 );
	Real BAHDtor_min( 0.01 ), BAHDtor_max( pi_2 - 0.01 ); Size BAHDtor_steps( 10 );
	//  Real cosBAH_min( MIN_xH ), cosBAH_max( MAX_xH ); Size cosBAH_steps( 5 );
	//  Real cosAHD_min( MIN_xD ), cosAHD_max( MAX_xD ); Size cosAHD_steps( 5 );
	//  Real other_tor_min( 0 ), other_tor_max( pi_2 ); Size other_tor_steps( 1 );
	//gimbal lock… <sigh>
	for ( Real AHdist = AHdist_min; AHdist <= AHdist_max;
			AHdist += (AHdist_max - AHdist_min)/static_cast<Real>(AHdist_steps-1) ) {
		for ( Real AHDang = AHDang_min - .0001; AHDang <= AHDang_max;
				AHDang += (AHDang_max - AHDang_min)/static_cast<Real>(AHDang_steps-1) ) {
			for ( Real BAHang = BAHang_min - 0.0001; BAHang <= BAHang_max;
					BAHang += (BAHang_max - BAHang_min)/static_cast<Real>(BAHang_steps-1) ) {
				//this loop should be longer than the  internal orientation loops
				for ( Real BAHDtor = BAHDtor_min; BAHDtor < BAHDtor_max;
						BAHDtor += (BAHDtor_max - BAHDtor_min)/static_cast<Real>(BAHDtor_steps-1) ) {
					for ( Real other_tor = other_tor_min; other_tor < other_tor_max;
							other_tor += (other_tor_max - other_tor_min)/static_cast<Real>(other_tor_steps-1) ) {

						//      Real cosBAH( std::cos( pi - BAHang ) );
						//      Real cosAHD( std::cos( pi - AHDang ) );
						//call hbonds::hbond_compute_energy( ... ) with these angles
						// and skip if is zero hb energy
						// this allows computing exactly what frac of actual HB phase space
						// can be filled w/ a favorable  molecule
						//      Real hb_energy( 0.0 );
						//      scoring::hbonds::hbond_compute_energy( *( hb_database_ ),
						//       scorefxn->energy_method_options().hbond_options(),
						//       hbe_type, AHdist, cosAHD, cosBAH, other_tor, hb_energy );
						//      TR << hb_energy << std::endl;
						//      if( hb_energy >= 0.0 ) continue; //was not actually an hbond

						++hb_states_tot;

						//reset chem  bond ftree
						pose.fold_tree( f_rot );

						pose.conformation().set_bond_angle( AtomID( batm, acc_pos ),
							AtomID( aatm, acc_pos  ),
							AtomID( hatm, don_pos  ),
							BAHang );

						pose.conformation().set_bond_angle( AtomID( aatm, acc_pos ),
							AtomID( hatm, don_pos  ),
							AtomID( datm, don_pos  ),
							AHDang );

						pose.conformation().set_torsion_angle( AtomID( batm, acc_pos ),
							AtomID( aatm, acc_pos  ),
							AtomID( hatm, don_pos  ),
							AtomID( datm, don_pos  ),
							BAHDtor );
						if ( new_rsd_is_acc ) {
							pose.conformation().set_torsion_angle( AtomID( bbatm, acc_pos ),
								AtomID( batm, acc_pos  ),
								AtomID( aatm, acc_pos  ),
								AtomID( hatm, don_pos  ),
								other_tor );
							/*
							pose.conformation().set_torsion_angle( AtomID( b2atm, acc_pos ),
							AtomID( batm, acc_pos  ),
							AtomID( aatm, acc_pos  ),
							AtomID( hatm, don_pos  ),
							other_tor );
							*/
						} else {
							pose.conformation().set_torsion_angle( AtomID( dbatm, don_pos ),
								AtomID( datm, don_pos  ),
								AtomID( hatm, don_pos  ),
								AtomID( aatm, acc_pos  ),
								other_tor );
						}

						pose.conformation().set_bond_length( AtomID( aatm, acc_pos ),
							AtomID( hatm, don_pos  ),
							AHdist );

						//do fast clash check, OH hbonds are only 0.8A!
						if ( fast_clash_check( pose, clash_check_atids, 0.8 ) ) continue;

						pose.fold_tree( f_jump );
						//minimize the new jump, too slow!
						//min_mover->apply( pose );
						scorefxn->score( pose );

						//get score
						Real new_rsd_score( pose.energies().residue_total_energies( new_seqpos ).dot( scorefxn->weights() ) );
						Real new_rsd_hbond_sc( pose.energies().residue_total_energies( new_seqpos )[ hbond_sc ] + pose.energies().residue_total_energies( new_seqpos )[ hbond_bb_sc ]  );
						Real new_rsd_fa_rep( pose.energies().residue_total_energies( new_seqpos )[ fa_rep ] );
						Size n_sc_hbonds( 0 ), n_bb_hbonds( 0 ), n_hbonds( 0 );
						int const lig_min_atom_sep( option[ hbscan::lig_min_atom_sep ] );
						Real const hbond_wtd_energy_cutoff( option[ hbscan::hbond_wtd_energy_cutoff ] );
						std::string hbond_pairs_str;
						get_n_hbonds_from_pair( pose, new_seqpos, seqpos, lig_min_atom_sep, hbond_wtd_energy_cutoff, n_sc_hbonds, n_bb_hbonds, n_hbonds, hbond_pairs_str );

						// skip if making less than 2 sidechain hbonds
						if ( n_sc_hbonds >= 2 && new_rsd_hbond_sc <= option[ hbscan::hbond_sc_cutoff ] && new_rsd_fa_rep <= option[ hbscan::fa_rep_cutoff ] ) {
							Real AHDangd( numeric::conversions::degrees( AHDang ) ), BAHangd( numeric::conversions::degrees( BAHang ) ),
								BAHDtord( numeric::conversions::degrees( BAHDtor ) ), other_tord( numeric::conversions::degrees( other_tor ) );
							std::string pdbname( pdbname_out_tag + ".hbscan_" + rsd.name1() + to_string( seqpos ) + utility::trim( rsd.atom_name( atomno ) ) + "-"
								+ new_rsd.name1() + to_string( new_seqpos ) + utility::trim( new_rsd.atom_name( new_atomno ) ) + "_"
								+ to_string( AHdist ) + "-" + to_string( Size( AHDangd ) ) + "-" + to_string( Size( BAHangd ) ) + "-"
								+ to_string( Size( BAHDtord ) ) + "-" + to_string( Size( other_tord ) ) + ".pdb" );
							TR << pdbname << " ";
							TR << hbond_pairs_str << " ";
							TR << "sc_hbonds: " << n_sc_hbonds << " ";
							TR << "new_rsd_score: " << new_rsd_score << " ";
							TR << "new_rsd_hbond_sc: " << new_rsd_hbond_sc << " ";
							TR << "new_rsd_fa_rep: " << new_rsd_fa_rep << " ";
							TR << "AHdist: " << AHdist << " AHDang: " << AHDangd << " BAHang: " << BAHangd << " BAHDtor: " << BAHDtord << " other_tor: " << other_tord << " ";
							TR << pose.energies().total_energies().weighted_string_of( scorefxn->weights() )
								+ " total_score: " + to_string( pose.energies().total_energies()[ total_score ] );
							TR << std::endl;
							++hb_states_good;
							pose.dump_pdb( pdbname );
							//now print cst file
							if ( option[ hbscan::write_cst_file ] ) {
								if ( new_rsd_is_acc ) { //if ligand is donor
									//must ask for first (unchanging) torsion
									Real tora( numeric::conversions::degrees( pose.conformation().torsion_angle( AtomID( dbatm, don_pos ),
										AtomID( datm, don_pos ),
										AtomID( hatm, don_pos ),
										AtomID( aatm, acc_pos ) ) ) );
									print_cst_file( pdbname + ".cst",
										rsd.name3(), new_rsd.name1(),
										rsd.atom_name( hatm ), rsd.atom_name( datm ), rsd.atom_name( dbatm ),
										new_rsd.atom_name( aatm ), new_rsd.atom_name( batm ), new_rsd.atom_name( bbatm ),
										AHdist,
										AHDangd, BAHangd,
										tora, BAHDtord, other_tord );
								} else { //else if ligand is acceptor
									//must ask for first (unchanging) torsion
									Real tora( numeric::conversions::degrees( pose.conformation().torsion_angle( AtomID( bbatm, acc_pos ),
										AtomID( batm, acc_pos ),
										AtomID( aatm, acc_pos ),
										AtomID( hatm, don_pos ) ) ) );
									print_cst_file( pdbname + ".cst",
										rsd.name3(), new_rsd.name1(),
										rsd.atom_name( aatm ), rsd.atom_name( batm ), rsd.atom_name( bbatm ),
										new_rsd.atom_name( hatm ), new_rsd.atom_name( datm ), new_rsd.atom_name( dbatm ),
										AHdist,
										BAHangd, AHDangd,
										tora, BAHDtord, other_tord );
								}
							}
						}
					}
				}
			}
		}
	}
	// return ( static_cast< Real >( hb_states_good ) /
	//  static_cast< Real >( hb_states_tot ) );
}



///////////////////////////////////////////////////////////////////////////////
void
byres_analysis(
	pose::Pose pose
)
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//tag is "." or "tag."
	std::string tag( option[ hbscan::tag ] );
	std::string const scname( tag + "sc" );

	std::string pdbname( core::pose::tag_from_pose( pose ) );

	//create a ScoreFunction from commandline options
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	//and add decompose bb_hbond energies option
	core::scoring::methods::EnergyMethodOptionsOP emopts(
		new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() )
	);
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
	scorefxn->set_energy_method_options( *emopts );


	//silent-type output
	io::silent::SilentFileData pose_silent_data;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//use resfile to note residues to analyze...
	Pose start_pose( pose );

	Size nres( pose.total_residue() );

	io::silent::SilentStructOP pose_ss( new io::silent::ScoreFileSilentStruct );
	core::scoring::hbonds::HBondSet hbset;
	// must score pose to update nbrlist first!
	scorefxn->score( pose );
	hbset.setup_for_residue_pair_energies( pose, false, false );

	chemical::ResidueTypeSetCOP rsd_set( pose.residue( 1 ).residue_type_set() );

	//user defined residue type to scan with
	std::string const restypes( option[ hbscan::restypes ] );
	for ( Size iaa = 0; iaa <= restypes.length() - 1; ++iaa ) {
		conformation::ResidueOP dock_rsd( conformation::ResidueFactory::create_residue(
			*( rsd_set->get_representative_type_aa( chemical::AA( chemical::aa_from_oneletter_code( restypes[ iaa ] ) ) ) ) ) );
		for ( Size seqpos = 1; seqpos <= nres; ++seqpos ) {
			//if option is set, only check the user defined sequence position
			if ( option[ hbscan::lig_seqpos ].user() && seqpos != Size( option[ hbscan::lig_seqpos ] ) ) continue;

			io::silent::SilentStructOP res_ss( new io::silent::ScoreFileSilentStruct );
			res_ss->decoy_tag( pdbname );

			//print byatom data
			io::silent::SilentStructOP atom_ss( new io::silent::ScoreFileSilentStruct );
			Residue const rsd( pose.residue( seqpos ) );

			for ( Size iatom = 1; iatom <= rsd.natoms(); ++iatom ) {
				for ( Size jatom = 1; jatom <= dock_rsd->natoms(); ++jatom ) {
					// skip if this is a backbone atoms
					if ( dock_rsd->atom_is_backbone( jatom ) ) continue;
					//skip if not an acceptor--polarH pair
					if ( !( (  rsd.atom_type( iatom ).is_acceptor() && dock_rsd->atom_type( jatom ).is_polar_hydrogen() )
							|| ( rsd.atom_type( iatom).is_polar_hydrogen() && dock_rsd->atom_type( jatom ).is_acceptor() ) ) ) continue;

					hb_database_ = core::scoring::hbonds::HBondDatabase::get_database( "sp2_elec_params" );
					scan_hbond_jumps( pose, scorefxn, seqpos, iatom, *dock_rsd, jatom );

					std::string atom_name( rsd.atom_name( iatom ) );
					atom_name.erase( 0, atom_name.find_first_not_of( " " ) );
					atom_name.erase( atom_name.find_last_not_of( " " ) + 1 );
					atom_ss->decoy_tag( rsd.name3() + "_" + string_of( seqpos ) +
						"_" + atom_name );
					atom_ss->add_string_value( "atom_type", rsd.atom_type( iatom ).name() );
					atom_ss->add_energy( "seqpos", static_cast< Real >( seqpos ) );
					atom_ss->add_string_value( "pdbname", pdbname );
					atom_ss->add_string_value( "pdb_chain", utility::string_split( pose.pdb_info()->pose2pdb( seqpos ) )[ 1 ] );
					atom_ss->add_string_value( "pdb_resnum", utility::string_split( pose.pdb_info()->pose2pdb( seqpos ) )[ 2 ] );
				}
			}
		}
	}
}

void
load_coords()
{
	using namespace core;

	import_pose::pose_stream::MetaPoseInputStream input =  import_pose::pose_stream::streams_from_cmd_line();
	core::chemical::ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	core::pose::Pose pose;
	while ( input.has_another_pose() ) {
		input.fill_pose( pose, *rsd_set );
		for ( Size i = 1; i<= option[ hbscan::repeat ]; ++i ) {
			byres_analysis( pose );
		}
	}
}

void*
my_main( void*)
{
	load_coords();
	exit(0);

}

int
main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace utility;

		option.add( hbscan::write_cst_file, "write a cst file for matching hbond docks" ).def( true );
		option.add( hbscan::tag, "nametag" ).def( "score" );
		option.add( hbscan::repeat, "rescore structure" ).def( 1 );
		option.add( hbscan::nloop_hbscan, "number of iter of hbscan per polar h or acc" ).def( 10 );
		option.add( hbscan::hbond_sc_cutoff, "max hbond energy" ).def( -0.0 );
		option.add( hbscan::fa_rep_cutoff, "max clash energy" ).def( 0.1 );
		option.add( hbscan::lig_seqpos, "sequence position of residue/ligand" ).def( 1 );
		option.add( hbscan::restypes, "aa1 string of scan res types" ).def( "DEKNQRSTWY" );
		option.add( hbscan::lig_min_atom_sep, "num atoms between ligand residue interactions" ).def( 3 );
		option.add( hbscan::hbond_wtd_energy_cutoff, "energy cutoff per hbond" ).def( -0.03 );

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


