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
 #include <protocols/frags/VallData.hh>
 #include <protocols/frags/TorsionFragment.hh>

 #include <core/scoring/dna/setup.hh>
 #include <core/scoring/dna/base_geometry.hh>
 #include <core/scoring/dna/BasePartner.hh>
 #include <core/scoring/GenBornPotential.hh>
 #include <core/scoring/LREnergyContainer.hh>
 #include <core/scoring/methods/Methods.hh>

 #include <protocols/simple_moves/BackboneMover.hh>
 #include <protocols/simple_moves/MinMover.hh>
 #include <protocols/moves/MonteCarlo.hh>
 #include <protocols/moves/Mover.hh>
 #include <protocols/moves/MoverContainer.hh>
 #include <protocols/moves/OutputMovers.hh>
 #include <protocols/rigid/RigidBodyMover.hh>
 // #include <protocols/moves/rigid_body_moves.hh>
 #include <protocols/moves/TrialMover.hh>
 #include <protocols/simple_moves/PackRotamersMover.hh>
 #include <protocols/simple_moves/RotamerTrialsMover.hh>
 #include <protocols/moves/RepeatMover.hh>

 #include <protocols/loops/ccd_closure.hh>
 #include <protocols/loops/loops_main.hh>

 #include <protocols/viewer/viewers.hh>

 #include <core/types.hh>

 #include <core/scoring/sasa.hh>

// #include <basic/prof.hh> // profiling
// #include <basic/CacheableData.hh> // profiling

 #include <core/id/SequenceMapping.hh>

 #include <core/chemical/AtomTypeSet.hh>
 #include <core/chemical/MMAtomTypeSet.hh>

 #include <core/chemical/AA.hh>
 #include <core/conformation/Residue.hh>
 #include <core/conformation/ResidueMatcher.hh>
 #include <core/pack/rotamer_set/RotamerCouplings.hh>
 #include <core/chemical/ResidueTypeSet.hh>
 #include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
 #include <core/chemical/VariantType.hh>

 #include <core/chemical/ChemicalManager.hh>

 #include <core/scoring/etable/Etable.hh>
 #include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
 #include <core/scoring/Ramachandran.hh>
 #include <core/pack/dunbrack/RotamerLibrary.hh>
 #include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
 #include <core/scoring/hbonds/HBondSet.hh>
 #include <core/scoring/hbonds/hbonds.hh>
 #include <core/scoring/etable/count_pair/CountPairFunction.hh>

 #include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/TaskOperation.hh>

 #include <core/kinematics/FoldTree.hh>
 #include <protocols/viewer/visualize.hh>
#include <core/kinematics/MoveMap.hh>
 #include <core/kinematics/util.hh>
 #include <core/id/AtomID_Map.hh>
 //#include <core/id/AtomID_Map.Pose.hh>

 #include <core/mm/MMTorsionLibrary.hh>
 #include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <basic/options/util.hh>//option.hh>
// #include <basic/options/after_opts.hh>

 #include <basic/basic.hh>

 #include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>

 #include <numeric/xyzVector.hh>
 #include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// //REMOVE LATER!
// #include <utility/io/izstream.hh>


// // C++ headers
 #include <cstdlib>
 #include <fstream>
 #include <iostream>
 #include <string>

//silly using/typedef
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pep_spec.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
//Auto Headers
#include <core/pose/util.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


 using basic::T;
 using basic::Error;
 using basic::Warning;


using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
//using namespace protocols;

using utility::vector1;
using std::string;


using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace pose;
using namespace chemical;
using namespace conformation;
using namespace scoring;
using namespace optimization;
using namespace kinematics;
using namespace protocols::moves;
using namespace id;
using namespace protocols::frags;
using utility::file::FileName;

void
dump_efactor_pdb(
	pose::Pose & pose,
	core::scoring::ScoreFunctionOP const scorefxn,
	std::string const & tag
)
{
//	id::AtomID_Mask const & mask;
//	id::initialize( mask, pose );

	( *scorefxn )( pose );

	char const *filename = tag.c_str();
	std::fstream out( filename, std::ios::out );


	Size number(0);

	static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );

	out << "MODEL     " << tag << "\n";
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		Real residue_total_energy( pose.energies().residue_total_energies( i ).dot( scorefxn->weights() ) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );

 //			if ( ! mask[ id::AtomID( j,i ) ] ) continue;

			//skip outputing virtual atom unless specified
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::output_virtual ]() &&
				(int)atom.type() == (int)rsd.atom_type_set().n_atomtypes() ) continue;

			++number;
			assert( rsd.chain() < int(chains.size()) ); // silly restriction
			char const chain( chains[ rsd.chain() ] );
			out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
				rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
				F(8,3,atom.xyz()(1)) <<
				F(8,3,atom.xyz()(2)) <<
				F(8,3,atom.xyz()(3)) <<
				F(6,2,1.0) << F(6,2, residue_total_energy ) << '\n';
		}
	}
	out << "ENDMDL\n";
}

Real
get_rmsd(
	pose::Pose pose,
	Size prot_begin,
	Size prot_end,
	Size pep_begin,
	Size pep_anchor,
	Size pep_end,
	pose::Pose ref_pose,
	Size ref_prot_begin,
	Size ref_prot_end,
	Size ref_pep_begin,
	Size ref_pep_anchor,
	Size ref_pep_end
)
{

	Size ref_prot_anchor( ref_prot_begin );

	Size pep_nterm( pep_anchor - pep_begin );
	Size ref_pep_nterm( ref_pep_anchor - ref_pep_begin );
	Size pep_cterm( pep_end - pep_anchor );
	Size ref_pep_cterm( ref_pep_end - ref_pep_anchor );

	Size nterm( pep_nterm );
	Size cterm( pep_cterm );
	if( pep_nterm < ref_pep_nterm ) ref_pep_begin += ref_pep_nterm - pep_nterm;
	else if( pep_nterm > ref_pep_nterm ){
		pep_begin += pep_nterm - ref_pep_nterm;
		nterm = ref_pep_nterm;
	}
	if( pep_cterm < ref_pep_cterm ) ref_pep_end -= ref_pep_cterm - pep_cterm;
	else if( pep_cterm > ref_pep_cterm ){
		pep_end -= pep_cterm - ref_pep_cterm;
		cterm = ref_pep_cterm;
	}
	//superpose if needed
	if( option[ pep_spec::ref_align ] ){
		id::AtomID_Map< id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, id::BOGUS_ATOM_ID );

		for ( Size i = prot_begin; i <= prot_end; ++i ) {
			id::AtomID const id1( pose.residue( i ).atom_index( "CA" ), i );
			id::AtomID const id2( ref_pose.residue( i + static_cast< int >( ref_prot_begin ) - static_cast< int >( prot_begin ) ).atom_index( "CA" ), i + static_cast< int >( ref_prot_begin ) - static_cast< int >( prot_begin ) );
			atom_map[ id1 ] = id2;
		}
		superimpose_pose( pose, ref_pose, atom_map );

	}
	Real total_rmsd( 0 );
	for( int i_seq = 0; i_seq <= nterm + cterm; ++i_seq ){
		Size ref_pep_seqpos = ref_pep_begin + i_seq;
		Size pep_seqpos = pep_begin + i_seq;
		Real sd( ref_pose.residue( ref_pep_seqpos ).xyz( "CA" ).distance_squared( pose.residue( pep_seqpos ).xyz( "CA" ) ) );
		total_rmsd += sd;
	}
	total_rmsd = std::sqrt( total_rmsd / ( nterm+cterm+1 ) );
	return total_rmsd;

}

std::string
pep_rmsd_analysis(
	pose::Pose pose,
	Size prot_begin,
	Size prot_end,
	Size pep_begin,
	Size pep_anchor,
	Size pep_end
)
{

	pose::Pose ref_pose;
	std::string ref_name( option[ pep_spec::ref_pose ] );
	core::import_pose::pose_from_file( ref_pose, ref_name , core::import_pose::PDB_file);

	Size ref_pep_anchor_in( option[ pep_spec::pep_anchor ] );
	if( option[ pep_spec::ref_pep_anchor ].user() ) ref_pep_anchor_in = option[ pep_spec::ref_pep_anchor ];
	std::string ref_pep_chain_in( option[ pep_spec::pep_chain ] );
	if( option[ pep_spec::ref_pep_chain ].user() ) ref_pep_chain_in = option[ pep_spec::ref_pep_chain ];
	Size ref_pep_anchor( ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_in ) );
	Size ref_pep_chain( ref_pose.chain( ref_pep_anchor ) );
	Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pep_chain ) );
	Size ref_pep_end( ref_pose.conformation().chain_end( ref_pep_chain ) );

	Size ref_prot_chain;
	for( Size i = 1; i <= ref_pose.conformation().num_chains(); ++i ){
		if( i == ref_pep_chain ) continue;
		if( !( ref_pose.residue( ref_pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			ref_prot_chain = i;
			break;
		}
	}
	Size ref_prot_begin( ref_pose.conformation().chain_begin( ref_prot_chain ) );
	Size ref_prot_end( ref_pose.conformation().chain_end( ref_prot_chain ) );
	Size ref_prot_anchor( ref_prot_begin );

	Size pep_nterm( pep_anchor - pep_begin );
	Size ref_pep_nterm( ref_pep_anchor - ref_pep_begin );
	Size pep_cterm( pep_end - pep_anchor );
	Size ref_pep_cterm( ref_pep_end - ref_pep_anchor );

	Size nterm( pep_nterm );
	Size cterm( pep_cterm );
	if( pep_nterm < ref_pep_nterm ) ref_pep_begin += ref_pep_nterm - pep_nterm;
	else if( pep_nterm > ref_pep_nterm ){
		pep_begin += pep_nterm - ref_pep_nterm;
		nterm = ref_pep_nterm;
	}
	if( pep_cterm < ref_pep_cterm ) ref_pep_end -= ref_pep_cterm - pep_cterm;
	else if( pep_cterm > ref_pep_cterm ){
		pep_end -= pep_cterm - ref_pep_cterm;
		cterm = ref_pep_cterm;
	}
	//superpose if needed
	if( option[ pep_spec::ref_align ] ){
		id::AtomID_Map< id::AtomID > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, id::BOGUS_ATOM_ID );

		for ( Size i = prot_begin; i <= prot_end; ++i ) {
			id::AtomID const id1( pose.residue( i ).atom_index( "CA" ), i );
			id::AtomID const id2( ref_pose.residue( i + static_cast< int >( ref_prot_begin ) - static_cast< int >( prot_begin ) ).atom_index( "CA" ), i + static_cast< int >( ref_prot_begin ) - static_cast< int >( prot_begin ) );
			atom_map[ id1 ] = id2;
		}
		superimpose_pose( pose, ref_pose, atom_map );

	}
	Real total_rmsd( 0 );
	std::string rmsd_analysis( "" );
	for( int i_seq = 0; i_seq <= nterm + cterm; ++i_seq ){
		Size ref_pep_seqpos = ref_pep_begin + i_seq;
		Size pep_seqpos = pep_begin + i_seq;
		Real sd( ref_pose.residue( ref_pep_seqpos ).xyz( "CA" ).distance_squared( pose.residue( pep_seqpos ).xyz( "CA" ) ) );
		total_rmsd += sd;
		rmsd_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_rmsd:\t" + string_of( std::sqrt( sd ) ) + "\t";
	}
	total_rmsd = std::sqrt( total_rmsd / ( nterm+cterm+1 ) );
	rmsd_analysis += "total_rmsd:\t" + string_of( total_rmsd ) + "\t";
	return rmsd_analysis;

}

std::string
pep_phipsi_analysis(
	pose::Pose pose,
	Size prot_begin,
	Size prot_end,
	Size pep_begin,
	Size pep_anchor,
	Size pep_end
)
{

	pose::Pose ref_pose;
	std::string ref_name( option[ pep_spec::ref_pose ] );
	core::import_pose::pose_from_file( ref_pose, ref_name , core::import_pose::PDB_file);

	Size ref_pep_anchor_in( option[ pep_spec::pep_anchor ] );
	if( option[ pep_spec::ref_pep_anchor ].user() ) ref_pep_anchor_in = option[ pep_spec::ref_pep_anchor ];
	std::string ref_pep_chain_in( option[ pep_spec::pep_chain ] );
	if( option[ pep_spec::ref_pep_chain ].user() ) ref_pep_chain_in = option[ pep_spec::ref_pep_chain ];
	Size ref_pep_anchor( ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_in ) );
	Size ref_pep_chain( ref_pose.chain( ref_pep_anchor ) );
	Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pep_chain ) );
	Size ref_pep_end( ref_pose.conformation().chain_end( ref_pep_chain ) );

	Size ref_prot_chain;
	for( Size i = 1; i <= ref_pose.conformation().num_chains(); ++i ){
		if( i == ref_pep_chain ) continue;
		if( !( ref_pose.residue( ref_pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			ref_prot_chain = i;
			break;
		}
	}
	Size ref_prot_begin( ref_pose.conformation().chain_begin( ref_prot_chain ) );
	Size ref_prot_end( ref_pose.conformation().chain_end( ref_prot_chain ) );
	Size ref_prot_anchor( ref_prot_begin );

	Size pep_nterm( pep_anchor - pep_begin );
	Size ref_pep_nterm( ref_pep_anchor - ref_pep_begin );
	Size pep_cterm( pep_end - pep_anchor );
	Size ref_pep_cterm( ref_pep_end - ref_pep_anchor );

	Size nterm( pep_nterm );
	Size cterm( pep_cterm );
	if( pep_nterm < ref_pep_nterm ) ref_pep_begin += ref_pep_nterm - pep_nterm;
	else if( pep_nterm > ref_pep_nterm ){
		pep_begin += pep_nterm - ref_pep_nterm;
		nterm = ref_pep_nterm;
	}
	if( pep_cterm < ref_pep_cterm ) ref_pep_end -= ref_pep_cterm - pep_cterm;
	else if( pep_cterm > ref_pep_cterm ){
		pep_end -= pep_cterm - ref_pep_cterm;
		cterm = ref_pep_cterm;
	}

	std::string phipsi_analysis( "" );
	for( int i_seq = 0; i_seq <= nterm + cterm; ++i_seq ){
		Size ref_pep_seqpos = ref_pep_begin + i_seq;
		Size pep_seqpos = pep_begin + i_seq;
		Real ramadev( 0 );

		//delta phi, psi, omega
		Real phi( std::abs( ref_pose.phi( ref_pep_seqpos ) - pose.phi( pep_seqpos ) ) );
		if( phi > 180 ) phi = std::abs( 360 - phi );
		ramadev += ( phi * phi );
//		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_phi:\t" + string_of( pose.phi( pep_seqpos ) ) + "\t";
		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_phidev:\t" + string_of( phi ) + "\t";

		Real psi( std::abs( ref_pose.psi( ref_pep_seqpos ) - pose.psi( pep_seqpos ) ) );
		if( psi > 180 ) psi = std::abs( 360 - psi );
		ramadev += ( psi * psi );
//		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_psi:\t" + string_of( pose.psi( pep_seqpos ) ) + "\t";
		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_psidev:\t" + string_of( psi ) + "\t";

		Real omega( std::abs( ref_pose.omega( ref_pep_seqpos ) - pose.omega( pep_seqpos ) ) );
		if( omega > 180 ) omega = std::abs( 360 - omega );
//		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_omega:\t" + string_of( pose.omega( pep_seqpos ) ) + "\t";
		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_omegadev:\t" + string_of( omega ) + "\t";

		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_ramadev:\t" + string_of( std::sqrt( ramadev / 2 ) ) + "\t";
	}
	return phipsi_analysis;

}

void
pep_energies_analysis(
	pose::Pose pose,
	core::scoring::ScoreFunctionOP full_scorefxn,
	std::string filename,
	core::scoring::EMapVector & emap_total_avg,
	vector1< core::scoring::EMapVector > & emap_res_avg,
	vector1< core::scoring::EMapVector > & emap_res_sd
)
{
	( *full_scorefxn )( pose );
	EMapVector emap_total( pose.energies().total_energies() );
	emap_total_avg += emap_total;

	for( Size i_seq = 1; i_seq <= pose.total_residue(); ++i_seq ){
		emap_res_avg[ i_seq ] += pose.energies().residue_total_energies( i_seq );
		char pep_aa( pose.residue( i_seq ).name1() );
		std::cout << filename << "\t" << pep_aa << "\t" << string_of( i_seq ) << "\ttotal_res\t" << pose.energies().residue_total_energies( i_seq ).weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< pose.energies().residue_total_energies( i_seq ).dot( full_scorefxn->weights() ) << "\n";
	}
	std::cout << std::endl;
}

void
pep_scan_analysis(
	pose::Pose pose,
	core::scoring::ScoreFunctionOP soft_scorefxn,
	core::scoring::ScoreFunctionOP full_scorefxn,
	std::string filename,
	Size pep_begin,
	Size pep_anchor,
	Size pep_end,
	vector1< core::scoring::EMapVector > & avg_res_avg_emap
)
{
	( *full_scorefxn )( pose );
	Pose start_pose( pose );

	vector1< bool > is_pep( pose.total_residue(), false );
	for( Size i = 1; i <= pose.total_residue(); ++i ){
		if( i >= pep_begin && i <= pep_end ) is_pep[ i ] = true;
	}

	if( option[ pep_spec::scan_seqpos ].user() ){
		pep_begin = option[ pep_spec::scan_seqpos ];
		pep_end = option[ pep_spec::scan_seqpos ];
	}
	for( Size pep_seqpos = pep_begin; pep_seqpos <= pep_end; ++pep_seqpos ){
		if( pep_seqpos == pep_anchor ) continue;
		pose = start_pose;

		char pep_aa( pose.residue( pep_seqpos ).name1() );

		float cutoff( option[ pep_spec::interface_cutoff ] );
		vector1< bool > is_mut_nbr( pose.total_residue(), false );
		Residue const & rsd1( pose.residue( pep_seqpos ) );
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			Residue const & rsd2( pose.residue(j) );
			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
					if ( rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) ) < cutoff*cutoff ) {
						is_mut_nbr[j] = true;
						break;
					}
				}
			}
		}

		{
			// the movable dof's
			kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
			mm_min->set_chi( is_pep );
			mm_min->set_chi( is_mut_nbr );
			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true );

			//define design task and repack task
			pack::task::RestrictResidueToRepackingOperationOP restrict_to_repack_taskop( new pack::task::RestrictResidueToRepackingOperation() );
			pack::task::PreventRepackingOperationOP prevent_repack_taskop( new pack::task::PreventRepackingOperation() );
			pack::task::PackerTaskOP rp_task( pack::task::TaskFactory::create_packer_task( pose ));
			rp_task->initialize_from_command_line().or_include_current( true );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( is_mut_nbr[i] || i == pep_seqpos ) {
					restrict_to_repack_taskop->include_residue( i );
					rp_task->nonconst_residue_task( i ).restrict_to_repacking();
				} else {
					rp_task->nonconst_residue_task( i ).prevent_repacking();
					prevent_repack_taskop->include_residue( i );
				}
			}
			pack::task::TaskFactoryOP rottrial_task_factory( new pack::task::TaskFactory );
			rottrial_task_factory->push_back( new pack::task::InitializeFromCommandlineOperation() );
			rottrial_task_factory->push_back( new pack::task::IncludeCurrentOperation() );
			rottrial_task_factory->push_back( restrict_to_repack_taskop );
			rottrial_task_factory->push_back( prevent_repack_taskop );

			protocols::simple_moves::PackRotamersMoverOP pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, rp_task, 1 ) );
			protocols::simple_moves::RotamerTrialsMoverOP rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, rottrial_task_factory ) );
			SequenceMoverOP seq = new SequenceMover;
			if( !option[ pep_spec::test_no_pack ] ){
				seq->add_mover( pack );
				seq->add_mover( rottrial );
			}
			if( !option[ pep_spec::test_no_min ] ){
				seq->add_mover( min_mover );
			}

			seq->apply( pose );
			( *full_scorefxn )( pose );

		}

		Pose mut_start_pose( pose );
		EMapVector start_emap( pose.energies().total_energies() );
		Real start_total( start_emap.dot( full_scorefxn->weights() ) );
		EMapVector res_avg_emap;
		Real res_avg_total( 0 );
		Size n_aas( 0 );
		for( Size resindex = 1; resindex <=20; ++resindex ){

			pose = mut_start_pose;
			char mut_aa( chemical::oneletter_code_from_aa( chemical::AA( resindex ) ) );

			if( option[ pep_spec::scan_restypes ].user() ){
				std::string scan_restypes( option[ pep_spec::scan_restypes ] );
				bool skip( true );
				for( Size ii = 0; ii <= scan_restypes.size() - 1; ++ii ){
					if( scan_restypes[ ii ] == mut_aa ){
						skip = false;
						break;
					}
				}
				if( skip ) continue;
			}
			else if( mut_aa == 'P' ) continue;

			chemical::make_sequence_change( pep_seqpos, chemical::AA(resindex), pose );
			++n_aas;

			// the movable dof's
			kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
			mm_min->set_chi( is_pep );
			mm_min->set_chi( is_mut_nbr );
			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true );

			//define design task and repack task
			pack::task::RestrictResidueToRepackingOperationOP restrict_to_repack_taskop( new pack::task::RestrictResidueToRepackingOperation() );
			pack::task::PreventRepackingOperationOP prevent_repack_taskop( new pack::task::PreventRepackingOperation() );
			pack::task::PackerTaskOP rp_task( pack::task::TaskFactory::create_packer_task( pose ));
			rp_task->initialize_from_command_line().or_include_current( true );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( is_mut_nbr[i] ) {
					restrict_to_repack_taskop->include_residue( i );
					rp_task->nonconst_residue_task( i ).restrict_to_repacking();
				} else {
					rp_task->nonconst_residue_task( i ).prevent_repacking();
					prevent_repack_taskop->include_residue( i );
				}
			}
			pack::task::TaskFactoryOP rottrial_task_factory( new pack::task::TaskFactory );
			rottrial_task_factory->push_back( new pack::task::InitializeFromCommandlineOperation() );
			rottrial_task_factory->push_back( new pack::task::IncludeCurrentOperation() );
			rottrial_task_factory->push_back( restrict_to_repack_taskop );
			rottrial_task_factory->push_back( prevent_repack_taskop );

			protocols::simple_moves::PackRotamersMoverOP pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, rp_task, 1 ) );
			protocols::simple_moves::RotamerTrialsMoverOP rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, rottrial_task_factory ) );
			SequenceMoverOP design_seq = new SequenceMover;
			if( !option[ pep_spec::test_no_pack ] ){
				design_seq->add_mover( pack );
				design_seq->add_mover( rottrial );
			}
			if( !option[ pep_spec::test_no_min ] ){
				design_seq->add_mover( min_mover );
			}

			design_seq->apply( pose );
			( *full_scorefxn )( pose );

			EMapVector diff_emap( pose.energies().total_energies() );
			diff_emap -= start_emap;
			Real diff_total( pose.energies().total_energies().dot( full_scorefxn->weights() ) - start_total );

			res_avg_emap += diff_emap;
			res_avg_total += diff_total;

			std::cout << pep_aa << "\t" << string_of( pep_seqpos ) << "\t" << mut_aa << "\tspec_res:\t" << diff_emap.weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< diff_total << "\n";

		}

		for( EMapVector::iterator itr = res_avg_emap.begin(); itr != res_avg_emap.end(); ++itr ){
			*itr *= ( 1.0 / n_aas );
		}
		res_avg_total = res_avg_total / n_aas;

		std::cout << pep_aa << "\t" << string_of( pep_seqpos ) << "\t" << "X" << "\tspec_res_avg:\t" << res_avg_emap.weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< res_avg_total << "\n";

		avg_res_avg_emap[ pep_seqpos ] += res_avg_emap;

	}

}


void
RunPepSpec()
{

	vector1< Pose > poses;
	vector1< double > nrgs;
	vector1< std::string > pdb_filenames;
	std::string pdb_list_filename( option[ pep_spec::pdb_list ] );
	Pose pose;
	{
		std::ifstream pdb_list_data( pdb_list_filename.c_str() );
		if ( !pdb_list_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + pdb_list_filename + '\n' );
		}
		std::string pdb_list_line;
		getline( pdb_list_data, pdb_list_line, '\t' );
		std::string filename( pdb_list_line );
		core::import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
	}

	Size pep_anchor_in( option[ pep_spec::pep_anchor ] );
	std::string const pep_chain_in( option[ pep_spec::pep_chain ] );
	Size pep_anchor( pose.pdb_info()->pdb2pose( pep_chain_in[0], pep_anchor_in ) );
	Size pep_chain( pose.chain( pep_anchor ) );
	Size pep_begin( pose.conformation().chain_begin( pep_chain ) );
	Size pep_end( pose.conformation().chain_end( pep_chain ) );

	Size prot_chain( 1 );
	for( Size i = 1; i <= pose.conformation().num_chains(); ++i ){
		if( i == pep_chain ) continue;
		if( !( pose.residue( pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			prot_chain = i;
			break;
		}
	}
	Size prot_begin( pose.conformation().chain_begin( prot_chain ) );
	Size prot_end( pose.conformation().chain_end( prot_chain ) );
	Size prot_anchor( prot_begin );

/*
	pose::Pose ref_pose;
	std::string ref_name( option[ pep_spec::ref_pose ] );
	core::import_pose::pose_from_file( ref_pose, ref_name , core::import_pose::PDB_file);

	Size ref_pep_anchor_in( option[ pep_spec::pep_anchor ] );
	if( option[ pep_spec::ref_pep_anchor ].user() ) ref_pep_anchor_in = option[ pep_spec::ref_pep_anchor ];
	std::string ref_pep_chain_in( option[ pep_spec::pep_chain ] );
	if( option[ pep_spec::ref_pep_chain ].user() ) ref_pep_chain_in = option[ pep_spec::ref_pep_chain ];
	Size ref_pep_anchor( ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_in ) );
	Size ref_pep_chain( ref_pose.chain( ref_pep_anchor ) );
	Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pep_chain ) );
	Size ref_pep_end( ref_pose.conformation().chain_end( ref_pep_chain ) );

	Size ref_prot_chain;
	for( Size i = 1; i <= ref_pose.conformation().num_chains(); ++i ){
		if( i == ref_pep_chain ) continue;
		if( !( ref_pose.residue( ref_pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			ref_prot_chain = i;
			break;
		}
	}
	Size ref_prot_begin( ref_pose.conformation().chain_begin( ref_prot_chain ) );
	Size ref_prot_end( ref_pose.conformation().chain_end( ref_prot_chain ) );
	Size ref_prot_anchor( ref_prot_begin );

	Size pep_nterm( pep_anchor - pep_begin );
	Size ref_pep_nterm( ref_pep_anchor - ref_pep_begin );
	Size pep_cterm( pep_end - pep_anchor );
	Size ref_pep_cterm( ref_pep_end - ref_pep_anchor );

	Size nterm( pep_nterm );
	Size cterm( pep_cterm );
	if( pep_nterm < ref_pep_nterm ) ref_pep_begin += ref_pep_nterm - pep_nterm;
	else if( pep_nterm > ref_pep_nterm ){
		pep_begin += pep_nterm - ref_pep_nterm;
		nterm = ref_pep_nterm;
	}
	if( pep_cterm < ref_pep_cterm ) ref_pep_end -= ref_pep_cterm - pep_cterm;
	else if( pep_cterm > ref_pep_cterm ){
		pep_end -= pep_cterm - ref_pep_cterm;
		cterm = ref_pep_cterm;
	}
	std::cout << "Nterm " << string_of( nterm ) << ", Cterm " << string_of( cterm ) << std::endl;
*/
	core::scoring::ScoreFunctionOP full_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::wts ] ) );
	core::scoring::ScoreFunctionOP soft_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::soft_wts ] ) );

	EMapVector emap_total_avg;
	vector1< EMapVector > emap_res_avg( pose.total_residue() );
	vector1< EMapVector > emap_res_sd( pose.total_residue() );

	vector1< EMapVector > emap_res_scan( pose.total_residue() );

	std::ifstream pdb_list_data( pdb_list_filename.c_str() );
	std::string pdb_list_line;
	Size n_structs( 0 );
	while( !getline( pdb_list_data, pdb_list_line, '\t' ).eof() ) {
		++n_structs;

		std::string filename( pdb_list_line );
		core::import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
		std::string data;
		getline( pdb_list_data, data );

/*
		if( option[ pep_spec::rescore_analysis ] ){
			( *full_scorefxn )( pose );
			vector1< bool > is_pep( pose.total_residue(), false );
			for ( Size i=1; i<= pose.total_residue(); ++i ) is_pep[i] = ( i >= pep_begin && i <= pep_end );
			vector1< bool > is_pep_nbr( pose.total_residue(), false );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const & rsd1( pose.residue(i) );
				if ( is_pep[i] ) continue;
				for ( Size j=1; j<= pose.total_residue(); ++j ) {
					Residue const & rsd2( pose.residue(j) );
					if ( !is_pep[j] ) continue;
					if( is_pep_nbr[i] ) break;
					for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
						if( is_pep_nbr[i] ) break;
						for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
							if ( rsd1.xyz(ii).distance( rsd2.xyz(jj) ) < cutoff ) {
								is_pep_nbr[i] = true;
								break;
							}
						}
					}
				}
			}
			kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
			mm_min->set_chi( is_pep );
			mm_min->set_chi( is_pep_nbr );
			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true );
		}
*/

		if( option[ pep_spec::coord_cluster_analysis ] ){
			Real cutoff( option[ pep_spec::coord_cluster_cutoff ] );
			std::ifstream pdb_list_data2( pdb_list_filename.c_str() );
			std::string pdb_list_line2;
			Pose pose2;
			bool has_twin( false );
			while( !getline( pdb_list_data2, pdb_list_line2, '\t' ).eof() ) {
				std::string filename2( pdb_list_line2 );
				core::import_pose::pose_from_file( pose2, filename2 , core::import_pose::PDB_file);
				std::string data2;
				getline( pdb_list_data2, data2 );
				Real rmsd( get_rmsd( pose, prot_begin, prot_end, pep_begin, pep_anchor, pep_end, pose2, prot_begin, prot_end, pep_begin, pep_anchor, pep_end ) );
				if( rmsd < cutoff ){
					has_twin = true;
					break;
				}
			}
			if( has_twin ) continue;
		}

		std::cout << filename << "\t" << data << "\t";
		if( option[ pep_spec::rmsd_analysis ] ) std::cout << pep_rmsd_analysis( pose, prot_begin, prot_end, pep_begin, pep_anchor, pep_end );
		if( option[ pep_spec::phipsi_analysis ] ) std::cout << pep_phipsi_analysis( pose, prot_begin, prot_end, pep_begin, pep_anchor, pep_end );
		std::cout << "\n";

		if( option[ pep_spec::energies_analysis ] ) pep_energies_analysis( pose, full_scorefxn, filename, emap_total_avg, emap_res_avg, emap_res_sd );
		if( option[ pep_spec::scan_analysis ] ) pep_scan_analysis( pose, soft_scorefxn, full_scorefxn, filename, pep_begin, pep_anchor, pep_end, emap_res_scan );

	}

	//energies avg over all poses
	if( option[ pep_spec::energies_analysis ] && n_structs > 1 ){
		for( Size i_seq = 1; i_seq <= pose.total_residue(); ++i_seq ){
			for( EMapVector::iterator itr = emap_res_avg[ i_seq ].begin(); itr != emap_res_avg[ i_seq ].end(); ++itr ){
				*itr *= ( 1.0 / n_structs );
			}
			char pep_aa( pose.residue( i_seq ).name1() );
			if( i_seq >= pep_begin && i_seq <= pep_end ) pep_aa = 'X';
			std::cout << "avg" << "\t" << pep_aa << "\t" << string_of( i_seq ) << "\ttotal_res\t" << emap_res_avg[ i_seq ].weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< emap_res_avg[ i_seq ].dot( full_scorefxn->weights() ) << "\n";
		}

		for( EMapVector::iterator itr = emap_total_avg.begin(); itr != emap_total_avg.end(); ++itr ){
			*itr *= ( 1.0 / n_structs );
		}
		std::cout << "total_avg\t" << emap_total_avg.weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< emap_total_avg.dot( full_scorefxn->weights() ) << "\n";
	}

	//scan avg over all poses
	if( option[ pep_spec::scan_analysis ] && n_structs > 1 ){
		if( option[ pep_spec::scan_seqpos ].user() ){
			pep_begin = option[ pep_spec::scan_seqpos ];
			pep_end = option[ pep_spec::scan_seqpos ];
		}
		for( Size i_seq = pep_begin; i_seq <= pep_end; ++i_seq ){
			if( i_seq == pep_anchor ) continue;
			for( EMapVector::iterator itr = emap_res_scan[ i_seq ].begin(); itr != emap_res_scan[ i_seq ].end(); ++itr ){
				*itr *= ( 1.0 / n_structs );
			}
			char pep_aa = 'X';
			if( option[ pep_spec::scan_seqpos ].user() ) pep_aa = pose.residue( i_seq ).name1();
			std::cout << "all_pdbs_avg\n";
			std::cout << pep_aa << "\t" << string_of( i_seq ) << "\tX\t" << "\ttotal_spec_res_avg\t" << emap_res_scan[ i_seq ].weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< emap_res_scan[ i_seq ].dot( full_scorefxn->weights() ) << "\n";
		}
	}

	pdb_list_data.close();
}

int main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init(argc, argv);

		RunPepSpec();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
