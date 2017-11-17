// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file sidechain_min.cc
/// @brief Minimizes non-disulfide-bonded side-chains in the input structure
/// @author Daniel J. Mandell

// Unit Headers
#include <protocols/loops/Loops.hh>
#include <protocols/viewer/viewers.hh>


// Rosetta Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <basic/options/option.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <devel/init.hh>


#include <basic/Tracer.hh> // tracer output

#include <utility/io/izstream.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;

static basic::Tracer TR( "procotols.looprelax" );

// some functions copied from loops_main (they were private)

//////////////////////////////////////////////////////////////////////////////////////
/// @details omega backbone torsion is always fixed. phi/psi backbone torsions within
/// the loop region are flexible. Depending on whether -fix_natsc flag, sidechain DOFs
/// of loop residues and/or their neighboring residues in the template will be set as
/// movable.
//////////////////////////////////////////////////////////////////////////////////////
void
loops_set_move_map(
				   protocols::loops::Loops const & loops,
				   bool const fix_template_sc,
				   core::kinematics::MoveMap & mm
				   )
{
	using namespace core::id;

	// allow chi to move
	mm.set_bb( false );
	mm.set_chi( !fix_template_sc );
	mm.set_jump( false );
	// allow phi/psi in loops to move
	for( protocols::loops::Loops::const_iterator it=loops.begin(), it_end=loops.end();
		 it != it_end; ++it ) {
		for( Size i=it->start(); i<=it->stop(); ++i ) {
			mm.set_bb(i, true);
			mm.set(TorsionID(i,BB,3), false); // omega is fixed
			mm.set_chi(i, true); // chi of loop residues
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details use TenANeighborGraph. As input, residue_positions[i] is true for residues to be counted.
/// As output, residue_position[i] is true for all neighbor residues including orginal input residues.
/// The function is used to find all neighboring residues of the loop residues in case they need to be
/// repacked or minimized in fullatom refinement.
void
get_tenA_neighbor_residues(
						   pose::Pose const & pose,
						   utility::vector1<bool> & residue_positions
						   )
{
	//make a local copy first because we will change content in residue_positions
	utility::vector1<bool> local_residue_positions = residue_positions;
	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );
	for ( int i=1; i <= int(local_residue_positions.size()); ++i ) {
		if ( ! local_residue_positions[i] ) continue;
		utility::graph::Node const * current_node( tenA_neighbor_graph.get_node(i)); // find neighbors for this node
		for ( utility::graph::Node::EdgeListConstIter it = current_node->const_edge_list_begin();
			  it != current_node->const_edge_list_end(); ++it ) {
			residue_positions[ (*it)->get_other_ind(i) ] = true;
		}
	}
}

void
fold_tree_from_loops(
					 Size const total_residue,
					 protocols::loops::Loops const & loops,
					 kinematics::FoldTree & f
					 )
{
	using namespace kinematics;

	f.clear();

	protocols::loops::Loops tmp_loops = loops;
	tmp_loops.sequential_order();

	Size jump_num = 0;
	for( protocols::loops::Loops::const_iterator it=tmp_loops.begin(), it_end=tmp_loops.end(),
		 it_next; it != it_end; ++it ) {
		it_next = it;
		it_next++;
		Size const jump_start =
			( it->start() == 1 ) ? it->start() : it->start() - 1;
		Size const jump_stop  =
			( it->stop() == total_residue ) ? it->stop() : it->stop() + 1;
		Size const jump_cut   = it->cut();
		Size const jump_next_start =
			( it_next == it_end ) ? total_residue : it_next->start()-1;
		if(  it->start() == 1 ){
			f.add_edge( jump_start, jump_stop, Edge::PEPTIDE );
			f.add_edge( jump_stop, jump_next_start, Edge::PEPTIDE );
			continue;
		}else if( it->stop() == total_residue ){
			f.add_edge( jump_start, jump_stop, Edge::PEPTIDE );
			continue;
		}
		jump_num++;
		f.add_edge( jump_start, jump_stop, jump_num );
		f.add_edge( jump_start, jump_cut,  Edge::PEPTIDE );
		f.add_edge( jump_cut+1, jump_stop, Edge::PEPTIDE );
		//		if ( jump_stop < jump_next_start )
		f.add_edge( jump_stop, jump_next_start, Edge::PEPTIDE );
	}
	Size const first_start =
		( tmp_loops.begin()->start() == 1 ) ? tmp_loops.begin()->start() : tmp_loops.begin()->start() - 1;
	//	if ( first_start != 1 )
	f.add_edge( 1, first_start, Edge::PEPTIDE );

	// reorder
	Size root;
	if( tmp_loops.begin()->start() == 1 &&
		tmp_loops.begin()->stop() != total_residue )
		root = tmp_loops.begin()->stop()+1;
	else root = 1;
	f.reorder(root);

}


//////////////////////////////////////////////////////////////////
/// @details read in secondary structure definition from a psipred_ss2 file and
/// set that to a Pose. A temporary setup since there is very limited ss support in
/// Pose right now.
///////////////////////////////////////////////////////////////////
bool
set_secstruct_from_psipred_ss2(
							   pose::Pose & pose
							   )
{
	using namespace basic::options;

	std::string filename( option[ OptionKeys::in::file::psipred_ss2 ]().name() );

	utility::io::izstream data( filename );
	if ( !data ) {
		TR.Error << "can not open psipred_ss2 file " << filename << std::endl;
		return false;
	}

	utility::vector1< char > secstructs;
	std::string line;
	Size count(0);
	while( getline( data, line ) ) {
		std::istringstream line_stream( line );
		Size pos;
		char aa, sec;
		line_stream >> pos >> aa >> sec;
		count++;
		if ( aa != oneletter_code_from_aa( pose.aa(count) ) ) {
			TR.Error << " sequence mismatch between pose and psipred_ss2 " << oneletter_code_from_aa( pose.aa(count) )
			<< " vs " << aa << " at seqpos " << count << std::endl;
			return false;
		}
		if ( sec != 'H' && sec != 'E' && sec != 'C' ) {
			TR.Error << "unrecognized secstruct char : " << sec << " at seqpos " << count << std::endl;
		}
		if ( sec == 'C' ) {
			secstructs.push_back( 'L' );
		} else {
			secstructs.push_back( sec );
		}
	}
	assert( secstructs.size() == pose.size() );
	for ( Size i = 1, iend=secstructs.size(); i <= iend; ++i ) {
		pose.set_secstruct( i, secstructs[i] );
	}

	return true;

}


/////////////////////////////////////////////////////////////////////////////
/// @details use 10A CB distance cutoff as neighboring residue defintion. The function
///is used for conveniently setting up sidechain movable residues in loop modeling.
/////////////////////////////////////////////////////////////////////////////////
void
select_loop_residues(
										 pose::Pose const & pose,
										 protocols::loops::Loops const & loops,
										 bool const include_neighbors,
										 utility::vector1<bool> & map
										 )
{
	for( protocols::loops::Loops::const_iterator it=loops.begin(), it_end=loops.end();
		 it != it_end; ++it ) {
		for( Size i=it->start(); i<=it->stop(); ++i ) {
			map[i] = true;
		}
	}

	if ( include_neighbors) get_tenA_neighbor_residues( pose, map );

	return;
}

/////////
///////// main function
/////////

void*
my_main( void* )
{
	using namespace basic::options;
	using namespace chemical;
	using namespace scoring;
	using namespace id;
	using namespace optimization;

	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, option[ OptionKeys::loops::input_pdb ]().name() , core::import_pose::PDB_file);
 	int const nres( pose.size() );

 	scoring::ScoreFunctionOP scorefxn( get_score_function() );
	//scorefxn->set_weight(scoring::mm_bend, 1.00);

	//bool const fa_input( option[OptionKeys::loops::fa_input] );

	// minimizer
	AtomTreeMinimizer minimizer;
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/, false /*deriv_check*/ );

	// set up the packer task
	pack::task::PackerTaskOP base_packer_task( pack::task::TaskFactory::create_packer_task( pose ));
	base_packer_task->initialize_from_command_line().restrict_to_repacking().or_include_current( false );
	base_packer_task->set_bump_check( true );

	// get the starting energy
	Real starting_total_energy = (*scorefxn)(pose);

	// do the repack and minimization
	pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
	utility::vector1<bool> allow_repacked( nres, true );

	// turn off side-chain repacking and minimization for disulfide-bonded cysteines
	// set up move map for minimizing non-disulfide bonded residues only
	kinematics::MoveMap mm_all_sc;
	mm_all_sc.set_bb( false );
	mm_all_sc.set_chi( false );
	for (core::Size i=1; i<= pose.size(); i++) {
		//if ( ( pose.residue(i).type().name() == "CYD" ) || ( pose.residue(i).type().name() == "CYX" ) ) {	// DJM: DEBUG
		// amw: changing to check variant types
		if ( pose.residue(i).type().is_disulfide_bonded() || pose.residue(i).has_variant_type( chemical::SIDECHAIN_CONJUGATION ) {	//
			allow_repacked[i] = false;
			TR << "Disabling side-chain optimization of disulfide bonded residue " << i << std::endl;
		}
		else if (! pose.residue(i).is_protein() ) {
			allow_repacked[i] = false;
		}
		else {
			mm_all_sc.set_chi( i, true );
		}
	}
	this_packer_task->restrict_to_residues( allow_repacked );

	// perform the minimization
	pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	pack::rotamer_trials( pose, *scorefxn, this_packer_task );
	minimizer.run( pose, mm_all_sc, *scorefxn, options );

	// output the minimized pose
	Real last_total_energy = (*scorefxn)(pose);
	std::string out_tag = option[ OptionKeys::out::output_tag ];
	std::string outname=option[ OptionKeys::out::path::path ]().name()+out_tag+"_min.pdb";
	std::ofstream out(outname.c_str(), std::ios::out | std::ios::binary);
	core::io::pdb::dump_pdb( pose, out );
	out << "total_energy_before_minimize: " << starting_total_energy << std::endl;
	out << "total_energy_after_minimize: " << last_total_energy << std::endl;

	return 0;
}

int
main( int argc, char * argv [] )
{
	try {

	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
  return 0;
}
