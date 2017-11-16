// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @author Eva-Maria Strauch ( evas01@u.washington.edu )

//Unite headers
#include <protocols/seeded_abinitio/SeededAbinitio_util.hh>
#include <protocols/seeded_abinitio/SeedFoldTree.hh>
#include <protocols/seeded_abinitio/SeedFoldTreeCreator.hh>
#include <protocols/seeded_abinitio/SeedFoldTree.fwd.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/xyzVector.hh>
#include <basic/datacache/DataMap.hh>
//Auto Headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/loops_main.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/protein_interface_design/util.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <set>
#include <string>
#include <utility>
#include <basic/Tracer.hh>
#include <protocols/simple_filters/AlaScan.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using namespace core::scoring;
using namespace protocols::seeded_abinitio;

static basic::Tracer TR( "protocols.seeded_abinitio.SeedFoldTree" );

namespace protocols {
namespace seeded_abinitio {

using namespace protocols::moves;
using namespace core;

// XRW TEMP std::string
// XRW TEMP SeedFoldTreeCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SeedFoldTree::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SeedFoldTreeCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SeedFoldTree );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SeedFoldTree::mover_name()
// XRW TEMP {
// XRW TEMP  return "SeedFoldTree";
// XRW TEMP }

SeedFoldTree::~SeedFoldTree() = default;

SeedFoldTree::SeedFoldTree() :
	protocols::moves::Mover( SeedFoldTree::mover_name() ),
	fold_tree_( /* NULL */ ),
	scorefxn_( /* NULL */ )
{}

SeedFoldTree::SeedFoldTree( core::kinematics::FoldTreeOP ft ) :
	protocols::moves::Mover( SeedFoldTree::mover_name() )
{
	fold_tree_ = ft;
	pdb_contains_target_ = false;
	anchor_specified_ = false;
	twochains_ = false;
	cut_points_.clear();
	all_seeds_.clear();
	set_jumps_manually = false;
	anchors_.clear();
	ddg_based_ = true;
	anchors_.clear();
	//manual_jump_pairs_.clear();
	folding_vertices_.clear();

}

protocols::moves::MoverOP
SeedFoldTree::clone() const {
	return( protocols::moves::MoverOP( new SeedFoldTree( *this ) ) );
}

protocols::moves::MoverOP
SeedFoldTree::fresh_instance() const {
	return protocols::moves::MoverOP( new SeedFoldTree );
}

void
SeedFoldTree::fold_tree( core::kinematics::FoldTreeOP ft ) {
	fold_tree_ = ft;
}

core::kinematics::FoldTreeOP
SeedFoldTree::fold_tree() const {
	return( fold_tree_ );
}

core::scoring::ScoreFunctionOP
SeedFoldTree::scorefxn() const{
	return scorefxn_;
}

bool
SeedFoldTree::ddg_based(){
	return ddg_based_;
}

void
SeedFoldTree::ddg_based( bool ddgb ){
	ddg_based_ = ddgb;
}


void
SeedFoldTree::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}


void
SeedFoldTree::set_anchor_res( utility::vector1< core::Size > anchor ){
	anchors_.clear();
	for ( Size i = 1; i <= anchor.size(); ++i ) {
		anchors_.push_back( anchor[i] );
	}
}

void
SeedFoldTree::anchor_specified( bool anchor_specified ){
	anchor_specified_ = anchor_specified;
}

bool
SeedFoldTree::anchor_specified(){
	return anchor_specified_;
}

Size define_cut_point_stochasticly (
	Size start_resi,
	Size stop_resi,
	std::string secondarystruct_seq,
	core::Size start_fold_pose
){
	TR <<"defining cut points stochasticly between the given two residues: "<<start_resi<<" and "<<stop_resi <<std::endl;
	utility::vector1<Size> loopy_regions;
	//char ss;

	TR<<"start and stop: " << start_resi << " " << stop_resi << "\nsecondary structure string between seeds: \n";
	//going through the regions beteween the seeds, to identify loops and assign a cut point at random
	for ( Size resi = start_resi + 1 /*+1 since strings count from 0 */ ; resi < stop_resi - 1/*as long as it is a set*/; ++resi ) {
		char ss = secondarystruct_seq[ resi - 1 ];
		TR << ss ;
		if ( ss == 'L' ) {
			loopy_regions.push_back( resi );//to adjust for string counting
		}
	}
	TR << std::endl;
	//picking one at random:
	if ( loopy_regions.size() < 1 ) {
		utility_exit_with_message("there are no loopy residues between the motifs, this is currently not supported.");
	}

	int low = 1;
	int high = loopy_regions.size();
	core::Size ran = numeric::random::rg().random_range( low, high ); // todo: bias more for center
	core::Size cutpoint = loopy_regions[ ran ] + start_fold_pose -1 ;

	TR.Debug <<"random number: "<< ran << ", number from loop container: " << loopy_regions[ran]<< ", adjusting by " << start_fold_pose - 1 <<", cutpoint: " << cutpoint << std::endl;
	TR<<"picked a cutpoint between "<<start_resi << " and " << stop_resi << " ( "<<cutpoint - start_fold_pose <<" ). Renumbering will be adjusted by "<< start_fold_pose -1 << std::endl;

	return cutpoint;
}


//this method identifies the closest pair of residues between the target and the seed
std::pair< Size, Size >
get_closest_residue_pair(
	Size seed_start,
	Size seed_stop,
	core::pose::PoseOP & target_seed_pose
){

	using namespace core::conformation;
	using namespace core::chemical;

	core::Size nearest_resi( 0 );
	core::Size nearest_resi_target( 0 );
	core::Real nearest_dist( 100000.0 );
	core::Real nearest_dist2 (100000.0 );
	std::pair<Size,Size> closest_pair;
	if ( target_seed_pose->conformation().num_chains() != 2 ) {
		utility_exit_with_message("only two chains as input supported" );
	}
	Size target_length = target_seed_pose->split_by_chain( 1 )->size();

	TR<<"iterating through each seed residue to find the closest target residue" <<std::endl;

	for ( Size seed_resi = seed_start; seed_resi <= seed_stop ; ++seed_resi ) {
		for ( Size target_resi = 1; target_resi <= target_length ; ++target_resi ) {

			//should be made fancier, using CB instead except for gly (todo)
			core::Real const distance( target_seed_pose->residue(target_resi).xyz( "CA" ).distance(target_seed_pose->residue(seed_resi).xyz( "CA" ) ) );

			if ( distance < nearest_dist ) {
				nearest_resi_target = target_resi;
				nearest_dist = distance;
			}

			if ( nearest_dist < nearest_dist2 ) {
				nearest_resi = seed_resi ;
				nearest_dist2 = nearest_dist;
			}
		}
	}
	runtime_assert( nearest_resi );
	runtime_assert( nearest_resi_target );
	TR<<"closest pair between seed and target: "<<nearest_resi <<" and: " <<nearest_resi_target <<", distance: "<<nearest_dist<<std::endl;
	closest_pair.second = nearest_resi;
	closest_pair.first = nearest_resi_target;

	return ( closest_pair  );
}

std::pair< Size, Size >
find_nearest_residue( Size anchor , pose::PoseOP & target_seed_pose ){
	using namespace core::conformation;
	using namespace core::chemical;

	Size nearest_resi( 0 );
	Real nearest_dist( 100000.0 );
	std::pair<Size,Size> closest_pair;
	if ( target_seed_pose->conformation().num_chains() != 2 ) {
		utility_exit_with_message("only two chains as input are currently supported" );
	}
	Size target_length = target_seed_pose->split_by_chain( 1 )->size();

	std::string anchor_atom = "CB";
	Residue res_anchor( target_seed_pose->residue( anchor ));
	if ( res_anchor.name3() == "GLY" ) {
		anchor_atom = "CA";
	}

	for ( Size target_resi = 1; target_resi <= target_length ; ++target_resi ) {

		std::string connect_atom = "CB";
		Residue res( target_seed_pose->residue( target_resi ));
		if ( res.name3() == "GLY" ) {
			connect_atom = "CA";
		}

		core::Real const distance( res.xyz( connect_atom ).distance( res_anchor.xyz( anchor_atom ) ) );

		if ( distance < nearest_dist ) {
			nearest_resi = target_resi;
			nearest_dist = distance;
		}
	}
	runtime_assert( nearest_resi );
	TR<<"closest residue to anchor residue "<< anchor <<" is "<<nearest_resi <<", distance: "<<nearest_dist<<std::endl;

	closest_pair.second = anchor;
	closest_pair.first = nearest_resi;

	return ( closest_pair  );

}

core::Size
SeedFoldTree::best_by_ala_scan( Size start, Size end, pose::PoseOP & ts_pose ){
	//runtime_assert( end <= ts_pose->size() );
	TR << "tspose: " << ts_pose->size() <<std::endl;
	TR <<"----------alanine scanning to identify the best jump atom in seed  ------------"<<std::endl;
	protocols::simple_filters::AlaScan ala_scan;
	ala_scan.repack( 0 );
	ala_scan.repeats( 1 );
	ala_scan.jump( 1 );
	ala_scan.scorefxn( scorefxn() );

	protocols::simple_filters::DdgFilter ddg_filter( 100/*ddg_threshold*/, scorefxn(), 1 /*jump*/, 1 /*repeats*/ );
	core::Real orig_dG(0.0);
	core::Real lowest_dG = 1000000;
	Size lowest_res = 0;
	pose::Pose pose( *ts_pose );
	orig_dG = ddg_filter.compute( pose );
	TR <<"\noriginal seed complex dG "<<orig_dG<<std::endl;

	TR.Debug <<"start: " << start << "stop : " << end << std::endl;

	for ( core::Size resi=start; resi<=end; ++resi ) {
		core::Real const ala_scan_dG( ala_scan.ddG_for_single_residue( pose, resi ) );
		core::Real const dG( ala_scan_dG - orig_dG );
		TR<<"dG for resi "<<pose.residue( resi ).name3()<<resi<<" is "<<dG<<std::endl;
		if ( dG < lowest_dG ) {
			lowest_dG = dG;
			lowest_res = resi;
		}
	}
	TR<<"seed residue with lowest dG based on ala scanning is : "<< lowest_res << " with dG of: " << lowest_dG << std::endl;
	return lowest_res;
}

core::kinematics::FoldTreeOP
SeedFoldTree::set_foldtree(
	pose::PoseOP & target_seed_pose,
	std::string secstr,
	protocols::loops::Loops & loops,
	bool protein_not_folded_yet ){

	using namespace core;
	using namespace kinematics;
	using namespace protocols::seeded_abinitio;

	fold_tree_ = core::kinematics::FoldTreeOP( new core::kinematics::FoldTree );
	fold_tree_->clear();
	Size seed_num = loops.size();

	if ( target_seed_pose->conformation().num_chains() == 2 ) {

		TR<<"two chains were were submitted for the seed pdb, reading target info"<< std::endl;

		target_chain_ = target_seed_pose->split_by_chain( 1 );
		TR<<"input pdb: "<< secstr.length() <<" target chain: " <<target_chain_->size() << std::endl;
		seeds_only_ = target_seed_pose->split_by_chain( 2 );

		Size rb_jump =1;
		Size target_length = target_chain_->size();
		/// this one needs better assertions....
		Size total_size_complex = secstr.length() + target_length;
		Size start_new_protein = target_length + 1 ;
		Size seed_start;// indeces for seeds in the complete pose
		Size seed_stop;
		Size pdb_start; // indeces for seeds in trunctated version
		Size pdb_stop;
		Size position_adjustment = 0;//1
		Size seed_res_counter = 1;
		std::set< core::Size > res_on_target;
		std::set< core::Size > res_on_design;
		utility::vector1 < std::pair < Size, Size > > jump_pair_collection ;

		//informing vertex container with the total length of the new protein
		folding_vertices_.insert( total_size_complex );
		folding_vertices_.insert( start_new_protein );

		if ( seed_num <= 0 ) {
			utility_exit_with_message( "NO SEEDS SPECIFIED!!!" );
		}

		//if there is more than one seed, an extra cutpoint is necessary
		if ( seed_num > 1 ) {
			//need to make sure that the cutpoints are not the same as the seed starts and stops below
			//will cause problem with the set container and the growing peptides mover
			TR<<"finding cutpoints..."<<std::endl;

			for ( Size seed_it = 2 ; seed_it <= seed_num ; ++seed_it ) {
				TR.Debug <<"loops[seed_it - 1].stop()"<< loops[seed_it - 1].stop()<<" and loops[seed_it - 1].stop()" <<loops[seed_it - 1].stop() << std::endl;
				Size end_seed = loops[seed_it - 1].stop();
				Size start_new_seed = loops[seed_it].start();
				TR<<"... between: "<< end_seed << " and " << start_new_seed  << std::endl;
				Size cut = define_cut_point_stochasticly( end_seed, start_new_seed, secstr , start_new_protein);
				cut_points_.push_back( cut );
				folding_vertices_.insert( cut );
				folding_vertices_.insert( cut + 1 );
				TR<<"vector: "<<cut_points_[seed_it - 1] << "method cut: " << cut << std::endl;
			}
		}//end additional cutpoints


		for ( Size seed_it = 1 ; seed_it <= seed_num ; ++seed_it ) {

			std::pair<Size,Size> jump_pair;

			if ( protein_not_folded_yet ) {

				TR << "assuming that the pose does not have full length yet" << std::endl;
				TR.Debug<<"seed_res_counter "<< seed_res_counter << std::endl;
				TR.Debug<<"seed start " << loops[seed_it].start() + target_length << std::endl;
				TR.Debug<<"seed stop " << loops[seed_it].stop() + target_length << std::endl;

				//values for seeds
				seed_start = loops[seed_it].start() + target_length;
				seed_stop  = loops[seed_it].stop() + target_length;

				//adjusting values for truncated version from the input pdb
				pdb_start = seed_res_counter + target_length;
				pdb_stop =  pdb_start + loops[seed_it].stop() - loops[seed_it].start(); //// - 1;

				TR.Debug<<"pdb_start for seed " << pdb_start << std::endl;
				TR.Debug<<"pdb_stop for seed: " << pdb_stop <<  std::endl;

				TR<<"seed_residue_counter: "<<seed_res_counter<<std::endl;
				//populating container for peptide growth
				folding_vertices_.insert(seed_start);
				folding_vertices_.insert(seed_stop);

				TR<<"numbering for seed(only)-target pose with target\n ----- SEED: "<< seed_it << " start: "<<seed_start<<", stop: "<<seed_stop<<" ---------" <<std::endl;
				TR<<"position adjustment of TRUNCATED seed motif by "<<position_adjustment<<std::endl;

				if ( anchor_specified() ) {
					if ( anchors_.size() < 1 ) {
						utility_exit_with_message("no anchor specified?!");
					}
					Size adjust_anchor = anchors_[seed_it] - loops[seed_it].start()  + seed_res_counter;
					TR.Debug << "anchor defined: "<< anchors_[seed_it]<< ", adjusted to " << adjust_anchor << std::endl;
					jump_pair = find_nearest_residue( adjust_anchor, target_seed_pose );
					jump_pair.second += loops[seed_it].start()  - seed_res_counter;
					TR<< "jump pairs: " << jump_pair.first << " " << jump_pair.second << std::endl;
				}

				if ( !anchor_specified() ) {
					if ( ddg_based_ ) {
						TR<<"computing dG for seed " << seed_it <<" to identify jump atom "<< std::endl;
						Size seed_jump_residue = best_by_ala_scan( pdb_start, pdb_stop, target_seed_pose );
						jump_pair = find_nearest_residue( seed_jump_residue, target_seed_pose );
					}

					if ( !ddg_based_ ) {
						//get cloesest pairs between target and seeds to set the jumps
						jump_pair = get_closest_residue_pair( pdb_start, pdb_stop, target_seed_pose );
					}
					TR.Debug<<"loops[seed_it].start(): "<<loops[seed_it].start() << " seed_res_counter " << seed_res_counter << std::endl;

					position_adjustment = loops[seed_it].start() - seed_res_counter; //-1
					jump_pair.second += position_adjustment;
				}

				seed_res_counter += loops[seed_it].stop() - loops[seed_it].start() + 1;
				TR<<"updating seed res counter: "<<seed_res_counter<<std::endl;

			} else { //end not folded yet
				//when the protein is already folded!
				TR << "assuming pose has its full length" <<std::endl;
				seed_start = loops[seed_it].start()+ target_length;
				seed_stop = loops[seed_it].stop() + target_length;
				TR<<"--------- SEED: " << seed_start <<" " << seed_stop << std::endl;
				//reset the total size of the complex:
				//total_size_complex = template_pose.size();
				TR<<"total size complex: " << total_size_complex << std::endl;
				//pose::PoseOP templ = new pose::Pose( template_pose );

				/// setting jump pairs:
				if ( anchor_specified_ ) {
					if ( anchors_[ seed_it ] == 0 ) {
						jump_pair = get_closest_residue_pair( seed_start, seed_stop, target_seed_pose );
					} else {
						jump_pair = find_nearest_residue( anchors_[seed_it], target_seed_pose );
					}
				} else {
					if ( ddg_based_ ) {
						TR<<"computing dG for seed " << seed_it <<" to identify jump atom "<< std::endl;
						Size seed_jump_residue = best_by_ala_scan( seed_start, seed_stop, target_seed_pose );
						jump_pair = find_nearest_residue( seed_jump_residue, target_seed_pose );
					} else {
						jump_pair = get_closest_residue_pair( seed_start, seed_stop, target_seed_pose );
					}
				}
			}

			rb_jump = seed_it;
			TR<<"finding closest opposing residues pair for seed starting IDs: " << seed_start <<" " <<seed_stop<<" as jump: "<<rb_jump<<std::endl;
			TR.Debug<<"after adjustment "<< jump_pair.second << std::endl;

			res_on_target.insert( jump_pair.first );
			res_on_design.insert( jump_pair.second );
			fold_tree_->add_edge( jump_pair.first, jump_pair.second, rb_jump );
			//cant connected the edges yet since it is going to be unordered. need to do that in a second loop :(
			TR<<"SEED: "<<seed_it <<", seed and target jump pairs: " <<jump_pair.first<<" and " <<jump_pair.second <<"\n";
			TR<<"registered jump pairs after full-length adjustments: " << position_adjustment <<std::endl;

		}//end jump loop

		TR<<"new SeedFoldTree jumps: " <<*fold_tree_ << std::endl;

		//now connect folding pieces/chunks

		///target first:
		Size target_head = 1;
		for ( core::Size const res : res_on_target ) {
			/// connect chain1 with no breaks
			fold_tree_->add_edge( target_head, res, Edge::PEPTIDE );
			target_head = res;
		}

		/// connect the last anchor residue on the target chain with the last residue on the target chain
		core::Size const target_lastjump( *res_on_target.rbegin() );
		fold_tree_->add_edge( target_lastjump, target_length , Edge::PEPTIDE );

		///connecting jump positions on fold pose with cutpoints
		Size cut_iter = 1;
		Size last_cut = 0;
		//Size last_jump = 0;

		for ( core::Size const jpos : res_on_design ) {
			TR<<"foldpose iterator: "<< jpos <<"and +1 "<< jpos+1 << std::endl;
			if ( last_cut != 0 ) {
				fold_tree_->add_edge( last_cut+1 , jpos , Edge::PEPTIDE ); //this way the cut is after the specified cut position, should it be before?
			}
			if ( jpos != *res_on_design.rbegin() ) {
				fold_tree_->add_edge( jpos, cut_points_[cut_iter], Edge::PEPTIDE );
				last_cut = cut_points_[cut_iter];
				++cut_iter;
			} else break;
		}

		/// refold chain2 from the first jump residue to the beginning of the chain and from the last key residue to its end
		core::Size const first( *res_on_design.begin() );
		core::Size const last(  *res_on_design.rbegin() );
		TR<<"first: " << *res_on_design.begin() << " and " << last << std::endl;

		//N-terminus of the new protein
		if ( first - 1 >= start_new_protein ) {
			fold_tree_->add_edge( first, start_new_protein , Edge::PEPTIDE );
		}
		//C-terminus of the new protein
		if ( last < total_size_complex ) {
			TR << "--- total_size_complex: " << total_size_complex << std::endl;
			fold_tree_->add_edge( last, total_size_complex, Edge::PEPTIDE );
		}

		TR <<"before deleting self edges: " << *fold_tree_ << std::endl;
		fold_tree_->delete_self_edges();
		TR<<"before reordering: " << *fold_tree_ <<std::endl;
		fold_tree_->reorder( 1 );

		TR<<"Fold tree:\n"<<*fold_tree_<<std::endl;

	} else if ( target_seed_pose->conformation().num_chains() == 1 ) { //end set_foldtree with target
		//////////////////////// foldtree set up without target chain ///////////////////////////////////
		/////this part was taken from bcorreia's code and hasnt been tested yet or integrated
		//this should be also activatable and general in case the usere wants to fold in the absence of the target

		TR<<"there is no target chain, either because you turned off the option, or it was not loaded" <<std::endl;

		if ( loops.size() > 1 ) {
			using namespace core;
			using namespace kinematics;
			//utility::vector1<Size> cut_points_;
			Size start_protein = 1;//silly type conversion

			TR<<"more than one seed is defined " << std::endl;
			for ( Size seed_it = 2 ; seed_it <= seed_num; ++seed_it ) {
				Size starting = loops[seed_it - 1].stop();
				Size ending = loops[seed_it].start();
				Size cut = define_cut_point_stochasticly( starting, ending, target_seed_pose->secstruct(), start_protein /*or 0? start fold pose */ );
				TR<<"adding cut: " << cut <<std::endl;
				cut_points_.push_back( cut );
				fold_tree_->add_edge( 1, target_seed_pose->size(), Edge::PEPTIDE );
			}

			for ( Size i=1;  i < loops.size() ; ++i ) {
				TR<<"seed "<<  i <<std::endl;
				TR<<"cut point"<< cut_points_[i]<<std::endl;
				Size cutpoint = cut_points_[i];
				core::pose::add_variant_type_to_pose_residue(*target_seed_pose, chemical::CUTPOINT_LOWER, cutpoint); // residue on the pose has to be assigned as a cut
				core::pose::add_variant_type_to_pose_residue(*target_seed_pose, chemical::CUTPOINT_UPPER, cutpoint+1);
				Size loop1_midpoint = ((loops[1].stop()-loops[1].start())/2) + loops[1].start();
				TR<<"loop 1 mid point"<< loop1_midpoint<<std::endl;
				Size variable_midpoint = ((loops[i+1].stop()-loops[i+1].start())/2) + loops[i+1].start();
				TR<<"Variable mid_point"<< variable_midpoint <<std::endl;
				fold_tree_->new_jump( loop1_midpoint, variable_midpoint, cutpoint );
			}
			TR << "Fold Tree for the scaffold " << *fold_tree_ << std::endl;

			return fold_tree_;

		}//end more than 1 seed
	} else { //end without target section
		TR<<"no special foldtree needed. There is no target chain addition and less then 2 or no seed defined"<<std::endl;
		fold_tree_->add_edge( 1, target_seed_pose->size(), Edge::PEPTIDE );
		TR << "Pose fold tree " << fold_tree_ << std::endl;
	}

	return fold_tree_;

}//end set foldtree method

/*
//needs some fixes... (todo)
core::kinematics::FoldTreeOP
SeedFoldTree::set_foldtree_manually(
utility::vector1 < std::pair < Size, Size > > & manual_jump_pairs,
core::pose::Pose & pose
){

using namespace core;
using namespace kinematics;
using namespace protocols::seeded_abinitio;

fold_tree_ = new core::kinematics::FoldTree;
fold_tree_->clear();

if( pose.conformation().num_chains() == 2 ){
//Size start_design = pose.conformation().chain_begin( 2 );

cut_points_.push_back( pose.conformation().chain_end( 1 ) );

//if there is more than one jump defined we need cutpoints
if( manual_jump_pairs.size() > 1 ){

for(Size it = 2 ; it <=manual_jump_pairs.size(); ++it){
if( manual_jump_pairs.size() > 1 ){
Size end_last_jump  = manual_jump_pairs[it-1].second;
Size start_new_jump = manual_jump_pairs[it].second;
Size cut = define_cut_point_stochasticly( end_last_jump, start_new_jump, pose, 0 );
cut_points_.push_back( cut );
}
}
}

fold_tree_->add_edge( 1 , pose.conformation().chain_end( 1 ), Edge::PEPTIDE );
fold_tree_->add_edge(pose.conformation().chain_begin( 2 ), pose.conformation().chain_end( 2 ), Edge::PEPTIDE );
TR<<"foldtree before adding extra jumps: "<<*fold_tree_ <<std::endl;

for( Size jpair_it = 1; jpair_it <= manual_jump_pairs.size(); ++jpair_it ){
fold_tree_->new_jump( manual_jump_pairs[jpair_it].first, manual_jump_pairs[jpair_it].second, cut_points_[jpair_it]);
TR<<"adding jump between residues: "<< manual_jump_pairs[jpair_it].first <<" and " << manual_jump_pairs[jpair_it].second <<" with cut " << cut_points_[jpair_it] << std::endl;
}

TR<<"before deleting self edges: " << *fold_tree_ << std::endl;
fold_tree_->delete_self_edges();
TR<<"before reordering: " << *fold_tree_ <<std::endl;
fold_tree_->reorder( 1 );

}

else {
utility_exit_with_message( "manual jumps are currently only set up for 2 chain targets, but super easy to set up for just 1" );
}

return fold_tree_;

}//end manual jump setfoldtree

*/

///////apply///////////
void
SeedFoldTree::apply( core::pose::Pose & pose )
{

	/* (todo)
	if (set_jumps_manually )
	fold_tree_ = set_foldtree_manually( manual_jump_pairs_, input_pose );
	*/

	bool protein_not_folded = true;
	Size chain_num = pose.conformation().num_chains();

	//if last chain and template pose have the same length, then the protein is at its full length
	if ( pose.split_by_chain( chain_num )->size() ==  template_pdb_->size() ) {
		protein_not_folded = false;
		TR<<"assuming pose has full size" << std::endl;
	}

	if ( chain_num <= 2 ) {
		TR<<"Previous fold tree: "<< pose.fold_tree()<<'\n';
		TR<<"reseting foldtree"<<std::endl;
		pose::PoseOP poseOP( new pose::Pose( pose ) );
		fold_tree_ = set_foldtree( poseOP, template_pdb_->secstruct(), all_seeds_ , protein_not_folded );
		//fold_tree_ = set_foldtree( pose, template_pdb_ , all_seeds_ , protein_not_folded );
	}

	if ( pose.conformation().num_chains() > 2 ) {
		utility_exit_with_message( "more than 2 chains as input are currently not supported" );
	}

	runtime_assert( fold_tree_ != nullptr );

	TR<<"Previous fold tree: "<< pose.fold_tree()<<'\n';
	pose.fold_tree( *fold_tree_ );
	TR<<"New fold tree: "<< pose.fold_tree()<<std::endl;
	protocols::loops::add_cutpoint_variants( pose );
	TR.flush();
}

utility::vector1 < core::Size >
SeedFoldTree::get_cutpoints(){ return cut_points_ ;}

std::set< core::Size >
SeedFoldTree::get_folding_vertices(){ return folding_vertices_;}

// XRW TEMP std::string
// XRW TEMP SeedFoldTree::get_name() const {
// XRW TEMP  return SeedFoldTree::mover_name();
// XRW TEMP }


void
SeedFoldTree::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	Movers_map const &,
	Pose const & /*input_pose*/){

	TR<<"SeedFoldTree has been invoked"<<std::endl;

	ddg_based_ = tag->getOption< bool >( "ddG_based", 0 );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );

	//parsing branch tags
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );

	for ( TagCOP const btag : branch_tags ) {
		/* this parsing option works, it is just not hooked in yet
		//in case anybody ever wanted to set them manually
		if( btag->getName() == "cut_points" ) {

		utility::vector1< Size > cut_points_;
		core::Size const resnum( protocols::rosetta_scripts::get_resnum( btag, input_pose ) );
		cut_points_.push_back( resnum );
		TR<< "adding cut point: "<< resnum << std::endl;
		}//end cut points
		*/

		anchor_specified_ = false;

		if ( btag->getName() == "Seeds" ) { //need an assertion for the presence of these or at least for the option file

			core::Size const begin( btag->getOption<core::Size>( "begin", 0 ) );
			core::Size const end( btag->getOption<core::Size>( "end", 0 ) );
			all_seeds_.add_loop( begin , end , 0, 0, false );
			if ( btag->hasOption( "anchor" ) ) {
				Size anchor_res = btag->getOption< core::Size >("anchor", 0 );
				TR<<"anchor residue: " << anchor_res << std::endl;
				anchors_.push_back( anchor_res );
				anchor_specified_ = true;
			} else {
				anchors_.push_back( 0 );
			}
		}//end seed tags


		//add option for manual setting of the jumps. this optioh should turn off the automatic foldtree
		if ( btag->getName() == "Jumps" ) { //need an assertion for the presence of these or at least for the option file
			set_jumps_manually = true;
			std::pair< Size, Size > jump_pair;
			jump_pair.first  = btag->getOption<core::Size>( "from", 0 ) ;
			jump_pair.second = ( btag->getOption<core::Size>( "to", 0 ) );
			if ( jump_pair.first > jump_pair.second ) {
				utility_exit_with_message("specifiied jumps need to be defined in sequence order" );
			}
			manual_jump_pairs_.push_back( jump_pair );

		}//end jump tags
	}//end b-tags

	std::string const template_pdb_fname( tag->getOption< std::string >( "template_pdb" ));
	template_pdb_ = core::pose::PoseOP( new core::pose::Pose ) ;
	core::import_pose::pose_from_file( *template_pdb_, template_pdb_fname , core::import_pose::PDB_file);

}//end parse my tag

std::string SeedFoldTree::get_name() const {
	return mover_name();
}

std::string SeedFoldTree::mover_name() {
	return "SeedFoldTree";
}

void SeedFoldTree::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("ddg_based", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::required_attribute("template_pdb", xs_string, "XRW TO DO");

	// Subelements
	AttributeList subelement_attributes;
	subelement_attributes
		+ XMLSchemaAttribute::required_attribute("begin", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute::required_attribute("end", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute::attribute_w_default("anchor", xsct_non_negative_integer, "XRW TO DO", "0");

	AttributeList jump_attributes;
	jump_attributes
		+ XMLSchemaAttribute::attribute_w_default("from", xsct_non_negative_integer, "XRW TO DO","0")
		+ XMLSchemaAttribute::attribute_w_default("to", xsct_non_negative_integer, "XRW TO DO","0");


	XMLSchemaSimpleSubelementList subelement_list;
	subelement_list.add_simple_subelement("Seeds", subelement_attributes, "XRW TO DO");
	subelement_list.add_simple_subelement("Jump", jump_attributes, "XRW TO DO");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, subelement_list );
}

std::string SeedFoldTreeCreator::keyname() const {
	return SeedFoldTree::mover_name();
}

protocols::moves::MoverOP
SeedFoldTreeCreator::create_mover() const {
	return protocols::moves::MoverOP( new SeedFoldTree );
}

void SeedFoldTreeCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SeedFoldTree::provide_xml_schema( xsd );
}

}//end seeded_abinitio
}//end protocols
