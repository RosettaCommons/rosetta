// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/public/scenarios/beta_strand_homodimer_design/exposed_strand_finder.cc
/// @brief  finds exposed strands with option to check rmsd to another strand

// Unit headers

// devel headers
#include <devel/init.hh>
// Utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/vector1.functions.hh>
#include <utility/vector1.hh>
// core headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
//#include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>
//#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/EtableEnergy.hh>
//#include <core/scoring/etable/CoarseEtableEnergy.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>

//Numeric headers
#include <numeric/xyzVector.hh>

//Basic headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

//Protocols headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/StructureRestrictor.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

//protocols JD2
#include <protocols/jd2/JobDistributor.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <utility>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <utility/vector0.hh>


static THREAD_LOCAL basic::Tracer TR( "apps.public.beta_strand_homodimer_design.exposed_strand_finder" );
static THREAD_LOCAL basic::Tracer TRout( "ExposedStrand" );

using namespace core;
using namespace utility;
using namespace protocols;
using namespace protocols::moves;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//  Index for location of tokens in a defined mutation
enum Index{
	chain_index = 1,
	start_index = 2,
	end_index = 3,
	max_neighbors = 16
};

// application specific options
basic::options::BooleanOptionKey const check_rmsd( "check_rmsd" );
basic::options::StringVectorOptionKey const strand_span( "strand_span" ); //format A 5 9
basic::options::IntegerOptionKey const beta_length( "beta_length" ); // min beta lenght to consider
basic::options::IntegerOptionKey const sat_allow( "sat_allow" ); // max allowed hbonds in i+2...
basic::options::RealOptionKey const max_E_allow( "max_E_allow" ); //max energy for bb clashes
basic::options::RealOptionKey const max_RMSD( "max_RMSD" ); //max RMSD allowed

basic::options::StringOptionKey const struct_file( "struct_file" ); //filename of structure restrictor

// mover deffinition
class ExposedStrandMover : public Mover {
public:

	ExposedStrandMover();

	virtual void apply( core::pose::Pose& pose );

	virtual void print_sheets(
		core::pose::Pose & pose,
		core::Size & start_sheet,
		core::Size & end_sheet,
		core::Size & nhbonds);

	virtual void print_sheets_extras(
		core::pose::Pose & pose,
		core::pose::Pose & native_pose,
		core::Size & start_sheet,
		core::Size & end_sheet,
		core::Size & nhbonds,
		Real & rmsd,
		Real & energy,
		Size & match_res);

	virtual void parse_strand_ids (core::pose::Pose & pose, utility::vector1<std::string> & strand_id );

	virtual Real bb_score(pose::Pose & pose, core::Size aligned_chain_num, core::scoring::ScoreFunctionOP & scorefxn);

	virtual pose::Pose move_superimpose(core::pose::Pose & pose1, core::pose::Pose & pose2,
		core::Size & start_res1, core::Size & start_res2,
		core::Size & end_res1, core::Size & end_res2 );


	virtual Real bb_rmsd( const core::pose::Pose & pose1,
		const core::pose::Pose & pose2 );

	virtual bool is_exposed( pose::Pose & pose, Size & resid, vector1< Real > sasa_values );

	virtual MoverOP clone() const {
		return MoverOP( new ExposedStrandMover( *this ) );
	}

	virtual std::string get_name() const{
		return "ExposedStrandMover";
	}

	virtual MoverOP fresh_instance() const {
		return clone();
	}

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	bool check_rmsd_;
	vector1< std::string > strand_def_vector_ ;
	vector1< core::Size > strand_pose_nums_ ; //contains start and end pose_ids
	core::pose::Pose native_pose_ ;
	Size beta_length_ ;
	Size sat_allow_ ;
	Size num_satisfied_;
	core::Real maxE_;
	core::Real maxRMSD_;
	char chain_char_;
	std::string struct_filename_;
	vector1<core::Real> best_rmsd_values_;
	utility::vector1< core::Size > res_to_loose_;
	vector1<pose::Pose> master_poses_;
	vector1<Real> full_scores_;
	vector1<Size> start_res_list_, end_res_list_;
};

ExposedStrandMover::ExposedStrandMover() {
	beta_length_ = option[ beta_length ];
	sat_allow_= option[ sat_allow ];
	check_rmsd_ = option[ check_rmsd ];
	//strand_def_vector_ = option[ strand_span ];
	scorefxn_ = core::scoring::get_score_function();
	// AMW: cppcheck wanted this to be initialized in the ctor
	chain_char_ = 'A';
	//struct_filename_ = option[ struct_file ];
}

void ExposedStrandMover::parse_strand_ids(core::pose::Pose & pose, utility::vector1<std::string> & strand_id ){
	chain_char_ = strand_id[chain_index][0] ;

	int start_num, end_num;
	std::istringstream ss_start( strand_id[ start_index ] );
	std::istringstream ss_end( strand_id[ end_index ] );
	ss_start >> start_num;
	ss_end >> end_num;
	strand_pose_nums_.push_back( pose.pdb_info()->pdb2pose(chain_char_, start_num ) );
	strand_pose_nums_.push_back( pose.pdb_info()->pdb2pose(chain_char_, end_num ) );

	TR << "Strand to match against is chain: " << chain_char_ << " poseid start: "
		<< strand_pose_nums_[1] << " poseid end: "<< strand_pose_nums_[2] <<std::endl;

}//end parse_strand_ids

//print out exposed sheets
void ExposedStrandMover::print_sheets(
	core::pose::Pose & pose,
	core::Size & start_sheet,
	core::Size & end_sheet,
	core::Size & nhbonds)
{
	using namespace std;
	utility::file::FileName filename (pose.pdb_info()->name());
	TRout << "FILE:  " << setw(15) << filename.base() << "   "
		<< "CHAIN:  " << setw(3) << pose.pdb_info()->chain(start_sheet) << "   "
		<< "START:  " << setw(4) << pose.pdb_info()->number(start_sheet) << "   "
		<< "END:  " << setw(4)  <<  pose.pdb_info()->number(end_sheet) << "   "
		<< "H_BONDS:  " << setw(2) << nhbonds << endl;
}// end print_sheets

//print out exposed sheets
void ExposedStrandMover::print_sheets_extras(
	core::pose::Pose & pose,
	core::pose::Pose & native_pose,
	core::Size & start_sheet,
	core::Size & end_sheet,
	core::Size & nhbonds,
	Real & rmsd,
	Real & energy,
	Size & match_res)
{
	using namespace std;
	utility::file::FileName filename (pose.pdb_info()->name());
	TRout << "FILE:  " << setw(15) << filename.base() << "   "
		<< "CHAIN:  " << setw(3) << pose.pdb_info()->chain(start_sheet) << "   "
		<< "START:  " << setw(4) << pose.pdb_info()->number(start_sheet) << "   "
		<< "END:  " << setw(4)  <<  pose.pdb_info()->number(end_sheet) << "   "
		<< "H_BONDS:  " << setw(2) << nhbonds << "   "
		<< "bb_RMSD:  " << setw(5) << setprecision(3) << fixed << rmsd << "   "
		<< "Chains_E:  "<< setw(8) << setprecision(2)<< fixed << energy<< "   "
		<< "MATCH_RES:   "  << setw(4) << native_pose.pdb_info()->number(match_res)
		<< native_pose.pdb_info()->chain(match_res)<< "   "
		<< endl;
}// end print_sheets


//////////////////////////////////////////////////////////////////

// these are the same size poses always
core::Real ExposedStrandMover::bb_rmsd(const pose::Pose & pose1 ,const pose::Pose & pose2){
	Real rms ( scoring::rmsd_with_super( pose1, pose2,  scoring::is_protein_backbone ) );
	return rms;
}

///////////////////////////////////////////////////////////////////////////////////
pose::Pose
ExposedStrandMover::move_superimpose(core::pose::Pose & pose1 /*input*/, core::pose::Pose & pose2 /*master*/,
	core::Size & start_res1, core::Size & start_res2,
	core::Size & end_res1, core::Size & end_res2){

	Size length(end_res1-start_res1);
	runtime_assert(end_res2-start_res2 == length);
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	//some initialization for each new run
	// maps every atomid to bogus atom
	atom_map.clear();
	core::pose::initialize_atomid_map( atom_map, pose1, core::id::BOGUS_ATOM_ID );

	for ( core::Size ii = 0; ii < length; ++ii ) {
		//Sloppy way of adding all bb atoms, must be a better way...
		core::id::AtomID const id1( pose1.residue( start_res1 + ii).atom_index("CA"), start_res1 + ii );
		core::id::AtomID const id2( pose2.residue( start_res2 + ii).atom_index("CA"), start_res2 + ii );
		atom_map[ id1 ] = id2;

		core::id::AtomID const id3( pose1.residue( start_res1 + ii).atom_index("C"),  start_res1 + ii );
		core::id::AtomID const id4( pose2.residue( start_res2 + ii).atom_index("C"),  start_res2 + ii );
		atom_map[ id3 ] = id4;

		core::id::AtomID const id5( pose1.residue( start_res1 + ii).atom_index("N"),  start_res1 + ii );
		core::id::AtomID const id6( pose2.residue( start_res2 + ii).atom_index("N"),  start_res2 + ii );
		atom_map[ id5 ] = id6;

		core::id::AtomID const id7( pose1.residue( start_res1 + ii).atom_index("O"),  start_res1 + ii );
		core::id::AtomID const id8( pose2.residue( start_res2 + ii).atom_index("O"),  start_res2 + ii );
		atom_map[ id7 ] = id8;
	}
	//do superimpose
	scoring::superimpose_pose(pose1, pose2, atom_map);

	//now add these two poses together
	pose2.append_residue_by_jump(pose1.residue( 1 ), pose2.total_residue(), "" , "",  true /*start new chain*/);
	for ( core::Size n = 2; n<= pose1.total_residue(); ++n ) {
		//add into a new chain if need be
		if ( pose1.residue(n).chain() != pose1.residue(n-1).chain() ) {
			pose2.append_residue_by_jump(pose1.residue( n ), pose2.total_residue(), "" , "",  true /*start new chain*/);
		} else { /*if (( pose1.residue(n-1).is_polymer() && !pose1.residue(n-1).is_upper_terminus() ) &&
			(    pose1.residue(n).is_polymer() &&    !pose1.residue(n).is_lower_terminus() ) )*/
			pose2.append_residue_by_bond( pose1.residue ( n ) );
		}
	} //end bringing poses together into pose2

	//return the master
	return pose2;

}//end rms_and_score

///////////////////////////////////////
// bb score
///////////////////////////////////////
core::Real ExposedStrandMover::bb_score(pose::Pose & pose, core::Size aligned_chain_num, core::scoring::ScoreFunctionOP & scorefxn){

	// score the bb-bb energy between chains
	// This part written by P.Doug Renfrew
	// make vectors of bb atoms in each chain individually
	// the master pose will always be chain 1.
	// need to make a vector of all atoms in the chain you are comparing too

	utility::vector1<core::conformation::Atom> chain1_bb_atoms;
	utility::vector1<core::conformation::Atom> chain2_bb_atoms;
	utility::vector1<core::conformation::Atom> all_bb_atoms;

	for ( Size j = 1; j <= pose.total_residue(); ++j ) {
		core::conformation::Residue const & res( pose.residue(j) );
		core::chemical::AtomIndices bb_ai( res.mainchain_atoms() );
		//assert( bb_ai.size() == 4 );
		core::Size chain_num( res.chain() );
		for ( Size jj = 1; jj <= bb_ai.size(); ++jj ) {
			if ( chain_num == 1 ) {
				chain1_bb_atoms.push_back( res.atom(jj) );
			} else if ( chain_num == aligned_chain_num ) {
				chain2_bb_atoms.push_back( res.atom(jj) );
			}
			//optional get all the atoms not in allinged chain,
			//only need to do if more than two chains in pose
			if ( pose.conformation().num_chains() >= 3 && chain_num != 1 ) {
				all_bb_atoms.push_back( res.atom(jj) );
			}
			//end optional
		}
	}

	//NOW SCORE!
	// get instance of etable energy method
	core::scoring::methods::EnergyMethodOptions const & emo(scorefxn->energy_method_options());
	core::scoring::etable::Etable const & et(*(core::scoring::ScoringManager::get_instance()->etable(emo).lock()));
	core::scoring::etable::TableLookupEtableEnergy ete( et, emo, false /*do_classic_intrares*/ );

	// iterate over both sets of atom and add into one emapvector
	//core::scoring::TwoBodyEMapVector tbemv;
	core::scoring::EMapVector tbemv;
	core::Real atr_wt( (*scorefxn).get_weight(core::scoring::fa_atr) );
	core::Real rep_wt( (*scorefxn).get_weight(core::scoring::fa_rep) );
	for ( Size ii = 1; ii <= chain1_bb_atoms.size(); ++ii ) {
		for ( Size jj = 1; jj <= chain2_bb_atoms.size(); ++jj ) {
			//calc distr squared
			Real d2( chain1_bb_atoms[ii].xyz().distance_squared( chain2_bb_atoms[jj].xyz() ) );
			ete.atom_pair_energy( chain1_bb_atoms[ii], chain2_bb_atoms[jj], 1, tbemv, d2 );
		}
	}
	core::Real bb_energy (rep_wt * tbemv[core::scoring::fa_rep] + atr_wt * tbemv[core::scoring::fa_atr] );

	// begin optional  ie skip if not needed
	core::Real all_energy;
	if ( pose.conformation().num_chains() >= 3 ) {
		//core::scoring::TwoBodyEMapVector tbemv_all;
		core::scoring::EMapVector tbemv_all;
		core::scoring::etable::TableLookupEtableEnergy ete_all( et, emo, false /*do_classic_intrares*/ );
		for ( Size ii = 1; ii <= chain1_bb_atoms.size(); ++ii ) {
			for ( Size jj = 1; jj <= all_bb_atoms.size(); ++jj ) {
				//calc distr squared
				Real d2_all( chain1_bb_atoms[ii].xyz().distance_squared(  all_bb_atoms[jj].xyz() ) );
				ete.atom_pair_energy( chain1_bb_atoms[ii], all_bb_atoms[jj], 1, tbemv_all, d2_all );
			}
		}
		all_energy = (rep_wt * tbemv_all[core::scoring::fa_rep] + atr_wt * tbemv_all[core::scoring::fa_atr] );
	} else { //end optional for many chains
		all_energy = bb_energy ;
	}
	TR<< "Number of chains: " <<pose.conformation().num_chains()
		<<"    Backbone-backbone score: " << all_energy << std::endl;
	return all_energy;
	//return bb_energy;
}//end bb_score

////////////////////////////////////////////////
// is this residue exposed function
///////////////////////////////////////////////
bool ExposedStrandMover::is_exposed( pose::Pose & pose, Size & resid, vector1< Real > sasa_values ){
	int neighbor_cutoff ( 16 );
	int neighbors_lots( 30 );
	Real sasa_per_atom_cut(2.0);
	bool exposed;
	core::scoring::TenANeighborGraph tang (pose.energies().tenA_neighbor_graph());
	int num_neighbors( tang.get_node(resid)->num_neighbors_counting_self() );
	if ( num_neighbors > neighbors_lots ) {
		exposed=false;
	} else if ( num_neighbors < neighbor_cutoff ) {
		exposed=true;
	} else if ( (sasa_values[resid] / pose.residue(resid).natoms())  > sasa_per_atom_cut ) {
		exposed=true;
	} else exposed=false;

#ifndef NDEBUG

	TR<<"Residue: " << resid  <<"   #neighbors: " << num_neighbors
		<< "  SASA:  "  << sasa_values[resid]
		<< "  SASA/atom:  "  << sasa_values[resid] / pose.residue(resid).natoms()
		<< "   is_exposed?:  "<< exposed << std::endl;
#endif


	return exposed;
}

////////////////////////////////////////////////
// Mover Apply
////////////////////////////////////////////////
void ExposedStrandMover::apply (core::pose::Pose & pose ) {

	//fill native pose if needed
	if ( basic::options::option[ in::file::native ].user() ) {
		core::import_pose::pose_from_file( native_pose_, basic::options::option[ in::file::native ](), core::import_pose::PDB_file);
		//parse_strand_ids( native_pose_, strand_def_vector_);
	}

	utility::file::FileName posename ( pose.pdb_info()->name() );
	TR << "Working on: " << posename.base() << std::endl;
	if ( basic::options::option[ struct_file ].active() ) {
		struct_filename_ = option[ struct_file ];
		TR << "Deleting according to structure file..."<< std::endl;
		protocols::moves::StructureRestrictorOP restrictor( new protocols::moves::StructureRestrictor(struct_filename_ ) );
		restrictor->apply(pose);
	}

	//score the pose
	(*scorefxn_)( pose );

	//core::scoring::TenANeighborGraph tang (pose.energies().tenA_neighbor_graph());

	//set up sasa calculator
	core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa_calc" ) ) {
		//TR << "sasa calculator already exists...continuing" << std::endl;
	} else {
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa_calc", sasa_calculator );
	}
	basic::MetricValue< vector1<core::Real> > mv_res_sasa;
	pose.metric( "sasa_calc", "residue_sasa", mv_res_sasa);
	vector1< Real > res_sasa( mv_res_sasa.value() );


	//init things needed later
	core::scoring::hbonds::HBondSet hbond_set;
	core::scoring::hbonds::fill_hbond_set(pose,
		false, /*calc deriv*/
		hbond_set,
		false /*bb only*/ );
	hbond_set.setup_for_residue_pair_energies(pose);

	//set up dssp info
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	// #ifndef NDEBUG
	//  //core::scoring::dssp::Dssp new_dssp( pose );
	//  //output to see if dssp info matches pose info
	//  for(Size ii = 1; ii<=pose.total_residue(); ++ii ){
	//   int num_neighbors_debug( tang.get_node(ii)->num_neighbors_counting_self() );
	//   TR << posename.base() <<"  residue: " << ii << "  PoseSS: " << pose.secstruct(ii) << "  num_neighbors: "
	//     << num_neighbors_debug << std::endl;
	//  }
	// #endif

	//set of what begining and end points of a strand are
	vector1< std::pair< Size,Size > > sheet_endpts;
	//where we left off last time through
	core::Size last_res_seen(0);
	//itterate through entire pose
	for (  core::Size ii=1; ii<=pose.total_residue(); ++ii ) {
		//skip ahead if we've already been through here
		if ( ii <= last_res_seen ) {
			continue;
		}

		//skip ahead if this residue is in a bb-bb hbond
		if ( hbond_set.acc_bbg_in_bb_bb_hbond( ii ) && hbond_set.don_bbg_in_bb_bb_hbond( ii ) ) {
			continue;
		}

		//try getting ss info
		char res_ss( pose.secstruct( ii ) ) ;
		if ( res_ss != 'E' ) {
			continue;
		} else { //ie if residue is a sheet
			//how long is the sheet?
			//int num_neighbors( tang.get_node(ii)->num_neighbors_counting_self() );
			Size current_res (ii), end_res(ii);
			Size sat_hbonds(0);


			////////////////////////////
			// itterate through an exposed strand
			////////////////////////////
			while ( pose.secstruct(current_res) == 'E' && sat_hbonds <= sat_allow_ && is_exposed(pose, current_res, res_sasa) ) {
				//are sheet hbonds exposed
				//only check every other residue...
				if ( (current_res - ii) % 2 == 0 ) {
					if ( hbond_set.acc_bbg_in_bb_bb_hbond(current_res) == true ) {
						++sat_hbonds;
					}
					if ( hbond_set.don_bbg_in_bb_bb_hbond(current_res) == true ) {
						++sat_hbonds;
					}
				}
				end_res=current_res;
				++current_res;
				//num_neighbors = tang.get_node(current_res)->num_neighbors_counting_self()
			}  //end while loop

			core::Size sheet_length(end_res - ii + 1);
			//long enough, exposed enough...
			if ( sheet_length >= beta_length_ && sat_hbonds <= sat_allow_ ) {
				last_res_seen = end_res;
				num_satisfied_ = sat_hbonds;
				//fill pair with begining and end of different sheets
				sheet_endpts.push_back( std::make_pair(ii, end_res) );
				//print this sheet info if we are not checking the RMSD
				if ( !check_rmsd_ ) {
					print_sheets(pose, ii, end_res, sat_hbonds);
				}
			} else {
				continue;
			}
		}//end if ss=E
	}//end itterate over all residues

	///////////////////////////////////////////
	//begin the stuff if  RMSD calc is desired
	///////////////////////////////////////////
	if ( check_rmsd_ ) {
		if ( !basic::options::option[ in::file::native ].active() ) {
			utility_exit_with_message_status( "Need to define -native to check RMSD stuff \n", 1 );
			if ( !option[strand_span].active() ) {
				utility_exit_with_message_status( "Need to define -strand_span to check RMSD to. \n", 1 );
			}
		}
		//read inputs
		strand_def_vector_ = option[ strand_span ];
		parse_strand_ids( native_pose_, strand_def_vector_);

		//fill variables
		utility::file::FileName native_name ( native_pose_.pdb_info()->name() );
		maxE_ = option[ max_E_allow ];
		maxRMSD_ = option[ max_RMSD ];

		//find chain to cut out of master_pose later
		//utility::vector1< core::Size >  res_to_loose;
		for ( core::Size ii =1; ii <= native_pose_.total_residue(); ++ii ) {
			if ( native_pose_.pdb_info()->chain(ii) == chain_char_ ) {
				res_to_loose_.push_back(ii);
			}
		}

		//now align and store
		TR << posename.base() <<" has " << sheet_endpts.size() << " sheets to check." << std::endl;

		//setup for itterating and printing
		Size native_start( strand_pose_nums_[1] ), native_end( strand_pose_nums_[2] );
		Size native_sheet_length( native_end - native_start + 1);
		Real rms_value (99.0);
		pose::Pose combined_pose;
		//itterate over all the sheets
		for ( Size ii =1; ii<= sheet_endpts.size(); ++ii ) {
			//start with known sheet between start_res & end_res
			Size pose_start(sheet_endpts[ii].first), pose_end(sheet_endpts[ii].second);
			Size pose_sheet_length(pose_end-pose_start + 1);
			if ( pose_sheet_length < beta_length_ || native_sheet_length < beta_length_ ) {
				utility_exit_with_message_status( "Lengths for pose input are not big enough \n", 1 );
			}

			// look for fragments that are beta_length in the given regions
			Size nfrags_pose   = Size( floor( (pose_sheet_length-beta_length_+2)/2 ) );
			Size nfrags_native =  Size( floor( (native_sheet_length-beta_length_+2)/2 ) );

			Size native_current (native_start), pose_current;
			for ( Size jj = 1; jj<=nfrags_native; ++jj ) {
				pose_current = pose_start;
				for ( Size kk = 1; kk <= nfrags_pose; ++kk ) {
					Size native_current_end (native_current + beta_length_ -1);
					Size pose_current_end (pose_current + beta_length_ -1);
					//reset things to original state
					pose::Pose native_copy (native_pose_), pose_copy(pose);

					////////////////////
					// dummy poses for RMSD calc
					////////////////////
					Pose temp_pose1;
					Pose temp_pose2; //tmp pose to put beta sheets in for rmsd calc
					temp_pose1.clear();
					temp_pose2.clear();
					for ( core::Size k = 0 ; k < beta_length_; ++k ) {
						temp_pose1.append_residue_by_bond(native_copy.residue( native_current + k ) );
						temp_pose2.append_residue_by_bond(pose_copy.residue( pose_current + k ) );
					}
					rms_value = bb_rmsd(temp_pose1, temp_pose2);
					/////////////////////////////////////////

					///Superposition and asorted fixes
					combined_pose = move_superimpose(pose_copy, native_copy,
						pose_current, native_current,
						pose_current_end, native_current_end);

					//now remove the peptide before scoring
					combined_pose.conformation().delete_residue_range_slow( min(res_to_loose_), max(res_to_loose_));
					TR << "Master pose size after peptide removal: " << combined_pose.total_residue() << std::endl;

					//score what's left
					Size new_sheet_start( pose_current + combined_pose.conformation().chain_end(1) ) ;
					Size aligned_chain_num(  combined_pose.residue( new_sheet_start ).chain() );
					TR << "Chain from input used in alignment: " << aligned_chain_num << std::endl;
					Real bb_clash_score ( bb_score(combined_pose, aligned_chain_num, scorefxn_) );


					//Now let's check our filters
					if ( bb_clash_score < maxE_ && rms_value < maxRMSD_ ) {

						std::string full_pose_name;
						std::stringstream blah3;
						blah3 << posename.base() << "_"<<  pose.pdb_info()->chain( pose_current )<<pose.pdb_info()->number( pose_current ) <<"_"
							<< native_name.base()<< "_"<< native_pose_.pdb_info()->chain( native_current )<<native_pose_.pdb_info()->number( native_current )
							<<".pdb";
						full_pose_name=blah3.str();

						//finally dump out info about sheets
						combined_pose.dump_pdb(full_pose_name);

						print_sheets_extras(pose,
							native_pose_,
							pose_current,
							pose_current_end,
							num_satisfied_,
							rms_value,
							bb_clash_score,
							native_current);
					}//check filters
					pose_current=pose_current+2;
				}//end pose fragment itteration
				native_current=native_current+2;
			}//end native fragment itteration

		}//itterate over all open sheets

	}//end ?check_rmsd?

}//end apply


////////////////////////////////////////////////
// Main
////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		//adding app specific options
		option.add( check_rmsd, "Check RMSD of this strand against another." ).def(false);
		option.add( strand_span, "The span which the strand to compare is.  Format: ChainChar startPDBnum endPDBnum" );
		option.add( beta_length, "The min lenght of a beta sheet to count as exposed." ).def( 3 );
		option.add( sat_allow, "The max number of satisfied hbonds to still count as exposed." ).def( 3 );
		option.add( max_E_allow, "The max allowed E between chains." ).def( 20.0 );
		option.add( max_RMSD, "Max RMSD between native sheet and searched sheet" ).def( 1.5 );
		option.add( struct_file, "file with info about chains and such");


		// init
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new ExposedStrandMover ) );

		std::cout << "Done! -------------------------------"<< std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
} //end main
