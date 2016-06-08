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
/// @author
/// @author Christopher Miles (cmiles@uw.edu)

#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <string>
#include <time.h>

//Auto Headers
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
#include <boost/algorithm/string.hpp>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

namespace basic {

utility::vector1< std::string > tag2string;

utility::vector1< clock_t > start_clock( n_prof_tags, 0 );
utility::vector1< double > total_clock( n_prof_tags, 0 );
utility::vector1< int > calls( n_prof_tags, 0 );
utility::vector1< int > bad_calls( n_prof_tags, 0 );

clock_t const SHRINK_FACTOR( 2 );

double const clock_factor( ( (double) SHRINK_FACTOR * 100.0 ) / CLOCKS_PER_SEC );

bool show_time_on_cerr( false );
void show_time( basic::Tracer& tr, std::string const& msg ) {
	using namespace std;
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	char formatted_time[51];
	std::strftime( formatted_time, 50, "%c", timeinfo ); // Will null terminate
	tr.Error << "TIME_STAMP: " << formatted_time << " " << msg << std::endl;
	if ( show_time_on_cerr ) std::cerr << tr.channel() << ": TIME_STAMP: " << formatted_time << " " << msg << std::endl;
}

void setup_tag2string() {
	tag2string.clear();
	for ( int i=1; i<= n_prof_tags; ++i ) {
		tag2string.push_back( ObjexxFCL::string_of(i) );
	}

	// now fill in
	tag2string[ TEST1 ] = "TEST1";
	tag2string[ TEST2 ] = "TEST2";
	tag2string[ TEST3 ] = "TEST3";
	tag2string[ TEST4 ] = "TEST4";
	tag2string[ ATOM_TREE_UPDATE_INTERNAL_COORDS ] = "ATOM_TREE_UPDATE_INTERNAL_COORDS";
	tag2string[ ATOM_TREE_UPDATE_XYZ_COORDS ] = "ATOM_TREE_UPDATE_XYZ_COORDS";
	tag2string[ ROTAMER_TRIALS ] = "ROTAMER_TRIALS";
	tag2string[ PACK_ROTAMERS ] = "PACK_ROTAMERS";
	tag2string[ UPDATE_RESIDUE_NEIGHBORS ] = "UPDATE_RESIDUE_NEIGHBORS";
	tag2string[ UPDATE_RESIDUE_TORSIONS ] = "UPDATE_RESIDUE_TORSIONS";
	tag2string[ UPDATE_RESIDUE_COORDINATES ] = "UPDATE_RESIDUE_COORDINATES";
	tag2string[ SETUP_NBLIST ] = "SETUP_NBLIST";

	tag2string[ SCORE ] = "SCORE";
	tag2string[ SCORE_BEGIN_NOTIFY ] = "SCORE_BEGIN_NOTIFY";
	tag2string[ SCORE_SETUP ] = "SCORE_SETUP";
	tag2string[ SCORE_FINALIZE ] = "SCORE_FINALIZE";
	tag2string[ SCORE_ONEBODY_ENERGIES ] = "SCORE_ONEBODY_ENERGIES";
	tag2string[ SCORE_NEIGHBOR_ENERGIES ] = "SCORE_NEIGHBOR_ENERGIES";
	tag2string[ SCORE_LONG_RANGE_ENERGIES ] = "SCORE_LONG_RANGE_ENERGIES";
	tag2string[ SCORE_DOT ] = "SCORE_DOT";
	tag2string[ SCORE_END_NOTIFY ] = "SCORE_END_NOTIFY";

	tag2string[ VDW_ENERGY ] = "VDW_ENERGY";
	tag2string[ ENERGY_ENVPAIR_POTENTIAL ] ="ENERGY_ENVPAIR_POTENTIAL";
	tag2string[ SECONDARY_STRUCTURE_ENERGY ] ="SECONDARY_STRUCTURE_ENERGY";
	tag2string[ SECONDARY_STRUCTURE_SSPAIR_ENERGY ] ="SECONDARY_STRUCTURE_SSPAIR_ENERGY";
	tag2string[ SECONDARY_STRUCTURE_HSPAIR_ENERGY ] ="SECONDARY_STRUCTURE_HSPAIR_ENERGY";
	tag2string[ SECONDARY_STRUCTURE_SHEETS_FROM_DIMERS_ENERGY ] ="SECONDARY_STRUCTURE_SHEETS_FROM_DIMERS_ENERGY";
	tag2string[ POSE_COPY ] = "POSE_COPY";

	tag2string[ ENERGY_GRAPH_COPY ] = "ENERGY_GRAPH_COPY";
	tag2string[ ENERGIES_COPY ] = "ENERGIES_COPY";

	tag2string[ CONFORMATION_DETECT_DISULF ] = "CONFORMATION_DETECT_DISULF";
	tag2string[ CONFORMATION_FIX_DISULF ] = "CONFORMATION_FIX_DISULF";
	tag2string[ CONFORMATION_COPY ] ="CONFORMATION_COPY";
	tag2string[ CHEMICAL_MAKE_POSE ] = "CHEMICAL_MAKE_POSE";

	tag2string[ CONSTRAINT_SCORE ] = "CONSTRAINT_SCORE";
	tag2string[ CONSTRAINT_SET_COPY ] = "CONSTRAINT_SET_COPY";
	tag2string[ CCD_CLOSE ] = "CCD_CLOSE";
	tag2string[ FUNC ] = "FUNC";
	tag2string[ DFUNC ] = "DFUNC";
	tag2string[ GET_ENERGIES ] = "GET_ENERGIES";
	tag2string[ SIMANNEALING ] = "SIMANNEALING";
	tag2string[ MC_ACCEPT ] = "MC_ACCEPT";
	tag2string[ INSERT_FRAGS ] = "INSERT_FRAGS";
	tag2string[ GB_GET_ALL_BORN_RADII ] = "GB_GET_ALL_BORN_RADII";
	tag2string[ GB_SETUP_FOR_PACKING ] = "GB_SETUP_FOR_PACKING";
	tag2string[ GEN_BORN_ROTAMER_PAIR_ENERGIES ] = "GEN_BORN_ROTAMER_PAIR_ENERGIES";
	tag2string[ GEN_BORN_ROTAMER_BACKGROUND_ENERGIES ] = "GEN_BORN_ROTAMER_BACKGROUND_ENERGIES";
	tag2string[ MULTIPOLE_SETUP ] = "MULTIPOLE_SETUP";
	tag2string[ MULTIPOLE_ENERGIES ] = "MULTIPOLE_ENERGIES";
	tag2string[ FACTS_GET_ALL_BORN_RADII ] = "FACTS_GET_ALL_BORN_RADII";
	tag2string[ FACTS_SETUP_FOR_PACKING ] = "FACTS_SETUP_FOR_PACKING";
	tag2string[ FACTS_ROTAMER_PAIR_ENERGIES ] = "FACTS_ROTAMER_PAIR_ENERGIES";
	tag2string[ FACTS_ROTAMER_BACKGROUND_ENERGIES ] = "FACTS_ROTAMER_BACKGROUND_ENERGIES";
	tag2string[ MINMOVER_APPLY ] = "MINMOVER_APPLY";
	tag2string[ BACKRUB_MOVER ] = "BACKRUB_MOVER";
	tag2string[ FIND_SUGAR_AND_SUITE_FRAGS_I ] = "FIND_SUGAR_AND_SUITE_FRAGS_I";
	tag2string[ FIND_SUGAR_AND_SUITE_FRAGS_II ] = "FIND_SUGAR_AND_SUITE_FRAGS_II";
	tag2string[ MAKE_BASE_PAIR_MOVE ] = "MAKE_BASE_PAIR_MOVE";
	tag2string[ MAKE_BASE_STEP_MOVE ] = "MAKE_BASE_STEP_MOVE";
	tag2string[ TOTAL ] = "TOTAL";

	// abinitio debugging tags
	tag2string[ ABINITIO ] = "ABINITIO";
	tag2string[ STAGE1 ] = "STAGE1";
	tag2string[ STAGE2 ] = "STAGE2";
	tag2string[ STAGE3 ] = "STAGE3";
	tag2string[ STAGE4 ] = "STAGE4";
	tag2string[ STAGE5 ] = "STAGE5";
	tag2string[ FRAGMENT_MOVER ] = "FRAGMENT_MOVER";
	tag2string[ RG ] = "RG";
	tag2string[ RG_LOCAL ] = "RG_LOCAL";
	tag2string[ SEQUENCE_COMPARISON ] = "SEQUENCE_COMPARISON";
	tag2string[ KDTREE_CONSTRUCT] = "KDTREE_CONSTRUCT";
	tag2string[ KDTREE_SEARCH] = "KDTREE_SEARCH";
	tag2string[ CONSTRUCT_DISTANCE_MATRIX] = "CONSTRUCT_DISTANCE_MATRIX";

	tag2string[ JD2 ] = "JD2";
	tag2string[ JD2_OUTPUT ] = "JD2_OUTPUT";
	tag2string[ JD2_SILENT_OUTPUTTER ] = "JD2_SILENT_OUTPUTTER";
	tag2string[ JD2_INIT_MOVER ] ="JD2_INIT_MOVER";
	tag2string[ ARCHIVE_SYNC_BATCHES ] = "ARCHIVE_SYNC_BATCHES";
	tag2string[ ARCHIVE_JOBSCOMPLETE ] = "ARCHIVE_JOBSCOMPLETE";
	tag2string[ ARCHIVE_CRITICAL_JOBSCOMPLETE ] = "ARCHIVE_CRITICAL_JOBSCOMPLETE";
	tag2string[ ARCHIVE_READ_DECOYS ] = "ARCHIVE_READ_DECOYS";
	tag2string[ ARCHIVE_GEN_BATCH ] = "ARCHIVE_GEN_BATCH";
	tag2string[ ARCHIVE_BLOCK_FILE ] = "ARCHIVE_BLOCK_FILE";
	tag2string[ ARCHIVE_FILL_POSE ] = "ARCHIVE_FILL_POSE";
	tag2string[ ARCHIVE_SCORE_POSE ] = "ARCHIVE_SCORE_POSE";
	tag2string[ ARCHIVE_EVALUATORS ] = "ARCHIVE_EVALUATORS";
	tag2string[ CA_RMSD_EVALUATION ] = "CA_RMSD_EVALUATION";
	tag2string[ TRUNCATED_SCORE_EVALUATOR ] ="TRUNCATED_SCORE_EVALUATOR";

	tag2string[ SAVE_ARCHIVE ] = "SAVE_ARCHIVE";

	tag2string[ ARCHIVE_EVAL_DECOYS ] = "ARCHIVE_EVAL_DECOYS";
	tag2string[ SILENT_READ_TAG_TEST ] = "SILENT_READ_TAG_TEST";
	tag2string[ MPI_FILE_BUF ] = "MPI_FILE_BUF";
	tag2string[ MPI_JD2_WAITS_FOR_ARCHIVE ] ="MPI_JD2_WAITS_FOR_ARCHIVE";
	tag2string[ MPI_NOTIFY_ARCHIVE ] = "MPI_NOTIFY_ARCHIVE";
	tag2string[ SAXS ] = "SAXS";

	// Fragmentpicker stuff
	tag2string[ FRAGMENTPICKING_CS_SCORE ] = "FRAGMENTPICKING_CS_SCORE";
	tag2string[ FRAGMENTPICKING_PROFILE_SCORE ] = "FRAGMENTPICKING_PROFILE_SCORE";
	tag2string[ FRAGMENTPICKING_PROFILE_CAHING ] = "FRAGMENTPICKING_PROFILE_CAHING";
	tag2string[ FRAGMENTPICKING_SECONDARY_SCORE ] = "FRAGMENTPICKING_SECONDARY_SCORE";
	tag2string[ FRAGMENTPICKING_READ_VALL ] = "FRAGMENTPICKING_READ_VALL";
	tag2string[ FRAGMENTPICKING ] = "FRAGMENTPICKING";
	tag2string[ FRAGMENTPICKING_CANDIDATES_COLLECTING ] = "FRAGMENTPICKING_CANDIDATES_COLLECTING";
	tag2string[ FRAGMENTPICKING_ATOMPAIR_SCORE ] = "FRAGMENTPICKING_ATOMPAIR_SCORE";
	tag2string[ FRAGMENTPICKING_PHIPSI_SCORE ] = "FRAGMENTPICKING_PHIPSI_SCORE";
	tag2string[ FRAGMENTPICKING_DIHEDRALCONSTR_SCORE ] = "FRAGMENTPICKING_DIHEDRALCONSTR_SCORE";

	tag2string[ MPICANONICALSAMPLING ] = "MPICANONICALSAMPLING";
	tag2string[ MPIPOOLCOMMUNICATION ] = "MPIPOOLCOMMUNICATION";
	tag2string[ MPICOMMCREATION ] =  "MPICOMMCREATION";
	tag2string[ MPIBARRIER ] =  "MPIBARRIER";
	tag2string[ MPIBARRIER_BEGIN ] = "MPIBARRIER_BEGIN";
	tag2string[ MPIBARRIER_END ] = "MPIBARRIER_END";
	tag2string[ MPI_GATHER_BARRIER ] = "MPI_GATHER_BARRIER";
	tag2string[ FARRAY_MANIPULATION ] =  "FARRAY_MANIPULATION";

	tag2string[ MPI_SLAVE_REPORT_NEW_COORDS ] = "MPI_SLAVE_REPORT_NEW_COORDS";
	tag2string[ MPI_SLAVE_REPORT_SIZES ] = "MPI_SLAVE_REPORT_SIZES";

	tag2string[ MPI_SEND_UPDATE ] = "MPI_SEND_UPDATE";
	tag2string[ MPI_SYNC_POOL_DIFF ] = "MPI_SYNC_POOL_DIFF";
	tag2string[ MPI_SEND_ACCEPTED ] = "MPI_SEND_ACCEPTED";

	tag2string[ POOL_RMSD_ADD_STRUCTURE ] = "POOL_RMSD_ADD_STRUCTURE";
	tag2string[ POOL_RMSD_EVALUATE ] = "POOL_RMSD_EVALUATE";
	tag2string[ POOL_RMSD_MASTER_EVALUATE ] = "POOL_RMSD_MASTER_EVALUATE";
	tag2string[ MPI_MASTER_BCAST_COORDS ] = "MPI_MASTER_BCAST_COORDS";
	tag2string[ MPI_MASTER_BCAST_WINNING_RANKS ] = "MPI_MASTER_BCAST_WINNING_RANKS";
	tag2string[ MPI_MASTER_BCAST_WINNING_STRUCTURES ] = "MPI_MASTER_BCAST_WINNING_STRUCTURES";
	tag2string[ MPI_MASTER_BCAST_NEW_COMM_SIZE ] = "MPI_MASTER_BCAST_NEW_COMM_SIZE";
	tag2string[ MPI_MASTER_BCAST_NEW_POOL_RANKS ] = "MPI_MASTER_BCAST_NEW_POOL_RANKS";
	tag2string[ MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD ] = "MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD";
	tag2string[ MPI_POOL_MASTER_THINKS ] = "MPI_POOL_MASTER_THINKS";
	tag2string[ MPI_POOL_SLAVE_THINKS ] = "MPI_POOL_SLAVE_THINKS";
	tag2string[ SIDECHAINMCMOVER ] = "SIDECHAINMCMOVER";
	tag2string[ SIMPLEINTGRAPH ] = "SIMPLEINTGRAPH";
	tag2string[ SIDECHAINMOVER ] = "SIDECHAINMOVER";
	tag2string[ INITIALIZE ] = "INITIALIZE";
	tag2string[ COMM_REDUCE_SIZE ] = "COMM_REDUCE_SIZE";
	tag2string[ CANONICALMOVER_WRITE_TO_FILE ] = "CANONICALMOVER_WRITE_TO_FILE";
	tag2string[ WRITE_TO_FILE ] = "WRITE_TO_FILE";
	tag2string[ CHECK_COMM_SIZE ] = "CHECK_COMM_SIZE";
	tag2string[ APPLY_MOVE ] = "APPLY_MOVE";
	tag2string[ DATA_STORAGE ] = "DATA_STORAGE";
	tag2string[ MASTER_PROCESS_NEW_STRUCTURES ] = "MASTER_PROCESS_NEW_STRUCTURES";
	tag2string[ COPY_COORDS ] = "COPY_COORDS";
	tag2string[ APPLY_SC_MOVE ] = "APPLY_SC_MOVE";
	tag2string[ APPLY_BB_MOVE ] = "APPLY_BB_MOVE";

	tag2string[ HIERARCHICAL_EVALUATE ] = "HIERARCHICAL_EVALUATE";
	tag2string[ HIERARCHICAL_ADD ] = "HIERARCHICAL_ADD";
	tag2string[ LOAD_HIERARCHY ] = "LOAD_HIERARCHY";
	tag2string[ HIERARCHY_SEND_COORDS ] = "HIERARCHY_SEND_COORDS";
	tag2string[ HIERARCHY_RECV_COORDS ] = "HIERARCHY_RECV_COORDS";
	tag2string[ INITIALIZE_HIERARCHY ] = "INITIALIZE_HIERARCHY";
	tag2string[ WRITE_DECOYS_TO_HIERARCHY ] = "WRITE_DECOYS_TO_HIERARCHY";
	tag2string[ HIERARCHY_GET_NEXT_CANDIDATE ] = "HIERARCHY_GET_NEXT_CANDIDATE";
	tag2string[ HIERARCHY_FIND_ADDRESS ] = "HIERARCHY_FIND_ADDRESS";
	tag2string[ MPIH_EVAL_CHECK_PROGRESS ] = "MPIH_EVAL_CHECK_PROGRESS";
	tag2string[ MPIH_EVAL_COMMUNICATE_NEW ] = "MPIH_EVAL_COMMUNICATE_NEW";
	tag2string[ MPIH_EVAL_AGAINST_NBR ] = "MPIH_EVAL_AGAINST_NBR";
	tag2string[ MPIH_PREPARE_WRITE_STRUCTURES ] = "MPIH_PREPARE_WRITE_STRUCTURES";
	tag2string[ MPIH_UPDATE_EVAL ] = "MPIH_UPDATE_EVAL";
	tag2string[ MPIH_ADD_FIRST_STRUCTURE ] = "MPIH_ADD_FIRST_STRUCTURE";
	tag2string[ MPIH_WRITE_STRUCT ] = "MPIH_WRITE_STRUCT";
	tag2string[ HIERARCHY_SETUP_TO_RECV ] = "HIERARCHY_SETUP_TO_RECV";
	tag2string[ FINALIZE ] = "FINALIZE";
	tag2string[ SORT_POOL ] = "SORT_POOL";
	tag2string[ HIERARCHICAL_FIND ] = "HIERARCHICAL_FIND";
	tag2string[ HIERARCHICAL_FIND_ADDRESS ] = "HIERARCHICAL_FIND_ADDRESS";
	tag2string[ HIERARCHICAL_SORT_ADDRESSES ] = "HIERARCHICAL_SORT_ADDRESSES";
	tag2string[ HIERARCHICAL_ROUND ] = "HIERARCHICAL_ROUND";
	tag2string[ HIERARCHICAL_POOL_SIZE ] = "HIERARCHICAL_POOL_SIZE";
	tag2string[ HIERARCHICAL_ADD_ELEM_TO_CACHE ] = "HIERARCHICAL_ADD_ELEM_TO_CACHE";
	tag2string[ LIB_FULL_PATH ] = "LIB_FULL_PATH";

	tag2string[ NOESY_ASSIGN_TOTAL ] = "NOESY_ASSIGN_TOTAL";
	tag2string[ NOESY_ASSIGN_INITIAL ] = "NOESY_ASSIGN_INITIAL";
	tag2string[ NOESY_ASSIGN_DIAGONAL ] = "NOESY_ASSIGN_DIAGONAL";
	tag2string[ NOESY_ASSIGN_CHEMSHIFT ] = "NOESY_ASSIGN_CHEMSHIFT";
	tag2string[ NOESY_ASSIGN_SYMMETRY ] = "NOESY_ASSIGN_SYMMETRY";
	tag2string[ NOESY_ASSIGN_DISTANCE ] = "NOESY_ASSIGN_DISTANCE";
	tag2string[ NOESY_ASSIGN_DECOY_COMP ] = "NOESY_ASSIGN_DECOY_COMP";
	tag2string[ NOESY_ASSIGN_NETWORK ] = "NOESY_ASSIGN_NETWORK";
	tag2string[ NOESY_ASSIGN_NETWORK_TOTAL ] = "NOESY_ASSIGN_NETWORK_TOTAL";
	tag2string[ NOESY_ASSIGN_NETWORK_FIND_RAW_ASSIGN ] = "NOESY_ASSIGN_NETWORK_FIND_RAW_ASSIGN";
	tag2string[ NOESY_ASSIGN_NETWORK_FILL_COV_GAMMA ] = "NOESY_ASSIGN_NETWORK_FILL_COV_GAMMA";
	tag2string[ NOESY_ASSIGN_NETWORK_RETRIEVE_ASSIGN ] = "NOESY_ASSIGN_NETWORK_RETRIEVE_ASSIGN";
	tag2string[ NOESY_ASSIGN_NETWORK_COMPUTE_NK ] = "NOESY_ASSIGN_NETWORK_COMPUTE_NK";
	tag2string[ NOESY_ASSIGN_NETWORK_PEAK_COUNT ] = "NOESY_ASSIGN_NETWORK_PEAK_COUNT";
	tag2string[ NOESY_ASSIGN_NETWORK_PEAK_COUNT_EVAL ] = "NOESY_ASSIGN_NETWORK_PEAK_COUNT_EVAL";
	tag2string[ NOESY_ASSIGN_NETWORK_INVALIDATE_SEQ_NOE ] = "NOESY_ASSIGN_NETWORK_INVALIDATE_SEQ_NOE";

	tag2string[ NOESY_ASSIGN_UPDATE_PEAK_VOL ] = "NOESY_ASSIGN_UPDATE_PEAK_VOL";
	tag2string[ NOESY_ASSIGN_CALIBRATE ] = "NOESY_ASSIGN_CALIBRATE";
	tag2string[ NOESY_ASSIGN_ELIMINATE ] = "NOESY_ASSIGN_ELIMINATE";
	tag2string[ NOESY_ASSIGN_GEN_CST ] = "NOESY_ASSIGN_GEN_CST";
	tag2string[ NOESY_ASSIGN_WRITE_CST ] = "NOESY_ASSIGN_WRITE_CST";
	tag2string[ NOESY_ASSIGN_MAP2CB ] = "NOESY_ASSIGN_MAP2CB";
	tag2string[ NOESY_ASSIGN_CP_GEN_CST ] = "NOESY_ASSIGN_CP_GEN_CST";
	tag2string[ NOESY_ASSIGN_PA_GEN_CST ] = "NOESY_ASSIGN_PA_GEN_CST";
	tag2string[ NOESY_ASSIGN_NMR_STRING ] = "NOESY_ASSIGN_NMR_STRING";
	tag2string[ NOESY_ASSIGN_REQUIRES_CB_MAPPING ] = "NOESY_ASSIGN_REQUIRES_CB_MAPPING";
	tag2string[ NOESY_ASSIGN_MAP2CB_NEW ] = "NOESY_ASSIGN_MAP2CB_NEW";
	tag2string[ NOESY_ASSIGN_WRITE_ASSIGNMENTS ] = "NOESY_ASSIGN_WRITE_ASSIGNMENTS";
	tag2string[ NOESY_ASSIGN_READ_INPUT ] = "NOESY_ASSIGN_READ_INPUT";
	tag2string[ NOESY_ASSIGN_DIST_INIT ] = "NOESY_ASSIGN_DIST_INIT";
	tag2string[ NOESY_ASSIGN_DIST_PREP_SCORE ] = "NOESY_ASSIGN_DIST_PREP_SCORE";
	tag2string[ NOESY_ASSIGN_DIST_SET_COMPABILITY_SCORE ] = "NOESY_ASSIGN_DIST_SET_COMPABILITY_SCORE";
	tag2string[ NOESY_ASSIGN_DIST_APPLY ] = "NOESY_ASSIGN_DIST_APPLY";
	tag2string[ NOESY_ASSIGN_DIST_MAKE_POSE ] = "NOESY_ASSIGN_DIST_MAKE_POSE";
	tag2string[ NOESY_ASSIGN_DIST_CST_EVAL ] = "NOESY_ASSIGN_DIST_CST_EVAL";
	tag2string[ NOESY_ASSIGN_DIST_CST_CAST ] = "NOESY_ASSIGN_DIST_CST_CAST";
	tag2string[ SILENT_FILL_POSE ] = "SILENT_FILL_POSE";
	tag2string[ SILENT_SET_POSE_COORDS ] = "SILENT_SET_POSE_COORDS";
	tag2string[ SILENT_FILL_STRUCT ] = "SILENT_FILL_STRUCT";
}


std::map<std::string, double> dynamic_prof_total;
std::map<std::string, int> dynamic_prof_calls;

DynamicProfileThis::DynamicProfileThis( std::string const& tag ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// don't profile unless instructed to via the option -run:profile
	if ( !option[basic::options::OptionKeys::run::profile] ) {
		return;
	}

	tag_ = tag;
	start_clock_ = clock() / SHRINK_FACTOR;
}

DynamicProfileThis::~DynamicProfileThis() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	// don't profile unless instructed to via the option -run:profile
	if ( !option[basic::options::OptionKeys::run::profile] ) {
		return;
	}

	clock_t const current( clock() / SHRINK_FACTOR );
	clock_t const start( start_clock_ );

	if ( current >= start ) {
		dynamic_prof_total[ tag_ ] += clock_factor * ( current - start );
		dynamic_prof_calls[ tag_ ] += 1;
	}
}

void prof_show() {
	using namespace ObjexxFCL;
	basic::Tracer tt( "core.util.prof", basic::t_info, true /*muted by default*/ );

	// clocks are shown in the unit of 0.01 sec
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// don't profile unless instructed to via the option -run:profile
	if ( !option[basic::options::OptionKeys::run::profile] ) {
		return;
	}

	static bool init( false );
	if ( !init ) {
		init = true;
		setup_tag2string();
	}

	PROF_STOP( basic::TOTAL );

	tt << "\n========================================\n";
	tt << "========================================\n";
	tt << "============ PROFILE INFO ==============\n";
	tt << "========================================\n";
	tt << "========================================\n";
	tt << A(12,"clock") << ' ' << A(9,"ncalls") << ' ' << A(9,"bad_calls") << ' ' << A(12,"c/call") << ' ' <<
		"tag" << '\n';
	for ( int i=1; i<= n_prof_tags; ++i ) {
		ProfTag const tag( static_cast< ProfTag >( i ) );
		double const t( total_clock[ tag ] );
		int const ncalls( calls[tag] );
		int const bcalls( bad_calls[tag] );
		double const clocks_per_call( ncalls != 0 ? t/ncalls : 0.0 );
		if ( ncalls ) {
			tt << F(12,2,t) << ' ' << I(9,ncalls) << ' ' << I(9,bcalls)
				<< ' ' << F(12,3, clocks_per_call ) << ' ' << tag2string[tag] << '\n';
		}
	}
	for ( std::map< std::string, double >::const_iterator it=dynamic_prof_total.begin(); it!=dynamic_prof_total.end(); ++it ) {
		std::string const& tag( it->first );
		double const t( it->second );
		int const ncalls( dynamic_prof_calls[tag] );
		double const clocks_per_call( ncalls != 0 ? t/ncalls : 0.0 );
		if ( ncalls ) {
			tt << F(12,2,t) << ' ' << I(9,ncalls) << ' ' << I(9,0)
				<< ' ' << F(12,3, clocks_per_call ) << ' ' << tag << '\n';
		}

	}
	tt << "========================================\n";
	tt << "========================================\n";
	tt << "========================================" << std::endl;
	PROF_START( TOTAL );
}

void prof_reset() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// don't profile unless instructed to via the option -run:profile
	if ( !option[basic::options::OptionKeys::run::profile] ) {
		return;
	}

	start_clock.clear();
	total_clock.clear();
	calls.clear();
	bad_calls.clear();

	start_clock.resize( n_prof_tags, 0 );
	total_clock.resize( n_prof_tags, 0 );
	calls.resize( n_prof_tags, 0 );
	bad_calls.resize( n_prof_tags, 0 );

	PROF_START( TOTAL );
}


} // basic
