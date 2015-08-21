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

#ifndef INCLUDED_basic_prof_hh
#define INCLUDED_basic_prof_hh

#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <string>

#include <platform/types.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <time.h>
#include <utility>
#include <vector>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

namespace basic {
/**
not intended for profiling inside tight loops
the clock() routine currently being used has fairly crappy
resolution and it will introduce some small overhead with the
function calls and the if-check even if not using -profile on
the command line

you can wrap it around large-ish chunks of code, like fullatom_energy
or rotamer_trials...

A simple setup for timing code fragments. Probably not optimal timing
functions -- I'm open to suggestions.

looks like (see eg fullatom_energy or scorefxn)

PROF_START( prof::TAG );

<function call>

PROF_STOP( prof::TAG );

where TAG is in the enum "Prof_tag" below (feel free to add new ones)
also add to tag2string if you want friendly output.

PROF_STOP checks the time and increments the total time assigned to TAG


2. later on, in your simulation code you can do:

prof_reset();

-- miscellaneous simulation --

prof_show();

The final call to prof::show() will display the time usage measured
by all the PROF_* calls between reset() and show()
**/

#ifdef NO_PROF
#define PROF_START(expr)
#define PROF_STOP(expr)
#else
#define PROF_START(expr) ( prof_start_function_body( expr) )
#define PROF_STOP(expr) ( prof_stop_function_body( expr ) )
#endif

enum ProfTag {
	TEST1 = 1,
	TEST2,
	TEST3,
	TEST4,
	TEST5,
	ATOM_TREE_UPDATE_INTERNAL_COORDS,
	ATOM_TREE_UPDATE_XYZ_COORDS,
	UPDATE_RESIDUE_TORSIONS,
	UPDATE_RESIDUE_COORDINATES,
	ROTAMER_TRIALS,
	PACK_ROTAMERS,
	UPDATE_RESIDUE_NEIGHBORS,
	SETUP_NBLIST,
	CONSTRAINT_SCORE,
	CONSTRAINT_SET_COPY,

	SCORE,
	SCORE_SETUP,
	SCORE_FINALIZE,
	SCORE_ONEBODY_ENERGIES,
	SCORE_NEIGHBOR_ENERGIES,
	SCORE_LONG_RANGE_ENERGIES,
	SCORE_DOT,
	SCORE_END_NOTIFY,
	SCORE_BEGIN_NOTIFY,
	VDW_ENERGY,
	ENERGY_ENVPAIR_POTENTIAL,
	SECONDARY_STRUCTURE_ENERGY,
	SECONDARY_STRUCTURE_SSPAIR_ENERGY,
	SECONDARY_STRUCTURE_HSPAIR_ENERGY,
	SECONDARY_STRUCTURE_SHEETS_FROM_DIMERS_ENERGY,

	POSE_COPY,
	ENERGY_GRAPH_COPY,

	ENERGIES_COPY,

	CONFORMATION_DETECT_DISULF,
	CONFORMATION_FIX_DISULF,
	CONFORMATION_COPY,
	CHEMICAL_MAKE_POSE,

	CCD_CLOSE,
	FUNC,
	DFUNC,
	GET_ENERGIES,
	SIMANNEALING,
	INSERT_FRAGS,
	MC_ACCEPT,
	GB_SETUP_FOR_PACKING,
	GB_GET_ALL_BORN_RADII,
	GEN_BORN_ROTAMER_PAIR_ENERGIES,
	GEN_BORN_ROTAMER_BACKGROUND_ENERGIES,
	MULTIPOLE_SETUP,
	MULTIPOLE_ENERGIES,
	FACTS_SETUP_FOR_PACKING,
	FACTS_GET_ALL_BORN_RADII,
	FACTS_ROTAMER_PAIR_ENERGIES,
	FACTS_ROTAMER_BACKGROUND_ENERGIES,
	ABINITIO,
	STAGE1,
	STAGE2,
	STAGE3,
	STAGE4,
	STAGE5,
	FRAGMENT_MOVER,
	MINMOVER_APPLY,
	BACKRUB_MOVER,
	FIND_SUGAR_AND_SUITE_FRAGS_I,
	FIND_SUGAR_AND_SUITE_FRAGS_II,
	MAKE_BASE_PAIR_MOVE,
	MAKE_BASE_STEP_MOVE,
	RG,
	RG_LOCAL,
	SEQUENCE_COMPARISON,
	KDTREE_CONSTRUCT,
	KDTREE_SEARCH,
	CONSTRUCT_DISTANCE_MATRIX,
	DALPHABALL,
	DALPHABALL_DERIV,
	MPI_FILE_BUF,
	JD2,
	JD2_OUTPUT,
	JD2_SILENT_OUTPUTTER,
	JD2_INIT_MOVER,
	ARCHIVE_SYNC_BATCHES,
	ARCHIVE_JOBSCOMPLETE,
	ARCHIVE_BLOCK_FILE,
	SAVE_ARCHIVE,
	ARCHIVE_CRITICAL_JOBSCOMPLETE,
	ARCHIVE_GEN_BATCH,
	ARCHIVE_READ_DECOYS,
	ARCHIVE_EVAL_DECOYS,
	ARCHIVE_FILL_POSE,
	ARCHIVE_SCORE_POSE,
	ARCHIVE_EVALUATORS,
	MPI_JD2_WAITS_FOR_ARCHIVE,
	MPI_NOTIFY_ARCHIVE,
	SILENT_READ_TAG_TEST,
	CA_RMSD_EVALUATION,
	CA_DME_EVALUATION,
	TRUNCATED_SCORE_EVALUATOR,

	SAXS,
	FRAGMENTPICKING_CS_SCORE,
	FRAGMENTPICKING_PROFILE_SCORE,
	FRAGMENTPICKING_PROFILE_CAHING,
	FRAGMENTPICKING_SECONDARY_SCORE,
	FRAGMENTPICKING_READ_VALL,
	FRAGMENTPICKING,
	FRAGMENTPICKING_CANDIDATES_COLLECTING,
	FRAGMENTPICKING_ATOMPAIR_SCORE,
	FRAGMENTPICKING_PHIPSI_SCORE,
	FRAGMENTPICKING_DIHEDRALCONSTR_SCORE,

	MPICANONICALSAMPLING,
	MPIPOOLCOMMUNICATION,
	MPICOMMCREATION,
	MPIBARRIER,
	MPIBARRIER_BEGIN,
	MPIBARRIER_END,
	MPI_GATHER_BARRIER,
	FARRAY_MANIPULATION,
	MPI_SLAVE_REPORT_NEW_COORDS,
	MPI_SLAVE_REPORT_SIZES,

	MPI_SEND_UPDATE,
	MPI_SYNC_POOL_DIFF,
	MPI_SEND_ACCEPTED,

	POOL_RMSD_ADD_STRUCTURE,
	POOL_RMSD_EVALUATE,
	POOL_RMSD_MASTER_EVALUATE,
	MPI_MASTER_BCAST_COORDS,
	MPI_MASTER_BCAST_WINNING_RANKS,
	MPI_MASTER_BCAST_WINNING_STRUCTURES,
	MPI_MASTER_BCAST_NEW_COMM_SIZE,
	MPI_MASTER_BCAST_NEW_POOL_RANKS,
	MPI_MASTER_BCAST_NUM_STRUCTURES_TO_ADD,
	MPI_POOL_MASTER_THINKS,
	MPI_POOL_SLAVE_THINKS,
	SIDECHAINMCMOVER,
	SIMPLEINTGRAPH,
	SIDECHAINMOVER,
	INITIALIZE,
	COMM_REDUCE_SIZE,
	CANONICALMOVER_WRITE_TO_FILE,
	WRITE_TO_FILE,
	CHECK_COMM_SIZE,
	APPLY_MOVE,
	DATA_STORAGE,
	MASTER_PROCESS_NEW_STRUCTURES,
	COPY_COORDS,
	APPLY_SC_MOVE,
	APPLY_BB_MOVE,

	HIERARCHICAL_EVALUATE,
	HIERARCHICAL_ADD,
	LOAD_HIERARCHY,
	HIERARCHY_SEND_COORDS,
	HIERARCHY_RECV_COORDS,
	INITIALIZE_HIERARCHY,
	HIERARCHY_SETUP_TO_RECV,
	WRITE_DECOYS_TO_HIERARCHY,
	HIERARCHY_GET_NEXT_CANDIDATE,
	HIERARCHY_FIND_ADDRESS,
	MPIH_EVAL_CHECK_PROGRESS,
	MPIH_EVAL_COMMUNICATE_NEW,
	MPIH_EVAL_AGAINST_NBR,
	MPIH_PREPARE_WRITE_STRUCTURES,
	MPIH_UPDATE_EVAL,
	MPIH_ADD_FIRST_STRUCTURE,
	MPIH_WRITE_STRUCT,
	FINALIZE,
	SORT_POOL,
	HIERARCHICAL_FIND,
	HIERARCHICAL_FIND_ADDRESS,
	HIERARCHICAL_SORT_ADDRESSES,
	HIERARCHICAL_ROUND,
	HIERARCHICAL_POOL_SIZE,
	HIERARCHICAL_ADD_ELEM_TO_CACHE,
	LIB_FULL_PATH,

	NOESY_ASSIGN_TOTAL,
	NOESY_ASSIGN_INITIAL,
	NOESY_ASSIGN_DIAGONAL,
	NOESY_ASSIGN_CHEMSHIFT,
	NOESY_ASSIGN_SYMMETRY,
	NOESY_ASSIGN_DISTANCE,
	NOESY_ASSIGN_DECOY_COMP,
	NOESY_ASSIGN_NETWORK_TOTAL,
	NOESY_ASSIGN_NETWORK,
	NOESY_ASSIGN_NETWORK_FIND_RAW_ASSIGN,
	NOESY_ASSIGN_NETWORK_FILL_COV_GAMMA,
	NOESY_ASSIGN_NETWORK_RETRIEVE_ASSIGN,
	NOESY_ASSIGN_NETWORK_COMPUTE_NK,
	NOESY_ASSIGN_NETWORK_PEAK_COUNT,
	NOESY_ASSIGN_NETWORK_PEAK_COUNT_EVAL,
	NOESY_ASSIGN_NETWORK_INVALIDATE_SEQ_NOE,
	NOESY_ASSIGN_UPDATE_PEAK_VOL,
	NOESY_ASSIGN_CALIBRATE,
	NOESY_ASSIGN_ELIMINATE,
	NOESY_ASSIGN_GEN_CST,
	NOESY_ASSIGN_WRITE_CST,
	NOESY_ASSIGN_MAP2CB,
	NOESY_ASSIGN_CP_GEN_CST,
	NOESY_ASSIGN_PA_GEN_CST,
	NOESY_ASSIGN_NMR_STRING,
	NOESY_ASSIGN_REQUIRES_CB_MAPPING,
	NOESY_ASSIGN_MAP2CB_NEW,
	NOESY_ASSIGN_WRITE_ASSIGNMENTS,
	NOESY_ASSIGN_READ_INPUT,
	NOESY_ASSIGN_DIST_INIT,
	NOESY_ASSIGN_DIST_PREP_SCORE,
	NOESY_ASSIGN_DIST_SET_COMPABILITY_SCORE,
	NOESY_ASSIGN_DIST_APPLY,
	NOESY_ASSIGN_DIST_MAKE_POSE,
	NOESY_ASSIGN_DIST_CST_EVAL,
	NOESY_ASSIGN_DIST_CST_CAST,
	// NOESY_ASSIGN_DIST_,
	SILENT_FILL_POSE,
	SILENT_SET_POSE_COORDS,
	SILENT_FILL_STRUCT,

	TOTAL, // keep these two last
	n_prof_tags = TOTAL
};


extern utility::vector1< std::string > tag2string;
extern utility::vector1< clock_t > start_clock;
extern utility::vector1< double > total_clock;
extern utility::vector1< int > calls;
extern utility::vector1< int > bad_calls;

extern double const clock_factor;
extern clock_t const SHRINK_FACTOR; // prevent overflow

inline
void
prof_start_function_body( ProfTag const tag )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// don't profile unless instructed to via the option -run:profile
	if ( !option[basic::options::OptionKeys::run::profile] ) {
		return;
	}

	start_clock[ tag ] = clock() / SHRINK_FACTOR;
	++calls[ tag ];
}


inline
void
prof_stop_function_body( ProfTag const tag )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// don't profile unless instructed to via the option -run:profile
	if ( !option[basic::options::OptionKeys::run::profile] ) {
		return;
	}

	clock_t const current( clock() / SHRINK_FACTOR );
	clock_t const start( start_clock[tag] );

	if ( current >= start ) {
		total_clock[ tag ] += clock_factor * ( current - start );
	} else {
		--calls[ tag ];
		++bad_calls[ tag ];
	}
}

void prof_reset();

void prof_show();

class ProfileThis {
public:
	ProfileThis( ProfTag tag ) : tag_( tag ) {
		PROF_START( tag_ );
	}
	~ProfileThis() {
		PROF_STOP( tag_ );
	}
private:
	ProfTag tag_;
};

class DynamicProfileThis {
public:
	DynamicProfileThis( std::string const& prof_tag );
	~DynamicProfileThis();
private:
	clock_t start_clock_;
	std::string tag_;
};

/// @brief print "TIME_STAMP: Www Mmm dd hh:mm:ss yyyy msg" on tr.Error and on std::cerr (if boolean is true)
extern bool show_time_on_cerr;
void show_time( basic::Tracer& tr, std::string const& msg );

} // basic

#endif
