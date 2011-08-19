// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/lh.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_lh_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_lh_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace lh { extern BooleanOptionKey const lh; }
namespace lh { extern StringOptionKey const db_prefix; }
namespace lh { extern IntegerVectorOptionKey const loopsizes; }
namespace lh { extern IntegerOptionKey const num_partitions; }
namespace lh { extern PathOptionKey const db_path; }
namespace lh { extern BooleanOptionKey const exclude_homo; }
namespace lh { extern StringOptionKey const refstruct; }
namespace lh { extern StringOptionKey const homo_file; }
namespace lh { extern RealVectorOptionKey const createdb_rms_cutoff; }
namespace lh { extern RealOptionKey const min_bbrms; }
namespace lh { extern RealOptionKey const max_bbrms; }
namespace lh { extern RealOptionKey const min_rms; }
namespace lh { extern RealOptionKey const max_rms; }
namespace lh { extern RealOptionKey const grid_space_multiplier; }
namespace lh { extern RealOptionKey const grid_angle_multiplier; }
namespace lh { extern IntegerOptionKey const skim_size; }
namespace lh { extern IntegerOptionKey const rounds; }
namespace lh { extern StringOptionKey const jobname; }
namespace lh { extern IntegerOptionKey const max_lib_size; }
namespace lh { extern IntegerOptionKey const max_emperor_lib_size; }
namespace lh { extern IntegerOptionKey const max_emperor_lib_round; }
namespace lh { extern IntegerOptionKey const library_expiry_time; }
namespace lh { extern StringOptionKey const objective_function; }
namespace lh { extern IntegerOptionKey const expire_after_rounds; }
namespace lh { extern StringOptionKey const mpi_resume; }
namespace lh { extern StringOptionKey const mpi_feedback; }
namespace lh { extern IntegerOptionKey const mpi_batch_relax_chunks; }
namespace lh { extern IntegerOptionKey const mpi_batch_relax_absolute_max; }
namespace lh { extern IntegerOptionKey const mpi_outbound_wu_buffer_size; }
namespace lh { extern IntegerOptionKey const mpi_loophash_split_size    ; }
namespace lh { extern RealOptionKey const mpi_metropolis_temp; }
namespace lh { extern IntegerOptionKey const mpi_save_state_interval; }
namespace lh { extern BooleanOptionKey const mpi_master_save_score_only; }
namespace lh { extern IntegerOptionKey const max_loophash_per_structure; }
namespace lh { extern RealOptionKey const rms_limit; }
namespace lh { extern BooleanOptionKey const centroid_only; }
namespace lh { extern BooleanOptionKey const write_centroid_structs; }
namespace lh { extern BooleanOptionKey const sandbox; }
namespace lh { extern BooleanOptionKey const create_db; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
