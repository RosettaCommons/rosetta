// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/mike.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_mike_OptionKeys_gen_hh
#define INCLUDED_basic_options_keys_mike_OptionKeys_gen_hh

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace mike { extern BooleanOptionKey const mike; }
namespace mike { extern BooleanOptionKey const lh_sandbox; }
namespace mike { extern BooleanOptionKey const lh_create_db; }
namespace mike { extern StringOptionKey const lh_db_prefix; }
namespace mike { extern IntegerVectorOptionKey const lh_loopsizes; }
namespace mike { extern IntegerOptionKey const lh_skim_size; }
namespace mike { extern StringOptionKey const lh_mpi_feedback; }
namespace mike { extern IntegerOptionKey const lh_mpi_batch_relax_chunks; }
namespace mike { extern IntegerOptionKey const lh_mpi_batch_relax_absolute_max; }
namespace mike { extern IntegerOptionKey const lh_mpi_outbound_wu_buffer_size; }
namespace mike { extern IntegerOptionKey const lh_mpi_loophash_split_size    ; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
