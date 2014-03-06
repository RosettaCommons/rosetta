// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/inout.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_inout_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_inout_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace inout { extern BooleanOptionKey const inout; }
namespace inout { extern BooleanOptionKey const fold_tree_io; }
namespace inout { extern BooleanOptionKey const dump_connect_info; }
namespace inout { namespace dbms { extern BooleanOptionKey const dbms; } }
namespace inout { namespace dbms { extern StringOptionKey const mode; } }
namespace inout { namespace dbms { extern StringOptionKey const database_name; } }
namespace inout { namespace dbms { extern StringOptionKey const pq_schema; } }
namespace inout { namespace dbms { extern StringOptionKey const host; } }
namespace inout { namespace dbms { extern StringOptionKey const user; } }
namespace inout { namespace dbms { extern StringOptionKey const password; } }
namespace inout { namespace dbms { extern IntegerOptionKey const port; } }
namespace inout { namespace dbms { extern BooleanOptionKey const readonly; } }
namespace inout { namespace dbms { extern BooleanOptionKey const separate_db_per_mpi_process; } }
namespace inout { namespace dbms { extern IntegerOptionKey const database_partition; } }
namespace inout { namespace dbms { extern BooleanOptionKey const use_compact_residue_schema; } }
namespace inout { namespace dbms { extern BooleanOptionKey const retry_failed_reads; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
