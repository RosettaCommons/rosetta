// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/option.cc
/// @brief  Program options global and initialization function
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Unit headers
#include <basic/options/option.hh>
#include <basic/options/option.cc.include.gen.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>

// option key includes

#include <basic/Tracer.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/options/option.cc.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
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
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
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
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/Tracer.fwd.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/DenovoProteinDesign.OptionKeys.gen.hh>
#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/helixAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/LoopModel.OptionKeys.gen.hh>
#include <basic/options/keys/MM.OptionKeys.gen.hh>
#include <basic/options/keys/MonteCarlo.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/PCS.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>
#include <basic/options/keys/SSrbrelax.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/abrelax.OptionKeys.gen.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/archive.OptionKeys.gen.hh>
#include <basic/options/keys/assembly.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/batch_relax.OptionKeys.gen.hh>
#include <basic/options/keys/bbg.OptionKeys.gen.hh>
#include <basic/options/keys/blivens.OptionKeys.gen.hh>
#include <basic/options/keys/boinc.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/options/keys/canonical_sampling.OptionKeys.gen.hh>
#include <basic/options/keys/casp.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/chrisk.OptionKeys.gen.hh>
#include <basic/options/keys/chunk.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/crossmatch.OptionKeys.gen.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/dwkulp.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/evolution.OptionKeys.gen.hh>
#include <basic/options/keys/fast_loops.OptionKeys.gen.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/keys/fldsgn.OptionKeys.gen.hh>
#include <basic/options/keys/flexPepDocking.OptionKeys.gen.hh>
#include <basic/options/keys/flexpack.OptionKeys.gen.hh>
#include <basic/options/keys/flxbb.OptionKeys.gen.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>
#include <basic/options/keys/fold_cst.OptionKeys.gen.hh>
#include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/gpu.OptionKeys.gen.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/krassk.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/loopfcst.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>
#include <basic/options/keys/ms.OptionKeys.gen.hh>
#include <basic/options/keys/murphp.OptionKeys.gen.hh>
#include <basic/options/keys/mysql.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <basic/options/keys/optE.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <basic/options/keys/orbitals.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/patterson.OptionKeys.gen.hh>
#include <basic/options/keys/pepspec.OptionKeys.gen.hh>
#include <basic/options/keys/phil.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/qsar.OptionKeys.gen.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/resample.OptionKeys.gen.hh>
#include <basic/options/keys/rescore.OptionKeys.gen.hh>
#include <basic/options/keys/residues.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/robert.OptionKeys.gen.hh>
#include <basic/options/keys/rot_anl.OptionKeys.gen.hh>
#include <basic/options/keys/rotamerdump.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/swa.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <basic/options/keys/threadsc.OptionKeys.gen.hh>
#include <basic/options/keys/ufv.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/cutoutdomain.OptionKeys.gen.hh>

namespace basic {
namespace options {

basic::Tracer TR("basic.options");

/// @details OptionCollection global
utility::options::OptionCollection option;


/// @brief Named verbosity levels
int const silent  ( 0 ); // No messages output
int const quiet   ( 1 );
int const standard( 2 );
int const inform  ( 4 );
int const chat    ( 6 );
int const yap     ( 7 );
int const gush    ( 8 );
int const verbose ( 9 ); // All messages output


/// @brief Initialize the options
OptionCollection &
initialize()
{

	using namespace utility::options;
	using namespace basic::options::OptionKeys;
	using utility::file::PathName;
#ifdef BOINC
	std::cerr << "Options::initialize()" << std::endl; std::cerr.flush();
#endif
	// Load the built-in options
	option.add_built_in_options();
#ifdef BOINC
	std::cerr << "Options::adding_options()" << std::endl; std::cerr.flush();
#endif
	// Include options generated by Python script

	add_all_rosetta_options( option );
#ifdef BOINC
	std::cerr << "Options::initialize() Check specs." << std::endl; std::cerr.flush();
#endif
	// Check for problems in the option specifications
	option.check_specs();
#ifdef BOINC
	std::cerr << "Options::initialize()  End reached" << std::endl; std::cerr.flush();
#endif

	return option;
}


/// @brief Process the specified options
/// note Do more complex value setup and checks here than the option system provides
OptionCollection &
process()
{
	using namespace utility::options;

	{ // Stupid binary fullatom options that really suck.
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if ( option[ in::file::fullatom ].user() ) {
			TR.Warning << "option[ in::file::fullatom ]() re-interpreted as setting "
				<< "option[ in::file::residue_type_set ]() to ";
			std::string type_set_name("fa_standard");
			if ( !option[ in::file::fullatom ]() ) {
				type_set_name = "centroid";
			}

			TR.Warning << type_set_name << std::endl;
			option[ in::file::residue_type_set ].value( type_set_name );
		}

		if ( option[ out::file::fullatom ].user() ) {
			TR.Warning << "option[ out::file::fullatom ]() re-interpreted as setting "
				<< "option[ out::file::residue_type_set ]() to ";
			std::string type_set_name("fa_standard");
			if ( !option[ out::file::fullatom ]() ) {
				type_set_name = "centroid";
			}

			TR.Warning << type_set_name << std::endl;
			option[ out::file::residue_type_set ].value( type_set_name );
		}
	}

	{ // Input paths
		using namespace basic::options::OptionKeys::in::path;

		option[ fragments ].default_to( option[ path ] );
		option[ pdb ].default_to( option[ path ] );
		option[ database ].default_to( option[ path ] );
	}

	{ // Output paths
		using namespace basic::options::OptionKeys::out::path;

		option[ pdb ].default_to( option[ path ] );
		option[ score ].default_to( option[ path ] );
		option[ movie ].default_to( option[ path ] );
	}

	{ // Packing options
		using namespace basic::options::OptionKeys::packing;

		if ( option[ solvate ] ) option[ explicit_h2o ].value( true ); // -solvate => -explicit_h2o
	}

	TR.flush();
	return option;
}


} // namespace options
} // namespace basic
