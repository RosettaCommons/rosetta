// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/OptionKeys.cc
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Unit headers
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.cc.include.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <utility/down_cast.hh>
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
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringVectorOption.fwd.hh>
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
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/DenovoProteinDesign.OptionKeys.gen.hh>
#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/helixAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/LoopModel.OptionKeys.gen.hh>
#include <basic/options/keys/MM.OptionKeys.gen.hh>
#include <basic/options/keys/MonteCarlo.OptionKeys.gen.hh>
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
#include <basic/options/keys/in.OptionKeys.gen.hh>
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
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
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
namespace OptionKeys {


// Include options generated by Python script
#include <basic/options/keys/OptionKeys.cc.gen2.hh>
/// @brief Lookup functors
//#include <utility/keys/KeyLookup.functors.ii>


} // namespace OptionKeys
} // namespace options
} // namespace basic
