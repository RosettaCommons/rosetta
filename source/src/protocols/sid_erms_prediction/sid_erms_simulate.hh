// -*- mode:c++;tab-width:2;incdent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/protocols/sid_erms_prediction/sid_erms_simulate
/// @brief Protocol for analyzing protein interfaces and simulating ERMS data
/// @author Robert Bolz (bolz.13@osu.edu)

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include <tuple>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/chains_util.hh>

#include <core/scoring/Energies.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <cmath>

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/pose_metric_calculators/SaltBridgeCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricContainer.hh>

#include <boost/config.hpp>
#include <vector>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <utility/vectorL.hh>
#include <utility/options/OptionCollection.hh>

#include <basic/options/keys/sid_erms_simulate.OptionKeys.gen.hh>

namespace protocols {
namespace sid_erms_prediction {

using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace basic::options;
using namespace basic::options::OptionKeys;

//calculate probability
core::Real calc_prob(core::Real x,core::Real a, core::Real b);

//read in the complex type: subunits and connectivities (nodes and edges)
std::tuple<core::Size, utility::vector1<char>, utility::vector1<utility::vector1<std::string>>> read_complex_type(std::string &complex_type_filename, core::Size &n_chains, utility::vector1<char> &nodes, utility::vector1<utility::vector1<std::string>> &edges);

//make sure intensities sum to 1
void check_intensities( const utility::vector1<utility::vector1<core::Real>> &ERMS_data );

//read in acceleration energies or full ERMS
std::tuple<utility::vector1<core::Real>, utility::vector1<utility::vector1<core::Real>>> read_ERMS(utility::vector1<core::Real> &ACE, utility::vector1<utility::vector1<core::Real>> &ERMS_read, const core::Size n_chains);

//read in B values if given via file input
utility::vector1<core::Real> read_B_vals(utility::vector1<core::Real> &B);

//check to make sure interfaces are symmetric
void check_interface_symmetry(const core::pose::Pose pose_check, const utility::vector1<std::string> &edges_check, core::Real &dSASA, core::Real &PRE);

//calculate B values from PDB (and check for disulfide bond if dimer to shift steepness)
utility::vector1<core::Real> calc_B_values(utility::vector1<core::Real>  &B, const utility::vector1<utility::vector1<std::string>> &edges, core::Real &A, const core::Size n_chains, const core::pose::Pose pose_ );

//function to simulate ERMS
utility::vector1<utility::vector1<core::Real>> simulate_ERMS( utility::vector1<utility::vector1<core::Real>> &ERMS_prediction, const utility::vector1<core::Real> &ACE, const utility::vector1<core::Real> &B, const core::Size n_chains, const utility::vector1<utility::vector1<std::string>> &edges, const core::Real A, const core::Real breakage_cut, const utility::vector1<char> &nodes);

//calculate RMSE
core::Real calc_RMSE(const utility::vector1<utility::vector1<core::Real>> & ERMS_prediction, const utility::vector1<utility::vector1<core::Real>> &ERMS, const core::Size n_chains, const core::Size n_ACE);

}
}
