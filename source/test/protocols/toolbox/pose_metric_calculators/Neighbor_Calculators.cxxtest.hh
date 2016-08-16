// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/pose_metric_calculators/Neighbor_Calculators.cxxtest.hh
/// @brief  test for NeighborsByDistanceCalculator and NeighborhoodByDistanceCalculator
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
//#include <test/UTracer.hh>

// Unit header
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/InterGroupNeighborsCalculator.hh>

// project headers
#include <core/types.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/Conformation.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
#include <basic/MetricValue.hh>

// option key includes
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/toolbox/pose_metric_calculators/InterGroupNeighborsCalculator.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.fwd.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

class Neighbor_CalculatorsTests : public CxxTest::TestSuite {

	core::pose::Pose pose; //core/conformation/dock_in.pdb

public:

	// --------------- Fixtures --------------- //
	Neighbor_CalculatorsTests() {
		core_init();

		core::import_pose::pose_from_file( pose, "core/conformation/dock_in.pdb" , core::import_pose::PDB_file);
		//213, 204, and 265 are appropriate residues from the dock_in pdb

		//set up the Neighbors calculators
		using core::pose::metrics::PoseMetricCalculatorOP;
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "nbdc_buried", PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator(213) ) );

		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "nbdc_surface", PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator(204) ) );

		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "nbdc_interface", PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator(265) ) );

		//set up the Neighborhood calculators (testing each constructor type)
		std::set< core::Size > crset;
		crset.insert(213); crset.insert(204), crset.insert(265);
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "nh_crset_calc", PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator(crset) ) );
	}

	virtual ~Neighbor_CalculatorsTests() {}

	static Neighbor_CalculatorsTests* createSuite() {
		return new Neighbor_CalculatorsTests();
	}

	static void destroySuite( Neighbor_CalculatorsTests *suite ) {
		delete suite;
	}

	void setUp() {}

	void tearDown() {}

	// ------------- Helper Functions ------------- //
	//runs through the list of legal metrics from a Neighbor calculator
	void check_nbrs(
		std::string const & calc,
		core::Size resid,
		std::string const & neighbors,
		core::Size num_neighbors,
		core::Real cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::neighbor_by_distance_cutoff]){

		basic::MetricValue< core::Size > central_residue_, num_neighbors_;
		basic::MetricValue< core::Real > cutoff_;

		pose.metric( calc, "central_residue", central_residue_);
		pose.metric( calc, "dist_cutoff",     cutoff_);
		pose.metric( calc, "num_neighbors",   num_neighbors_);

		TS_ASSERT( resid == central_residue_.value() );
		TS_ASSERT( cutoff == cutoff_.value() );
		TS_ASSERT( num_neighbors == num_neighbors_.value() );
		TS_ASSERT( neighbors == pose.print_metric(calc, "neighbors") );
		//std::cout << "X" << pose.print_metric(calc, "neighbors") << "X" << std::endl;
		return;
	}

	//runs through the list of legal metrics from a Neighborhood calculator
	void check_hood(
		std::string const & calc,
		std::string const & resids,
		std::string const & neighbors,
		core::Size num_neighbors,
		core::Real cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::neighbor_by_distance_cutoff]){

		basic::MetricValue< core::Size > num_neighbors_;
		basic::MetricValue< core::Real > cutoff_;

		pose.metric( calc, "num_neighbors", num_neighbors_);
		pose.metric( calc, "dist_cutoff",     cutoff_);

		TS_ASSERT( cutoff == cutoff_.value() );
		TS_ASSERT( resids == pose.print_metric(calc, "central_residues") );
		TS_ASSERT( num_neighbors == num_neighbors_.value() );
		TS_ASSERT( neighbors == pose.print_metric(calc, "neighbors") );
		return;
	}

	//runs through the list of legal metrics from an InterGroupNeighbors calculator
	void check_ignc(
		std::string const & calc,
		std::string const & groups,
		core::Size const num_neighbors,
		std::string const & neighbors,
		core::Real cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::neighbor_by_distance_cutoff]){

		basic::MetricValue< core::Size > num_neighbors_;
		basic::MetricValue< core::Real > cutoff_;


		pose.metric( calc, "dist_cutoff",     cutoff_);
		pose.metric( calc, "num_neighbors",   num_neighbors_);

		TS_ASSERT( groups == pose.print_metric(calc, "groups"));
		TS_ASSERT( cutoff == cutoff_.value() );
		TS_ASSERT( num_neighbors == num_neighbors_.value() );
		TS_ASSERT( neighbors == pose.print_metric(calc, "neighbors") );
		//std::cout << "X" << pose.print_metric(calc, "neighbors") << "X" << std::endl;
		return;
	}

	// --------------- Test Cases --------------- //

	/// @brief
	void test_NeighborsByDistanceCalculators() {
		//the strings and magic numbers are the return values from calling the calculator
		std::string const bur_nbr("42 43 44 54 55 58 102 138 139 140 142 183 190 191 193 194 195 196 197 198 199 211 212 213 214 215 216 226 227 228 229 262 263 264 ");
		std::string const sur_nbr("122 202 203 204 205 206 208 ");
		std::string const int_nbr("33 39 40 41 42 43 58 141 142 143 149 150 151 191 192 193 194 195 263 264 265 266 267 277 278 299 ");

		check_nbrs("nbdc_buried", 213, bur_nbr, 34);
		check_nbrs("nbdc_surface", 204, sur_nbr, 7);
		check_nbrs("nbdc_interface", 265, int_nbr, 26);
		return;
	}

	void test_NeighborHoodDistanceCalculators() {
		//the strings and magic numbers are the return values from calling the calculator
		std::string const resids("204 213 265 ");
		core::Size const n_nbrs(57);
		std::string const neighborhood("33 39 40 41 42 43 44 54 55 58 102 122 138 139 140 141 142 143 149 150 151 183 190 191 192 193 194 195 196 197 198 199 202 203 204 205 206 208 211 212 213 214 215 216 226 227 228 229 262 263 264 265 266 267 277 278 299 ");

		check_hood("nh_crset_calc", resids, neighborhood, n_nbrs);

	}//end test_NeighborsByDistanceCalculator

	void test_InterGroupNeighborsCalculator() {

		//OK, so dock_in.pdb is not a great example protein for this calculator - it's intended for doing domain-domain interfaces or interfaces between chemically linked proteins (like ubiquitin).  We're going to pretend the two chains in this protein are domains and set it up that way.
		std::set< core::Size > domain1, domain2;
		core::Size const ch1end(pose.conformation().chain_end(1)), nres(pose.total_residue());
		for ( core::Size i(1); i<=ch1end; ++i ) domain1.insert(i);
		for ( core::Size i(ch1end+1); i<=nres; ++i ) domain2.insert(i);

		//set up the pair object
		std::pair< std::set<core::Size> , std::set<core::Size> > pair_of_domains;
		pair_of_domains = std::make_pair(domain1, domain2);

		//set up the complex container
		utility::vector1< std::pair< std::set<core::Size>, std::set<core::Size> > > domains;
		domains.push_back(pair_of_domains);

		//set up the calculator
		using core::pose::metrics::PoseMetricCalculatorOP;
		std::string const calc_d("IGNC_d");
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc_d, PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::InterGroupNeighborsCalculator(domains) ) );

		//   std::cout << pose.print_metric("IGNC_d", "groups") << std::endl;
		//   std::cout << pose.print_metric("IGNC_d", "dist_cutoff") << std::endl;
		//   std::cout << pose.print_metric("IGNC_d", "num_neighbors") << std::endl;
		//   std::cout << pose.print_metric("IGNC_d", "neighbors") << std::endl;
		check_ignc(
			calc_d,
			"{ ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245) ; ( 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301) }",
			82,
			"33 39 40 41 42 43 55 56 57 58 59 95 96 97 98 99 100 102 141 142 143 146 147 148 149 150 151 172 174 175 189 190 191 192 193 194 195 196 197 212 213 214 215 216 217 218 219 220 221 224 226 227 228 254 255 256 257 258 259 260 261 262 263 264 265 266 267 274 275 276 277 278 280 281 283 284 285 287 288 299 300 301 ");


		//as a further example with more domains:
		//set up an appropriate group
		std::set< core::Size > group1, group2, group3;
		//we are arbitrarily defining three contiguous groups of residues as three domains - this is not realistic for the structure but will serve for our purposes here.  Normally you'd have actual protein domains.

		for ( core::Size i(01); i<=10; ++i ) group1.insert(i);
		for ( core::Size i(51); i<=60; ++i ) group2.insert(i);
		for ( core::Size i(81); i<=90; ++i ) group3.insert(i);

		//make pair objects - here we try 1/2 and 2/3 but do not bother with 1/3
		std::pair< std::set<core::Size> , std::set<core::Size> > pair12, pair23;
		pair12 = std::make_pair(group1, group2);
		pair23 = std::make_pair(group2, group3);

		//push them into the vector
		utility::vector1< std::pair< std::set<core::Size>, std::set<core::Size> > > groups;
		groups.push_back(pair12);
		groups.push_back(pair23);

		//make the calculator
		std::string const calc_g("IGNC_g");
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc_g, PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::InterGroupNeighborsCalculator(groups) ) );

		//   std::cout << pose.print_metric("IGNC_g", "groups") << std::endl;
		//   std::cout << pose.print_metric("IGNC_g", "dist_cutoff") << std::endl;
		//   std::cout << pose.print_metric("IGNC_g", "num_neighbors") << std::endl;
		//   std::cout << pose.print_metric("IGNC_g", "neighbors") << std::endl;
		check_ignc(
			calc_g,
			"{ ( 1 2 3 4 5 6 7 8 9 10) ; ( 51 52 53 54 55 56 57 58 59 60) }{ ( 51 52 53 54 55 56 57 58 59 60) ; ( 81 82 83 84 85 86 87 88 89 90) }",
			16,
			"51 52 53 54 55 56 57 58 59 60 83 85 87 88 89 90 ");

	}


};//end class
