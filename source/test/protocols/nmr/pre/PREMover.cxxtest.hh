// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pre/PREMover.cc
/// @brief   Mover that converts NMR PRE rates into CB-CB atom pair distances and
///          assigns them as atom pair constraints to the pose. In case of degenerate spins
///          (e.g. degenerate protons or symmetric spins) an AmbiguousNMRDistanceConstraint
///          is created.
///          As constraints function a SplineFunc is fitted to a PRE distance histogram.
///          By default, the histogram file is read from the spinlabel database and the generated
///          spline potential converts the HN PRE distances into CB-CB distance constraints.
///          Alternatively, the user can provide a different histogram file to generate CB-CB distance
///          constraints from PRE rates collected for a different type of PRE nucleus e.g. 15N or 1HA.
/// @details last Modified: 12/16/16
///          Note that the calculation of PRE distances is simplified and takes into account only
///          the dipolar and in part Curie relaxation so far. Cross-correlated relaxation is neglected.
///          This approach is reasonable for radicals or paramagnetic ions with a nearly isotropic g-tensor
///          (e.g. nitroxide, Mn2+, Cu2+) which are commonly used in NMR PRE experiments. However, for ions
///          with an anisotropic g-tensor (e.g. lanthanides) Curie and cross-correlated relaxation become
///          more prominent and a direct conversion of PRE rates to distances is not easily possible any longer.
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/rosettascripts.hh>

// Unit headers
#include <protocols/nmr/pre/PREMover.hh>
#include <protocols/nmr/pre/PREMoverCreator.hh>

// Package headers
#include <core/scoring/nmr/pre/PREData.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/MoverFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/SplineFunc.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <sstream>
#include <utility>

static basic::Tracer TR("protocols.nmr.pre.PREMover.cxxtest");

class PREMoverTests : public CxxTest::TestSuite {
private:
	core::pose::Pose gb1_;

public:

	/// @brief Setup Test
	void setUp() {
		// Initialize core & options system
		core_init();
		// Load pose from pdb
		core::import_pose::pose_from_file(gb1_, "protocols/nmr/pre/gb1_.pdb", core::import_pose::PDB_file);
	}

	void tearDown() {
		gb1_.clear();
	}

	/// @brief Test PREMover instantiation through MoverFactory
	void test_PREMover_instantiation_through_factory() {
		using namespace protocols::nmr::pre;

		protocols::moves::MoverOP base_mover_op = protocols::moves::MoverFactory::get_instance()->newMover("PREMover");
		PREMoverOP pm_op = utility::pointer::dynamic_pointer_cast< PREMover >(base_mover_op);
		TS_ASSERT(pm_op);
		TS_ASSERT_EQUALS(pm_op->get_name(), "PREMover");
	}

	void test_PREMover_instantiation_through_ctor() {
		using namespace protocols::nmr::pre;

		PREMoverOP pm_op = PREMoverOP( new PREMover("protocols/nmr/pre/pre_data_input.txt", gb1_) );
		TS_ASSERT_EQUALS(pm_op->get_name(), "PREMover");
	}

	/// @brief Test PREMover and PREDatasets setters and getters
	void test_PREMover_setters_and_getters() {
		using namespace protocols::nmr::pre;
		using namespace core::scoring::nmr::pre;
		using namespace core::scoring;

		PREMoverOP pm_op = PREMoverOP( new PREMover );
		PREDataOP data = PREDataOP( new PREData("protocols/nmr/pre/pre_data_input.txt", gb1_) );
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function("ref2015");

		pm_op->set_pre_data(data);
		pm_op->set_scorefunction(sfxn);
		pm_op->set_weighted_average(true);

		TS_ASSERT_EQUALS(pm_op->get_pre_data()->get_number_spinlabel_sites(), 1);
		TS_ASSERT_EQUALS(pm_op->get_scorefunction()->get_name(), "ref2015");
		TS_ASSERT(pm_op->weighted_average());
	}

	/// @brief Test parsing XML options to the PREMover
	void test_PREMover_parse_my_tag_function() {
		using namespace protocols::nmr::pre;
		using namespace utility::tag;

		PREMoverOP pm_op(new PREMover);
		basic::datacache::DataMap data;
		protocols::moves::Movers_map m_map;
		protocols::filters::Filters_map f_map;
		core::scoring::ScoreFunctionOP sfxn(new core::scoring::ScoreFunction);
		sfxn->set_weight(core::scoring::atom_pair_constraint,1.0);
		data.add("scorefxns","stage1", sfxn);

		// Options for scoring function, datafile, constraints weight, deviation and function
		std::string xml1 = "<PREMover scorefxn=\"stage1\" />";
		std::string xml2 = "<PREMover pre_input_file=\"protocols/nmr/pre/pre_data_input.txt\" />";
		std::string xml3 = "<PREMover minimize_w_csts=\"true\" />";
		std::string xml4 = "<PREMover weighted_average=\"true\" />";
		std::string xml5 =
			"<PREMover>"
			"  <Histograms spinlabel=\"R1A\" histogram_file=\"protocols/nmr/pre/r1a_pre_n_distance_potential.histogram\" bin_size=\"0.5\"/>"
			"</PREMover>";
		TagOP tag1 = Tag::create(xml1);
		TagOP tag2 = Tag::create(xml2);
		TagOP tag3 = Tag::create(xml3);
		TagOP tag4 = Tag::create(xml4);
		TagOP tag5 = Tag::create(xml5);

		// Has scorefunction been parsed?
		pm_op->parse_my_tag(tag1, data, f_map, m_map, gb1_);
		TS_ASSERT_EQUALS(pm_op->get_scorefunction(), sfxn);

		// Have PRE data been parsed?
		pm_op->parse_my_tag(tag2, data, f_map, m_map, gb1_);
		TS_ASSERT_EQUALS(pm_op->get_pre_data()->get_number_spinlabel_sites(), 1);

		// Have boolean options been parsed?
		pm_op->parse_my_tag(tag3, data, f_map, m_map, gb1_);
		TS_ASSERT(pm_op->minimize_w_pre_csts());
		pm_op->parse_my_tag(tag4, data, f_map, m_map, gb1_);
		TS_ASSERT(pm_op->weighted_average());

		pm_op->parse_my_tag(tag5, data, f_map, m_map, gb1_);

		//pm_op->show(TR.Debug);
	}

	/// @brief Test that PREMover assigns constraints to the pose
	void test_PREMover_apply_function() {
		using namespace protocols::nmr::pre;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;

		PREMoverOP pm_op = PREMoverOP( new PREMover("protocols/nmr/pre/pre_data_input.txt", gb1_) );
		core::scoring::ScoreFunctionOP sfxn = get_score_function();
		sfxn->set_weight(atom_pair_constraint,10.0);
		pm_op->set_scorefunction(sfxn);
		pm_op->set_weighted_average(true);

		ConstraintSetCOP csts = gb1_.constraint_set();
		TS_ASSERT(csts->is_empty());

		// Has a constraint set been created?
		pm_op->apply(gb1_);
		csts = gb1_.constraint_set();
		TS_ASSERT(!csts->is_empty());

		// Pull out the first atom pair constraint from second pre dataset which is constraint number 138
		utility::vector1< ConstraintCOP >csts_vec = csts->get_all_constraints();
		AtomPairConstraintOP atompair_csts_op = utility::pointer::dynamic_pointer_cast< AtomPairConstraint >(utility::pointer::const_pointer_cast< Constraint > (csts_vec[1]));
		SplineFunc & spline = dynamic_cast< SplineFunc& >(const_cast< Func& >(atompair_csts_op->get_func()));
		TS_ASSERT(atompair_csts_op);
		TS_ASSERT_EQUALS(atompair_csts_op->atom1().rsd(), 13);
		TS_ASSERT_EQUALS(atompair_csts_op->atom2().rsd(),  2);
		core::Real x(atompair_csts_op->dist(gb1_));
		TS_ASSERT_DELTA(x, 15.9772, 1.0e-3);
		TS_ASSERT_DELTA(spline.get_exp_val(), 18.0251, 1.0e-3);
		TS_ASSERT_DELTA(spline.get_weight(), 3.7109, 1.0e-3);
		TS_ASSERT_DELTA(spline.func(x), -3.49862, 1.0e-3);
	}

};
