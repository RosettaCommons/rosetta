// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/PoseMetricCalculators/RotamerBoltzCalculator.hh
/// @brief Calculates Rotamer occupancy of each rotameric state in a given set of residues.
/// @author Hetu Kamisetty
/// @author Tom Linsky (tlinsky at uw dot edu)


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_RotamerBoltzCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_RotamerBoltzCalculator_hh
#include <protocols/toolbox/pose_metric_calculators/RotamerBoltzCalculator.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/toolbox/EnergyLandscapeEvaluator.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>
//#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <utility/vector0.fwd.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

typedef core::Real RotamerProbability;
typedef std::map< core::Size, RotamerProbability > RotamerProbabilities;

class RotamerBoltzCalculator: public core::pose::metrics::StructureDependentCalculator {
public:

	RotamerBoltzCalculator(core::scoring::ScoreFunctionOP scorefxn, core::Real temp, core::Real repacking_radius=6.0);
	//computeBoltzWeight(utility::vector0<int> moltres_to_pack, utility::vector0<int> rotid_in_moltres, core::Real temp);//gives rotamer id into rotamerset for each position of interest.
	//computeBoltzWeight(utility::vector0<int> rot_to_pack);//gives rot_to_pack where entry is rotameric id into rotamersets
	void computeAllBoltz( core::pose::Pose const & pose );
	RotamerProbability computeBoltzWeight(core::pose::Pose& pose, core::Size resi);
	virtual core::pose::metrics::PoseMetricCalculatorOP clone() const;

	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	void set_energy_landscape_evaluator( protocols::toolbox::EnergyLandscapeEvaluatorCOP evaluator );

	void set_lazy( bool const lazy );

	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );
	core::Real computeBoltzWeight(core::pose::Pose& pose, core::Size resi, protocols::simple_moves::MinMoverOP min_mover, core::pack::task::PackerTaskOP task);

	core::scoring::ScoreFunctionOP scorefxn() const {return scorefxn_;};
	core::Real computeBoltzSum(core::Real init_score, utility::vector1<core::Real> scores);

	protocols::simple_moves::MinMoverOP init_minmover(core::pose::Pose& pose, core::Size resi, bool unbound, core::pack::task::PackerTaskOP  task);

	core::pack::task::PackerTaskOP
	init_task( core::pose::Pose const & pose, core::Size const resi );

	/// @note This function returns an uninitialized filter; this should probably not be used.
	//protocols::simple_filters::ScoreTypeFilter stf(){
	//	return stf_;
	//}
	void temperature(core::Real temp){
		temperature_= temp;
	}
	void repacking_radius(core::Real rad) {
		repacking_radius_ = rad;
	}
	core::Real repacking_radius() const{
		return repacking_radius_;
	}
	core::Real temperature() const{
		return temperature_;
	}

private:
	//core::pose:: Pose & pose_;
	core::Real rb_jump(){
		return rb_jump_;
	}

	core::Real
	compute_boltz_weight_packrotamers(
		core::pose::Pose & pose,
		core::Size const resi,
		protocols::simple_moves::MinMoverOP min_mover,
		core::pack::task::PackerTaskOP task ) const;

	core::Real
	computeBoltzWeight_lazy(core::pose::Pose& pose, Size resi,  protocols::simple_moves::MinMoverOP min_mover, core::pack::task::PackerTaskOP task);

	core::pack::rotamer_set::RotamerSetOP
	create_rotamer_set( core::pose::Pose const & pose, core::Size const resi ) const;

private:
	/// @brief selects residues to compute rotamer boltz on
	core::select::residue_selector::ResidueSelectorCOP selector_;
	/// @brief computes probability from scores of rotamers
	protocols::toolbox::EnergyLandscapeEvaluatorCOP evaluator_;
	bool lazy_;
	core::Real rb_jump_;
	core::Real repacking_radius_;
	//utility::vector0<int> init_rot_to_pack(core::pack::rotamer_set::RotamerSetsCOP rotamer_sets, core::Size rot_to_miss);
	utility::vector0<int> init_rot_to_pack(core::pack::rotamer_set::RotamerSetsCOP rotamer_sets, core::Size moltenres, core::Size rot_to_fix);
	core::scoring::ScoreFunctionOP scorefxn_;
	//core::kinematics::MoveMapOP mm_;
	//protocols::simple_moves::MinMover min_mover_;
	core::Real temperature_;
	//protocols::simple_filters::ScoreTypeFilter const stf_;
	utility::vector1< core::Real > all_boltz_;
	RotamerProbabilities probabilities_;
	core::pack::rotamer_set::RotamerSetOP rotset_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	RotamerBoltzCalculator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_RotamerBoltzCalculator )
#endif // SERIALIZATION


#endif
