// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/ddG.hh
/// @brief calculate interface ddg
/// @author Sarel Fleishman


#ifndef INCLUDED_protocols_simple_moves_ddG_hh
#define INCLUDED_protocols_simple_moves_ddG_hh

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/energy_methods/PoissonBoltzmannEnergy.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/types.hh>
#include <protocols/simple_ddg/ddG.fwd.hh>
#include <string>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <set>

#include <core/pack/task/PackerTask.fwd.hh> // AUTO IWYU For PackerTaskOP
#include <core/scoring/ScoreFunction.fwd.hh> // AUTO IWYU For ScoreFunctionCOP, ScoreFunctionOP
#include <map> // AUTO IWYU For map

namespace protocols {
namespace simple_ddg {

// This used to be a subclass of calc_taskop_movers::DesignRepackMover, but drifted over time so that it's now completely independant.
class ddG : public protocols::moves::Mover
{
public:

	typedef core::Real Real;
	typedef core::scoring::ScoreType ScoreType;
	typedef core::pose::Pose Pose;

public :
	ddG();
	ddG( core::scoring::ScoreFunctionCOP scorefxn_in, core::Size const jump=1 );
	ddG( core::scoring::ScoreFunctionCOP scorefxn_in, core::Size const jump/*=1*/, utility::vector1<core::Size> const & chain_ids );
	void apply (Pose & pose) override;
	void clean_pose( Pose & pose_copy );
	void compute_rmsd_with_super( Pose const & pose, Real & input_rmsd_super, Real & input_rmsd );
	virtual void calculate( Pose const & pose_in );
	virtual void report_ddG( std::ostream & out ) const;
	virtual Real sum_bound() const;
	virtual Real sum_unbound() const;
	virtual Real sum_ddG() const;

	core::Size rb_jump() const { return rb_jump_; }
	void rb_jump( core::Size j ) { rb_jump_ = j; }
	core::Size repeats() const { return repeats_; }
	void repeats( core::Size repeats ) { repeats_ = repeats; }
	core::Real translate_by() const { return translate_by_; }
	void translate_by( core::Real translate_by ) { translate_by_ = translate_by; }
	///@brief Get the hard-set chain ids. Note these are not necessarily the only chains to be moved.
	utility::vector1<core::Size> chain_ids() const { return chain_ids_; }
	std::set< core::Size > get_movable_jumps( core::pose::Pose const & pose ) const;
	~ddG() override;
	protocols::moves::MoverOP fresh_instance() const override { return (protocols::moves::MoverOP) utility::pointer::make_shared< ddG >(); }
	protocols::moves::MoverOP clone() const override;
	void parse_my_tag(  utility::tag::TagCOP, basic::datacache::DataMap & ) override;

	protocols::moves::MoverOP relax_mover() const{ return relax_mover_; }
	void relax_mover( protocols::moves::MoverOP m ){ relax_mover_ = m; }
	protocols::filters::FilterOP filter() const;
	void filter( protocols::filters::FilterOP f );

	void task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
	core::pack::task::TaskFactoryOP task_factory() const { return task_factory_; }
	void use_custom_task( bool uct ) { use_custom_task_ = uct; }
	bool use_custom_task() const { return use_custom_task_; }

	inline void repack_unbound( bool rpu ) { repack_unbound_ = rpu; }
	inline bool repack_unbound() const { return repack_unbound_; }

	void repack_bound( bool rpb ) { repack_bound_ = rpb; }
	bool repack_bound() const { return repack_bound_; }

	void relax_bound( bool rlb ) { relax_bound_ = rlb; }
	bool relax_bound() const { return relax_bound_; }

	inline void relax_unbound( bool const rlu ) { relax_unbound_ = rlu; }
	inline bool relax_unbound() const { return relax_unbound_; }

	inline void compute_rmsd( bool const crms ) { compute_rmsd_ = crms; }
	inline bool compute_rmsd() const { return compute_rmsd_; }

	virtual void scorefxn( core::scoring::ScoreFunctionCOP scorefxn_in );

	void dump_pdbs( bool dp ) { dump_pdbs_ = dp; }
	bool dump_pdbs() const { return dump_pdbs_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	utility::tag::XMLSchemaComplexTypeGeneratorOP
	define_ddG_schema();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private :

	/// @brief Helper method to appropriately form unbound complex. Returns false if monomer.
	bool unbind(Pose & pose) const;
	void setup_task(Pose const & pose);

	//fd: a few functions specific to solvated mode
	void duplicate_waters_across_jump( Pose & pose, core::Size jumpnum ) const;
	void setup_solvated_task( Pose const & pose, bool bound );
	void do_minimize( Pose & pose ) const;

	std::map< ScoreType, Real > bound_energies_;
	std::map< ScoreType, Real > unbound_energies_;

	std::map< core::Size, Real > bound_per_residue_energies_;
	std::map< core::Size, Real > unbound_per_residue_energies_;

	//Real bound_total_energy_;
	//Real unbound_total_energy_;
	core::Size repeats_;
	core::scoring::ScoreFunctionOP scorefxn_;
	void fill_energy_vector( Pose const & pose, std::map<ScoreType, Real > & energy_map );
	void fill_per_residue_energy_vector(Pose const & pose, std::map<core::Size,Real> & energy_map);
	core::Size rb_jump_;
	utility::vector1<core::Size> chain_ids_;
	utility::vector1<std::string> chain_names_;
	bool per_residue_ddg_;
	bool repack_unbound_;
	protocols::moves::MoverOP relax_mover_; //dflt NULL
	// Instead of score function delta, use the filter.
	protocols::filters::FilterOP filter_; //dflt NULL

	core::pack::task::TaskFactoryOP task_factory_;
	core::pack::task::PackerTaskOP task_;
	bool use_custom_task_;
	bool repack_bound_;
	bool relax_bound_;
	bool relax_unbound_;
	bool solvate_, solvate_unbound_, solvate_rbmin_, min_water_jump_; //fd solvate the bound and unbound states

	/// info carrier for poisson-boltzmann potential energy computation
	core::energy_methods::PBLifetimeCacheOP pb_cached_data_;

	/// true when PB potential is part of scorefxn
	bool pb_enabled_;
	/// distance in A to separate moledules
	core::Real translate_by_; //dflt set to 1000. Default resets to 100 for RosettaScripts with a scorefxn containing PB_elec.

	// native pose read if provided via command line
	core::pose::PoseOP native_;

	// number of waters in final bound and unbound state
	core::Size bound_HOH_, bound_HOH_V_, unbound_HOH_, unbound_HOH_V_;
	bool compute_rmsd_;
	core::Real bound_rmsd_, bound_rmsd_super_;

	// dump pdbs for debugging
	bool dump_pdbs_; //dflt false

};

} // movers
} // protocols
#endif /*INCLUDED_protocols_protein_interface_design_movers_ddG_HH*/
