// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mean_field/GenMeanFieldMover.hh
/// @brief Reports to Tracer a probability distribution per rotamer, per repackable/designable res, based on a mean-field method
/// can also report entropy per residue or per all repackable/designable res
/// if some residues are designable, will report logo plot as well
/// @author Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_GenMeanFieldMover_hh
#define INCLUDED_protocols_mean_field_GenMeanFieldMover_hh

// Unit headers
#include <protocols/mean_field/GenMeanFieldMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/mean_field/MeanField.fwd.hh>
#include <protocols/mean_field/AAMatrix.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <utility/vector1.hh>
#include <protocols/mean_field/jagged_array.hh>

namespace protocols {
namespace mean_field {

class GenMeanFieldMover : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	/// @brief default ctor
	GenMeanFieldMover();
	~GenMeanFieldMover() override;

	/// @brief performs main work of Mover
	void apply( core::pose::Pose & pose ) override;

	/// @brief reports the specificity profile (AAMatrix) and backbone boltzmann probs, as appropriate
	void report_aa_matrix() const;

	/// @brief reports rotamer probabilities (RotMatrix)
	void report_rot_prob() const;

	/// @brief reports distances between predicted AAMatrix and am
	void report_dist( AAMatrix const & am ) const;

	/// @brief reads in a AAMatrix (generally the experimental gold standard with which to compare the predicted)
	void read_aa_matrix() ;

	/// @brief reads in the input pdbs
	void read_input_pdbs();

	/// @brief prepares the task and unbinds the poses if necessary
	void prepare_task_poses( core::pose::Pose const & pose );

	/// @brief instantiates the appropriate MF class and calls processing methods
	void calc_mean_field();

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	//parsing functions
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	std::string get_name() const override;

	static std::string mover_name();

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//Accessors/mutators
	inline
	core::pack::task::TaskFactoryOP task_factory() const
	{
		return task_factory_;
	}

	inline
	void task_factory( core::pack::task::TaskFactoryOP task_factory )
	{
		task_factory_ = task_factory;
	}

	inline
	core::scoring::ScoreFunctionOP scorefxn() const{
		return scorefxn_;
	}

	inline
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
		scorefxn_ = scorefxn;
	}

	inline
	core::Real lambda_memory() const{
		return lambda_memory_;
	}

	inline
	void lambda_memory( core::Real const lm )
	{
		lambda_memory_ = lm;
	}

	inline
	core::Real tolerance() const{
		return tolerance_;
	}

	inline
	void tolerance( core::Real const tol )
	{
		tolerance_ = tol;
	}

	inline
	core::Real temperature() const{
		return temperature_;
	}

	inline
	void temperature( core::Real const temp )
	{
		temperature_ = temp;
	}


	inline
	core::Real threshold() const{
		return threshold_;
	}

	inline
	void threshold( core::Real const threshold )
	{
		threshold_ = threshold;
	}

	inline
	core::Size init_option() const{
		return init_option_;
	}

	inline
	void init_option( core::Size const iopt )
	{
		init_option_ = iopt;
	}

	inline
	bool unbound() const
	{
		return unbound_;
	}

	inline
	void unbound( bool const ub ) {
		unbound_ = ub;
	}

	inline
	std::string type() const
	{
		return type_;
	}

	inline
	void type(std::string const & s)
	{
		type_ = s;
	}

	//default value methods
	static inline
	core::Real default_value_for_lambda_memory()
	{
		return core::Real( 0.5 );
	}

	static inline
	core::Real default_value_for_tolerance()
	{
		return core::Real( 0.0001 );
	}

	static inline
	core::Real default_value_for_temperature()
	{
		return core::Real( 0.8 );
	}

	static inline
	core::Size default_value_for_init_option()
	{
		return core::Size( 1 );
	}

	static inline
	core::Real default_value_for_threshold()
	{
		return core::Real( 10.0 );
	}

	static inline
	bool default_value_for_unbound()
	{
		return bool( false );
	}

private:
	//options
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real lambda_memory_;
	core::Real tolerance_;
	core::Real temperature_; //this is actually kT rather than just T
	core::Size init_option_; //init rot_matrix to 1/nrot_per_res
	core::Real threshold_; //threshold to use in calculating energies matrix
	bool unbound_; //calculate mean-field matrix in unbound state
	std::string type_;

	//member variables for storing data used for scoring
	//no changes made to pose through mutable variables
	core::pack::task::PackerTaskOP task_;
	utility::vector1 < core::pack::task::PackerTaskOP > tasks_;
	core::pose::PoseOPs poses_;
	utility::vector1 < std::string > pdb_list_;

	protocols::mean_field::MeanFieldOP mean_field_;
	protocols::mean_field::AAMatrix aa_matrix_;
	protocols::mean_field::AAMatrix pred_aa_matrix_;

};

} // mean_field
} // protocols

#endif //INCLUDED_protocols_mean_field_GenMeanFieldMover_HH_

