// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief LoopRelaxMover.hh
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_LoopRelaxMover_hh
#define INCLUDED_protocols_comparative_modeling_LoopRelaxMover_hh

#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/LoopsFileIO.fwd.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <string>

#include <utility/vector1.hh>

#include <protocols/loops/Loops.hh>


namespace protocols {
namespace comparative_modeling {

class LoopRelaxMover : public moves::Mover {

public:

	// constructor
	LoopRelaxMover();

	/// @brief Alternative constructor.
	/// BE WARNED: THIS CONSTRUCTOR DOES NOT CALL SET_DEFAULTS().
	/// AS A RESULT, THE SCORE FUNCTIONS (AMONG OTHER THINGS) WILL
	/// NOT BE INITIALIZED.
	LoopRelaxMover(
		std::string remodel,
		std::string intermedrelax,
		std::string refine,
		std::string relax,
		loops::Loops const & loops
	);

	/// @brief Unresolved-loop-index constructor.
	/// BE WARNED: THIS CONSTRUCTOR DOES NOT CALL SET_DEFAULTS().
	/// AS A RESULT, THE SCORE FUNCTIONS (AMONG OTHER THINGS) WILL
	/// NOT BE INITIALIZED.
	LoopRelaxMover(
		std::string remodel,
		std::string intermedrelax,
		std::string refine,
		std::string relax,
		loops::LoopsFileData const & loops_from_file
	);

	/// @brief GuardedLoopsOP constructor.  This constructor copies the pointer
	/// to a GuardedLoops object; Loop indices must be resolved by the time
	/// that the guarded_loops_ object is accessed, but, they need not
	/// have been resolved at the time of this object's construction.
	/// BE WARNED: THIS CONSTRUCTOR DOES NOT CALL SET_DEFAULTS().
	/// AS A RESULT, THE SCORE FUNCTIONS (AMONG OTHER THINGS) WILL
	/// NOT BE INITIALIZED.
	LoopRelaxMover(
		std::string remodel,
		std::string intermedrelax,
		std::string refine,
		std::string relax,
		loops::GuardedLoopsFromFileOP guarded_loops
	);

	/// @brief Copy-ctor; Shallow copy of all data.
	LoopRelaxMover( LoopRelaxMover const & src );

	/// @brief assignment operator; Shallow copy of all data.
	LoopRelaxMover & operator = ( LoopRelaxMover const & rhs );


	//destructor
	~LoopRelaxMover() override;

	/// @brief Apply the loop-relax protocol to a Pose.  At the call to this function, the
	/// loop indices originating from the loop file (which may list indices with PDB identifiers
	/// i.e. res+insertioncode+chain ) will be resolved. Until this point, the Loops object
	/// in the LoopRelaxMover cannot be used.
	void apply( core::pose::Pose & pose ) override;

	// getters and setters
	void scorefxns( core::scoring::ScoreFunctionOP cen_scorefxn, core::scoring::ScoreFunctionOP fa_scorefxn );
	void fa_scorefxn( core::scoring::ScoreFunctionOP fa_scorefxn );
	void cen_scorefxn( core::scoring::ScoreFunctionOP cen_scorefxn );

	void cmd_line_csts( bool setting ) {
		cmd_line_csts_ = setting;
	}

	void copy_sidechains( bool setting ) {
		copy_sidechains_ = setting;
	}

	void rebuild_filter( core::Real setting ) {
		rebuild_filter_ = setting;
	}

	void n_rebuild_tries( core::Size setting ) {
		n_rebuild_tries_ = setting;
	}

	void remodel( std::string val ) {
		remodel_ = val;
	}

	void intermedrelax( std::string val ) {
		intermedrelax_ = val;
	}

	void refine( std::string val ) {
		refine_ = val;
	}

	void relax( std::string val ) {
		relax_ = val;
	}

	/// @brief Must be called before the Loops data can be read from.
	void resolve_loopfile_indices( core::pose::Pose const & pose );

	/// @brief set the Loops object (with resolved indices) directly.  This would override
	/// the unresolved-indices held in the LoopsFileData if that data were set in the constructor
	void loops( protocols::loops::LoopsOP const val );

	/// @brief Set the loop file data.  This will require that a Pose be presented to
	/// the LoopRelax object before get_loops() can be called.
	void loops_file_data( loops::LoopsFileData const & loopfiledata );

	void frag_libs( utility::vector1< core::fragment::FragSetOP > new_libs );

	bool cmd_line_csts() const {
		return cmd_line_csts_;
	}

	bool copy_sidechains() const {
		return copy_sidechains_;
	}

	core::Size n_rebuild_tries() const {
		return n_rebuild_tries_;
	}

	core::Real rebuild_filter() const {
		return rebuild_filter_;
	}

	/// @brief Loops accessor.  May only be retrieved after the loop indices have been resolved
	/// in a call to apply.
	protocols::loops::LoopsOP get_loops();

	/// @brief Loops accessor.  May only be retrieved after the loop indices have been resolved
	/// in a call to apply.
	protocols::loops::LoopsCOP get_loops() const;

	std::string remodel () const {
		return remodel_;
	}

	std::string relax() const {
		return relax_;
	}

	std::string intermedrelax() const {
		return intermedrelax_;
	}

	std::string refine() const {
		return refine_;
	}

	utility::vector1< core::fragment::FragSetOP > frag_libs() const;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void compute_rmsd( bool const c ){ compute_rmsd_ = c; }
	bool compute_rmsd() const { return compute_rmsd_; }

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	// void task_factory( core::pack::task::TaskFactoryOP tf ); /// currently taskfactory is not supported
	// core::pack::task::TaskFactoryOP task_factory() const;

private:
	void set_defaults_();

	bool cmd_line_csts_;
	bool copy_sidechains_;

	core::Size n_rebuild_tries_;
	core::Real rebuild_filter_;

	std::string remodel_;
	std::string intermedrelax_;
	std::string refine_;
	std::string relax_;

	loops::GuardedLoopsFromFileOP guarded_loops_;

	core::scoring::ScoreFunctionOP cen_scorefxn_;
	core::scoring::ScoreFunctionOP fa_scorefxn_;

	utility::vector1< core::fragment::FragSetOP > frag_libs_;
	bool compute_rmsd_; //dflt true; but set to false if you change the pose length between start and rmsd measurements;
	// core::pack::task::TaskFactoryOP task_factory_; //dflt NULL; a way to restrict packing from outside the mover
};

} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_LoopRelaxMover_HH

