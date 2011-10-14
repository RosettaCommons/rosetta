// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/loophash/FastGapMover.hh
/// @brief protocols for closing gaps
/// @author Ken Jung <kenjung@uw.ed>


#ifndef INCLUDED_protocols_loophash_FastGapMover_hh
#define INCLUDED_protocols_loophash_FastGapMover_hh


namespace protocols {
namespace idealize {


/// @brief Mover class for transforming a Pose to ideal bonds of the given Pose.
/// The idea is that this Mover stochastically picks a move-able position, forces
/// that position into ideal geometry, and tries to use minimization to bring
/// the coordinates back to very near their starting points.

class FastGapMover : public moves::Mover {

public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

public:
	FastGapMover():
		Mover("FastGapMover"),
		atom_pair_constraint_weight_( 0.0 ),
		coordinate_constraint_weight_( 0.01 ),
		fast_( false ),
		chainbreaks_( false ),
		report_CA_rmsd_(true)
	{}

	/// @brief clone has to be overridden only if clone invocation is expected.
	virtual moves::MoverOP clone() const {
		return new FastGapMover( *this );
	}

	virtual moves::MoverOP fresh_instance() const {
		return new FastGapMover;
	}

	FastGapMover &
	atom_pair_constraint_weight( core::Real const setting )
	{
		atom_pair_constraint_weight_ = setting;
		return *this;
	}


	FastGapMover &
	coordinate_constraint_weight( core::Real const setting )
	{
		coordinate_constraint_weight_ = setting;
		return *this;
	}

	FastGapMover &
	fast( bool const setting )
	{
		fast_ = setting;
		return *this;
	}

	FastGapMover &
	chainbreaks( bool const setting )
	{
		chainbreaks_ = setting;
		return *this;
	}

	FastGapMover &
	report_CA_rmsd( bool const setting )
	{
		report_CA_rmsd_ = setting;
		return *this;
	}

	void
	apply( Pose & pose );

	virtual std::string get_name() const;

	/// @brief Sets the list of residue positions to idealize.
	void set_pos_list( utility::vector1< Size > pos_list ) {
		pos_list_ = pos_list;
	}

	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

private:
	// methods
void find_next_gap( Pose & pose, Size & idx, Real & gap_distance );

std::map< std::string, core::pose::Pose >
	poses_from_cmd_line(utility::vector1< std::string > const & fn_list);

private:
	// data
	utility::vector1< Size > pos_list_;

	Real atom_pair_constraint_weight_;
	Real coordinate_constraint_weight_;
	bool fast_;
	bool chainbreaks_;
	bool report_CA_rmsd_;
};

} // idealize
} // protocols

#endif
