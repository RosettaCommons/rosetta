// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/ddG.hh
/// @brief calculate interface ddg
/// @author Sarel Fleishman


#ifndef INCLUDED_protocols_simple_moves_ddG_hh
#define INCLUDED_protocols_simple_moves_ddG_hh

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <core/types.hh>
#include <protocols/simple_moves/ddG.fwd.hh>
#include <string>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



namespace protocols {
namespace simple_moves {

class ddG : public simple_moves::DesignRepackMover
{
public:
	typedef core::Real Real;
	typedef core::scoring::ScoreType ScoreType;
	typedef core::pose::Pose Pose;
public :
	ddG();
	ddG( core::scoring::ScoreFunctionCOP scorefxn_in, core::Size const jump=1, bool const symmetry=false );
	virtual void apply (Pose & pose);
	void calculate( Pose const & pose );
	void symm_ddG( core::pose::Pose const & pose_in );
	void no_repack_ddG(core::pose::Pose const & pose_in);
	void report_ddG( std::ostream & out ) const;
	Real sum_ddG() const;
	core::Size rb_jump() const { return rb_jump_; }
	void rb_jump( core::Size j ) { rb_jump_ = j; }
	virtual ~ddG();
	protocols::moves::MoverOP fresh_instance() const { return (protocols::moves::MoverOP) new ddG; }
	protocols::moves::MoverOP clone() const;
	void parse_my_tag(  utility::tag::TagPtr const, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const& );

	virtual std::string get_name() const;
	protocols::moves::MoverOP relax_mover() const{ return relax_mover_; }
	void relax_mover( protocols::moves::MoverOP m ){ relax_mover_ = m; }

private :
	std::map< ScoreType, Real > bound_energies_;
	std::map< ScoreType, Real > unbound_energies_;

	std::map< Size, Real > bound_per_residue_energies_;
	std::map< Size, Real > unbound_per_residue_energies_;

	Real bound_total_energy_;
	Real unbound_total_energy_;
	Size repeats_;
	core::scoring::ScoreFunctionOP scorefxn_;
	void fill_energy_vector( Pose const & pose, std::map<ScoreType, Real > & energy_map );
	void fill_per_residue_energy_vector(Pose const & pose, std::map<Size,Real> & energy_map);
	core::Size rb_jump_;
	utility::vector1<core::Size> chain_ids_;
	bool symmetry_;
	bool per_residue_ddg_;
	bool repack_;
	protocols::moves::MoverOP relax_mover_; //dflt NULL; in the unbound state, relax before taking the energy
};

} // movers
} // protocols
#endif /*INCLUDED_protocols_protein_interface_design_movers_ddG_HH*/

