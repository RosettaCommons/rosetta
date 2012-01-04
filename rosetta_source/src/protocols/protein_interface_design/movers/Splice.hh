// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/Splice.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_Splice_hh
#define INCLUDED_protocols_protein_interface_design_movers_Splice_hh
#include <protocols/protein_interface_design/movers/Splice.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

//@brief lightweight class containing bb torsions
class BBDofs : public utility::pointer::ReferenceCount
{
	public:
		BBDofs() : resid_( 0 ), phi_( 0.0 ), psi_( 0.0 ), omega_( 0.0 ), resn_( "" ) {}
		BBDofs( core::Size const resid, core::Real const phi, core::Real const psi, core::Real const omega, std::string const resn ) : resid_( resid ), phi_( phi ), psi_( psi ), omega_( omega ), resn_( resn ) {}
		core::Size resid() const{ return resid_; }
		core::Real phi() const{ return phi_; }
		core::Real psi() const{ return psi_; }
		core::Real omega() const{ return omega_; }
		std::string resn() const {return resn_; }
		void resid( core::Size const r ){ resid_ = r; }
		void phi( core::Real const p ){ phi_ = p; }
		void psi( core::Real const p ){ psi_ = p; }
		void omega( core::Real const o ){ omega_ = o; }
		void resn( std::string const r ){ resn_ = r; }
		~BBDofs();
	private:
		core::Size resid_;
		core::Real phi_, psi_, omega_;
		std::string resn_;
};

//@ container for BBDofs
class ResidueBBDofs : public utility::pointer::ReferenceCount
{
	public:
		typedef utility::vector1< BBDofs > bbdof_list;
		typedef bbdof_list::iterator iterator;
		typedef bbdof_list::const_iterator const_iterator;

		ResidueBBDofs() : cut_site_( 0 ){ clear(); }
		~ResidueBBDofs();
		void cut_site( core::Size const c ){ cut_site_ = c; }
		core::Size cut_site() const { return cut_site_; }
		void clear() { bbdofs_.clear(); }
		void push_back( BBDofs const b ){ bbdofs_.push_back( b ); }
		const_iterator begin() const{ return bbdofs_.begin(); }
		const_iterator end() const{ return bbdofs_.end(); }
		iterator begin(){ return bbdofs_.begin(); }
		iterator end(){ return bbdofs_.end(); }
		core::Size size() const{ return bbdofs_.size(); }
		BBDofs & operator[]( int const i ) { return bbdofs_[ i ]; }
	private:
		core::Size cut_site_;
		bbdof_list bbdofs_;
};


/// @brief designs alanine residues in place of the residue identities at the interface. Retains interface glycines and prolines.
class Splice : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	Splice();
	void apply( Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new Splice ); }
		void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~Splice();
	void from_res( core::Size const f ){ from_res_ = f; }
	core::Size from_res() const { return from_res_; }
	void to_res( core::Size const t ){ to_res_ = t; }
	core::Size to_res() const { return to_res_; }
	std::string source_pdb() const { return source_pdb_; }
	void source_pdb( std::string const s ){ source_pdb_ = s; }
	void ccd( bool const c ){ ccd_ = c;}
	bool ccd() const { return ccd_; }
	void scorefxn( core::scoring::ScoreFunctionOP sf );
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::Real rms_cutoff() const{ return rms_cutoff_; }
	void rms_cutoff( core::Real const r ){ rms_cutoff_ = r; }
	void res_move( core::Size const r ){ res_move_ = r; }
	core::Size res_move() const{ return res_move_; }
	void randomize_cut( bool const r ){ randomize_cut_ = r; }
	bool randomize_cut() const{ return randomize_cut_; }
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );
	std::string torsion_database_fname() const{ return torsion_database_fname_; }
	void torsion_database_fname( std::string const d ){ torsion_database_fname_ = d; }
	core::Size database_entry()const {return database_entry_; }
	void database_entry( core::Size const d ){ database_entry_ = d; }
	void read_torsion_database();
	utility::vector1< ResidueBBDofs > torsion_database() const{ return torsion_database_; }
	void torsion_database( utility::vector1< ResidueBBDofs > const d ){ torsion_database_ = d; }

private:
	core::Size from_res_, to_res_;
	std::string source_pdb_;
	bool ccd_;//dflt true; do ccd?
	core::scoring::ScoreFunctionOP scorefxn_; //dflt score12 with reweighted sheet weight
	core::Real rms_cutoff_; //dflt 99999; after splicing, checks the average displacement of Ca atoms in the source and target segments. Failure leads to mover failure and no output
	core::Size res_move_; //dflt 4; how many residues to allow to move during ccd
	bool randomize_cut_; //dflt false; true: place cut in a randomly chosen loop residue, if available. false: place cut at loop's end
	core::pack::task::TaskFactoryOP task_factory_; // dflt NULL; Another access point to setting which residues to splice. This works at present only with one segment, so you set designable residues and Splice will then determine the first and last residues among these and splice that section out.
	std::string torsion_database_fname_; //dflt ""; set to true in order to read directly from a torsion database
	core::Size database_entry_; //dflt 0; in which case tests a random entry in each apply
	utility::vector1< ResidueBBDofs > torsion_database_;
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_Splice_HH*/
