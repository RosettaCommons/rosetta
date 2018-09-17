// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/Splice.hh
/// @author Gideon Lapidoth (glapidoth@gmail.com)

#ifndef INCLUDED_protocols_splice_SpliceOut_hh
#define INCLUDED_protocols_splice_SpliceOut_hh

#include <protocols/splice/SpliceOut.fwd.hh>
#include <protocols/splice/SpliceSegment.fwd.hh>
#include <protocols/splice/SpliceManager.hh>
#include <protocols/splice/Splice.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <core/sequence/SequenceProfile.hh>
#include <protocols/task_operations/SeqprofConsensusOperation.fwd.hh>
#include <basic/database/open.hh>
#include <core/pack/task/PackerTask.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <protocols/moves/Mover.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>



namespace protocols {
namespace splice {



class SpliceOut : virtual public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
	typedef utility::vector1< ResidueBBDofs >::const_iterator dbase_const_iterator;
	//This data structure is always updated in Splice and holds all the information that Splice helper function will use
public:
	SpliceOut();

	virtual void apply( Pose & pose ) override;
	virtual std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new SpliceOut ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & /* pose*/ ) override;
	virtual ~SpliceOut();
	std::string source_pdb() const { return source_pdb_; }
	void source_pdb( std::string const & s ){ source_pdb_ = s; }
	void scorefxn( core::scoring::ScoreFunctionOP sf );
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::Real rms_cutoff() const{ return rms_cutoff_; }
	void rms_cutoff( core::Real const r ){ rms_cutoff_ = r; }
	core::Real rms_cutoff_loop() const{ return rms_cutoff_loop_; }
	void rms_cutoff_loop( core::Real const r ){ rms_cutoff_loop_ = r; }
	void cut_secondarystruc( bool const r){ cut_secondarystruc_ =r; }
	bool cut_secondarystruc() const{ return cut_secondarystruc_; }
	//void set_fold_tree(core::pose::Pose & pose, core::Size const vl_vh_cut) ;
	bool thread_original_sequence() const { return thread_original_sequence_; }
	void thread_original_sequence( bool const s ){ thread_original_sequence_ = s; }
	std::string dbase_file_name() const;
	void dbase_file_name( std::string const & f );
	protocols::filters::FilterOP splice_filter() const;
	void splice_filter( protocols::filters::FilterOP f );
	void ccd_mover( core::pose::Pose & pose,core::kinematics::MoveMapOP mm);
	void min_mover( core::pose::Pose & pose,core::kinematics::MoveMapOP mm);
	void tail_mover( core::pose::Pose & pose,core::kinematics::MoveMapOP mm);
	int SpliceOutFilter(core::pose::Pose * pose) ;
	int SpliceOutRMSDFilter(core::pose::Pose * pose);

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static
	std::string
	mover_name( );


	/// sequence profiles
	/// Splice changes the backbone of the current pose and the following methods deal with dynamically constructing a
	/// sequence profile for the current backbone choices.
	//void read_splice_segments( std::string const & segment_type, std::string const & segment_name, std::string const & file_name );
	void load_pdb_segments_from_pose_comments( core::pose::Pose  const  & p); // get the segment names for those segments that are constant in this splice function
	//void modify_pdb_segments_with_current_segment( std::string const & pdb_name ); // set the current segment name
	//void add_sequence_constraints( core::pose::Pose & pose ); // add SequenceProfileConstraints based on the sequence profile
	/// @brief add dihedral constraint to grafted loop according to source pdb dihedral angles
	virtual void superimpose_source_on_pose( core::pose::Pose const & pose, core::pose::Pose & source_pose );
	bool rb_sensitive() const{ return rb_sensitive_; }
	void rb_sensitive( bool const r ){ rb_sensitive_ = r;}
	void delete_hairpin( bool const b ){ delete_hairpin_ = b; }
	bool delete_hairpin() const{ return delete_hairpin_; }
	core::Size delete_hairpin_n() const{ return delete_hairpin_n_; }
	void delete_hairpin_n( core::Size const c ){ delete_hairpin_n_ = c;}
	core::Size delete_hairpin_c() const{ return delete_hairpin_c_; }
	void delete_hairpin_c( core::Size const c ){ delete_hairpin_c_ = c;}
	void remove_hairpin( core::pose::Pose & pose ) const; // the function that actually removes the hairpin
	void check_source_segment() ;
	void copy_dofs_from_source();
	virtual void place_cut_site_in_segment(core::pose::Pose const & pose);
	void superimposed(bool const b){superimposed_=b;}
	bool superimposed(){return superimposed_;}
	//source_from_res, source_to_res
	void source_from_res(core::Size const c){source_from_res_=c;}
	core::Size source_from_res(){return source_from_res_;}
	void source_to_res(core::Size const c){source_to_res_=c;}
	core::Size source_to_res(){return source_to_res_;}
	void handle_mover_tag(TagCOP const tag,protocols::moves::Movers_map const & movers);
	virtual void set_fold_tree_nodes(core::pose::Pose const & pose);
	virtual void set_source_from_to_res();
	virtual core::Size set_anchor_res();
	virtual void build_ideal_segment(core::pose::Pose & pose);
	virtual void set_loop_length_change(protocols::protein_interface_design::movers::LoopLengthChange & llc);
	virtual void write_database_to_file(core::pose::Pose const & pose);
	void parse_SpliceOut_tags(TagCOP const tag,protocols::moves::Movers_map const & movers,protocols::filters::Filters_map const & filters);
	virtual void abstract_parse_tag(TagCOP const tag);
	void randomize_cut(bool const b){randomize_cut_=b;}
	bool randomize_cut(){return randomize_cut_;}
	void minimize_segment(core::pose::Pose & pose);
	virtual std::string name_for_filter();


protected:
	void save_values(); // call at beginning of apply. Used to keep the from_res/to_res values, which might be changed by apply during a run
	void retrieve_values(); // call at end of apply
	//void copy_stretch( core::pose::Pose & target, core::pose::Pose const & source, core::Size const from_res, core::Size const to_res );
	core::Size find_non_active_site_cut_site(core::pose::Pose const & pose);
	//void chainbreak_check( core::pose::Pose const & pose , core::Real const tolerance , bool fail_retry_if_found , bool crash_if_found );
	void check_pose_pssm_match(core::pose::PoseOP pose, utility::vector1<core::sequence::SequenceProfileOP> profile_vector);

protected: // Should be private!

	// is done by user defined order
	protocols::moves::MoverOP submover_;
	core::Real design_shell_ = 6.0;
	core::Real repack_shell_ = 8.0;
	core::scoring::ScoreFunctionOP scorefxn_; //dflt score12 with reweighted sheet weight
	// start - the starting pose for replacing the torsions at the start
	core::pose::PoseOP source_pose_;
	std::map < std::string, std::string> protein_family_to_database_;
	std::map< std::string/*1AHW*/, std::string/*L1.1*/ > pdb_to_H3_seq_map_; /* This object stores the H3 seqeunces of all PDBs in the database. The logic for this is that the H3, we build pssm from sequnce and blosum on the fly*/
	bool ignore_chain_break_ = false; // if we want to ignore the checking if the source PDB has a chainbreak
	core::Real tolerance_ = 0.23;
	std::vector < std::string> mover_type_ = {"MinMover","LoopMover_Refine_CCD","TailSegmentMover"};
	utility::vector1<std::array<int, 3>> fold_tree_nodes_;
	bool write_to_database_ = true;

	SpliceManager splicemanager;

	void (SpliceOut::*call_mover)(core::pose::Pose & pose,core::kinematics::MoveMapOP mm);
	void (SpliceOut::*call_mover_tail_)(core::pose::Pose & pose,core::kinematics::MoveMapOP mm);

private:

	core::Size saved_from_res_ = 0, saved_to_res_ = 0;
	std::string source_pdb_ = "";
	// after splicing, checks the average displacement of Ca atoms in the source and target segments.
	// Failure leads to mover failure and no output
	core::Real rms_cutoff_ = 999999;
	core::Real rms_cutoff_loop_;
	/// true: place cut in a randomly chosen loop residue, if available. false: place cut at loop's end
	bool randomize_cut_ = false;
	bool cut_secondarystruc_ = false; /// true: allows placing the cut within secondary structures
	core::pack::task::TaskFactoryOP task_factory_;
	core::kinematics::FoldTreeOP saved_fold_tree_;
	bool thread_original_sequence_ = false;//To force the original segment seqeunce on the pose segment (Gideon Lapdioth, 170814)
	//bool rtmin_ = true;//whether or not to let splice do rtmin following design (Ask assaf alon)
	std::string dbase_file_name_ = ""; // a file name into which the loop database is dumped
	// dflt NULL; to communicate the current Splice mover's loop origin to the GenericMC
	utility::pointer::shared_ptr< basic::datacache::DataMapObj< std::string > > mover_tag_;
	protocols::filters::FilterOP splice_filter_;
	// A map form protein family name to the order of the segments (eg. <"antibodies",<L1_L2,L3,H1_H2,H3>>)
	std::map< std::string, utility::vector1< std::string > > order_segments_;
	// dflt false; should we impose the RB dof of the current pose on the template before finding aligned residues.
	// (Ask Christoffer)
	bool rb_sensitive_ = false;
	std::map < std::string, std::string> database_segment_map_;//map between antibody segment and database file, e.g. <L1_L2,"l1_l2.db">
	bool CG_const_ = false;// if set to true then We aplly CG constraints from source pose onto pose according to PSSM rules,Gideon Aug14
	bool delete_hairpin_ = false; // if true cut out the top of the hairpin leaving two stumps
	core::Size delete_hairpin_n_ = 0, delete_hairpin_c_ = 0; // the n and c points for cutting; n_ is counted from the cysteine and c_ is counted from the next cutsite
	core::Size source_from_res_ = 0, source_to_res_ = 0;//the stem positions on the source pdb
	bool superimposed_ = false; //Is the source protein aligned to the template?,Gideon 20Nov16

};

} //splice
} //protocols

#endif //INCLUDED_protocols_splice_Splice_hh
