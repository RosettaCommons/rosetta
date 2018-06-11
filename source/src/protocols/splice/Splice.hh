// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/splice/Splice.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// Gidoen Lapidoth (glapidoth@gmail.com)
/// @brief Not used anymore. samples BB conformation from a PDB and pallies it to the pose

#ifndef INCLUDED_protocols_splice_Splice_hh
#define INCLUDED_protocols_splice_Splice_hh

#include <protocols/splice/Splice.fwd.hh>
#include <protocols/splice/SpliceSegment.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <protocols/task_operations/SeqprofConsensusOperation.fwd.hh>
#include <basic/database/open.hh>
#include <core/pack/task/PackerTask.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <protocols/splice/SpliceManager.hh>


namespace protocols {
namespace splice {


//@brief used for checking backbone segment compatibility
class BBMatching : public utility::pointer::ReferenceCount
{
public:
	BBMatching();
	void read_partners_from_file(std::string const & pdb_name,std::string const & file_name);
	virtual ~BBMatching();
	core::Size size(){return size_;}
	std::string show(){return pdb_names_;}
	std::map< std::string,  std::vector<std::string>  > pdb_to_next_allowed_partners(){return pdb_to_next_allowed_partners_;};
	std::map< std::string,  std::vector<std::string>  > pdb_to_prev_allowed_partners(){return pdb_to_prev_allowed_partners_;};
private:
	std::map< std::string,  std::vector<std::string>  > pdb_to_next_allowed_partners_;
	std::map< std::string, std::vector<std::string>  > pdb_to_prev_allowed_partners_;
	core::Size size_;
	std::string pdb_names_;
};

/*class SpliceData : public utility::pointer::ReferenceCount
{
public:
SpliceData();
virtual ~SpliceData();
core::Size from_res;
core::Size to_res;
core::Size chain_num;
};*/
/// @brief container for BBDofs, providing a convenient operator [], size, other methods and iterators that allow splice to treat
/// ResidueBBDofs as a simple vector (even though it contains other elements as well)


class Splice : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
	typedef utility::vector1< ResidueBBDofs >::const_iterator dbase_const_iterator;
	typedef utility::pointer::shared_ptr< basic::datacache::DataMapObj< bool > > DataccacheBoolDataOP;
public:
	Splice();
	void apply( Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new Splice ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	~Splice() override;

	void from_res( core::Size const f ){ from_res_ = f; } //stem residue of segment
	core::Size from_res() const { return from_res_; }
	void to_res( core::Size const t ){ to_res_ = t; }//stem residue of segment
	core::Size to_res() const { return to_res_; }
	std::string source_pdb() const { return source_pdb_; } //PDB to copy source conformation from
	void source_pdb( std::string const & s ){ source_pdb_ = s; }
	void ccd( bool const c ){ ccd_ = c;} //whether or not to apply CCD
	void design_shell(core::Real const c) {design_shell_ = c;} // Design shell radius around segment
	void repack_shell(core::Real const c) {repack_shell_ = c;} // packing shell radius around segment
	bool ccd() const { return ccd_; }
	core::Real dihedral_const() const { return dihedral_const_; } //getter for dihedral_const option
	core::Real coor_const() const {return coor_const_ ;} //coordinate constraints on pose
	core::Real design_shell() const {return design_shell_ ;}
	core::Real repack_shell() const {return repack_shell_ ;}
	void scorefxn( core::scoring::ScoreFunctionOP sf );
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::Real rms_cutoff() const{ return rms_cutoff_; } //rmsd cut off between designed segment and source segment
	void rms_cutoff( core::Real const r ){ rms_cutoff_ = r; }
	core::Real rms_cutoff_loop() const{ return rms_cutoff_loop_; } //rmsd cut off between designed segment and source segment - only for loop SS
	void rms_cutoff_loop( core::Real const r ){ rms_cutoff_loop_ = r; }
	void res_move( core::Size const r ){ res_move_ = r; } //how many residues to move during CCD - depracted
	core::Size res_move() const{ return res_move_; }
	void randomize_cut( bool const r ){ randomize_cut_ = r; } // where to insert cut in segment before CCD - aplly at random
	bool randomize_cut() const{ return randomize_cut_; }
	void set_cut_res( core::Size const r ){ set_cut_res_ = r; } // where to insert cut in segment before CCD
	core::Size set_cut_res() const{ return set_cut_res_; }
	void cut_secondarystruc( bool const r){ cut_secondarystruc_ =r; } //whether or not to cut in SS (non loop)
	bool cut_secondarystruc() const{ return cut_secondarystruc_; }
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );
	core::pack::task::TaskFactoryOP design_task_factory() const;
	void design_task_factory( core::pack::task::TaskFactoryOP tf );
	void set_segment_names_ordered (utility::vector1< std::string > SegNameOrder){segment_names_ordered_ = SegNameOrder;} //setter for segment_name_ordered
	utility::vector1< std::string > get_segment_names_ordered() const {return segment_names_ordered_;} //getter for segment_name_ordered
	void set_dofs_pdb_name (std::string dofsPDBname) {dofs_pdb_name=dofsPDBname;}
	std::string get_dofs_pdb_name() const {return dofs_pdb_name;}
	std::string torsion_database_fname() const{ return torsion_database_fname_; } //File where BB conformtion of segment is stored
	void torsion_database_fname( std::string const & d ){ torsion_database_fname_ = d; }
	core::Size database_entry()const {return database_entry_; }
	void database_entry( core::Size const d ){ database_entry_ = d; }
	void read_torsion_database();
	utility::vector1< ResidueBBDofs > torsion_database() const{ return torsion_database_; }
	void torsion_database( utility::vector1< ResidueBBDofs > const & d ){ torsion_database_ = d; }
	std::string template_file() const{ return template_file_; } //PDB file that is used as reference for determining for_res and to_res
	void template_file( std::string const & s ){ template_file_ = s; }
	void poly_ala( bool const p ){ poly_ala_ = p; } //whether to use poly ala in designed segment
	bool poly_ala() const{ return poly_ala_; }
	void equal_length( bool const e ){ equal_length_ = e; }
	bool equal_length() const{ return equal_length_; }
	void fold_tree( core::pose::Pose & pose, core::Size const start, core::Size const stop, core::Size const cut ) const;
	void tail_fold_tree(core::pose::Pose & pose, core::Size const vl_vh_cut,core::Size chain_break,std::string segment_name) const ;
	void set_fold_tree(core::pose::Pose & pose, core::Size const vl_vh_cut) ;
	bool design() const{ return design_; }
	void design( bool const d ) { design_ = d; }
	void delta_lengths( utility::vector1< int > const & dl ){ delta_lengths_ = dl; } //how much should the designed segment change relative to
	utility::vector1< int > delta_lengths() { return delta_lengths_; }
	bool dbase_iterate() const { return dbase_iterate_; }
	void dbase_iterate( bool const d ){ dbase_iterate_ = d; }
	utility::vector1< core::Size >::const_iterator dbase_begin() const;
	utility::vector1< core::Size >::const_iterator dbase_end() const;
	core::Size find_dbase_entry( core::pose::Pose const & pose ); // returns a dbase entry
	core::Size locked_res() const;
	void locked_res( core::Size const r );
	void locked_res_id( char const c ); //in casse we want to use with inverse rotamer, not used
	char locked_res_id() const;
	std::string checkpointing_file() const;
	void checkpointing_file( std::string const & cf );
	bool rtmin() const { return rtmin_; }
	void rtmin( bool const d ){ rtmin_ = d; }
	bool allow_all_aa() const { return allow_all_aa_; }
	void allow_all_aa( bool const d ){ allow_all_aa_ = d; }
	bool thread_original_sequence() const { return thread_original_sequence_; }
	void thread_original_sequence( bool const s ){ thread_original_sequence_ = s; }
	void load_from_checkpoint(); // load relevant internal data during a checkpoint recovery
	void save_to_checkpoint() const; // save relevant data for future checkpoint recovery
	std::string loop_dbase_file_name() const;
	void loop_dbase_file_name( std::string const & f );
	void loop_pdb_source( std::string const & l ); //deprected
	std::string loop_pdb_source() const;
	protocols::filters::FilterOP splice_filter() const;
	void splice_filter( protocols::filters::FilterOP f ); //which filter to use to check mover
	void database_pdb_entry( std::string const & s ){ database_pdb_entry_ = s; }
	std::string database_pdb_entry() const { return database_pdb_entry_; }
	//void minimize_segment(core::pose::Pose & pose);

	/// sequence profiles
	/// Splice changes the backbone of the current pose and the following methods deal with dynamically constructing a
	/// sequence profile for the current backbone choices.
	void read_splice_segments( std::string const & segment_type, std::string const & segment_name, std::string const & file_name );
	core::sequence::SequenceProfileOP generate_sequence_profile(core::pose::Pose & pose);
	void load_pdb_segments_from_pose_comments( core::pose::Pose  const  & p); // get the segment names for those segments that are constant in this splice function
	void modify_pdb_segments_with_current_segment( std::string const & pdb_name ); // set the current segment name
	void add_sequence_constraints( core::pose::Pose & pose ); // add SequenceProfileConstraints based on the sequence profile
	/// @brief add dihedral constraint to grafted loop according to source pdb dihedral angles
	void add_dihedral_constraints( core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size nearest_to_from,core::Size nearest_to_to );
	void add_coordinate_constraints( core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size nearest_to_from,core::Size nearest_to_to, core::Size anchor,std::string atom_type="CA",core::pack::task::PackerTaskOP task=nullptr);

	void profile_weight_away_from_interface( core::Real const p ); //change weight of PSSM as cutoff
	core::Real profile_weight_away_from_interface() const;
	bool restrict_to_repacking_chain2() const{ return restrict_to_repacking_chain2_; }
	void restrict_to_repacking_chain2( bool const r ){ restrict_to_repacking_chain2_ = r; }
	core::Size get_current_seg() { return current_segment_pos; } ; // getter for the current segemnt number that is being designed
	//bool check_aa(std::string s, utility::vector1<core::Real > profRow);
	bool add_sequence_constraints_only() const{ return add_sequence_constraints_only_; }
	void add_sequence_constraints_only( bool const a ){ add_sequence_constraints_only_ = a; }

	bool rb_sensitive() const{ return rb_sensitive_; } //if RB sensitive then modify the VL/VH orientation of the protein (only relevant for antibodies)
	void rb_sensitive( bool const r ){ rb_sensitive_ = r;}

	void rb_adjust_template( core::pose::Pose const & pose ) const; // adjust the template_pose_ according the rb state of the current pose (if rb_sensitive is on)
	void protein_family( std::string const & s) {protein_family_=s;}; //setter of the protein_family tag option
	core::Size chain_num()const {return chain_num_;};
	void chain_num( core::Size const c ) {chain_num_ = c;};

	bool skip_alignment() const{ return skip_alignment_; }
	void skip_alignment( bool const b ){ skip_alignment_ = b; }

	bool superimposed() const{ return superimposed_; } //is the pose aligned to the template PDB
	void superimposed( bool const b ){ superimposed_ = b; }

	core::Size source_pdb_from_res() const{ return source_pdb_from_res_; }
	void source_pdb_from_res( core::Size const s ){ source_pdb_from_res_ = s; }

	core::Size source_pdb_to_res() const{ return source_pdb_to_res_; }
	void source_pdb_to_res( core::Size const s ){ source_pdb_to_res_ = s; }
	void delete_hairpin( bool const b ){ delete_hairpin_ = b; } //delete hairpin residues
	bool delete_hairpin() const{ return delete_hairpin_; }
	core::Size delete_hairpin_n() const{ return delete_hairpin_n_; }
	void delete_hairpin_n( core::Size const c ){ delete_hairpin_n_ = c;}
	core::Size delete_hairpin_c() const{ return delete_hairpin_c_; }
	void delete_hairpin_c( core::Size const c ){ delete_hairpin_c_ = c;}
	void remove_hairpin( core::pose::Pose & pose ) const; // the function that actually removes the hairpin

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	void save_values(); // call at beginning of apply. Used to keep the from_res/to_res values, which might be changed by apply during a run
	void retrieve_values(); // call at end of apply
	std::string parse_pdb_code(std::string pdb_file_name);
	void copy_stretch( core::pose::Pose & target, core::pose::Pose const & source, core::Size const from_res, core::Size const to_res );
	core::Size find_non_active_site_cut_site(core::pose::Pose const & pose);
	void superimpose_source_on_pose( core::pose::Pose const &, core::pose::Pose & );
	void chainbreak_check( core::pose::Pose const & pose , core::Real const tolerance , bool fail_retry_if_found , bool crash_if_found );
	void correct_bb_O_placements(core::pose::Pose Pose,core::Size cut_site);
	// This vector will hold the segment names by order so when the segments are concatenated into a single profile it
	// is done by user defined order
	utility::vector1< std::string > segment_names_ordered_;

	std::string dofs_pdb_name; //This variable hold the name of the pdb in the torsion db
	core::Size from_res_, to_res_, saved_from_res_, saved_to_res_;
	std::string source_pdb_;
	bool ccd_; //dflt true; do ccd?
	core::Real dihedral_const_; //dflt 1; gideonla
	core::Real coor_const_; //dflt 1; gideonla
	core::Real design_shell_; //dflt 6.0 gideonla
	core::Real repack_shell_; //dflt 8.0 gideonla
	core::scoring::ScoreFunctionOP scorefxn_; //dflt score12 with reweighted sheet weight

	// dflt 99999; after splicing, checks the average displacement of Ca atoms in the source and target segments.
	// Failure leads to mover failure and no output
	core::Real rms_cutoff_;
	core::Real rms_cutoff_loop_;

	core::Size res_move_; //dflt 4; how many residues to allow to move during ccd

	// dflt false; true: place cut in a randomly chosen loop residue, if available. false: place cut at loop's end
	bool randomize_cut_;

	bool cut_secondarystruc_; //dflt false; true: allows placing the cut within secondary structures
	core::Size set_cut_res_; //dftl 0, mainly for debugging. Allows the user to set which residue will be used as the cut site.
	// dflt NULL; Another access point to setting which residues to splice.
	// This works at present only with one segment, so you set designable residues and Splice will then determine the
	// first and last residues among these and splice that section out.
	core::pack::task::TaskFactoryOP task_factory_;

	// dflt NULL; a task_factory used to restrict design during splicing.
	// A 'good' idea for this is to define the aligned segments through RestrictToAlignedSegments and send those to
	// this task_factory.  During splicing, this task_factory will be used to restrict the design operations in
	// addition to what DesignAroundOperation determines as the designable residues.  So, by applying the user-defined
	// RestrictToAlignedSegments as well as dao, you get design on the spliced segment + its vicinity in other aligned
	// segments, and repack in a slightly larger shell.
	core::pack::task::TaskFactoryOP design_task_factory_;

	std::string torsion_database_fname_; //dflt ""; set to true in order to read directly from a torsion database
	core::Size database_entry_; //dflt 0; in which case tests a random entry in each apply
	std::string database_pdb_entry_; // dflt ""; e.g., "1yihl" specify this only if you want just one loop to be spliced
	utility::vector1< ResidueBBDofs > torsion_database_;
	utility::vector1< ResidueBBDofs > tail_torsion_database_;


	// dflt ""; which source file to use as the template to determine what from_res() and to_res() refer to.
	// The input structure may change during a trajectory and so from_res() and to_res() might lose their sense.
	// If this is "", the input file is taken to be template
	std::string template_file_;

	// dflt true; thread ala residues in each position other than Gly/Pro or conserved in the source pdb.
	// If false, keeps the input sequence (except Gly/Pro, which are replaced)
	bool poly_ala_;

	bool equal_length_; // dflt false; restrict threading to loops equal in length to the original

	// template - relative to what is the torsion dbase computed (1x9q)
	// start - the starting pose for replacing the torsions at the start
	core::pose::PoseOP template_pose_, start_pose_, source_pose_;
	core::kinematics::FoldTreeOP saved_fold_tree_;

	bool design_; //dflt false; design all non-pro/gly residues in template
	utility::vector1< int > delta_lengths_; // dflt empty; change loop length by how much? 0 is always assumed
	bool dbase_iterate_; //dflt false;
	bool allow_all_aa_;//to allow all amino acids (Ask assaf alon)
	bool thread_original_sequence_;//To force the original segment seqeunce on the pose segment (Gideon Lapdioth, 170814)
	bool rtmin_;//whether or not to let splice do rtmin following design (Ask assaf alon)
	bool first_pass_; // dflt true;

	// indices to the subset of the dbase library over which multiple calls iterate
	utility::vector1< core::Size > dbase_subset_;

	// used if multiple calls to splice are made to iterate through the list
	utility::vector1< core::Size >::const_iterator current_dbase_entry_;

	// dflt false; this is a weird construct to allow placing the variable on the basic::datacache::DataMap
	DataccacheBoolDataOP end_dbase_subset_;

	// dflt NULL; a residue that serves as the root for a fold tree jump to the other chain.
	// This residue is expected to be within the loop span, and allows the loop to be refined while keeping the rigid
	// body jump between the two chains; it's a only ostensibly a vector, as it has to be compatible with placestub,
	// but it only looks at the first element of that vector
	utility::pointer::shared_ptr< basic::datacache::DataMapObj < utility::vector1< core::Size > > > locked_res_;

	char locked_res_id_; // dflt ''; the one-letter code for the locked residue

	// dflt ""; a file that contains checkpointing information to recover from job termination when iterating over a
	// loop database
	std::string checkpointing_file_;

	std::string loop_dbase_file_name_; //dflt ""; a file name into which the loop database is dumped

	// dflt ""; what is the source pdb from which the loop came?
	// This is used in writing the loop to the loop dbase, and helps keep track of where loops come from during design.
	std::string loop_pdb_source_;

	// dflt NULL; to communicate the current Splice mover's loop origin to the GenericMC
	utility::pointer::shared_ptr< basic::datacache::DataMapObj< std::string > > mover_tag_;

	protocols::filters::FilterOP splice_filter_;
	std::string Pdb4LetName_;
	std::string protein_family_;

	// A map form protein family name to the order of the segments (eg. <"antibodies",<L1_L2,L3,H1_H2,H3>>)
	std::map< std::string, utility::vector1< std::string > > order_segments_;
	core::Size chain_num_;

	//sequence profiles
	bool use_sequence_profiles_; // dflt false; set internally only, by whether or not the Segments are defined
	std::string segment_type_; //dflt ""; what segment is this? Used to decide which profiles to use. examples, L1,L2,L3
	core::Size current_segment_pos; // save the position of the segment_type_ in the segment name ordered vector.

	// stores sequence profiles for all possible segments (this doesn't change during a run), e.g., L1, ...; L2, ....
	std::map< std::string, SpliceSegmentOP > splice_segments_;

	// which pdb file did each segment in the current pose come from (used to build the current profile).
	// This uses the pose comment structure to retain the information through successive applies
	std::map< std::string/*which segment (L1,L2...)*/, std::string/*pdb name*/ > pdb_segments_;

	// dflt 1.0; you can define a different weight outside an 8A shell around the partner protein.
	// This should typically be set higher than 1.0, implying that the sequence profile carries a larger weight away
	// from the functional site
	core::Real profile_weight_away_from_interface_;

	bool restrict_to_repacking_chain2_; // dflt true; if false, does two-sided design during splice

	// dflt false; if true, only add constraints and return, don't do any splicing. (ask Assaf)
	bool add_sequence_constraints_only_;
	// dflt false; if true, only set fold tree no design no nothing.
	bool set_fold_tree_only_;

	// dflt false; should we impose the RB dof of the current pose on the template before finding aligned residues.
	// (Ask Christoffer)
	bool rb_sensitive_;
	std::map < std::string, std::string> protein_family_to_database_;
	bool skip_alignment_; // dflt false; use with care, ask Sarel. If you're not making any loop length changes and the alignment is perfect (you're simply splicing a same-length segment) drop the ambiguities of finding the alignment.

	std::map < std::string, std::string> database_segment_map_;//map between antibody segment and database file, e.g. <L1_L2,"l1_l2.db">
	core::Size vl_vh_cut;//keeps track of the residue number between the vl vh of the pose. this postion should change with loop length changes.
	std::string tail_segment_; //dflt ""; User should insert either "C"/"N" designating that that the N-termini tail of
	// the inserted segment should also be inserted and modeled. For Example if the user inserted
	// The L1_L2 segment from a source antibody and used "N" option than the residues that
	// are at the N terminal end of the segment all the way up to the first residue will also be added to the pose from the template
	//PDB
	bool min_seg_;//dflt false. if set to true then we perform minimization on bb and chi of segement after splice in,Gideon Aug14
	bool CG_const_;//dflt false. if set to true then We aplly CG constraints from source pose onto pose according to PSSM rules,Gideon Aug14
	bool debug_;//dflt false if set to true then the all sorts of dump pdb options are activated
	std::string mover_name_;//for debugging puposes, gets the mover name to add to dumped pdb during debugging.
	std::map< std::string/*1AHW*/, std::string/*L1.1*/ > pdb_to_H3_seq_map_; /* This object stores the H3 seqeunces of all PDBs in the database. The logic for this is that the H3, we build pssm from sequnce and blosum on the fly*/
	bool superimposed_; //dflt true; when doing splice out, is the source pdb superimposed on the template pdb? This is usually the case. If it isn't we expect the user to specify the stem sites on the source and template exactly. Only dihedral restraints will be used in this case.
	core::Size source_pdb_from_res_, source_pdb_to_res_; //dflt 0,0; see above superimposed.
	bool ignore_chain_break_;//dflt false, if we want to ignore the checking if the source PDB has a chainbreak
	core::Real tolerance_; //dflt 0.23
	bool delete_hairpin_; // dflt false; if true cut out the top of the hairpin leaving two stumps
	core::Size delete_hairpin_n_, delete_hairpin_c_; // the n and c points for cutting; n_ is counted from the cysteine and c_ is counted from the next cutsite
	core::Size allowed_cuts_; //dflt=1. if the number of chain_breaks in the psoe is more than the allowed number then fail
	bool enzdes_;//dflt 0, if set to true add Enzdes constraints before doing packing and RT_min,Gideon 2May15
	std::map< std::string, BBMatchingOP > bb_comp_db_;//this is a map between segment name and an object holding all possible conformation pairs (e.g. does pdb 1h04 from blade1 goes with pdb 2x3b from blade 2)
};

} //splice
} //protocols

#endif //INCLUDED_protocols_splice_Splice_hh
