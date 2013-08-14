// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/splice/Splice.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_devel_splice_Splice_hh
#define INCLUDED_devel_splice_Splice_hh

#include <devel/splice/Splice.fwd.hh>
#include <devel/splice/SpliceSegment.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/DataMapObj.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/moves/DataMapObj.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <protocols/toolbox/task_operations/SeqprofConsensusOperation.fwd.hh>

namespace devel {
namespace splice {

//@brief lightweight class containing bb torsions and residue identities
class BBDofs : public utility::pointer::ReferenceCount
{
	public:
		BBDofs() : resid_( 0 ), phi_( 0.0 ), psi_( 0.0 ), omega_( 0.0 ), resn_( "" ){}
		BBDofs( core::Size const resid, core::Real const phi, core::Real const psi, core::Real const omega, std::string const resn ) : resid_( resid ), phi_( phi ), psi_( psi ), omega_( omega ), resn_( resn ){}
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
		virtual ~BBDofs();
	private:
		core::Size resid_; /// this is currently not used in splice
		core::Real phi_, psi_, omega_;
		std::string resn_;
};

///@brief container for BBDofs, providing a convenient operator [], size, other methods and iterators that allow splice to treat
/// ResidueBBDofs as a simple vector (even though it contains other elements as well)
class ResidueBBDofs : public utility::pointer::ReferenceCount
{
	public:
		typedef utility::vector1< BBDofs > bbdof_list;
		typedef bbdof_list::iterator iterator;
		typedef bbdof_list::const_iterator const_iterator;

		ResidueBBDofs() : cut_site_( 0 ), start_loop_( 0 ), stop_loop_( 0 ), source_pdb_("") { clear(); }
		virtual ~ResidueBBDofs();
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
		core::Size start_loop() const{ return start_loop_; }
		void start_loop( core::Size const s ){ start_loop_ = s; }
		core::Size stop_loop() const{ return stop_loop_; }
		void stop_loop( core::Size const s ){ stop_loop_ = s; }
		std::string source_pdb() const{ return source_pdb_; }
		void source_pdb( std::string const s ){ source_pdb_ = s; }
	private:
		core::Size cut_site_, start_loop_, stop_loop_;
		bbdof_list bbdofs_;
		std::string source_pdb_; // the source pdb from which the loop is taken
};


/// @brief designs alanine residues in place of the residue identities at the interface. Retains interface glycines and prolines.
class Splice : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
	typedef utility::vector1< ResidueBBDofs >::const_iterator dbase_const_iterator;
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
	void design_shell(core::Real const c) {design_shell_ = c;}
	void repack_shell(core::Real const c) {repack_shell_ = c;}
	bool ccd() const { return ccd_; }
	core::Real dihedral_const() const { return dihedral_const_; } //getter for dihedral_const option
	core::Real coor_const() const {return coor_const_ ;}
	core::Real design_shell() const {return design_shell_ ;}
	core::Real repack_shell() const {return repack_shell_ ;}
	void scorefxn( core::scoring::ScoreFunctionOP sf );
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::Real rms_cutoff() const{ return rms_cutoff_; }
	void rms_cutoff( core::Real const r ){ rms_cutoff_ = r; }
	void res_move( core::Size const r ){ res_move_ = r; }
	core::Size res_move() const{ return res_move_; }
	void randomize_cut( bool const r ){ randomize_cut_ = r; }
	bool randomize_cut() const{ return randomize_cut_; }
	void cut_secondarystruc( bool const r){ cut_secondarystruc_ =r; }
	bool cut_secondarystruc() const{ return cut_secondarystruc_; }
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );
	core::pack::task::TaskFactoryOP design_task_factory() const;
	void design_task_factory( core::pack::task::TaskFactoryOP tf );
	void set_segment_names_ordered (utility::vector1< std::string > SegNameOrder){segment_names_ordered_ = SegNameOrder;} //setter for segment_name_ordered
	utility::vector1< std::string > get_segment_names_ordered() const {return segment_names_ordered_;} //getter for segment_name_ordered
	void set_dofs_pdb_name (std::string dofsPDBname) {dofs_pdb_name=dofsPDBname;}
	std::string get_dofs_pdb_name() const {return dofs_pdb_name;}
	std::string torsion_database_fname() const{ return torsion_database_fname_; }
	void torsion_database_fname( std::string const d ){ torsion_database_fname_ = d; }
	core::Size database_entry()const {return database_entry_; }
	void database_entry( core::Size const d ){ database_entry_ = d; }
	void read_torsion_database();
	utility::vector1< ResidueBBDofs > torsion_database() const{ return torsion_database_; }
	void torsion_database( utility::vector1< ResidueBBDofs > const d ){ torsion_database_ = d; }
	std::string template_file() const{ return template_file_; }
	void template_file( std::string const s ){ template_file_ = s; }
	void poly_ala( bool const p ){ poly_ala_ = p; }
	bool poly_ala() const{ return poly_ala_; }
	void equal_length( bool const e ){ equal_length_ = e; }
	bool equal_length() const{ return equal_length_; }
	void fold_tree( core::pose::Pose & pose, core::Size const start, core::Size const stop, core::Size const cut ) const;
	bool design() const{ return design_; }
	void design( bool const d ) { design_ = d; }
	void delta_lengths( utility::vector1< int > const dl ){ delta_lengths_ = dl; }
	utility::vector1< int > delta_lengths() { return delta_lengths_; }
	bool dbase_iterate() const { return dbase_iterate_; }
	void dbase_iterate( bool const d ){ dbase_iterate_ = d; }
	utility::vector1< core::Size >::const_iterator dbase_begin() const;
	utility::vector1< core::Size >::const_iterator dbase_end() const;
	core::Size find_dbase_entry( core::pose::Pose const & pose ); // returns a dbase entry
	core::Size locked_res() const;
	void locked_res( core::Size const r );
	void locked_res_id( char const c );
	char locked_res_id() const;
	std::string checkpointing_file() const;
	void checkpointing_file( std::string const cf );

	void load_from_checkpoint(); // load relevant internal data during a checkpoint recovery
	void save_to_checkpoint() const; // save relevant data for future checkpoint recovery

	std::string loop_dbase_file_name() const;
	void loop_dbase_file_name( std::string const f );
	void loop_pdb_source( std::string const l );
	std::string loop_pdb_source() const;
	protocols::filters::FilterOP splice_filter() const;
	void splice_filter( protocols::filters::FilterOP f );
	void database_pdb_entry( std::string const s ){ database_pdb_entry_ = s; }
	std::string database_pdb_entry() const { return database_pdb_entry_; }

/// sequence profiles
/// Splice changes the backbone of the current pose and the following methods deal with dynamically constructing a
/// sequence profile for the current backbone choices.
	void read_splice_segments( std::string const segment_type, std::string const segment_name, std::string const file_name );
	core::sequence::SequenceProfileOP generate_sequence_profile(core::pose::Pose & pose);
	void load_pdb_segments_from_pose_comments( core::pose::Pose  const  & p); // get the segment names for those segments that are constant in this splice function
	void modify_pdb_segments_with_current_segment( std::string const pdb_name ); // set the current segment name
	void add_sequence_constraints( core::pose::Pose & pose ); // add SequenceProfileConstraints based on the sequence profile
	/// @brief add dihedral constraint to grafted loop according to source pdb dihedral angles
	void add_dihedral_constraints( core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size nearest_to_from,core::Size nearest_to_to, core::Size cut_site  );
	void add_coordinate_constraints( core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size nearest_to_from,core::Size nearest_to_to);

	void profile_weight_away_from_interface( core::Real const p );
	core::Real profile_weight_away_from_interface() const;
	bool restrict_to_repacking_chain2() const{ return restrict_to_repacking_chain2_; }
	void restrict_to_repacking_chain2( bool const r ){ restrict_to_repacking_chain2_ = r; }
	core::Size get_current_seg() { return current_segment_pos; } ; // getter for the current segemnt number that is being designed
	bool check_aa(std::string s, utility::vector1<core::Real > profRow);
  bool add_sequence_constraints_only() const{ return add_sequence_constraints_only_; }
	void add_sequence_constraints_only( bool const a ){ add_sequence_constraints_only_ = a; }

private:
	void save_values(); // call at beginning of apply. Used to keep the from_res/to_res values, which might be changed by apply during a run
	void retrieve_values(); // call at end of apply
	utility::vector1< std::string > segment_names_ordered_;//This vector will hold the segment names by order so when the segments are concatented into a single profile it is done by user defined order
	std::string dofs_pdb_name; //This variable hold the name of the pdb in the torsion db
	core::Size from_res_, to_res_, saved_from_res_, saved_to_res_;
	std::string source_pdb_;
	bool ccd_;//dflt true; do ccd?
	core::Real dihedral_const_;//dflt 1; gideonla
	core::Real coor_const_;//dflt 1; gideonla
	core::Real design_shell_;//dflt 6.0 gideonla
	core::Real repack_shell_;//dflt 8.0 gideonla
	core::scoring::ScoreFunctionOP scorefxn_; //dflt score12 with reweighted sheet weight
	core::Real rms_cutoff_; //dflt 99999; after splicing, checks the average displacement of Ca atoms in the source and target segments. Failure leads to mover failure and no output
	core::Size res_move_; //dflt 4; how many residues to allow to move during ccd
	bool randomize_cut_; //dflt false; true: place cut in a randomly chosen loop residue, if available. false: place cut at loop's end
	bool cut_secondarystruc_; //dflt false; true: allows placing the cut within secondary structures
	core::pack::task::TaskFactoryOP task_factory_; // dflt NULL; Another access point to setting which residues to splice. This works at present only with one segment, so you set designable residues and Splice will then determine the first and last residues among these and splice that section out.
	core::pack::task::TaskFactoryOP design_task_factory_; // dflt NULL; a task_factory used to restrict design during splicing. A 'good' idea for this is to define the aligned segments through RestrictToAlignedSegments and send those to this task_factory. During splicing, this task_factory will be used to restrict the design operations in addition to what DesignAroundOperation determines as the designable residues. So, by applying the user-defined RestrictToAlignedSegments as well as dao, you get design on the spliced segment + its vicinity in other aligned segments, and repack in a slightly larger shell.
	std::string torsion_database_fname_; //dflt ""; set to true in order to read directly from a torsion database
	core::Size database_entry_; //dflt 0; in which case tests a random entry in each apply
	std::string database_pdb_entry_; // dflt ""; e.g., "1yihl" specify this only if you want just one loop to be spliced
	utility::vector1< ResidueBBDofs > torsion_database_;
	std::string template_file_; //dflt ""; which source file to use as the template to determine what from_res() and to_res() refer to. The input structure may change during a trajectory and so from_res() and to_res() might lose their sense. If this is "", the input file is taken to be template
	bool poly_ala_; /// dflt true; thread ala residues in each position other than Gly/Pro or conserved in the source pdb. If false, keeps the input sequence (except Gly/Pro, which are replaced)
	bool equal_length_; // dflt false; restrict threading to loops equal in length to the original
	core::pose::PoseOP template_pose_, start_pose_; // template - relative to what is the torsion dbase computed (1x9q); start - the starting pose for replacing the torsions at the start
	core::kinematics::FoldTreeOP saved_fold_tree_;
	bool design_; //dflt false; design all non-pro/gly residues in template
	utility::vector1< int > delta_lengths_; // dflt empty; change loop length by how much? 0 is always assumed
	bool dbase_iterate_; //dflt false;
	bool first_pass_; // dflt true;
	utility::vector1< core::Size > dbase_subset_; // indices to the subset of the dbase library over which multiple calls iterate
	utility::vector1< core::Size >::const_iterator current_dbase_entry_; // used if multiple calls to splice are made to iterate through the list
	utility::pointer::owning_ptr< protocols::moves::DataMapObj< bool > > end_dbase_subset_; // dflt false; this is a weird construct to allow placing the variable on the DataMap
	utility::pointer::owning_ptr< protocols::moves::DataMapObj < utility::vector1< core::Size > > > locked_res_; // dflt NULL; a residue that serves as the root for a fold tree jump to the other chain. This residue is expected to be within the loop span, and allows the loop to be refined while keeping the rigid body jump between the two chains; it's a only ostensibly a vector, as it has to be compatible with placestub, but it only looks at the first element of that vector
	char locked_res_id_; // dflt ''; the one-letter code for the locked residue
	std::string checkpointing_file_; // dflt ""; a file that contains checkpointing information to recover from job termination when iterating over a loop database
	std::string loop_dbase_file_name_; //dflt ""; a file name into which the loop database is dumped
	std::string loop_pdb_source_; //dflt ""; what is the source pdb from which the loop came? This is used in writing the loop to the loop dbase, and helps keep track of where loops come from during design.
	utility::pointer::owning_ptr< protocols::moves::DataMapObj< std::string > > mover_tag_; /// dflt NULL; to communicate the current Splice mover's loop origin to the GenericMC
	protocols::filters::FilterOP splice_filter_;
	std::string Pdb4LetName_;


///sequence profiles
	bool use_sequence_profiles_; // dflt false; set internally only, by whether or not the Segments are defined
	std::string segment_type_; //dflt ""; what segment is this? Used to decide which profiles to use. examples, L1,L2,L3
	core::Size current_segment_pos; // save the position of the segment_type_ in the segment name ordered vector.
	std::map< std::string, SpliceSegmentOP > splice_segments_; // stores sequence profiles for all possible segments (this doesn't change during a run), e.g., L1, ...; L2, ....
	std::map< std::string/*which segment (L1,L2...)*/, std::string/*pdb name*/ > pdb_segments_; // which pdb file did each segment in the current pose come from (used to build the current profile). This uses the pose comment structure to retain the information through successive applies
	core::Real profile_weight_away_from_interface_; //dflt 1.0; you can define a different weight outside an 8A shell around the partner protein. This should typically be set higher than 1.0, implying that the sequence profile carries a larger weight away from the functional site
	bool restrict_to_repacking_chain2_; // dflt true; if false, does two-sided design during splice
	bool add_sequence_constraints_only_; // dflt false; if true, only add constraints and return, don't do any splicing. (ask Assaf)
};

} //splice
} //devel

#endif //INCLUDED_devel_splice_Splice_hh
