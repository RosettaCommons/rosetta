// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/SpliceManager.hh
/// @author Gideon Lapidoth

#ifndef INCLUDED_protocols_splice_SpliceManager_hh
#define INCLUDED_protocols_splice_SpliceManager_hh

#include <protocols/splice/SpliceManager.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <map>
#include <protocols/splice/SpliceSegment.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/id/SequenceMapping.hh>
#include <utility/tag/Tag.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/datacache/DataMapObj.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/rosetta_scripts/util.hh>




namespace protocols {
namespace splice {
typedef utility::tag::TagCOP TagCOP;
typedef utility::pointer::shared_ptr< basic::datacache::DataMapObj< bool > > DataccacheBoolDataOP;


//@brief lightweight class containing bb torsions and residue identities
class BBDofs : public utility::pointer::ReferenceCount
{
public:
	BBDofs() : resid_( 0 ), phi_( 0.0 ), psi_( 0.0 ), omega_( 0.0 ), resn_( "" ){}
	BBDofs( core::Size const resid, core::Real const phi, core::Real const psi, core::Real const omega, std::string const & resn ) : resid_( resid ), phi_( phi ), psi_( psi ), omega_( omega ), resn_( resn ){}
	core::Size resid() const{ return resid_; }
	core::Real phi() const{ return phi_; }
	core::Real psi() const{ return psi_; }
	core::Real omega() const{ return omega_; }
	std::string resn() const {return resn_; }
	void resid( core::Size const r ){ resid_ = r; }
	void phi( core::Real const p ){ phi_ = p; }
	void psi( core::Real const p ){ psi_ = p; }
	void omega( core::Real const o ){ omega_ = o; }
	void resn( std::string const& r ){ resn_ = r; }
	~BBDofs() override;
private:
	core::Size resid_; /// this is currently not used in splice
	core::Real phi_, psi_, omega_;
	std::string resn_;
};

class ResidueBBDofs : public utility::pointer::ReferenceCount
{
public:
	typedef utility::vector1< BBDofs > bbdof_list;
	typedef bbdof_list::iterator iterator;
	typedef bbdof_list::const_iterator const_iterator;

	ResidueBBDofs() : cut_site_( 0 ), start_loop_( 0 ), stop_loop_( 0 ), source_pdb_(""), aa_sequence_(""), dssp_("") { clear(); }
	~ResidueBBDofs() override;
	void cut_site( core::Size const c ){ cut_site_ = c; }
	core::Size cut_site() const { return cut_site_; }
	void clear() { bbdofs_.clear(); }
	void push_back( BBDofs const & b ){ bbdofs_.push_back( b ); }
	const_iterator begin() const{ return bbdofs_.begin(); }
	const_iterator end() const{ return bbdofs_.end(); }
	iterator begin(){ return bbdofs_.begin(); }
	iterator end(){ return bbdofs_.end(); }
	core::Size size() const{ return bbdofs_.size(); }// -1 is because the size of bbdofs also includes the last entry which is the name of the source pdb. I want only the number of residues
	BBDofs & operator[]( int const i ) { return bbdofs_[ i ]; }
	core::Size start_loop() const{ return start_loop_; }
	void start_loop( core::Size const s ){ start_loop_ = s; }
	core::Size stop_loop() const{ return stop_loop_; }
	void stop_loop( core::Size const s ){ stop_loop_ = s; }
	std::string source_pdb() const{ return source_pdb_; }
	void source_pdb( std::string const & s ){ source_pdb_ = s; }
	std::string tail_segment() const{ return tail_segment_; }
	void tail_segment( std::string const & s ){ tail_segment_ = s; }
	core::Size disulfide() const{ return disulfide_; }
	void disulfide( core::Size const s ){ disulfide_ = s; }
	void aa_sequence( std::string const & s ){ aa_sequence_ = s; }
	std::string aa_sequence() const{ return aa_sequence_; }
	void dssp( std::string const & s ){ dssp_ = s; }
	std::string dssp() const{ return dssp_; }
	bbdof_list bbdofs(){return bbdofs_;}



private:
	core::Size cut_site_, start_loop_, stop_loop_, disulfide_/*what is the disulfide on the template*/;
	bbdof_list bbdofs_;
	std::string source_pdb_, tail_segment_/*either n or c*/; // the source pdb from which the loop is taken
	std::string aa_sequence_, dssp_;


};


class SpliceManager : public utility::pointer::ReferenceCount
{

public:
	SpliceManager();
	virtual ~SpliceManager();
	void fold_tree( core::pose::Pose & pose,utility::vector1<std::array<int, 3>> positions ) const;

	//void superimpose_source_on_pose( core::pose::Pose const & pose, core::pose::Pose & source_pose);
	void set_BB_dofs(core::pose::Pose & pose);
	core::sequence::SequenceProfileOP generate_sequence_profile(core::pose::Pose & pose);
	void generate_profile_from_seq(std::map< std::string/*1AHW*/, std::string/*L1.1*/ >  segment_seq_map, std::string segmentName, std::string sourcename,SpliceSegmentOP SpliceSegment);
	void add_sequence_constraints(core::pose::Pose & pose);
	void rb_adjust_template(core::pose::Pose const & pose) const ;
	void adjust_stem_positions_by_template(core::pose::Pose const & pose);
	core::Size template_from_res(){return template_from_res_;}
	//can only be set once
	void template_from_res(core::Size const i){ template_from_res_=i;}
	//can only be set once
	core::Size template_to_res(){return template_to_res_;}
	void template_to_res(core::Size const i){template_to_res_=i;}
	core::Size pose_from_res(){return pose_from_res_;}
	void pose_from_res(core::Size const i){pose_from_res_=i;}
	core::Size pose_to_res(){return pose_to_res_;}
	void pose_to_res(core::Size const i){pose_to_res_=i;}
	void chain_num(core::Size const i){chain_num_=i;}
	core::Size chain_num() const {return chain_num_;}
	void use_sequence_profiles(bool const b){use_sequence_profiles_=b;}
	core::Size use_sequence_profiles(){return use_sequence_profiles_;}
	void add_sequence_constraints_only(bool const b){add_sequence_constraints_only_=b;}
	bool add_sequence_constraints_only(){return add_sequence_constraints_only_;}
	void segment_type(std::string const & s){segment_type_=s;}
	std::string segment_type(){return segment_type_;}
	void source_pdb_name(std::string const & s){source_pdb_name_=s;}
	std::string source_pdb_name(){return source_pdb_name_;}
	void pdb_segments(std::map< std::string, std::string > const & m){pdb_segments_=m;}
	std::map< std::string, std::string >& pdb_segments(){return pdb_segments_;}
	void segment_names_ordered (utility::vector1< std::string > SegNameOrder){segment_names_ordered_ = SegNameOrder;} //setter for segment_name_ordered
	utility::vector1< std::string >& segment_names_ordered() {return segment_names_ordered_;} //getter for segment_name_ordered
	std::map< std::string, SpliceSegmentOP > & splice_segments(){return splice_segments_;}
	void splice_segments(std::map< std::string, SpliceSegmentOP > const & m){splice_segments_=m;}
	core::pose::PoseOP template_pose(){return template_pose_;}
	void template_pose(core::pose::PoseOP const pose){template_pose_=pose;}
	void rb_sensitive(bool const b){rb_sensitive_=b;}
	bool rb_sensitive()const {return rb_sensitive_;}
	void dofs(ResidueBBDofs const & d){dofs_=d;}
	ResidueBBDofs& dofs() {return dofs_;}
	core::Real profile_weight_away_from_interface() const { return profile_weight_away_from_interface_;}
	void profile_weight_away_from_interface(core::Real const p) {profile_weight_away_from_interface_ = p;}
	void scorefxn(core::scoring::ScoreFunctionOP sf) {scorefxn_ = sf;}
	core::scoring::ScoreFunctionOP scorefxn() const {return scorefxn_;}
	void cut_site(core::Size const c){cut_site_=c;}
	core::Size cut_site(){return cut_site_;}
	void residue_diff(core::Size const c){residue_diff_=c;}
	int residue_diff(){return residue_diff_;}
	std::string tail_seg(){return tail_seg_;}
	void tail_seg(std::string const & s){tail_seg_=s;}
	core::kinematics::MoveMapOP mm(){return mm_;}
	void mm(core::kinematics::MoveMapOP const m){mm_=m;}
	void update_pose_stem_positions();
	void check_sequence_profile(core::pose::Pose & pose, core::id::SequenceMappingOP smap, core::sequence::SequenceProfileOP seqprof);
	void parse_tags(TagCOP const tag,basic::datacache::DataMap &data);
	core::pack::task::TaskFactoryOP task_factory() const {return task_factory_;}
	void task_factory(core::pack::task::TaskFactoryOP tf) {task_factory_ = tf;}
	void dbase_file_name(std::string const & s) {dbase_file_name_ = s;}
	std::string dbase_file_name() const {return dbase_file_name_;}
	void design_shell(core::Real const c) {design_shell_ = c;}
	void repack_shell(core::Real const c) {repack_shell_ = c;}
	core::Real design_shell() const {return design_shell_ ;}
	core::Real repack_shell() const {return repack_shell_ ;}
	std::string template_file() const{ return template_file_; }
	void template_file( std::string const & s ){ template_file_ = s; }
	bool thread_original_sequence() const { return thread_original_sequence_; }
	void thread_original_sequence( bool const s ){ thread_original_sequence_ = s; }
	bool debug() const { return debug_; }
	void debug( bool const s ){ debug_ = s; }
	std::string mover_name() const{ return mover_name_; }
	void mover_name( std::string const & s ){ mover_name_ = s; }
	void modify_pdb_segments_with_current_segment(std::string const & pdb_name);
	void chainbreak_check( core::pose::Pose const & pose , core::Real const tolerance , bool fail_retry_if_found , bool crash_if_found);
	core::Size allowed_cuts() const { return allowed_cuts_; }
	void allowed_cuts( core::Size const s ){ allowed_cuts_ = s; }
	void add_dihedral_constraints( core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size nearest_to_from,core::Size nearest_to_to );
	void add_coordinate_constraints( core::pose::Pose & pose, core::pose::Pose const & source_pose,core::Size nearest_to_from,core::Size nearest_to_to, core::Size anchor,std::string atom_type="CA",core::pack::task::PackerTaskOP task=NULL);
	void parse_segments(utility::vector1<TagCOP> const & sub_tags,TagCOP const tag, basic::datacache::DataMap &data);
	core::Size anchor_res(){return anchor_res_;}
	void anchor_res(core::Size const c){anchor_res_=c;}
	bool adjust_stem_positions_by_template_;
private:
	core::Size template_from_res_;//Should remain constant
	core::Size template_to_res_;//Should remain constant
	core::Size pose_from_res_;
	core::Size pose_to_res_;
	std::string  tail_seg_;
	int residue_diff_;
	core::Size chain_num_;
	ResidueBBDofs dofs_;
	bool use_sequence_profiles_; // dflt false;
	std::string segment_type_; //dflt ""; what segment is this? Used to decide which profiles to use. examples, L1,L2,L3
	bool add_sequence_constraints_only_;
	std::string source_pdb_name_;
	std::map< std::string/*which segment (L1,L2...)*/, std::string/*pdb name*/ > pdb_segments_;
	std::map< std::string, SpliceSegmentOP > splice_segments_;
	core::Real profile_weight_away_from_interface_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< std::string > segment_names_ordered_;
	core::pose::PoseOP template_pose_;
	bool rb_sensitive_;
	core::Size cut_site_;//position of chain break on segment
	core::kinematics::MoveMapOP mm_;
	core::pack::task::TaskFactoryOP task_factory_;
	std::string dbase_file_name_; //dflt ""; a file name into which the loop database is dumped
	core::Real design_shell_; //dflt 6.0 gideonla
	core::Real repack_shell_; //dflt 8.0 gideonla
	std::string mover_name_;//for debugging puposes, gets the mover name to add to dumped pdb during debugging.
	bool debug_;
	std::string template_file_;
	bool thread_original_sequence_;//To force the original segment seqeunce on the pose segment (Gideon Lapdioth, 170814)
	core::Size allowed_cuts_; //dflt=1. if the number of chain_breaks in the psoe is more than the allowed number then fail
	core::Size anchor_res_;


};//end SpliceManager definition


//@brief lightweight class containing bb torsions and residue identities

} //splice
} //protocols

#endif //INCLUDED_protocols_splice_SpliceManager_hh
