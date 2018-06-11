// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RotLibOut.hh
/// @brief generate a database of rotamers from a list of pdbs in a text file. this is used
/// @author Gideon Lapidoth (glapidoth@gmail.com)

#ifndef INCLUDED_protocols_splice_RotLibOut_hh
#define INCLUDED_protocols_splice_RotLibOut_hh

#include <protocols/splice/RotLibOut.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <algorithm>
#include <core/kinematics/Jump.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// C++ Headers
namespace protocols {
namespace splice {

class RotLibOut : public protocols::moves::Mover {
public:
	RotLibOut();
	~RotLibOut();

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void set_min_dist( core::Real );//setter for min_dist_ variable
	core::Real get_min_dist();
	void Rotamer_db_file( std::string s ) {Rotamer_db_file_=s;}//setter for min_dist_ variable
	std::string Rotamer_db_file() {return Rotamer_db_file_;}
	void seq_aln_file( std::string s ) {seq_aln_file_=s;}//setter for min_dist_ variable
	std::string seq_aln_file(){ return seq_aln_file_;}
	void dump_pdb( bool b ) {dump_pdb_=b;}//setter for min_dist_ variable
	bool dump_pdb() {return dump_pdb_ ;}
	void min_frag_length( core::Size i ) {min_frag_length_=i;}//setter for min_dist_ variable
	core::Size min_frag_length() {return min_frag_length_;}

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ,
		protocols::filters::Filters_map const & ,
		protocols::moves::Movers_map const & ,
		core::pose::Pose const &  ) override;
	void find_matching_res( core::pose::Pose & pose, core::pose::Pose & hit_pose);
	void print();
	void print_lib_to_file();
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static std::string mover_name();


	//  std::string template_pdb_fname() const {return template_pdb_fname_; }
	//  void template_pdb_fname( std::string const s ){ template_pdb_fname_ = s; }

	//  void jump_dbase_fname( std::string const s ) { jump_dbase_fname_ = s ; }
	//  std::string jump_dbase_fname() const{ return jump_dbase_fname_; }
	//
	// bool jump_from_foldtree() const{ return jump_from_foldtree_; }
	// void jump_from_foldtree( bool const jfft ){ jump_from_foldtree_ = jfft;}

private:
	core::Real min_dist_;// user defined cutoff of distance between atoms (defines what is close)
	utility::vector1 <utility::vector1< core::conformation::ResidueOP > > resi_vec; //a 2D matrix holding the aligned residues for each template position.
	std::map< std::string, utility::vector1< std::string > > AA_vec; //a 2D matrix holding the aligned residues for each template position.
	utility::vector1< std::string > pdb_name_list;//holds the names of the input pdbs
	std::string Rotamer_db_file_;
	std::string seq_aln_file_;
	bool dump_pdb_;
	core::Size min_frag_length_;//only use fragments equal or greater than set length
};

} // splice
} // protocols

#endif
