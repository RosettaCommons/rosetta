// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/splice/SampleRotamersFromPDB.cc
/// @brief  Given source PDBs limits the rotamer sampling to those found in equivalent positions in the PDB
/// @detail At initialization the taskoperation reads a database with allowed rotamers for each position
/// and changes the Rotmer set accordingly
/// @author Gideon Lapidoth ( glapidoth@gmail.com )



#ifndef INCLUDED_protocols_splice_SampleRotamersFromPDB_hh
#define INCLUDED_protocols_splice_SampleRotamersFromPDB_hh


// Unit Headers
#include <protocols/splice/SampleRotamersFromPDB.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <map>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace splice {
typedef utility::vector1 <utility::vector1< core::conformation::ResidueOP > > res_matrix;



//typedef utility::vector1 <utility::vector1< ROT > > rot_matrix;

//@brief used for checking backbone segment compatibilty
class RotLibdb : public utility::pointer::ReferenceCount
{

public:
	RotLibdb();
	void read_profile(std::string const & file_name, std::string const & segment_name );//read profile file
	void add_segment_fname(std::string const & file_name, std::string const & segment_name ); //add segment name file name pair to map
	std::string get_segment_fname(std::string const & segment_name); // getter for  add_segment_fname
	rot_matrix fill_rotamer_matrix_from_db_file(std::string fname); //generate matrix object from files
	void read_pdb_profile( std::string const & file_name );// read the pdb-profile match from a disk file
	rot_matrix pdb_profile( std::string const & pdb_name ); //return matrix object from pdb_name

private:
	//std::map< std::string /*blade3*/,  res_matrix  > rot_matrix_;
	std::map< std::string, std::string  > pdb_to_profile_map_;
	//core::Size size_;
	std::string pdb_names_;
	std::map< std::string/*L3.10.1*/, rot_matrix > rot_matrix_profile_;
	std::map< std::string, std::string > ROTdb_segments_fileName_;

};

class SampleRotamersFromPDB_RotamerSetOperation : public core::pack::rotamer_set::RotamerSetOperation
{
public:


	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::rotamer_set::RotamerSet RotamerSet;
	typedef core::pack::rotamer_set::Rotamers Rotamers;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef utility::graph::GraphCOP GraphCOP;
	typedef core::pack::rotamer_set::RotamerSetOperationOP RotamerSetOperationOP;


public:


	SampleRotamersFromPDB_RotamerSetOperation();
	SampleRotamersFromPDB_RotamerSetOperation(bool add_rot, utility::vector1< core::Size > SampleAtAlignedpositions,bool d,bool ccd,std::string file_name);


	virtual ~SampleRotamersFromPDB_RotamerSetOperation();

	virtual
	RotamerSetOperationOP
	clone() const;

	virtual
	void
	alter_rotamer_set(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & /*sfxn*/,
		core::pack::task::PackerTask const & /*ptask*/,
		utility::graph::GraphCOP /*packer_neighbor_graph*/,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);
	//change the rotamer set from pdb file
	void alter_rotamer_set_from_pdb(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & ptask,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);
	//change the rotamer set from db file
	void alter_rotamer_set_from_db(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & ptask,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);
	void add_rotamer_to_rotamer_set(
		core::pose::Pose const & pose,
		core::pack::rotamer_set::RotamerSet & rotamer_set,
		core::conformation::ResidueOP cur_rot
	);
	//@brief rather than filtering rotamers that match those in the db file I can rescore them to be favoured by the packer

	void add_rotamer_constraints(
		core::pose::Pose & pose,
		core::Size seqnum,
		core::conformation::ResidueOP cur_rot
	);
	bool is_residue_allowed (
		core::chemical::ResidueType const & restype,
		core::pack::task::ResidueLevelTask const & rtask
	);
	core::conformation::ResidueOP ROT2res(ROT rot,core::pose::Pose const & pose);

public:

	void add_pose(core::pose::PoseCOP pose);
	void initialize_from_command_line();
	void fill_rotamer_matrix_from_db_file();
	void copy_rotamer_matrix(rot_matrix const & rm);
private: // data
	utility::vector1< core::pose::PoseCOP > poses_;
	bool add_rotamer_;//dflt true, if set to true add the pdb rotamer to the rotamer vector
	utility::vector1< core::Size > SampleAtAlignedpositions_;// specifies positions on the template pdb from which rotmaers will be sampled. similar to  RestrictIdentitiesAtAlignedPositions TO
	bool debug_; //dflt false. Setting to true spits out bunch of text to the tracer.
	bool ccd_;
	std::string db_file_;// the file name for the db text file
	rot_matrix resi_vec_;// hold the rotamer vectors given from db file


};


//////////////////////////////////////////////////////////////////////////////////////////////
class SampleRotamersFromPDB : public core::pack::task::operation::TaskOperation {
public:


	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;


public:


	/// @brief default constructor
	SampleRotamersFromPDB();
	/// @brief destructor
	virtual ~SampleRotamersFromPDB();

	/// @brief make clone
	virtual TaskOperationOP clone() const override;


public:


	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task) const override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "SampleRotamersFromPDB"; }

public:

	utility::vector1< core::Size > SampleAtAlignedpositions() const{ return SampleAtAlignedpositions_; }
	void SampleAtAlignedpositions( utility::vector1< core::Size > const s ){ SampleAtAlignedpositions_ = s; }
	void parse_tag(TagCOP tag , DataMap & ) override;
	rot_matrix combine_rot_dbs(core::pose::Pose const & pose)const ;
	void concatenate_rot_matrix(rot_matrix & a, rot_matrix const & b) const;

private: // data
	bool add_rotamer_;//dflt false, if set to true add the pdb rotamer to the rotamer vector
	utility::vector1< core::Size > SampleAtAlignedpositions_;// specifies positions on the template pdb from which rotmaers will be sampled. similar to  RestrictIdentitiesAtAlignedPositions TO
	bool debug_; //dflt false. Setting to true
	bool ccd_;

	std::string db_fname_;// the file name for the db text file

	std::map< std::string, RotLibdbOP > ROTdb_segments_;
	utility::vector1< std::string > segment_names_ordered_;
	//res_matrix profile_vector_;

};


} // protocols
} // splice


#endif
