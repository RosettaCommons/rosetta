// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/TaskOperations.hh
/// @brief  An operation to perform on a packer task --
///         usually, by a PackerTaskFactory right after the task's construction
///         This is an incomplete list.  Freely define your own TaskOperation and hand it
///         to a PackerTaskFactory.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_operation_TaskOperations_hh
#define INCLUDED_core_pack_task_operation_TaskOperations_hh

// Unit Headers
#include <core/pack/task/operation/TaskOperation.hh> // abstract base classes
#include <core/pack/task/operation/TaskOperations.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerLinks.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

class RestrictToRepacking : public TaskOperation {
public:
	typedef TaskOperation parent;

public:
	virtual ~RestrictToRepacking();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const & pose, PackerTask & task ) const;

	virtual void parse_tag( TagCOP, DataMap & );
	virtual void parse_def( utility::lua::LuaObject const & def);

};

/// @brief RestrictResidueToRepacking
class RestrictResidueToRepacking : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	virtual ~RestrictResidueToRepacking();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

	void include_residue( core::Size resid );
	void clear();

	virtual void parse_tag( TagCOP, DataMap & );

private:
	utility::vector1< core::Size > residues_to_restrict_to_repacking_;
};

// @brief  RestrictAbsentCanonicalAAS
class RestrictAbsentCanonicalAAS : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	RestrictAbsentCanonicalAAS();
	RestrictAbsentCanonicalAAS( core::Size resid, utility::vector1< bool > keep_aas );
	virtual ~RestrictAbsentCanonicalAAS();

	virtual TaskOperationOP clone() const;

	/// @brief if resid_ is zero, will apply RestrictAbsentCanonicalAAs to all residues
	/// in the pose. O/w only to resid_
	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

	/// @brief a human-readible string-based mutator
	void keep_aas( std::string const keep_aas );
	/// @brief direct vector1-basd mutator.  If an amino acid is not present (false) in the boolean vector, then do not allow it at this position.  The boolean vector is a 20-length vector in alphabetical order by one-letter code.
	void keep_aas( utility::vector1< bool > keep_aas );
	virtual void parse_tag( TagCOP, DataMap & );
	void include_residue( core::Size const resid );
private:
	core::Size resid_;
	utility::vector1< bool > keep_aas_;
};


/// @brief DisallowIfNonnative allows you to define what residues are NOT allowed in packing unless that residue is present in the input.  Behaves like RestrictAbsentCanonicalAAS and NOTAA except will allow  a resitricted residue at a position if it is there to begin with at the time of Task creation.  Will do all residues unless otherwise defined by selection syntax below
class DisallowIfNonnative:  public TaskOperation{

public:
	typedef TaskOperation parent;

public:
	DisallowIfNonnative();
	DisallowIfNonnative( utility::vector1< bool > disallowed_aas );
	DisallowIfNonnative( utility::vector1< bool > disallowed_aas, utility::vector1<core::Size> res_selection);
	virtual ~DisallowIfNonnative();
	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;
	//helper functions to define desired AAs
	void clear();
	//define as true which residues are NOT allowed.
	// The functions are not additive, will not remember previous disallowed
	void disallow_aas( utility::vector1< bool > const & cannonical_disallowed );
	void disallow_aas( std::string const & aa_string );
	// allow restriction of residues, additive
	void restrict_to_residue( core::Size const & resid);
	void restrict_to_residue( utility::vector1< core::Size > const & residues);
	virtual void parse_tag( TagCOP, DataMap & );
	virtual void parse_def( utility::lua::LuaObject const & def);


private:
	utility::vector1< bool > invert_vector( utility::vector1< bool > disallowed_aas);
	utility::vector1< core::Size > residue_selection_;
	utility::vector1< bool > disallowed_aas_;
	utility::vector1< bool > allowed_aas_;
};

/// @brief rotamer explosion for a residue
class RotamerExplosion : public TaskOperation
{
public:
	typedef TaskOperation parent;
public:
	RotamerExplosion();
	RotamerExplosion( core::Size const resid, ExtraRotSample const sample_level, core::Size const chi );
	virtual ~RotamerExplosion();
	virtual TaskOperationOP clone() const;
	virtual void apply( pose::Pose const &, PackerTask & ) const;
	virtual void parse_tag( TagCOP, DataMap & );
	void resid( core::Size const r );
	void chi( core::Size const c );
	void sample_level( ExtraRotSample const s );
private:
	core::Size resid_, chi_;
	ExtraRotSample sample_level_;
};

class InitializeFromCommandline : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	virtual ~InitializeFromCommandline();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

	virtual void parse_tag( TagCOP, DataMap & ); //parses nothing
	virtual void parse_def( utility::lua::LuaObject const & def);
};

class InitializeExtraRotsFromCommandline : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	virtual ~InitializeExtraRotsFromCommandline();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

};

class IncludeCurrent : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	virtual ~IncludeCurrent();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;
	virtual void parse_def( utility::lua::LuaObject const & def);

};

///@details control the extra chi rotamers for all residues
class ExtraRotamersGeneric : public TaskOperation {
public:
	typedef TaskOperation parent;

public:
  ExtraRotamersGeneric();
	virtual ~ExtraRotamersGeneric();

	virtual TaskOperationOP clone() const;

	virtual void parse_tag( TagCOP, DataMap & );

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

private:

	bool ex1_;
	bool ex2_;
	bool ex3_;
	bool ex4_;
	bool ex1aro_;
	bool ex2aro_;
	bool ex1aro_exposed_;
	bool ex2aro_exposed_;

	ExtraRotSample ex1_sample_level_;
	ExtraRotSample ex2_sample_level_;
	ExtraRotSample ex3_sample_level_;
	ExtraRotSample ex4_sample_level_;
	ExtraRotSample ex1aro_sample_level_;
	ExtraRotSample ex2aro_sample_level_;
	ExtraRotSample ex1aro_exposed_sample_level_;
	ExtraRotSample ex2aro_exposed_sample_level_;
	ExtraRotSample exdna_sample_level_;
	Size extrachi_cutoff_;
};

class ReadResfile : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	ReadResfile();
	ReadResfile( std::string const & );
	virtual ~ReadResfile();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

	void filename( std::string const & filename );

	///@brief queries options system for resfile name
	void default_filename();

	std::string const & filename() const;

	virtual void parse_tag( TagCOP, DataMap & );
	virtual void parse_def( utility::lua::LuaObject const & def);

private:
	std::string resfile_filename_;
};

/// @brief written by flo, feb 2011
/// class that can apply a resfile to a pose
/// that had its length changed at some point
/// in a protocol. A LengthEventCollector must
/// be set in the pose's observer cache for this
/// to work properly
class ReadResfileAndObeyLengthEvents : public ReadResfile
{

public:
	typedef ReadResfile parent;

	ReadResfileAndObeyLengthEvents();
	ReadResfileAndObeyLengthEvents( std::string const &);

	~ReadResfileAndObeyLengthEvents();

	virtual TaskOperationOP clone() const;

	void
	apply( pose::Pose const & pose, PackerTask & ptask ) const;

	virtual void parse_tag( TagCOP, DataMap & );

	std::list< ResfileCommandCOP > const &
	resfile_commands(
		Size const resfile_seqpos,
		ResfileContents const & contents,
		PackerTask const & ptask ) const;

private:
	bool apply_default_commands_to_inserts_;

};


class SetRotamerCouplings : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	SetRotamerCouplings();
	SetRotamerCouplings(SetRotamerCouplings const & );
	virtual ~SetRotamerCouplings();
	SetRotamerCouplings const & operator = ( SetRotamerCouplings const & );

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

	void
	set_couplings( rotamer_set::RotamerCouplingsOP couplings );


private:
	rotamer_set::RotamerCouplingsCOP rotamer_couplings_;
};


class SetRotamerLinks : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	SetRotamerLinks();
	SetRotamerLinks(SetRotamerLinks const & );
	virtual ~SetRotamerLinks();
	SetRotamerLinks const & operator = ( SetRotamerLinks const & );

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

	void
	set_links( rotamer_set::RotamerLinksOP links );


private:
	rotamer_set::RotamerLinksCOP rotamer_links_;
};


///@brief when a PackerTask is created by the Factory, the RotamerOperation will be given to it
class AppendRotamer : public TaskOperation
{
public:
	typedef TaskOperation parent;
public:
	AppendRotamer();
	virtual ~AppendRotamer();
	AppendRotamer( rotamer_set::RotamerOperationOP rotamer_operation );
	AppendRotamer( AppendRotamer const & );
	virtual TaskOperationOP clone() const;
	virtual void apply( pose::Pose const &, PackerTask & ) const;
	void set_rotamer_operation( rotamer_set::RotamerOperationOP rotamer_operation );


private:
	rotamer_set::RotamerOperationOP rotamer_operation_;
};


///@brief when a PackerTask is created by the Factory, the RotamerSetOperation will be given to it
///@author Barak Raveh (Dec 2008)
class AppendRotamerSet : public TaskOperation
{
public:
	typedef TaskOperation parent;
public:
	AppendRotamerSet();
	virtual ~AppendRotamerSet();
	AppendRotamerSet( rotamer_set::RotamerSetOperationOP rotamer_set_operation );
	AppendRotamerSet( AppendRotamerSet const & );
	virtual TaskOperationOP clone() const;
	virtual void apply( pose::Pose const &, PackerTask & ) const;
	void set_rotamer_set_operation( rotamer_set::RotamerSetOperationOP rotamer_operation );

private:
	rotamer_set::RotamerSetOperationOP rotamer_set_operation_;
};

///@brief Apply rotamerSetOperation to only the rotamerSet for the given residue
///@author Tim Jacobs (2011)
class AppendResidueRotamerSet : public TaskOperation
{
public:
	typedef TaskOperation parent;
public:
	AppendResidueRotamerSet();
	virtual ~AppendResidueRotamerSet();
	AppendResidueRotamerSet( core::Size resnum, rotamer_set::RotamerSetOperationOP rotamer_set_operation );
	AppendResidueRotamerSet( AppendResidueRotamerSet const & );
	virtual TaskOperationOP clone() const;
	virtual void apply( pose::Pose const &, PackerTask & ) const;
	void set_resnum( core::Size resnum );
	void set_rotamer_set_operation( rotamer_set::RotamerSetOperationOP rotamer_operation );

private:
	core::Size resnum_;
	rotamer_set::RotamerSetOperationOP rotamer_set_operation_;
};

class PreserveCBeta : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	virtual ~PreserveCBeta();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

};


/// @brief PreventRepacking allows you to prevent repacking (NATRO behavior) through the Factory.  Useful if you do not know the residue numbers when the resfile is created.  Note that this is unlike RestrictToRepacking; you have to specify which residues.  If PreventRepacking worked on the entire Task you'd have a do-nothing task.
class PreventRepacking : public TaskOperation
{
public:
	typedef TaskOperation parent;

public:
	virtual ~PreventRepacking();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( pose::Pose const &, PackerTask & ) const;

	void include_residue( core::Size resid );
	void clear();

	virtual void parse_tag( TagCOP, DataMap & );

private:
	utility::vector1< core::Size > residues_to_prevent_;
	std::string residue_selection_;
};


/// @brief RestrictYSDesign restricts positions to a binary Tyr/Ser alphabet
class RestrictYSDesign : public TaskOperation {

public:
	typedef TaskOperation parent;

public:
	RestrictYSDesign();
	RestrictYSDesign( RestrictYSDesign const & src );
	RestrictYSDesign( utility::vector1< core::Size > const & resids );
	virtual ~RestrictYSDesign();

	virtual TaskOperationOP clone() const;
	virtual void apply( pose::Pose const &, PackerTask & task) const;
	void include_resid( core::Size const resid );
	void include_gly( bool const gly );
	// void clear();  // Not defined, commenting out to make Python binding compile

private:
	utility::vector1< core::Size > YSresids_;
	bool gly_switch_;
};

/// @brief ExtraRotamer for a residue. You know, -ex1, -ex2, all that.
class ExtraRotamers : public TaskOperation
{
public:
	typedef TaskOperation parent;
public:
	ExtraRotamers();
	ExtraRotamers( core::Size const resid, core::Size const chi, core::Size level = 1 );
	virtual ~ExtraRotamers();
	virtual TaskOperationOP clone() const;
	virtual void apply( pose::Pose const &, PackerTask & ) const;
	virtual void parse_tag( TagCOP, DataMap & );
private:
	core::Size resid_, chi_, level_;
};

/// @brief ExtraChiCutoff (correponding to flag "-extrachi_cutoff <float>" )
class ExtraChiCutoff : public TaskOperation
{
public:
	typedef TaskOperation parent;
public:
	ExtraChiCutoff();
	ExtraChiCutoff( core::Size const resid, core::Size const extrachi_cutoff );
	virtual ~ExtraChiCutoff();
	virtual TaskOperationOP clone() const;
	virtual void apply( pose::Pose const &, PackerTask & ) const;
	virtual void parse_tag( TagCOP, DataMap & );
private:
	core::Size resid_, extrachi_cutoff_;
};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
