// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResFilters.hh
/// @brief  core-level (very general) classes that take a pose and a residue index, and return true or false
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResFilters_hh
#define INCLUDED_core_pack_task_operation_ResFilters_hh

// Unit Headers
#include <core/pack/task/operation/ResFilters.fwd.hh>

#include <core/pack/task/operation/ResFilter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

#include <string>
#include <set>


namespace core {
namespace pack {
namespace task {
namespace operation {

class ResFilterComposition : public ResFilter
{
	public:
		typedef ResFilter parent;

	public:
		ResFilterComposition();
		ResFilterComposition(utility::vector1<ResFilterCOP> const &);

		virtual void parse_tag( TagCOP );

	protected:
		utility::vector1<ResFilterCOP> sub_filters_;
	  void parse_sub_filters_tag( TagCOP );
};

class AnyResFilter : public ResFilterComposition
{
public:
	typedef ResFilterComposition parent;

	AnyResFilter();
	AnyResFilter(utility::vector1<ResFilterCOP> const &);
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
};

class AllResFilter : public ResFilterComposition
{
public:
	typedef ResFilterComposition parent;

	AllResFilter();
	AllResFilter(utility::vector1<ResFilterCOP> const &);
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
};

class NoResFilter : public ResFilterComposition
{
public:
	typedef ResFilterComposition parent;

	NoResFilter();
	NoResFilter(utility::vector1<ResFilterCOP> const &);
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
};

// @brief Convenience filter to filter by residue type (polar, apolar, aromatic, charged)
class ResidueTypeFilter : public ResFilter
{
public:
	typedef ResFilter parent;

	ResidueTypeFilter();
	ResidueTypeFilter(bool polar, bool apolar, bool aromatic, bool charged);

	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
	virtual void parse_tag( TagCOP );

private:
	bool polar_, apolar_, aromatic_, charged_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// NOTE: in most cases, each 'Inst' class inherits largely from its corresponding 'Is' class

class ResidueHasProperty : public ResFilter {
public:
	typedef ResFilter parent;
public:
	ResidueHasProperty();
	ResidueHasProperty( std::string const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
	virtual void parse_tag( TagCOP );
	virtual std::string const & property() const { return property_; }
private:
	std::string property_;
};

class ResidueLacksProperty : public ResidueHasProperty {
public:
	typedef ResidueHasProperty parent;
public:
	ResidueLacksProperty();
	ResidueLacksProperty( std::string const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
};

class ResiduePDB_InfoHasLabel : public ResFilter {
public:
  typedef ResFilter parent;
public:
  ResiduePDB_InfoHasLabel();
  ResiduePDB_InfoHasLabel( std::string const & );
  virtual bool operator() ( Pose const &, Size ) const;
  virtual ResFilterOP clone() const;
  virtual void parse_tag( TagCOP );
  virtual std::string const & property() const { return property_; }
private:
  std::string property_;
};

class ResiduePDB_InfoLacksLabel : public ResiduePDB_InfoHasLabel {
public:
  typedef ResiduePDB_InfoHasLabel parent;
public:
  ResiduePDB_InfoLacksLabel();
  ResiduePDB_InfoLacksLabel( std::string const & );
  virtual bool operator() ( Pose const &, Size ) const;
  virtual ResFilterOP clone() const;
};

class ResidueName3Is : public ResFilter {
public:
	typedef ResFilter parent;
public:
	ResidueName3Is();
	ResidueName3Is( std::string const & );
	ResidueName3Is( std::set<std::string> const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
	virtual void parse_tag( TagCOP );
private:
	std::set<std::string> name3_set;
};

class ResidueName3Isnt : public ResidueName3Is {
public:
	typedef ResidueName3Is parent;
public:
	ResidueName3Isnt();
	ResidueName3Isnt( std::string const & );
	ResidueName3Isnt( std::set<std::string> const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
};

class ResidueIndexIs : public ResFilter {
public:
	typedef ResFilter parent;

public:
	ResidueIndexIs();
	ResidueIndexIs( Size );
	ResidueIndexIs( utility::vector1< Size > const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
	virtual void parse_tag( TagCOP );
	virtual utility::vector1< Size > const & indices() const;

private:
	utility::vector1< Size > indices_;
};

class ResidueIndexIsnt : public ResidueIndexIs {
public:
	typedef ResidueIndexIs parent;

public:
	ResidueIndexIsnt();
	ResidueIndexIsnt( Size );
	ResidueIndexIsnt( utility::vector1< Size > const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
};

class ResiduePDBIndexIs : public ResFilter {
public:
	typedef ResFilter parent;

	struct ChainPos { // for (optional) pdb indexing
		char chain_; int pos_;
		ChainPos( char chain, int pos ) : chain_(chain), pos_(pos) {}
		bool operator == ( ChainPos const & other ) const
			{ return( chain_ == other.chain_ && pos_ == other.pos_ ); }
	};

public:
	ResiduePDBIndexIs();
	ResiduePDBIndexIs( char, int );
	ResiduePDBIndexIs( utility::vector1< ChainPos > const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
	virtual void parse_tag( TagCOP );
	virtual utility::vector1< ChainPos > const & indices() const;

private:
	utility::vector1< ChainPos > indices_;
};

class ResiduePDBIndexIsnt : public ResiduePDBIndexIs {
public:
	typedef ResiduePDBIndexIs parent;

public:
	ResiduePDBIndexIsnt();
	ResiduePDBIndexIsnt( char, int );
	ResiduePDBIndexIsnt( utility::vector1< ChainPos > const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
};

class ChainIs : public ResFilter {
public:
	typedef ResFilter parent;
public:
	ChainIs();
	ChainIs( char const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
	virtual void parse_tag( TagCOP );
	virtual char const & chain() const { return chain_; }
private:
	char chain_;
};

class ChainIsnt : public ChainIs {
public:
	typedef ChainIs parent;
public:
	ChainIsnt();
	ChainIsnt( char const & );
	virtual bool operator() ( Pose const &, Size ) const;
	virtual ResFilterOP clone() const;
};

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
