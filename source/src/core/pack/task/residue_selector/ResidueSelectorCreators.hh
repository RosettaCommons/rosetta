// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueSelectorCreators.hh
/// @brief  Class declarations for the ResidueSelectorCreators for a set of simple ResidueSelectors
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_task_residue_selector_ResidueSelectorCreators_HH
#define INCLUDED_core_pack_task_residue_selector_ResidueSelectorCreators_HH

// Package headers
#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreator.hh>

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class AndResidueSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class ChainSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class ClashBasedRepackShellSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class InterGroupInterfaceByVectorSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class NotResidueSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class ResidueIndexSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class ResidueNameSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class NeighborhoodResidueSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class NumNeighborsSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class OrResidueSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class JumpUpstreamSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class JumpDownstreamSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class SecondaryStructureSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

class TrueResidueSelectorCreator : public ResidueSelectorCreator {
public:
  virtual ResidueSelectorOP create_residue_selector() const;
  virtual std::string keyname() const;
};

} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#endif
