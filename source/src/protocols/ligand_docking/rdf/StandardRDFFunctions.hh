// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/rdf/StandardRDFFunctions.hh
///
/// @brief headers for standard RDF Functions used by RosettaHTS
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_ligand_docking_rdf_StandardRDFFunctions_hh
#define INCLUDED_protocols_ligand_docking_rdf_StandardRDFFunctions_hh

#include <protocols/ligand_docking/rdf/StandardRDFFunctions.fwd.hh>
#include <protocols/ligand_docking/rdf/RDFFunctionCreator.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/ligand_docking/rdf/RDFBase.hh>
#include <core/scoring/etable/coulomb/Coulomb.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/orbitals/OrbitalsScore.fwd.hh>

namespace protocols {
namespace ligand_docking {
namespace rdf {

/// @brief Creator to geneate a new RDFEtableFunction
class RDFEtableCreator : public RDFFunctionCreator
{
public:
	
	/// @brief return a pointer to a newly created RDFEtableFunction
	virtual RDFBaseOP create_rdf_function() const;
	/// @brief return the name of the RDFEtableFunction
	virtual std::string type_name() const;
};
	
/// @brief RDFEtableFunction computes fa_sol,fa_rep,fa_atr for a pair of atoms
class RDFEtableFunction : public RDFBase
{
public:
	RDFEtableFunction();
	
	virtual ~RDFEtableFunction();
	
	/// @brief parse tags for RDFEtableFunction tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map);
	
	/// @brief compute atr, rep and solvation energy for atom pair
	virtual RDFResultList operator()(AtomPairData const & atom_data );

	
private:
	core::scoring::etable::AnalyticEtableEvaluatorOP etable_evaluator_;
	
};
	
/// @brief Creator to geneate a new RDFElecFunction
class RDFElecCreator : public RDFFunctionCreator
{
public:
	
	/// @brief return a pointer to a newly created RDFElecFunction
	virtual RDFBaseOP create_rdf_function() const;
	/// @brief return the name of the RDFEtableFunction
	virtual std::string type_name() const;
};

/// @brief RDFElecFunction computes fa_elec for a pair of atoms
class RDFElecFunction : public RDFBase
{
public:
	RDFElecFunction();
	
	virtual ~RDFElecFunction();
	
	/// @brief parse tags for RDFElecFunction tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map);
	
	/// @brief compute fa_elec for atom pair
	virtual RDFResultList operator()(AtomPairData const & atom_data );
	
	
private:
	core::scoring::etable::coulomb::CoulombOP coloumb_;
};
	
	
/// @brief Creator to geneate a new RDFChargeFunction
class RDFChargeCreator : public RDFFunctionCreator
{
public:
	
	/// @brief return a pointer to a newly created RDFElecFunction
	virtual RDFBaseOP create_rdf_function() const;
	/// @brief return the name of the RDFEtableFunction
	virtual std::string type_name() const;
};

/// @brief RDFChargeFunction computes fa_elec for a pair of atoms
class RDFChargeFunction : public RDFBase
{
public:
	RDFChargeFunction();
	
	virtual ~RDFChargeFunction();
	
	/// @brief parse tags for RDFElecFunction tag
	virtual void parse_my_tag(
							  utility::tag::TagCOP tag,
							  basic::datacache::DataMap & );
	
	/// @brief compute fa_elec for atom pair
	virtual RDFResultList operator()(AtomPairData const & atom_data );
	
	
private:
	FunctionSign function_sign_;
	std::string function_name_;

};

/// @brief Creator to geneate a new RDFHbondFunction
class RDFHbondCreator : public RDFFunctionCreator
{
public:
	
	/// @brief return a pointer to a newly created RDFHbondFunction
	virtual RDFBaseOP create_rdf_function() const;
	/// @brief return the name of the RDFHbondFunction
	virtual std::string type_name() const;
};
	
/// @brief RDFHbondFunction computes h-bonding energy for a pair of atoms.
class RDFHbondFunction : public RDFBase
{
public:
	RDFHbondFunction();
	
	virtual ~RDFHbondFunction();
	
	/// @brief parse tags for RDFHbondFunction tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map);
	
	/// @brief compute hbond energy for atom pair
	virtual RDFResultList operator()(AtomPairData const & atom_data );
	
	/// @brief setup hbond database for each pose
	virtual void preamble(core::pose::Pose & pose);
	
private:
	core::scoring::hbonds::HBondSetOP hbond_set_;
	FunctionSign function_sign_;
	std::string function_name_;
	
};
	
/// @brief Creator to geneate a new RDFBinaryHbondFunction
class RDFBinaryHbondCreator : public RDFFunctionCreator
{
public:
	
	/// @brief return a pointer to a newly created RDFBinaryHbondFunction
	virtual RDFBaseOP create_rdf_function() const;
	/// @brief return the name of the RDFBinaryHbondFunction
	virtual std::string type_name() const;
};
	
/// @brief RDFBinaryHbondFunction returns 1.0 if a pair of atoms are donor and acceptor, 0.0 otherwise.
class RDFBinaryHbondFunction : public RDFBase
{
public:
	RDFBinaryHbondFunction();
	
	virtual ~RDFBinaryHbondFunction();
	
	/// @brief parse tags for RDFBinaryHbondFunction tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map);
	
	/// @brief compute binary hbond status for atom pair
	virtual RDFResultList operator()(AtomPairData const & atom_data );
	
	
private:
	FunctionSign function_sign_;
	std::string function_name_;
};

/// @brief Creator to geneate a new RDFOrbitalFunction
class RDFOrbitalFunctionCreator : public RDFFunctionCreator
{
public:
	
	/// @brief return a pointer to a newly created RDFOrbitalFunction
	virtual RDFBaseOP create_rdf_function() const;
	/// @brief return the name of the RDFOrbitalFunction
	virtual std::string type_name() const;
};
	
/// @brief RDFOrbitalFunction computes the orbital score energies of a pair of atoms
class RDFOrbitalFunction : public RDFBase
{
public:
	RDFOrbitalFunction();
	
	virtual ~RDFOrbitalFunction();
	
	/// @brief parse tags for RDFOrbitalFunction tag
	virtual void parse_my_tag(
							  utility::tag::TagCOP tag,
							  basic::datacache::DataMap & data_map);
	
	/// @brief compute binary hbond status for atom pair
	virtual RDFResultList operator()(AtomPairData const & atom_data );
	
	/// @brief update residue neighbors and cache current pose
	virtual void preamble(core::pose::Pose & pose);
	
	
private:
	core::pose::PoseAP pose_;
	core::scoring::orbitals::OrbitalsScoreOP orbital_score_;
};
	
/// @brief Creator to geneate a new RDFBinaryOrbitalFunction
class RDFBinaryOrbitalFunctionCreator : public RDFFunctionCreator
{
public:
	
	/// @brief return a pointer to a newly created RDFBinaryOrbitalFunction
	virtual RDFBaseOP create_rdf_function() const;
	/// @brief return the name of the RDFBinaryOrbitalFunction
	virtual std::string type_name() const;
};
	
/// @brief RDFBinaryOrbitalFunction returns 1 for various orbital pairs and 0 otherwise
class RDFBinaryOrbitalFunction : public RDFBase
{
public:
	RDFBinaryOrbitalFunction();
	
	virtual ~RDFBinaryOrbitalFunction();
	
	/// @brief parse tags for RDFOrbitalFunction tag
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map);
	
	/// @brief compute binary hbond status for atom pair
	virtual RDFResultList operator()(AtomPairData const & atom_data );
	
	/// @brief update residue neighbors and cache current pose
	virtual void preamble(core::pose::Pose & pose);
	
	
private:
	core::pose::PoseAP pose_;
};
	
}
}
}



#endif
