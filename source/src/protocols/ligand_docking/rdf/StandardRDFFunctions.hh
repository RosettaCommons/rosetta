// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/rdf/StandardRDFFunctions.hh
///
/// @brief headers for standard RDF Functions used by RosettaHTS
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_ligand_docking_rdf_StandardRDFFunctions_hh
#define INCLUDED_protocols_ligand_docking_rdf_StandardRDFFunctions_hh

#include <protocols/ligand_docking/rdf/StandardRDFFunctions.fwd.hh>
#include <protocols/ligand_docking/rdf/RDFFunctionCreator.hh>
#include <protocols/ligand_docking/rdf/RDFBase.hh>
#include <core/scoring/etable/coulomb/Coulomb.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/orbitals/OrbitalsScore.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace ligand_docking {
namespace rdf {

/// @brief Creator to geneate a new RDFEtableFunction
class RDFEtableCreator : public RDFFunctionCreator
{
public:

	/// @brief return a pointer to a newly created RDFEtableFunction
	RDFBaseOP create_rdf_function() const override;
	/// @brief return the name of the RDFEtableFunction
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

/// @brief RDFEtableFunction computes fa_sol,fa_rep,fa_atr for a pair of atoms
class RDFEtableFunction : public RDFBase
{
public:
	RDFEtableFunction();

	~RDFEtableFunction() override;

	/// @brief parse tags for RDFEtableFunction tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map) override;

	/// @brief compute atr, rep and solvation energy for atom pair
	RDFResultList operator()(AtomPairData const & atom_data ) override;
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::scoring::etable::AnalyticEtableEvaluatorOP etable_evaluator_;

};

/// @brief Creator to geneate a new RDFElecFunction
class RDFElecCreator : public RDFFunctionCreator
{
public:

	/// @brief return a pointer to a newly created RDFElecFunction
	RDFBaseOP create_rdf_function() const override;
	/// @brief return the name of the RDFEtableFunction
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

/// @brief RDFElecFunction computes fa_elec for a pair of atoms
class RDFElecFunction : public RDFBase
{
public:
	RDFElecFunction();

	~RDFElecFunction() override;

	/// @brief parse tags for RDFElecFunction tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map) override;

	/// @brief compute fa_elec for atom pair
	RDFResultList operator()(AtomPairData const & atom_data ) override;
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::scoring::etable::coulomb::CoulombOP coloumb_;
};


/// @brief Creator to geneate a new RDFChargeFunction
class RDFChargeCreator : public RDFFunctionCreator
{
public:

	/// @brief return a pointer to a newly created RDFElecFunction
	RDFBaseOP create_rdf_function() const override;
	/// @brief return the name of the RDFEtableFunction
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

/// @brief RDFChargeFunction computes fa_elec for a pair of atoms
class RDFChargeFunction : public RDFBase
{
public:
	RDFChargeFunction();

	~RDFChargeFunction() override;

	/// @brief parse tags for RDFElecFunction tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;

	/// @brief compute fa_elec for atom pair
	RDFResultList operator()(AtomPairData const & atom_data ) override;
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	FunctionSign function_sign_;
	std::string function_name_;

};

/// @brief Creator to geneate a new RDFHbondFunction
class RDFHbondCreator : public RDFFunctionCreator
{
public:

	/// @brief return a pointer to a newly created RDFHbondFunction
	RDFBaseOP create_rdf_function() const override;
	/// @brief return the name of the RDFHbondFunction
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

/// @brief RDFHbondFunction computes h-bonding energy for a pair of atoms.
class RDFHbondFunction : public RDFBase
{
public:
	RDFHbondFunction();

	~RDFHbondFunction() override;

	/// @brief parse tags for RDFHbondFunction tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map) override;

	/// @brief compute hbond energy for atom pair
	RDFResultList operator()(AtomPairData const & atom_data ) override;

	/// @brief setup hbond database for each pose
	void preamble(core::pose::Pose & pose) override;
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

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
	RDFBaseOP create_rdf_function() const override;
	/// @brief return the name of the RDFBinaryHbondFunction
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

/// @brief RDFBinaryHbondFunction returns 1.0 if a pair of atoms are donor and acceptor, 0.0 otherwise.
class RDFBinaryHbondFunction : public RDFBase
{
public:
	RDFBinaryHbondFunction();

	~RDFBinaryHbondFunction() override;

	/// @brief parse tags for RDFBinaryHbondFunction tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map) override;

	/// @brief compute binary hbond status for atom pair
	RDFResultList operator()(AtomPairData const & atom_data ) override;
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	FunctionSign function_sign_;
	std::string function_name_;
};

/// @brief Creator to geneate a new RDFOrbitalFunction
class RDFOrbitalFunctionCreator : public RDFFunctionCreator
{
public:

	/// @brief return a pointer to a newly created RDFOrbitalFunction
	RDFBaseOP create_rdf_function() const override;
	/// @brief return the name of the RDFOrbitalFunction
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd) const override;
};

/// @brief RDFOrbitalFunction computes the orbital score energies of a pair of atoms
class RDFOrbitalFunction : public RDFBase
{
public:
	RDFOrbitalFunction();

	~RDFOrbitalFunction() override;

	/// @brief parse tags for RDFOrbitalFunction tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map) override;

	/// @brief compute binary hbond status for atom pair
	RDFResultList operator()(AtomPairData const & atom_data ) override;

	/// @brief update residue neighbors and cache current pose
	void preamble(core::pose::Pose & pose) override;
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::pose::PoseAP pose_;
	core::scoring::orbitals::OrbitalsScoreOP orbital_score_;
};

/// @brief Creator to geneate a new RDFBinaryOrbitalFunction
class RDFBinaryOrbitalFunctionCreator : public RDFFunctionCreator
{
public:

	/// @brief return a pointer to a newly created RDFBinaryOrbitalFunction
	RDFBaseOP create_rdf_function() const override;
	/// @brief return the name of the RDFBinaryOrbitalFunction
	std::string type_name() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd) const override;
};

/// @brief RDFBinaryOrbitalFunction returns 1 for various orbital pairs and 0 otherwise
class RDFBinaryOrbitalFunction : public RDFBase
{
public:
	RDFBinaryOrbitalFunction();

	~RDFBinaryOrbitalFunction() override;

	/// @brief parse tags for RDFOrbitalFunction tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map) override;

	/// @brief compute binary hbond status for atom pair
	RDFResultList operator()(AtomPairData const & atom_data ) override;

	/// @brief update residue neighbors and cache current pose
	void preamble(core::pose::Pose & pose) override;
	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::pose::PoseAP pose_;
};

}
}
}


#endif
