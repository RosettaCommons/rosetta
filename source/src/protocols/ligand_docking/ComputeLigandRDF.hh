// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/ComputeLigandRDF.hh
///
/// @brief header file for ComputeLigandRDF mover
/// @author Sam DeLuca (sam@decarboxy.com)


#ifndef INCLUDED_protocols_ligand_docking_ComputeLigandRDF_hh
#define INCLUDED_protocols_ligand_docking_ComputeLigandRDF_hh

#include <protocols/ligand_docking/ComputeLigandRDF.fwd.hh>
#include <protocols/ligand_docking/rdf/RDFBase.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/chemical/AtomType.hh>

#include <utility/tag/Tag.fwd.hh>

#include <map>
#include <string>

namespace protocols {
namespace ligand_docking {

	
class ComputeLigandRDF : public protocols::moves::Mover{

public:
	ComputeLigandRDF();
	virtual ~ComputeLigandRDF();
	ComputeLigandRDF(ComputeLigandRDF const & that);
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag
	(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
	);

private:
	
	///@brief compute the RDF between pairs of protein and ligand atoms
	std::map<std::string, utility::vector1<core::Real> > ligand_protein_rdf(core::pose::Pose & pose);
	///@brief compute the RDF between pairs of protein protein atoms
	std::map<std::string, utility::vector1<core::Real> > protein_protein_rdf(core::pose::Pose & pose);
	///@brief compute the RDF given a set of atom-atom pairs.  does most of the work.
	std::map<std::string, utility::vector1<core::Real> >  compute_rdf(
		core::pose::Pose & pose,
		utility::vector1<std::pair<core::id::AtomID, core::id::AtomID> > const & atom_pairs );


private:
	std::string mode_;
	std::string ligand_chain_;
	core::Size bin_count_;
	core::Real smoothing_factor_;
	core::Real range_squared_;
	utility::vector1<rdf::RDFBaseOP> rdf_functions_;

}; // class ComputeLigandRDF


} // namespace ligand_docking
} // namespace protocols

#endif
