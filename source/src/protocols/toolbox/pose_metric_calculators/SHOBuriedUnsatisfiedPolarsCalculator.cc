// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Calculator for buried unsatisfied polar atoms under the SHO solvation model
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#include <protocols/toolbox/pose_metric_calculators/SHOBuriedUnsatisfiedPolarsCalculator.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AtomType.hh>
#include <basic/MetricValue.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <map>
#include <fstream>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

using core::Size;
using core::Real;
using core::id::AtomID;
using core::scoring::geometric_solvation::create_ExactSHOEnergy_from_cmdline;


static basic::Tracer TR("protocols.toolbox.pose_metric_calculators.SHOBuriedUnsatisfiedPolarsCalculator");


/// @brief constructs the calculator for a target atom in residues of a target
///  amino acid type.
///
/// @params[in] sho_cutoff maximum SHO energy value for a polar group to be
///  considered as exposed (i.e., not buried). It must be expressed in original SHO units.
/// @params[in] tgt_amino target amino acid type
/// @params[in] tgt_atom target atom name
/// @params[in] sfxn score function used to (previously) score the pose
/// @details The SHO method is initialized as in
///  ExactOccludedHbondSolEnergyCreator::create_energy_method() for consistency
///  with any scoring functions using SHO.
///
SHOBuriedUnsatisfiedPolarsCalculator::SHOBuriedUnsatisfiedPolarsCalculator(
	core::Real sho_cutoff, std::string tgt_amino, std::string tgt_atom,
	core::scoring::ScoreFunctionCOP sfxn) :

	sho_meth_(create_ExactSHOEnergy_from_cmdline(sfxn->energy_method_options())),
	sho_cutoff_( sho_cutoff ),
	tgt_amino_(tgt_amino),
	tgt_atom_(tgt_atom),
	sfxn_(sfxn) {

	search_typ_ = (tgt_amino_ == "any") ? ANY_AA__TGT_ATOM : TGT_AA__TGT_ATOM;
}


/// @brief constructs the calculator for a target residue set.
///
/// @params[in] sho_cutoff maximum SHO energy value for a polar group to be
///  considered as exposed (i.e., not buried). It must be expressed in original SHO units.
/// @params[in] tgt_res_idxs set of target residue indexes. An empty set means
///  "all residues".
/// @params[in] sfxn score function used to (previously) score the pose
///
/// @details The SHO method is initialized as in
///   ExactOccludedHbondSolEnergyCreator::create_energy_method() for consistency
///   with any scoring functions using SHO.
///
SHOBuriedUnsatisfiedPolarsCalculator::SHOBuriedUnsatisfiedPolarsCalculator(
	Real sho_cutoff, utility::vector1<Size> const& tgt_res_idxs,
	core::scoring::ScoreFunctionCOP sfxn) :

	sho_meth_(create_ExactSHOEnergy_from_cmdline(sfxn->energy_method_options())),
	sho_cutoff_( sho_cutoff ) ,
	tgt_res_idxs_(tgt_res_idxs),
	tgt_amino_("none"),
	tgt_atom_("none"),
	sfxn_(sfxn) {
	search_typ_ = tgt_res_idxs_.size() ? TGT_RES : ALL_ATOMS;
}


/// @brief clones this calculator
///
/// @details members sho_meth_ and sfxn_ of the clone will point, respectively,
///   to the same objects as members sho_meth_ and sfxn_ of this calculator.
///
core::pose::metrics::PoseMetricCalculatorOP
SHOBuriedUnsatisfiedPolarsCalculator::clone() const {

	return core::pose::metrics::PoseMetricCalculatorOP(
		new SHOBuriedUnsatisfiedPolarsCalculator(*this));
}


/// @brief outputs a quantity's value that is cached in the calculator
///
/// @param[in] key identifier of the quantity
/// @param[out] valptr storage for the value
///
void SHOBuriedUnsatisfiedPolarsCalculator::lookup(std::string const &key,
	basic::MetricValueBase* valptr) const {

	if ( key == "num_buns" ) {
		basic::check_cast(valptr, &num_burunsat_atoms_,
			"the number of SHO buried unsatisfied atoms must be of type Size");
		(static_cast<basic::MetricValue<Size> *>(valptr))->
			set(num_burunsat_atoms_);
	} else {
		basic::Error() <<
			"SHOBuriedUnsatisfiedPolarsCalculator cannot compute metric" << key <<
			std::endl;
		utility_exit();
	}
}


/// @brief returns the string representation of a quantity's value that is
///  cached in the calculator
///
/// @param[in] key identifier of the quantity
///
std::string SHOBuriedUnsatisfiedPolarsCalculator::print(std::string const& key) const {

	if ( key == "num_buns" ) {
		return utility::to_string(num_burunsat_atoms_);
	} else {
		basic::Error() << "SHOBuriedUnsatisfiedPolarsCalculator cannot print metric"
			<< key << std::endl;
		utility_exit();
	}
}


/// @brief finds buried unsatisfied atoms by building all needed data structures
///  from scratch.
///
/// @param[in] ps pose for which the data are to be computed
///
void SHOBuriedUnsatisfiedPolarsCalculator::recompute(
	core::pose::Pose const& ps) {

	// compute SHO energy for each polar atom
	core::Size N = ps.size();
	for ( Size i=1; i<=N; ++i ) {

		core::conformation::Residue const& rsd = ps.residue(i);

		// donors
		for ( core::chemical::AtomIndices::const_iterator
				hnum = rsd.Hpos_polar().begin(),
				hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {

			Size const h( *hnum );
			sho_energies_[AtomID(h, i)] =
				sho_meth_->compute_sho_donor_atom_energy(rsd, i, h, ps);
		}

		// acceptors
		for ( core::chemical::AtomIndices::const_iterator
				anum = rsd.accpt_pos().begin(),
				anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {

			Size const acc( *anum );
			sho_energies_[AtomID(acc, i)] =
				sho_meth_->compute_sho_acceptor_atom_energy(rsd, i, acc, ps);
		}
	}

	partition(ps);

	num_burunsat_atoms_ = burunsat_atoms_.size();

	num_other_atoms_ = other_atoms_.size();
}


/// @brief prints on screen the SHO energies of all atoms in a pose
///
/// @param[in]: ps the pose
///
/// @details The output can be seen as a sequence of N blocks, where N is the
///  number of residues in the pose. Block i contains the energies of residue
///  i (i=1,...N). Within block i, line j contains the energy of atom j
///  (j=1,...,M(i), where M(i) is the number of atoms in residue i).
///  Each energy is preceded by the atom's identifier.
///
void SHOBuriedUnsatisfiedPolarsCalculator::print_sho_energies(
	core::pose::Pose const& ps) {

	Size const N = ps.size();
	for ( Size i=1; i<=N; ++i ) {
		core::conformation::Residue const& rsd = ps.residue(i);
		Size const M = rsd.natoms();
		for ( Size j=1; j<=M; ++j ) {
			AtomID aid(j,i);
			print_atom_info(aid, ps);
			std::cout << ": " << sho_energies_[AtomID(j, i)] << std::endl;
		}
	}
}


/// @brief partitions a target residue's polar atoms into "buried unsatisfied"
///  and "other"
///
/// @params[in] ps the pose at issue
/// @params[in] tgtidx pose index of the target residue
///
/// @details At the end of this function, burunsat_atoms_[B+i] will contain the
///  ith polar atom in the target residue for which functions is_buried() and
///  is_unsat() are both true (i=1,...N, where N is the number of such
///  atoms; B is the size of burunsat_atoms_ before this function is called).
///
/// @details At the end of this function, other_atoms_[P+i] will contain the
///  ith polar atom in the target residue for which either function is_buried()
///  or function is_unsat() are false (i=1,...M, where M is the number of such
///  atoms; P is the size of other_atoms_ before this function is called).
///
/// @details As to the ordering of atoms, donors are assumed to precede
///  acceptors, donors follow the order specified by Hpos_polar(), and
///  acceptors follow the order specified by accpt_pos().
///
void SHOBuriedUnsatisfiedPolarsCalculator::residue_partition(
	core::pose::Pose const& ps, Size tgtidx) {

	core::conformation::Residue const& rsd = ps.residue(tgtidx);

	// donors
	for ( core::chemical::AtomIndices::const_iterator
			hnum = rsd.Hpos_polar().begin(),
			hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {

		Size const h( *hnum );
		assign_atom(AtomID(h, tgtidx));
	}

	// acceptors
	for ( core::chemical::AtomIndices::const_iterator
			anum = rsd.accpt_pos().begin(),
			anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {

		Size const acc( *anum );
		assign_atom(AtomID(acc, tgtidx));
	}
}


/// @brief partitions target atoms into "buried unsatisfied" and "other".
///
/// @params[in] ps the pose at issue
///
/// @details At the end of this function burunsat_atoms_[i] will contain
/// the identifier of the ith target atom for which functions is_buried()
/// and is_unsat() are both true (i=1,...N, where N is the number of such
/// atoms).
///
/// @details At the end of this function other_atoms_[i] will contain
/// the identifier of the ith target atom for which either function
/// is_buried() or function is_unsat() are false (i=1,...,M, where M is
/// the number of such atoms)
///
/// @details In the above description, target atoms are assumed to be
///  ordered by pose residue number under options ALL_ATOMS, TGT_AA__TGT_ATOM,
///  and ANY_AA__TGT_ATOM, and by index in tgt_res_idxs_ under option
///  TGT_RES.
///
void SHOBuriedUnsatisfiedPolarsCalculator::partition(
	core::pose::Pose const& ps) {

	Size const N = ps.size();

	switch(search_typ_) {

	// select all polar atoms in the pose
	case ALL_ATOMS : {
		for ( Size i=1; i<=N; ++i ) {
			residue_partition(ps, i);
		}
		break;
	}

		// select all polar atoms in the target residue set
	case TGT_RES : {
		const Size T = tgt_res_idxs_.size();
		for ( Size i=1; i<=T; ++i ) {
			residue_partition(ps, tgt_res_idxs_[i]);
		}
		break;
	}

		// select the target atom from every residue of the target amino acid type
	case TGT_AA__TGT_ATOM : {

		using namespace core::chemical;
		AA TGTAA = aa_from_oneletter_code(tgt_amino_[0]);
		for ( Size i=1; i<=N; ++i ) {
			core::conformation::Residue const& rsd = ps.residue(i);
			if ( rsd.aa() == TGTAA ) {
				assign_atom(AtomID(rsd.atom_index(tgt_atom_), i));
			}
		}
		break;
	}

		// select the target atom from every residue that has it
	case ANY_AA__TGT_ATOM : {

		for ( Size i=1; i<=N; ++i ) {
			core::conformation::Residue const& rsd = ps.residue(i);
			if ( rsd.has(tgt_atom_) ) {
				assign_atom(AtomID(rsd.atom_index(tgt_atom_), i));
			}
		}
		break;
	}
	}
}


/// @brief assigns an atom to either the "buried unsatisfied" class or the
///  "other" class
///
/// @param[in] aid the atom
///
void SHOBuriedUnsatisfiedPolarsCalculator::assign_atom(
	AtomID aid) {

	if ( is_buried(aid) && is_unsat(aid) ) {
		burunsat_atoms_.push_back(aid);
	} else {
		other_atoms_.push_back(aid);
	}
}


/// @brief is atom aid buried?
///
/// @param[in] aid
///
bool SHOBuriedUnsatisfiedPolarsCalculator::is_buried(AtomID aid) {

	return sho_energies_[aid] > sho_cutoff_;
}


/// @brief is atom aid unsatisfied?
///
/// @param[in] aid
///
bool SHOBuriedUnsatisfiedPolarsCalculator::is_unsat(AtomID aid) const {

	return hbond_set_->nhbonds(aid) == 0;
}


/// @brief prints on screen info about the SHO energy and the hydrogen bonds of
///  each atom in a subset (i.e., either the "buried unsatisfied" atoms or the
///  "other" atoms).
///
/// @param[in] subset the subset of atoms
/// @param[in] ps the pose containing the subset of atoms
///
/// @details The ith output line contains the identifier of the ith atom in the
///  subset followed by its SHO energy and hydrogen bond info (i=1,...,N,
///  where N is the number of atoms in the subset).
///
void SHOBuriedUnsatisfiedPolarsCalculator::print_atom_subset(
	utility::vector1<AtomID> const& subset,
	core::pose::Pose const& ps) const {

	Size N = subset.size();
	for ( Size i=1; i<=N; ++i ) {
		AtomID aid = subset[i];
		print_atom_info(aid, ps);

		// primary info
		Real shoe = sho_energies_.find(aid)->second;
		std::cout << " " <<
			"Esho:" << shoe << " " <<
			"#hb:" << hbond_set_->nhbonds(aid);

		// secondar,y energy info
		std::cout << "  (";

		Real hbe = hbond_energy(aid, ps);

		Real wtote = sfxn_->get_weight(core::scoring::occ_sol_exact) *
			// core::scoring::geometric_solvation::LK_MATCHING_WEIGHT_EXACT * temporary fix to allow suppression of LK_MATCHING_WEIGHT_EXACT in header file
			shoe + hbe;

		std::cout <<
			"wEhb:" << hbe << " " <<
			"wEsho+wEhb: " << wtote << ')' << std::endl;
	}
}


/// @brief returns the total H-bond energy of an atom.
///
/// @param[in] aid AtomID of the atom
/// @param[in] ps pose to which the atom belongs
///
/// @details the total H-bond energy of the atom is defined as half the sum
///  of the weighted energies of the H-bonds that it forms.
///
Real SHOBuriedUnsatisfiedPolarsCalculator::hbond_energy(AtomID aid,
	core::pose::Pose const& ps) const {

	Real tot=0;
	utility::vector1<core::scoring::hbonds::HBondCOP> hbonds =
		hbond_set_->atom_hbonds(aid);

	core::Size N = hbonds.size();
	for ( Size i=1; i<=N; ++i ) {

		// get weight for current H-bond's energy
		core::scoring::hbonds::HBondCOP const& bond = hbonds[i];

		Size dhatm = bond->don_hatm();
		Size dresi = bond->don_res();
		core::conformation::Residue const& dres = ps.residue(dresi);
		Size datm = dres.atom_base(dhatm);

		Size aatm = bond->acc_atm();
		Size aresi = bond->acc_res();
		core::conformation::Residue const& ares = ps.residue(aresi);

		core::scoring::hbonds::HBEvalTuple hbe_type(datm, dres, aatm, ares);
		core::scoring::hbonds::HBondWeightType const hbwt =
			core::scoring::hbonds::get_hbond_weight_type(hbe_type.eval_type());

		Real hbw = 0;
		switch(hbwt){

		case core::scoring::hbonds::hbw_NONE:
		case core::scoring::hbonds::hbw_SR_BB :
			hbw = sfxn_->get_weight(core::scoring::hbond_sr_bb);
			break;

		case core::scoring::hbonds::hbw_LR_BB :
			hbw = sfxn_->get_weight(core::scoring::hbond_lr_bb);
			break;

		case core::scoring::hbonds::hbw_SR_BB_SC:
		case core::scoring::hbonds::hbw_LR_BB_SC :
			hbw = sfxn_->get_weight(core::scoring::hbond_bb_sc);
			break;

		case core::scoring::hbonds::hbw_SC :
			hbw = sfxn_->get_weight(core::scoring::hbond_sc);
			break;

		default :
			break;
		}

		// add weighted energy
		tot += 0.5*hbw*bond->energy()*bond->weight();
	}

	return tot;
}


/// @brief prints an atom's identifier on screen
///
/// @param[in] aid AtomID of the atom
/// @param[in] ps pose containing the atom
///
void SHOBuriedUnsatisfiedPolarsCalculator::print_atom_info(
	AtomID aid, core::pose::Pose const& ps) const {

	std::cout <<
		ps.pdb_info()->chain(aid.rsd()) <<
		ps.pdb_info()->number(aid.rsd()) <<
		"(" << aid.rsd() << ')' <<
		ps.residue(aid.rsd()).atom_name(aid.atomno());
}


/// @brief prints all information in this calculator
///
/// @prams[in] ps the pose at issue
///
void SHOBuriedUnsatisfiedPolarsCalculator::print_all_info(
	core::pose::Pose const& ps) const {

	std::cout << "----------------------------------------" << std::endl;
	std::cout << "target amino acid type: " << tgt_amino_ << std::endl;
	std::cout << "target atom name: " << tgt_atom_ << std::endl;

	std::cout << "target residue indexes in pose: ";
	Size const T = tgt_res_idxs_.size();
	if ( T ) {
		// ith output field contains ith target residue index (i=1,...,T)
		for ( Size i=1; i<=T; ++i ) {
			std::cout << tgt_res_idxs_[i] << ' ';
		}
	} else {
		if ( search_typ_ == ALL_ATOMS ) {
			std::cout << "all";
		} else {
			std::cout << "none";
		}
	}
	std::cout << std::endl;

	std::cout << "SHO cutoff: " << sho_cutoff_ << std::endl;
	std::cout << "----------------------------------------" << std::endl;

	std::cout << std::endl;
	std::cout << "NOTE: The energies output here are SHO energies even for H-bonded atoms" << std::endl;

	std::cout << std::endl << "BURIED UNSATISFIED ATOMS:" << std::endl;
	print_atom_subset(burunsat_atoms_, ps);

	std::cout << std::endl << "OTHER ATOMS:" << std::endl;
	print_atom_subset(other_atoms_, ps);

	std::cout << std::endl <<
		"number of buried unsatisfied atoms: " << burunsat_atoms_.size();

	std::cout << std::endl <<
		"number of other atoms: " << other_atoms_.size() <<
		std::endl;
}


/// @brief recomputes all data from scratch and prints it to screen
///
/// @param[in] ps pose for which the data are to be computed
///
/// @details the pose must already have been scored by the score function
///  pointed to by sfxn_.
///
void SHOBuriedUnsatisfiedPolarsCalculator::recompute_and_print(
	core::pose::Pose const& ps) {

	sho_meth_->init_hbond_data(ps);
	hbond_set_ = sho_meth_->hbond_set();

	recompute(ps);

	print_all_info(ps);
}


/// @brief extracts the pose indexes of a selected subset of residues
///
/// @param[in] setf path to a file specifying the residue subset according
///  to the following format:
///
///  C1 R1 I1\n
///  ...
///  CN RN IN\n
///
///   Here, Ci, Ri, and Ii indicate the chain identifier, residue index, and
///  insertion code (as specified in the pose's input PDB file) of the ith
///  residue in the subset (i=1,...,N; N>=1).
/// @param[out] rset vector to hold the residue indexes. The vector must be
///  passed empty.
/// @param[in] ps the pose.
///
/// @details: after this function has been called, rset[i] is the pose index
///  of the residue specified by the ith input line (i=1,...,N).
///
/// @details: blank chain identifiers and insertion codes must be specified
///  with the '_' character.
///
void residue_subset(std::string setf, utility::vector1<Size>& rset,
	core::pose::Pose &ps) {

	std::ifstream setfs(setf.c_str());
	if ( !setfs ) {
		TR << "can't open " << setf << std::endl;
		exit(0);
	}

	char cid;
	int idx;
	char ico;
	while ( setfs >> cid >> idx >> ico ) {

		if ( cid == '_' ) {
			cid = ' ';
		}

		if ( ico == '_' ) {
			ico = ' ';
		}

		rset.push_back(ps.pdb_info()->pdb2pose(cid, idx, ico));
	}
}

} // pose_metric_calculators
} // toolbox
} // protocols

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculator::SHOBuriedUnsatisfiedPolarsCalculator() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( search_typ_ ) ); // enum protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculator::SearchTyp
	// Don't serialize the energy method; instead recreate it from the command line;
	// since the only way it ever gets created is from the command line
	// arc( CEREAL_NVP( sho_meth_ ) ); // core::scoring::geometric_solvation::ExactOccludedHbondSolEnergyOP
	// EXEMPT sho_meth_
	arc( CEREAL_NVP( sho_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( tgt_res_idxs_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( tgt_amino_ ) ); // std::string
	arc( CEREAL_NVP( tgt_atom_ ) ); // std::string
	arc( CEREAL_NVP( hbond_set_ ) ); // core::scoring::hbonds::HBondSet
	arc( CEREAL_NVP( sho_energies_ ) ); // std::map<AtomID, Real>
	arc( CEREAL_NVP( burunsat_atoms_ ) ); // utility::vector1<AtomID>
	arc( CEREAL_NVP( num_burunsat_atoms_ ) ); // Size
	arc( CEREAL_NVP( other_atoms_ ) ); // utility::vector1<AtomID>
	arc( CEREAL_NVP( num_other_atoms_ ) ); // Size
	arc( CEREAL_NVP( sfxn_ ) ); // core::scoring::ScoreFunctionCOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( search_typ_ ); // enum protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculator::SearchTyp
	// arc( sho_meth_ ); // core::scoring::geometric_solvation::ExactOccludedHbondSolEnergyOP

	arc( sho_cutoff_ ); // core::Real
	arc( tgt_res_idxs_ ); // utility::vector1<Size>
	arc( tgt_amino_ ); // std::string
	arc( tgt_atom_ ); // std::string
	std::shared_ptr< core::scoring::hbonds::HBondSet > local_hbs;
	arc( local_hbs ); // core::scoring::hbonds::HBondSet
	hbond_set_ = local_hbs;
	arc( sho_energies_ ); // std::map<AtomID, Real>
	arc( burunsat_atoms_ ); // utility::vector1<AtomID>
	arc( num_burunsat_atoms_ ); // Size
	arc( other_atoms_ ); // utility::vector1<AtomID>
	arc( num_other_atoms_ ); // Size
	std::shared_ptr< core::scoring::ScoreFunction > local_sfxn;
	arc( local_sfxn ); // core::scoring::ScoreFunctionCOP
	sfxn_ = local_sfxn; // copy the non-const pointer(s) into the const pointer(s)

	// instead, recreate the sho_meth_ from the command line
	sho_meth_ = create_ExactSHOEnergy_from_cmdline(sfxn_->energy_method_options());
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculator );
CEREAL_REGISTER_TYPE( protocols::toolbox::pose_metric_calculators::SHOBuriedUnsatisfiedPolarsCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_SHOBuriedUnsatisfiedPolarsCalculator )
#endif // SERIALIZATION
