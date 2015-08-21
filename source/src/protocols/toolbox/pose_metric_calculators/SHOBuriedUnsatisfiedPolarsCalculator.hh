// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Calculator for buried unsatisfied polar atoms under the SHO solvation model
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#ifndef INCLUDED_SHOBuriedUnsatisfiedPolarsCalculator_hh
#define INCLUDED_SHOBuriedUnsatisfiedPolarsCalculator_hh

#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>
#include <string>
#include <map>

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

using core::Size;
using core::Real;
using core::id::AtomID;

class SHOBuriedUnsatisfiedPolarsCalculator :
	public core::pose::metrics::StructureDependentCalculator {

public:

	/// @brief constructs the calculator for a target atom in residues of a target
	///  amino acid type
	SHOBuriedUnsatisfiedPolarsCalculator(Real sho_cutoff, std::string tgt_amino,
		std::string tgt_atom, core::scoring::ScoreFunctionCOP sfxn);

	/// @brief constructs the calculator for a target residue
	SHOBuriedUnsatisfiedPolarsCalculator(Real sho_cutoff,
		utility::vector1<Size> const& tgt_res_idxs,
		core::scoring::ScoreFunctionCOP sfxn);

	/// @brief clones this calculator
	core::pose::metrics::PoseMetricCalculatorOP clone() const;

	/// @brief returns the set of target atoms classified as "buried unsatisfied"
	utility::vector1<AtomID> const& burunsat_atoms() {return burunsat_atoms_;}

	/// @brief returns the set of target atoms classified as "other"
	utility::vector1<AtomID> const& other_atoms() {return other_atoms_;}

	/// @brief prints all information in this calculator
	void print_all_info(core::pose::Pose const& ps) const;

	/// @brief returns the total H-bond energy of an atom
	Real hbond_energy(AtomID aid, core::pose::Pose const& ps) const;

	/// @brief recomputes all data from scratch and prints it to screen
	void recompute_and_print(core::pose::Pose const& ps);

protected:

	/// @brief outputs a quantity's value that is cached in the calculator
	void lookup(std::string const &key, basic::MetricValueBase* valptr) const;

	/// @brief returns the string representation of a quantity's value that is
	///  cached in the calculator
	std::string print(std::string const& key) const;

	/// @brief finds buried unsatisfied atoms by building all needed data structures
	///  from scratch
	void recompute(core::pose::Pose const& ps);

private:

	/// @brief identifiers for the allowed types of search
	enum SearchTyp {ANY_AA__TGT_ATOM, TGT_AA__TGT_ATOM, ALL_ATOMS, TGT_RES};

	/// @brief partitions target atoms into "buried unsatisfied" and "other"
	void partition(core::pose::Pose const& ps);

	/// @brief partitions a target residue's atoms into "buried unsatisfied"
	///  and "other"
	void residue_partition(core::pose::Pose const& ps, Size tgtidx);

	/// @brief prints on screen the SHO energies of all atoms in a pose
	void print_sho_energies(core::pose::Pose const& ps);

	/// @brief assigns an atom to either the "buried unsatisfied" class or the
	///  "other" class
	void assign_atom(AtomID aid);

	/// @brief is aid buried?
	bool is_buried(AtomID aid);

	/// @brief is aid unsatisfied?
	bool is_unsat(AtomID aid) const;

	/// @brief prints info on a subset of atoms
	void print_atom_subset(
		utility::vector1<AtomID> const& subset,
		core::pose::Pose const& ps) const;

	/// @brief returns the average SHO energy of a subset of atoms
	Real get_average_energy(
		utility::vector1<AtomID> const& subset,
		std::map<AtomID, Real> & sho_energies) const;

	/// @brief prints an atom's identifier on screen
	void print_atom_info(AtomID aid, core::pose::Pose const& ps) const;

	/// @brief selected type of search
	SearchTyp search_typ_;

	/// @brief SHO energy calculator
	core::scoring::geometric_solvation::ExactOccludedHbondSolEnergyOP sho_meth_;

	/// @brief maximum SHO energy value for a polar group to be considered as
	///  exposed (i.e., not buried)
	core::Real sho_cutoff_;

	/// @brief target residue indexes
	utility::vector1<Size> tgt_res_idxs_;

	/// @brief target amino acid type
	std::string tgt_amino_;

	/// @brief target atom type
	std::string tgt_atom_;

	/// @brief stores H-bond info for the entire pose
	core::scoring::hbonds::HBondSet hbond_set_;

	/// @brief stores the SHO energies of all polar atoms in the pose
	std::map<AtomID, Real> sho_energies_;

	/// @brief set of buried unsatisfied atoms
	utility::vector1<AtomID> burunsat_atoms_;

	/// @brief number of buried unsatisfied atoms
	Size num_burunsat_atoms_;

	/// @brief set of atoms that are NOT buried unsatisfied
	utility::vector1<AtomID> other_atoms_;

	/// @brief number of atoms that are NOT buried unsatisfied
	Size num_other_atoms_;

	/// @brief score function used to (previously) score the pose
	core::scoring::ScoreFunctionCOP sfxn_;
};

/// @brief extracts the pose indexes of a selected subset of residues
void residue_subset(std::string setf, utility::vector1<Size>& rset,
	core::pose::Pose &ps);

} // pose_metric_calculators
} // toolbox
} // protocols

#endif
