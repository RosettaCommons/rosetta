/// @file
/// @brief


#include <devel/init.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/comparative_modeling/hybridize/TMalign.hh>
#include <protocols/comparative_modeling/hybridize/util.hh>
#include <utility/vector1.hh>


void
superimpose_tmalign( core::pose::Pose const & ref_pose, core::pose::Pose & pose) {
	protocols::comparative_modeling::hybridize::TMalign tm_align;
	int reval = tm_align.apply(pose, ref_pose);
	if (reval != 0) {
		std::cerr << "ERROR!" << std::endl;
		return;
	}

	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID );
	std::string seq_pose, seq_ref, aligned;
	core::Size n_mapped_residues = 0;

	tm_align.alignment2AtomMap(pose, ref_pose, n_mapped_residues, atom_map);
	tm_align.alignment2strings(seq_pose, seq_ref, aligned);

	if (n_mapped_residues >= 6) {
		utility::vector1< core::Real > aln_cutoffs;
		aln_cutoffs.push_back(6);
		aln_cutoffs.push_back(4);
		aln_cutoffs.push_back(3);
		aln_cutoffs.push_back(2);
		aln_cutoffs.push_back(1.5);
		aln_cutoffs.push_back(1);
		core::Real min_coverage = 0.2;
		protocols::comparative_modeling::hybridize::partial_align(pose, ref_pose, atom_map, true, aln_cutoffs, min_coverage); // iterate_convergence = true
	}
}

void
multimodel_gdt( core::pose::Pose const & ref_pose, utility::vector1<core::pose::Pose> const &models ) {
	core::Size nmodels = models.size();
	core::Size nres = ref_pose.total_residue();

	utility::vector1<core::Real> best_align(nres, (core::Real)999.0);

	for (int j=1; j<=nmodels; ++j) {
		core::Real r1=0,r2=0,r3=0,r4=0,r7=0;

		for (int i=1; i<=nres; ++i) {
			core::Size n = ref_pose.pdb_info()->number(i);
			numeric::xyzVector<core::Real> x_n = ref_pose.residue(i).xyz(2);

			core::Size n_idx = models[j].pdb_info()->pdb2pose( 'A', n );
			if (n_idx == 0) 
				n_idx = models[j].pdb_info()->pdb2pose( ' ', n );

			if (n_idx != 0) {
				numeric::xyzVector<core::Real> y_n = models[j].residue(n_idx).xyz(2);
				best_align[i] = std::min( best_align[i], (y_n-x_n).length());
			}

			if (best_align[i]<=1.0) r1+=1.0;
			if (best_align[i]<=2.0) r2+=1.0;
			if (best_align[i]<=3.0) r3+=1.0;
			if (best_align[i]<=4.0) r4+=1.0;
			if (best_align[i]<=7.0) r7+=1.0;
		}

		r1/=nres; r2/=nres; r3/=nres; r4/=nres; r7/=nres;
		std::cout << j << " models:  maxsub1=" << r1 << "  maxsub2=" << r2 << "  maxsub3=" << r3 << "  maxsub4=" << r4 << "  maxsub7=" << r7 
		  << "  gdtmm=" << (r1+r2+r3+r4+r7)/5.0 << std::endl;
	}
}


int
main( int argc, char * argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	devel::init(argc, argv);
	
	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	utility::vector1<core::pose::Pose> models = core::import_pose::poses_from_pdbs( *residue_set, option[OptionKeys::in::file::s]());
	core::pose::PoseOP native = core::import_pose::pose_from_pdb(*residue_set, option[OptionKeys::in::file::native]());

	for (int i=1; i<=models.size(); ++i) {
		superimpose_tmalign ( *native , models[i]);
	}
	multimodel_gdt( *native, models );
}

