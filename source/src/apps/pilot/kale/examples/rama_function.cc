// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[]) {
	using core::Real;
	using core::Size;
	using core::scoring::Ramachandran;
	using core::scoring::ScoringManager;
	using core::chemical::aa_gly;

	devel::init(argc, argv);

	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	for (int phi = 0; phi < 360; phi += 1) {
		for (int psi = 0; psi < 360; psi += 1) {
			//Real value; rama.probability_of_sampling_phipsi(aa_gly, phi, psi, value);
			Real value = rama.eval_rama_score_residue(aa_gly, phi, psi);
			//bool value = rama.phipsi_in_allowed_rama(aa_gly, phi, psi);

			cout << setw(4) << phi << " " << setw(4) << psi << " ";
			cout << fixed << setprecision(10) << setw(15) << value << endl;
		}
	}
}	

// vim:fo-=a
