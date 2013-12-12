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

	Real phi, psi;
	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	for (int i = 1; i <= 1000000; i++) {
		rama.random_phipsi_from_rama(aa_gly, phi, psi);
		//rama.uniform_phipsi_from_allowed_rama(aa_gly, phi, psi);

		if (phi < 0) phi += 360;
		if (psi < 0) psi += 360;

		cout << fixed << setprecision(10) << setw(15) << phi << " ";
		cout << fixed << setprecision(10) << setw(15) << psi << endl;
	}
}	
