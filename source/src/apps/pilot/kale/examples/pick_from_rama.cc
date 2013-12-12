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
	using core::chemical::aa_ala;

	devel::init(argc, argv);

	Real phi, psi, prob;
	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	for (int i = 0; i < 1000000; i++) {
		rama.random_phipsi_from_rama(aa_ala, phi, psi, prob);
		cout << fixed << setw(8) << setprecision(3) << phi << "  ";
		cout << fixed << setw(8) << setprecision(3) << psi << "  ";
		cout << fixed << setw(8) << setprecision(5) << prob << endl;
	}
}	
