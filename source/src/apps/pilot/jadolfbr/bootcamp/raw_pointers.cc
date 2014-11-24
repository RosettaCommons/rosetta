

#include <vector>
#include <iostream>

///@brief A basic example class for bootcamp
///
///@details Within Rosetta, would be expanded into a 'Mover' type of class to perturb a specific subset 
///  of residues in Cartesian Space.
///
class CartPerturber {
	
public:

	//Constructors
	CartPerturber():
		rounds_(10)
	{}
	CartPerturber(int rounds):
		rounds_(rounds)
	{};
	
	CartPerturber(int rounds, std::vector<bool> const & residues):
		rounds_(rounds),
		residues_(residues)
	{};
	
	
	~CartPerturber(){};
	
	///@brief Print info out - would apply some change to a pose. Will apply to all residues if no residues were given.
	void
	apply(/*a pose would be here*/) {
		std::cout << "Applying CartPerturber for " << rounds_ << " rounds "<<std::endl;
		if (residues_.size() == 0){
			std::cout << "Applying to all residues" << std::endl;
		}
	}
	
	
	///////////////////////////////////////////////////////////////////
	// Setters
	//
	//
	void set_rounds(int rounds) {
		rounds_ = rounds;
	}
	
	void set_residues(std::vector<bool> const & residues){
		residues_ = residues;
	};
	
	
	/////////////////////////////////////////////////////////////////
	// Getters
	//
	//
	void add_rounds(int rounds) {
		rounds_ = rounds_ + rounds;
	}
	
	
	int get_rounds() const {
		return rounds_;
	}
	
	std::vector<bool> get_residues() const {
		return residues_;
	}
	
	

	
private:
	int rounds_;
	std::vector<bool> residues_;
	
};

//////////////////////////////////// Tasks //////////////////////////////////////////////////////////////////////////////////
//* allocating and deallocating CartPerturber
//* calling a function of CartPerturber with the -> operator
//* passing a pointer to CartPerturber into a function and having that function change CartPerturber
//* having a CartPerturber shared between two other classes





//////////////////////////////////// Example Completion/test (Remove before send) //////////////////////////////////

class Perturbers1 {

public:
	Perturbers1(){}
	~Perturbers1(){}
	
	void
	add_perturber(CartPerturber * perturber){
		
	}

private:
	std::vector<CartPerturber*> perturbers_;
	
};

class Perturbers2 {

public:
	Perturbers2(){}
	~Perturbers2(){}
	
	void
	add_perturber(CartPerturber * perturber){
		
	}

private:
	std::vector<CartPerturber*> perturbers_;
	
};


int main (){
		
	CartPerturber *perturber = new CartPerturber();
	perturber->set_rounds(5);
	perturber->add_rounds(10);
	
	std::vector<bool> residues(30, false);
	for (int i = 10; i <= 20; ++i){
		residues[i] = true;
	}
	
	perturber->apply();
	perturber->set_residues(residues);
	perturber->apply();
	
	
	
	Perturbers1 pert1 =  Perturbers1();
	Perturbers2 pert2 = Perturbers2();
	
	pert1.add_perturber(perturber);
	pert2.add_perturber(perturber);
			
	delete perturber;
	return(0);
	
}
