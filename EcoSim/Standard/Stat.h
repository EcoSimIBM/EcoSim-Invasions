#ifndef STAT_H_
#define STAT_H_

#include "Species.h"
#include "FCMPrey.h"
#include "Manipulation.h"

class Ecosystem;

using namespace std;

class Stat {
private:

	int nbMalePrey;
	int nbFemalePrey;
	int nbMalePred;
	int nbFemalePred;
	
	int nbInvMalePrey;  //inv equivalents
	int nbInvFemalePrey;
	int nbInvMalePred;
	int nbInvFemalePred;
	
	int nbSensPrey;
	int nbConceptsPrey;
	int nbMoteursDepPrey;
	int nbMoteursFixPrey;
	int nbSensPred;
	int nbConceptsPred;
	int nbMoteursDepPred;
	int nbMoteursFixPred;
	//MRE RandomGoodGene
	
	double energyAvgMalePrey;					//-- Average energy of all prey
	double strengthAvgMalePrey;					//-- Average strength of all prey
	double energyAvgMalePred;					//-- Average energy of all predators
	double strengthAvgMalePred;					//-- Average strength of all pred
	int nbBirthMalePrey;						//-- Number of new prey born
	int nbBirthMalePred;						//-- Number of new predators born
	double ageAvgDeadMalePrey;					//-- Average age of the prey that died
	double ageAvgDeadMalePred;					//-- Average age of the predators that died
	double avgDistEvMalePrey;					//-- Average distance evolved for the prey
	double avgDistEvMalePred;					//-- Average distance evolved for the predators
	
	//inv equivalents
	double energyAvgInvMalePrey;					//-- Average energy of all prey
	double strengthAvgInvMalePrey;					//-- Average strength of all prey
	double energyAvgInvMalePred;					//-- Average energy of all predators
	double strengthAvgInvMalePred;					//-- Average strength of all pred
	int nbBirthInvMalePrey;						//-- Number of new prey born
	int nbBirthInvMalePred;						//-- Number of new predators born
	double ageAvgDeadInvMalePrey;					//-- Average age of the prey that died
	double ageAvgDeadInvMalePred;					//-- Average age of the predators that died
	double avgDistEvInvMalePrey;					//-- Average distance evolved for the prey
	double avgDistEvInvMalePred;					//-- Average distance evolved for the predators
	
	
	double avgEnergySpentMalePrey;
	double avgEnergySpentFemalePrey;
	double avgEnergySpentMalePred;
	double avgEnergySpentFemalePred;
	//inv
	double avgEnergySpentInvMalePrey;
	double avgEnergySpentInvFemalePrey;
	double avgEnergySpentInvMalePred;
	double avgEnergySpentInvFemalePred;
	
	double avgEnergyTransferredMalePrey;
	double avgEnergyTransferredFemalePrey;
	double avgEnergyTransferredMalePred;
	double avgEnergyTransferredFemalePred;
	//inv
	double avgEnergyTransferredInvMalePrey;
	double avgEnergyTransferredInvFemalePrey;
	double avgEnergyTransferredInvMalePred;
	double avgEnergyTransferredInvFemalePred;

	double energyAvgFemalePrey;					//-- Average energy of all prey
	double strengthAvgFemalePrey;					//-- Average strength of all prey
	double energyAvgFemalePred;					//-- Average energy of all predators
	double strengthAvgFemalePred;					//-- Average strength of all pred
	int nbBirthFemalePrey;						//-- Number of new prey born
	int nbBirthFemalePred;						//-- Number of new predators born
	double ageAvgDeadFemalePrey;					//-- Average age of the prey that died
	double ageAvgDeadFemalePred;					//-- Average age of the predators that died
	double avgDistEvFemalePrey;					//-- Average distance evolved for the prey
	double avgDistEvFemalePred;					//-- Average distance evolved for the predators
	//inv
	double energyAvgInvFemalePrey;					//-- Average energy of all prey
	double strengthAvgInvFemalePrey;					//-- Average strength of all prey
	double energyAvgInvFemalePred;					//-- Average energy of all predators
	double strengthAvgInvFemalePred;					//-- Average strength of all pred
	int nbBirthInvFemalePrey;						//-- Number of new prey born
	int nbBirthInvFemalePred;						//-- Number of new predators born
	double ageAvgDeadInvFemalePrey;					//-- Average age of the prey that died
	double ageAvgDeadInvFemalePred;					//-- Average age of the predators that died
	double avgDistEvInvFemalePrey;					//-- Average distance evolved for the prey
	double avgDistEvInvFemalePred;					//-- Average distance evolved for the predators

	double avgDistMatingPrey;
	double avgDistMatingPred;
	//inv
	double avgDistMatingInvPrey;
	double avgDistMatingInvPred;

	vector <double> sumPreyMaleGenome;					//-- keep the sum of genetome M.M
	vector <double> sumPredMaleGenome;					//-- keep the sum of genetome M.M
	vector <double> sumPreyFemaleGenome;					//-- keep the sum of genetome M.M
	vector <double> sumPredFemaleGenome;					//-- keep the sum of genetome M.M
	//inv
	vector <double> sumInvPreyMaleGenome;					//-- keep the sum of genetome M.M
	vector <double> sumInvPredMaleGenome;					//-- keep the sum of genetome M.M
	vector <double> sumInvPreyFemaleGenome;					//-- keep the sum of genetome M.M
	vector <double> sumInvPredFemaleGenome;					//-- keep the sum of genetome M.M

	int nbDeadMalePreyA;						//-- Number of prey that died because of age
	int nbDeadMalePreyE;						//-- Number of prey that died because energy <= 0
	int nbDeadMalePreyK;						//-- Number of prey that died because they were eaten by predator
	int nbDeadMalePredA;						//-- Number of pred that died because of age
	int nbDeadMalePredE;						//-- Number of pred that died because energy <= 0
	int nbDeadMalePredF;						//-- Number of predator that died because of fighting with preys  M.M

	//inv
	int nbDeadInvMalePreyA;						//-- Number of prey that died because of age
	int nbDeadInvMalePreyE;						//-- Number of prey that died because energy <= 0
	int nbDeadInvMalePreyK;						//-- Number of prey that died because they were eaten by predator
	int nbDeadInvMalePredA;						//-- Number of pred that died because of age
	int nbDeadInvMalePredE;						//-- Number of pred that died because energy <= 0
	int nbDeadInvMalePredF;						//-- Number of predator that died because of fighting with preys  M.M

	int nbDeadFemalePreyA;						//-- Number of prey that died because of age
	int nbDeadFemalePreyE;						//-- Number of prey that died because energy <= 0
	int nbDeadFemalePreyK;						//-- Number of prey that died because they were eaten by predator
	int nbDeadFemalePredA;						//-- Number of pred that died because of age
	int nbDeadFemalePredE;						//-- Number of pred that died because energy <= 0
	int nbDeadFemalePredF;						//-- Number of predator that died because of fighting with preys  M.M
	//inv
	int nbDeadInvFemalePreyA;						//-- Number of prey that died because of age
	int nbDeadInvFemalePreyE;						//-- Number of prey that died because energy <= 0
	int nbDeadInvFemalePreyK;						//-- Number of prey that died because they were eaten by predator
	int nbDeadInvFemalePredA;						//-- Number of pred that died because of age
	int nbDeadInvFemalePredE;						//-- Number of pred that died because energy <= 0
	int nbDeadInvFemalePredF;						//-- Number of predator that died because of fighting with preys  M.M

	double stateAvgMalePrey;						//
	double stateAvgMalePred;						//
	double stateAvgFemalePrey;						//
	double stateAvgFemalePred;						//
	//inv
	double stateAvgInvMalePrey;						//
	double stateAvgInvMalePred;						//
	double stateAvgInvFemalePrey;						//
	double stateAvgInvFemalePred;						//
	
	int nbNativeSpeciesPrey;
	int nbNativeSpeciesPred;
	int nbInvasiveSpeciesPrey;
	int nbInvasiveSpeciesPred;
	
	double persuasionTotalMalePrey;
	double persuasionTotalMalePred;
	double persuasionTotalFemalePrey;
	double persuasionTotalFemalePred;
	
	//inv
	double persuasionTotalInvMalePrey;
	double persuasionTotalInvMalePred;
	double persuasionTotalInvFemalePrey;
	double persuasionTotalInvFemalePred;

	double speedAvgMalePrey;						//-- Average speed of the prey
	double speedAvgMalePred;						//-- Average speed of the predators
	double speedAvgFemalePrey;						//-- Average speed of the prey
	double speedAvgFemalePred;						//-- Average speed of the predators
	//inv
	double speedAvgInvMalePrey;						//-- Average speed of the prey
	double speedAvgInvMalePred;						//-- Average speed of the predators
	double speedAvgInvFemalePrey;						//-- Average speed of the prey
	double speedAvgInvFemalePred;						//-- Average speed of the predators

	int ageMalePrey;							//
	int ageMalePred;							//
	int ageFemalePrey;							//
	int ageFemalePred;							//
	//inv
	int ageInvMalePrey;							//
	int ageInvMalePred;							//
	int ageInvFemalePrey;							//
	int ageInvFemalePred;							//

	int nbPreyByCase;						//
	int nbPredByCase;						//
	int nbCasePrey;							//
	int nbCasePred;							//
	//inv
	int nbInvPreyByCase;						//
	int nbInvPredByCase;						//
	int nbCaseInvPrey;							//
	int nbCaseInvPred;							//

	int nbArcAvgMalePrey;						//-- Average number of arcs for the MalePrey
	int nbArcAvgMalePred;						//-- Average number of arcs for the MalePredators
	//inv
	int nbArcAvgInvMalePrey;						//-- Average number of arcs for the MalePrey
	int nbArcAvgInvMalePred;						//-- Average number of arcs for the MalePredators

	int nbArcAvgFemalePrey;						//-- Average number of arcs for the FemalePrey
	int nbArcAvgFemalePred;						//-- Average number of arcs for the FemalePredators
	//inv
	int nbArcAvgInvFemalePrey;						//-- Average number of arcs for the FemalePrey
	int nbArcAvgInvFemalePred;						//-- Average number of arcs for the FemalePredators

	long double nbTotalGrass;						//-- Total grass in the world
#ifdef TWO_RESOURCES
	double nbTotalGrass2;					//-- Armin
#endif
	long double nbTotalMeat;						//-- Total meat in the world

	//vector<vector<int> > distribAge;		//-- Distribution of ages

	vector<int> nbActionPreyMale;				//-- Actions of the prey
	vector<int> nbActionPreyFemale;				//-- Actions of the prey
	vector<int> nbActionPredMale;				//-- Actions of the predators
	vector<int> nbActionPredFemale;				//-- Actions of the predators
	vector<float> activations_prey_male;			//-- Prey's concepts
	vector<float> activations_prey_female;			//-- Prey's concepts
	vector<float> activations_pred_male;			//-- Predator's concepts
	vector<float> activations_pred_female;			//-- Predator's concepts
	//inv
	vector<int> nbActionInvPreyMale;				//-- Actions of the prey
	vector<int> nbActionInvPreyFemale;				//-- Actions of the prey
	vector<int> nbActionInvPredMale;				//-- Actions of the predators
	vector<int> nbActionInvPredFemale;				//-- Actions of the predators
	vector<float> activations_Invprey_male;			//-- Prey's concepts
	vector<float> activations_Invprey_female;			//-- Prey's concepts
	vector<float> activations_Invpred_male;			//-- Predator's concepts
	vector<float> activations_Invpred_female;			//-- Predator's concepts

	double gEntropy;						    //-- Marwa: Global entropy for prey population
	double maxEntropy;						//-- Marwa: Maximum global entropy for prey population
	//inv
	double gInvEntropy;						    //-- Marwa: Global entropy for prey population
	double maxInvEntropy;						//-- Marwa: Maximum global entropy for prey population

	float SpeciesRatioPrey;					//-- Meisam
	float SpeciesRatioPred;					//-- Meisam
	float ExtinctionRatioPrey;				//-- Meisam
	float ExtinctionRatioPred;				//-- Meisam
	//inv
	float SpeciesRatioInvPrey;					//-- Meisam
	float SpeciesRatioInvPred;					//-- Meisam
	float ExtinctionRatioInvPrey;				//-- Meisam
	float ExtinctionRatioInvPred;				//-- Meisam

	Manipulation Manip;							//-- Class Instances: Manipulating the input

public:

	Stat();
	virtual ~Stat();

	Stat(int, int, int, int, int, int, int, int, Ecosystem *);

	void setSpeciesRatioPrey(float s)		{ SpeciesRatioPrey = s; };	//-- Meisam
	void setSpeciesRatioPred(float s)		{ SpeciesRatioPred = s; };	//-- Meisam
	void setExtinctionRatioPrey(float s)  	{ ExtinctionRatioPrey = s; };	//-- Meisam
	void setExtinctionRatioPred(float s)	{ ExtinctionRatioPred = s; };	//-- Meisam
	//Inv
	void setSpeciesRatioInvPrey(float s)		{ SpeciesRatioInvPrey = s; };	//-- Meisam
	void setSpeciesRatioInvPred(float s)		{ SpeciesRatioInvPred = s; };	//-- Meisam
	void setExtinctionRatioInvPrey(float s)  	{ ExtinctionRatioInvPrey = s; };	//-- Meisam
	void setExtinctionRatioInvPred(float s)	{ ExtinctionRatioInvPred = s; };	//-- Meisam

	void incNbMalePrey()					{ nbMalePrey++; };	//-- RS
	void incNbFemalePrey()					{ nbFemalePrey++; };	//-- RS
	void incNbMalePred()  					{ nbMalePred++; };	//-- RS
	void incNbFemalePred()					{ nbFemalePred++; };	//-- RS
	//inv
	void incNbInvMalePrey()					{ nbInvMalePrey++; };	//-- RS
	void incNbInvFemalePrey()				{ nbInvFemalePrey++; };	//-- RS
	void incNbInvMalePred()  				{ nbInvMalePred++; };	//-- RS
	void incNbInvFemalePred()				{ nbInvFemalePred++; };	//-- RS

	void reset();
	
	void incNbArcAvgMalePrey(int a)				{ nbArcAvgMalePrey += a; };
	void incNbArcAvgMalePred(int a) 			{ nbArcAvgMalePred += a; };
	//inv
	void incNbArcAvgInvMalePrey(int a)			{ nbArcAvgInvMalePrey += a; };
	void incNbArcAvgInvMalePred(int a) 			{ nbArcAvgInvMalePred += a; };

	void incNbPreyByCase(int n)  			{ nbPreyByCase += n; };
	void incNbPredByCase(int n)    			{ nbPredByCase += n; };
	void incNbCasePrey()   					{ nbCasePrey++; };
	void incNbCasePred()  					{ nbCasePred++; };
	//inv
	void incNbInvPreyByCase(int n)  			{ nbInvPreyByCase += n; };
	void incNbInvPredByCase(int n)    			{ nbInvPredByCase += n; };
	void incNbCaseInvPrey()   					{ nbCaseInvPrey++; };
	void incNbCaseInvPred()  					{ nbCaseInvPred++; };
	
	//
	int getNbPreyByCase()  			{ return nbPreyByCase; };
	int getNbPredByCase()    			{ return nbPredByCase; };
	int getNbCasePrey()   					{ return nbCasePrey; };
	int getNbCasePred()  					{ return nbCasePred; };
	//inv
	int getNbInvPreyByCase()  			{ return nbInvPreyByCase; };
	int getNbInvPredByCase()    			{ return nbInvPredByCase; };
	int getNbCaseInvPrey()   					{ return nbCaseInvPrey; };
	int getNbCaseInvPred()  					{ return nbCaseInvPred; };
	//

	void incPersuasionMalePrey(float a)				{ persuasionTotalMalePrey += a; };
	void incPersuasionMalePred(float a)				{ persuasionTotalMalePred += a; };
	void incPersuasionFemalePrey(float a)				{ persuasionTotalFemalePrey += a; };
	void incPersuasionFemalePred(float a)				{ persuasionTotalFemalePred += a; };
	//inv
	void incPersuasionInvMalePrey(float a)				{ persuasionTotalInvMalePrey += a; };
	void incPersuasionInvMalePred(float a)				{ persuasionTotalInvMalePred += a; };
	void incPersuasionInvFemalePrey(float a)				{ persuasionTotalInvFemalePrey += a; };
	void incPersuasionInvFemalePred(float a)				{ persuasionTotalInvFemalePred += a; };
	
	void incNbArcAvgFemalePrey(int a)				{ nbArcAvgFemalePrey += a; };
	void incNbArcAvgFemalePred(int a) 			{ nbArcAvgFemalePred += a; };
	//inv
	void incNbArcAvgInvFemalePrey(int a)				{ nbArcAvgInvFemalePrey += a; };
	void incNbArcAvgInvFemalePred(int a) 			{ nbArcAvgInvFemalePred += a; };
	
	void incEnergySpentMalePrey(float a)				{ avgEnergySpentMalePrey += a; };
	void incEnergySpentFemalePrey(float a) 			{ avgEnergySpentFemalePrey += a; };
	void incEnergySpentMalePred(float a)				{ avgEnergySpentMalePred += a; };
	void incEnergySpentFemalePred(float a) 			{ avgEnergySpentFemalePred += a; };
	//inv
	void incEnergySpentInvMalePrey(float a)				{ avgEnergySpentInvMalePrey += a; };
	void incEnergySpentInvFemalePrey(float a) 			{ avgEnergySpentInvFemalePrey += a; };
	void incEnergySpentInvMalePred(float a)				{ avgEnergySpentInvMalePred += a; };
	void incEnergySpentInvFemalePred(float a) 			{ avgEnergySpentInvFemalePred += a; };
	
	void incAvgEnergyTransferredMalePrey(float a)				{ avgEnergyTransferredMalePrey += a; };
	void incAvgEnergyTransferredFemalePrey(float a) 			{ avgEnergyTransferredFemalePrey += a; };
	void incAvgEnergyTransferredMalePred(float a)				{ avgEnergyTransferredMalePred += a; };
	void incAvgEnergyTransferredFemalePred(float a) 			{ avgEnergyTransferredFemalePred += a; };
	//inv
	void incAvgEnergyTransferredInvMalePrey(float a)				{ avgEnergyTransferredInvMalePrey += a; };
	void incAvgEnergyTransferredInvFemalePrey(float a) 			{ avgEnergyTransferredInvFemalePrey += a; };
	void incAvgEnergyTransferredInvMalePred(float a)				{ avgEnergyTransferredInvMalePred += a; };
	void incAvgEnergyTransferredInvFemalePred(float a) 			{ avgEnergyTransferredInvFemalePred += a; };

	void incNbTotalGrass(float n) 			{ nbTotalGrass += (double) n; };
#ifdef TWO_RESOURCES
	void incNbTotalGrass2(float n) 			{ nbTotalGrass2 += (double) n; };
#endif
	void incNbTotalMeat(float n) 			{ nbTotalMeat += n; };
	
	void incStateAvgMalePred(float e) 			{ stateAvgMalePred += e; };
	void incStateAvgMalePrey(float e) 			{ stateAvgMalePrey += e; };
	void incStateAvgFemalePred(float e) 			{ stateAvgFemalePred += e; };
	void incStateAvgFemalePrey(float e) 			{ stateAvgFemalePrey += e; };
	void incNbActionPreyMale(int i)  			{ nbActionPreyMale.at(i)++; };
	void incNbActionPreyFemale(int i)  			{ nbActionPreyFemale.at(i)++; };
	void incNbActionPredMale(int i)  			{ nbActionPredMale.at(i)++; };
	void incNbActionPredFemale(int i)  			{ nbActionPredFemale.at(i)++; };
	//inv
	void incStateAvgInvMalePred(float e) 			{ stateAvgInvMalePred += e; };
	void incStateAvgInvMalePrey(float e) 			{ stateAvgInvMalePrey += e; };
	void incStateAvgInvFemalePred(float e) 			{ stateAvgInvFemalePred += e; };
	void incStateAvgInvFemalePrey(float e) 			{ stateAvgInvFemalePrey += e; };
	void incNbActionInvPreyMale(int i)  			{ nbActionInvPreyMale.at(i)++; };
	void incNbActionInvPreyFemale(int i)  			{ nbActionInvPreyFemale.at(i)++; };
	void incNbActionInvPredMale(int i)  			{ nbActionInvPredMale.at(i)++; };
	void incNbActionInvPredFemale(int i)  			{ nbActionInvPredFemale.at(i)++; };

	void incNbDeadMalePreyA()  					{ nbDeadMalePreyA++; };
	void incNbDeadMalePreyE()  					{ nbDeadMalePreyE++; };
	void incNbDeadMalePreyK()  					{ nbDeadMalePreyK++; };
	
	void incNbDeadMalePredA()  					{ nbDeadMalePredA++; };
	void incNbDeadMalePredE()  					{ nbDeadMalePredE++; };
	void incnbDeadMalePredF() 					{ nbDeadMalePredF++; };  		// M.M
	
	void incNbDeadFemalePreyA()  					{ nbDeadFemalePreyA++; };
	void incNbDeadFemalePreyE()  					{ nbDeadFemalePreyE++; };
	void incNbDeadFemalePreyK()  					{ nbDeadFemalePreyK++; };
	
	void incNbDeadFemalePredA()  					{ nbDeadFemalePredA++; };
	void incNbDeadFemalePredE()  					{ nbDeadFemalePredE++; };
	void incnbDeadFemalePredF() 					{ nbDeadFemalePredF++; };  		// M.M
	
	//inv
	void incNbDeadInvMalePreyA()  					{ nbDeadInvMalePreyA++; };
	void incNbDeadInvMalePreyE()  					{ nbDeadInvMalePreyE++; };
	void incNbDeadInvMalePreyK()  					{ nbDeadInvMalePreyK++; };
	
	void incNbDeadInvMalePredA()  					{ nbDeadInvMalePredA++; };
	void incNbDeadInvMalePredE()  					{ nbDeadInvMalePredE++; };
	void incnbDeadInvMalePredF() 					{ nbDeadInvMalePredF++; };  		// M.M
	
	void incNbDeadInvFemalePreyA()  				{ nbDeadInvFemalePreyA++; };
	void incNbDeadInvFemalePreyE()  				{ nbDeadInvFemalePreyE++; };
	void incNbDeadInvFemalePreyK()  				{ nbDeadInvFemalePreyK++; };
	
	void incNbDeadInvFemalePredA()  				{ nbDeadInvFemalePredA++; };
	void incNbDeadInvFemalePredE()  				{ nbDeadInvFemalePredE++; };
	void incnbDeadInvFemalePredF() 					{ nbDeadInvFemalePredF++; };  		// M.M

	void inc_activation_prey_male(int i, float c)	{ activations_prey_male.at(i) += c; };
	void inc_activation_prey_female(int i, float c)	{ activations_prey_female.at(i) += c; };
	void inc_activation_pred_male(int i, float c)	{ activations_pred_male.at(i) += c; };
	void inc_activation_pred_female(int i, float c)	{ activations_pred_female.at(i) += c; };
	void incEnergyAvgMalePrey(float e)  			{ energyAvgMalePrey += e; };
	void incStrengthAvgMalePrey(float s)	 		{ strengthAvgMalePrey += s; };  // M.M
	void setEnergyAvgMalePrey(float e) 				{ energyAvgMalePrey = e; };
	void setStrengthAvgMalePrey(float s)	 		{ strengthAvgMalePrey = s; };   // M.M
	void incEnergyAvgMalePred(float e) 				{ energyAvgMalePred += e; };
	void incStrengthAvgMalePred(float s)	 		{ strengthAvgMalePred += s; };  // M.M
	void setEnergyAvgMalePred(float e) 				{ energyAvgMalePred = e; };
	void setStrengthAvgMalePred(float s) 			{ strengthAvgMalePred = s; };   // M.M
	void incNbBirthMalePrey() 						{ nbBirthMalePrey++; };
	void incNbBirthMalePred() 						{ nbBirthMalePred++; };
	//inv
	void inc_activation_Invprey_male(int i, float c)	{ activations_Invprey_male.at(i) += c; };
	void inc_activation_Invprey_female(int i, float c)	{ activations_Invprey_female.at(i) += c; };
	void inc_activation_Invpred_male(int i, float c)	{ activations_Invpred_male.at(i) += c; };
	void inc_activation_Invpred_female(int i, float c)	{ activations_Invpred_female.at(i) += c; };
	void incEnergyAvgInvMalePrey(float e)  			{ energyAvgInvMalePrey += e; };
	void incStrengthAvgInvMalePrey(float s)	 		{ strengthAvgInvMalePrey += s; };  // M.M
	void setEnergyAvgInvMalePrey(float e) 				{ energyAvgInvMalePrey = e; };
	void setStrengthAvgInvMalePrey(float s)	 		{ strengthAvgInvMalePrey = s; };   // M.M
	void incEnergyAvgInvMalePred(float e) 				{ energyAvgInvMalePred += e; };
	void incStrengthAvgInvMalePred(float s)	 		{ strengthAvgInvMalePred += s; };  // M.M
	void setEnergyAvgInvMalePred(float e) 				{ energyAvgInvMalePred = e; };
	void setStrengthAvgInvMalePred(float s) 			{ strengthAvgInvMalePred = s; };   // M.M
	void incNbBirthInvMalePrey() 						{ nbBirthInvMalePrey++; };
	void incNbBirthInvMalePred() 						{ nbBirthInvMalePred++; };

	void incEnergyAvgFemalePrey(float e)  			{ energyAvgFemalePrey += e; };
	void incStrengthAvgFemalePrey(float s) 			{ strengthAvgFemalePrey += s; };  // M.M
	void setEnergyAvgFemalePrey(float e) 			{ energyAvgFemalePrey = e; };
	void setStrengthAvgFemalePrey(float s) 			{ strengthAvgFemalePrey = s; };   // M.M
	void incEnergyAvgFemalePred(float e) 			{ energyAvgFemalePred += e; };
	void incStrengthAvgFemalePred(float s) 			{ strengthAvgFemalePred += s; };  // M.M
	void setEnergyAvgFemalePred(float e) 			{ energyAvgFemalePred = e; };
	void setStrengthAvgFemalePred(float s) 			{ strengthAvgFemalePred = s; };   // M.M
	void incNbBirthFemalePrey() 					{ nbBirthFemalePrey++; };
	void incNbBirthFemalePred() 					{ nbBirthFemalePred++; };
	void incSpeedAvgMalePrey(float e) 				{ speedAvgMalePrey += e; };//MRE RandomGoodGene
	void setSpeedAvgFemalePrey(float e)  			{ speedAvgFemalePrey = e; };//MRE RandomGoodGene
	void incSpeedAvgMalePred(float e) 				{ speedAvgMalePred += e; };//MRE RandomGoodGene
	void setSpeedAvgFemalePred(float e) 			{ speedAvgFemalePred = e; };//MRE RandomGoodGene
	void incSpeedAvgFemalePrey(float e) 			{ speedAvgFemalePrey += e; };//MRE RandomGoodGene
	void setSpeedAvgMalePrey(float e)  				{ speedAvgMalePrey = e; };//MRE RandomGoodGene
	void incSpeedAvgFemalePred(float e) 			{ speedAvgFemalePred += e; };//MRE RandomGoodGene
	void setSpeedAvgMalePred(float e) 				{ speedAvgMalePred = e; };//MRE RandomGoodGene
	void incAgeAvgDeadMalePrey(int a) 				{ ageAvgDeadMalePrey += a; };
	void incAgeAvgDeadMalePred(int a) 				{ ageAvgDeadMalePred += a; };
	void incAgeMalePrey(int a) 						{ ageMalePrey += a; };
	void incAgeMalePred(int a) 						{ ageMalePred += a; };
	void incAvgDistEvMalePrey(float g) 				{ avgDistEvMalePrey += g; };
	void incAvgDistEvMalePred(float g) 				{ avgDistEvMalePred += g; };
	void incAgeAvgDeadFemalePrey(int a) 			{ ageAvgDeadFemalePrey += a; };
	void incAgeAvgDeadFemalePred(int a) 			{ ageAvgDeadFemalePred += a; };
	void incAgeFemalePrey(int a) 					{ ageFemalePrey += a; };
	void incAgeFemalePred(int a) 					{ ageFemalePred += a; };
	void incAvgDistEvFemalePrey(float g) 			{ avgDistEvFemalePrey += g; };
	void incAvgDistEvFemalePred(float g) 			{ avgDistEvFemalePred += g; };
	void incAvgDistMatingPrey(float g) 				{ avgDistMatingPrey += g; };	    //-- Meisma: Distance Parents
	void incAvgDistMatingPred(float g) 				{ avgDistMatingPred += g; };		//-- Meisma: Distance Parents

	//inv
	void incEnergyAvgInvFemalePrey(float e)  			{ energyAvgInvFemalePrey += e; };
	void incStrengthAvgInvFemalePrey(float s) 			{ strengthAvgInvFemalePrey += s; };  // M.M
	void setEnergyAvgInvFemalePrey(float e) 			{ energyAvgInvFemalePrey = e; };
	void setStrengthAvgInvFemalePrey(float s) 			{ strengthAvgInvFemalePrey = s; };   // M.M
	void incEnergyAvgInvFemalePred(float e) 			{ energyAvgInvFemalePred += e; };
	void incStrengthAvgInvFemalePred(float s) 			{ strengthAvgInvFemalePred += s; };  // M.M
	void setEnergyAvgInvFemalePred(float e) 			{ energyAvgInvFemalePred = e; };
	void setStrengthAvgInvFemalePred(float s) 			{ strengthAvgInvFemalePred = s; };   // M.M
	void incNbBirthInvFemalePrey() 					{ nbBirthInvFemalePrey++; };
	void incNbBirthInvFemalePred() 					{ nbBirthInvFemalePred++; };
	void incSpeedAvgInvMalePrey(float e) 				{ speedAvgInvMalePrey += e; };//MRE RandomGoodGene
	void setSpeedAvgInvFemalePrey(float e)  			{ speedAvgInvFemalePrey = e; };//MRE RandomGoodGene
	void incSpeedAvgInvMalePred(float e) 				{ speedAvgInvMalePred += e; };//MRE RandomGoodGene
	void setSpeedAvgInvFemalePred(float e) 			{ speedAvgInvFemalePred = e; };//MRE RandomGoodGene
	void incSpeedAvgInvFemalePrey(float e) 			{ speedAvgInvFemalePrey += e; };//MRE RandomGoodGene
	void setSpeedAvgInvMalePrey(float e)  				{ speedAvgInvMalePrey = e; };//MRE RandomGoodGene
	void incSpeedAvgInvFemalePred(float e) 			{ speedAvgInvFemalePred += e; };//MRE RandomGoodGene
	void setSpeedAvgInvMalePred(float e) 				{ speedAvgInvMalePred = e; };//MRE RandomGoodGene
	void incAgeAvgDeadInvMalePrey(int a) 				{ ageAvgDeadInvMalePrey += a; };
	void incAgeAvgDeadInvMalePred(int a) 				{ ageAvgDeadInvMalePred += a; };
	void incAgeInvMalePrey(int a) 						{ ageInvMalePrey += a; };
	void incAgeInvMalePred(int a) 						{ ageInvMalePred += a; };
	void incAvgDistEvInvMalePrey(float g) 				{ avgDistEvInvMalePrey += g; };
	void incAvgDistEvInvMalePred(float g) 				{ avgDistEvInvMalePred += g; };
	void incAgeAvgDeadInvFemalePrey(int a) 			{ ageAvgDeadInvFemalePrey += a; };
	void incAgeAvgDeadInvFemalePred(int a) 			{ ageAvgDeadInvFemalePred += a; };
	void incAgeInvFemalePrey(int a) 					{ ageInvFemalePrey += a; };
	void incAgeInvFemalePred(int a) 					{ ageInvFemalePred += a; };
	void incAvgDistEvInvFemalePrey(float g) 			{ avgDistEvInvFemalePrey += g; };
	void incAvgDistEvInvFemalePred(float g) 			{ avgDistEvInvFemalePred += g; };
	void incAvgDistMatingInvPrey(float g) 				{ avgDistMatingInvPrey += g; };	    //-- Meisma: Distance Parents
	void incAvgDistMatingInvPred(float g) 				{ avgDistMatingInvPred += g; };		//-- Meisma: Distance Parents

	void writeHeader(Ecosystem *);
	void writeHeaderInv(Ecosystem *);
	void writeStat(Ecosystem *);
	void catchupStatInv(Ecosystem *); 
	void writeStatInv(Ecosystem *);
	void RemoveEnd(Ecosystem *);		//-- Meisam: Remove additional generations after Restoring

	void globalEntropy(Ecosystem *);    //-- Marwa: added for calculating global entropy  and maximum entropy for prey population
	void globalInvEntropy(Ecosystem *);
};

#endif /* STAT_H_ */
