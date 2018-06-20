/** MonteCarlo.cpp
*
*  Anna Sinelnikova
*  Uppsala, Sweden 2017
*/

#include "PCMC/MonteCarlo.h"
#include "PCMC/ReadWriteFiles/ParamFileReader.h"
#include "PCMC/PCAmacros.h"
#include "PCMC/Timer.h"
#include "PCMC/PolymerObservable.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

namespace PCA{

MonteCarlo::MonteCarlo(const Dictionary& dictionary)
{
    std::string etalon;
    int number;
    
    etalon = "THREADS_PER_CORE";
    number =dictionary.searchAndCheck(etalon);
    if(number>=0)
	loopsPerCore = (int)dictionary.value(number);
    else if (number==-1)
	loopsPerCore = 1;
    else{
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once.\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MAX_LOG_T";
    maxLogT = dictionary[etalon];
    
    etalon = "MIN_LOG_T";
    minLogT = dictionary[etalon];
    
    etalon = "LOG_T_STEP";
    logTstep = dictionary[etalon];
    
    etalon = "SWEEPS_PER_STEP";
    sweepsPerStep = (int)dictionary[etalon];

    etalon = "CORES";
    number =dictionary.searchAndCheck(etalon);
    if(number>=0)
	cores = (int)dictionary.value(number);
    else if (number==-1)
	cores = 0;
    else{
	printf("Error in format of parameter file .pcap:\n");
	printf("Parameter '%s' defines more than once.\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "REGIME";
    regime = (Regime)dictionary[etalon];
    stepsPerLoop = (int)((maxLogT-minLogT)/logTstep)+1;
    
}

MonteCarlo::MonteCarlo(const char* fileName, PolymerMC* polymer_in,
			const Hamiltonian* hamiltonian_in,
			const Interaction* interaction_in)
{
    readFromParamFile(fileName);
    stepsPerLoop = (int)((maxLogT-minLogT)/logTstep)+1;
    
    polymerEtalon = polymer_in;
    hamiltonian = hamiltonian_in;
    interaction = interaction_in;
    
    if(regime==Regime::twoChains)
	_PCA_ERROR("MonteCarlo constructor:\n\tFor 2 chains you should use another constructor.");
}

MonteCarlo::MonteCarlo(const char* fileName, PolymerMC* polymer_in1,
			PolymerMC* polymer_in2,
			const Hamiltonian* hamiltonian_in,
			const Interaction* interaction_in,
			double minDist_in)
{
    readFromParamFile(fileName);
    stepsPerLoop = (int)((maxLogT-minLogT)/logTstep)+1;
    
    polymerEtalon = polymer_in1;
    polymerEtalon2 = polymer_in2;
    hamiltonian = hamiltonian_in;
    interaction = interaction_in;
    regime = Regime::twoChains;
    minDist = minDist_in;
}

MonteCarlo::~MonteCarlo(){};

void MonteCarlo::run(PolymerMC* polymer_in,
			const Hamiltonian* hamiltonian_in,
			const Interaction* interaction_in,
			int myCoreNumber, int totalCoreNumber,
			double minDist_in)
{	
    polymerEtalon = polymer_in;
    polymerEtalon2 = nullptr;
    if(regime==Regime::twoChains)
	polymerEtalon2 = polymer_in;
    hamiltonian = hamiltonian_in;
    interaction = interaction_in;
    minDist = minDist_in;
    cores = totalCoreNumber;

    int fakeCoreNumber; //effective core number
    char* fname1;
    FILE* ktfp, *ktfp2, *logfp, *tfp, *checkfp, *conffp, *conffp2,
	*rgyrFinalfp, *rgyrCurrentfp,*eFinalfp, *eCurrentfp;
    PolymerMC *polymer;
    PolymerMC *polymer2 = nullptr;
    if(myCoreNumber==0){
	logfp = fopen("logs/monteCarlo.log", "w");
	_PCMC_WRITE_RUNTIME_CONTEXT(logfp);
	fprintf(logfp,"********\n* PCMC * \n********\n");
	fprintf(logfp,"CPU MULTY CORE\n\n");
	fprintf(logfp,"Cores %i\n", totalCoreNumber);
	fprintf(logfp,"fakeCores %i\n", loopsPerCore);
	fprintf(logfp,"Steps/fakeCore: %i\n\n", stepsPerLoop);
	fprintf(logfp,"Regime: %s\n", getRegimeStr().c_str());
	checkfp = fopen("results/paramFileLog.pcmc","w");
	std::string info = PCMC_COMMENTED_RUNTIME_CONTEXT_STRING;
	fprintf(checkfp,"%s\n",info.c_str());
	fprintf(checkfp, "# The file was automaticaly generated by MonteCarlo::run(.).\n");
	polymerEtalon->writeInParamFile(checkfp);
	if(hamiltonian!=nullptr)
	    hamiltonian->writeInParamFile(checkfp);
	if(interaction!=nullptr)
	    interaction->writeInParamFile(checkfp);
	this->writeInParamFile(checkfp);
	fclose(checkfp);
    }
    double temperature;
    for(int k=0; k<loopsPerCore; k++){
	
	if(myCoreNumber==0)
	    Timer::tick(k);
	
	polymer = new PolymerMC(*polymerEtalon);
//	if(k>0)
//	    polymer->initWithRandomTaus();
	    
	if(regime==Regime::twoChains){
	    polymer2 = new PolymerMC(*polymerEtalon2);
	//    if(k>0)
	    Vector r0_2(4,0,0);
	    polymer2->initWithRandomTaus(r0_2);
	    
	}
	
	fakeCoreNumber = myCoreNumber+k*totalCoreNumber;
	
	if(fakeCoreNumber==0){
	    tfp = fopen("results/TemperatureMap.dat","w");
	    rgyrCurrentfp = fopen("results/RadiusOfGyrationCurrent.dat","w");
	    rgyrFinalfp = fopen("results/RadiusOfGyrationFinal.dat","w");
	    eCurrentfp = fopen("results/EnergyCurrent.dat","w");
	    eFinalfp = fopen("results/EnergyFinal.dat","w");
	}
	
	fname1 = new char[100];
	sprintf(fname1,"results/configurations/%iconf.dat",fakeCoreNumber);
	ktfp = fopen(fname1, "w");
	delete [] fname1;
	
	fname1 = new char[100];
	sprintf(fname1,"results/configurations/%iconf.pca",fakeCoreNumber);
	conffp = fopen(fname1, "w");
	delete [] fname1;
	
	polymer->writeRadiusVectorsInFile(conffp);
	if(regime==Regime::twoChains){
	    fname1 = new char[100];
	    sprintf(fname1,"results/configurations/%iconf2.dat",fakeCoreNumber);
	    ktfp2 = fopen(fname1, "w");
	    delete [] fname1;
	    
	    fname1 = new char[100];
	    sprintf(fname1,"results/configurations/%iconf2.pca",fakeCoreNumber);
	    conffp2 = fopen(fname1, "w");
	    delete [] fname1;
	
	    polymer2->writeRadiusVectorsInFile(conffp2);
	}

	for(double t=maxLogT; t>=minLogT; t-=logTstep){
		
	    if(_PCA_IS_EQUAL(t,0.0))
		t = 0.0;
		
	    temperature = pow(10,t);
	    
	    if(myCoreNumber==0){
		fprintf(logfp,"fakeCore: %i;\tlogT: %g\n", k, t);
		fflush(logfp);
	    }
	    
	    switch(regime){
	    /* Thermalization */
		case Regime::normal:
		    fprintf(rgyrCurrentfp,"\t\tlogT: %g\n",t);
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%5==0){
//			    fprintf(logfp,"%g\t%i\n",t, i);fflush(logfp);
//			    polymer->writeRadiusVectorsInFile(conffp);
			    fprintf(rgyrCurrentfp,"%.6le\n",PolymerObservable::radiusOfGyration(*polymer));
			}
			polymer->updateAllSites(temperature, *hamiltonian, *interaction);
		    }
		break;
		case Regime::withoutH:
		    fprintf(rgyrCurrentfp,"\t\tlogT: %g\n",t);
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%5==0)
			    fprintf(rgyrCurrentfp,"%.6le\n",PolymerObservable::radiusOfGyration(*polymer));
			
			polymer->updateAllSitesWithoutH(temperature, *interaction);
		    }
		break;
		case Regime::withSA:
		    fprintf(rgyrCurrentfp,"\t\tlogT: %g\n",t);
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%5==0)
			    fprintf(rgyrCurrentfp,"%.6le\n",PolymerObservable::radiusOfGyration(*polymer));
			polymer->updateAllSitesWithOnlySA(temperature, *hamiltonian);
		    }
		break;
		case Regime::twoChains:
		    fprintf(rgyrCurrentfp,"\t\tlogT: %g\n",t);
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%5==0)
			    fprintf(rgyrCurrentfp,"%.6le\n",PolymerObservable::radiusOfGyration(*polymer));
			polymer->updateAllSitesWithOnlySA(temperature, *hamiltonian);
			polymer->updateR02chains(temperature, *interaction, *polymer2);
			polymer->updateAllSites2chains(temperature, *hamiltonian, *interaction, *polymer2, minDist);
			polymer2->updateAllSites2chainsBW(temperature, *hamiltonian, *interaction, *polymer, minDist);
			
			polymer->updateAllSites2chainsBW(temperature, *hamiltonian, *interaction, *polymer2, minDist);
			polymer2->updateR02chains(temperature, *interaction, *polymer);
			polymer2->updateAllSites2chains(temperature, *hamiltonian, *interaction, *polymer, minDist);
//if(i%50==0)
//polymer->writeRadiusVectorsInFile(conffp);
			
			
			
		    }
		break;
		
		case Regime::backwards:
		    fprintf(rgyrCurrentfp,"\t\tlogT: %g\n",t);
		    for(int i=0; i<sweepsPerStep;i++){
			if(i%5==0)
			    fprintf(rgyrCurrentfp,"%.6le\n",PolymerObservable::radiusOfGyration(*polymer));
			polymer->updateAllSites(temperature, *hamiltonian, *interaction);
			polymer->updateAllSitesBW(temperature, *hamiltonian, *interaction);
			
		    }
		break;
		default:
		    _PCA_ERROR("MonteCarlo.run: unknown regime");
	    }
//polymer->printAcceptNumberR();
//polymer2->printAcceptNumberR();

	    polymer->writeKappaTauInFile(ktfp);
	    polymer->writeRadiusVectorsInFile(conffp);
	    
	    if(regime==Regime::twoChains){
		polymer2->writeKappaTauInFile(ktfp2);
		polymer2->writeRadiusVectorsInFile(conffp2);
//		if(!polymer2->selfAvoidingCondition(0, 3.8))
//	        printf("Finally 2 !SA\n");

//		else
//		printf("Finally 2 SA is Ok\n");
	    }
	    
//if(!polymer->selfAvoidingCondition(0, 3.8))
//printf("Finally !SA\n");

//else
//printf("Finally SA is Ok\n");
	    if(fakeCoreNumber==0){
		fprintf(tfp,"%g\t%.15le\n", t, temperature);
		fflush(tfp);
		
		fprintf(rgyrFinalfp, "%.6le\n", PolymerObservable::radiusOfGyration(*polymer));
		fflush(rgyrFinalfp);
		
		if(regime==Regime::twoChains){
		    const Vector* r1=polymer->getRadiusVectors();
		    double e1 = hamiltonian->energyAllSites(*polymer);
		    fprintf(eFinalfp, "%.15le\t", e1);
		    fflush(eFinalfp);
		
		    const Vector* r2=polymer2->getRadiusVectors();
		    double e2 = hamiltonian->energyAllSites(*polymer2);
		    
		    fprintf(eFinalfp, "%.15le\t", e2);
		    fflush(eFinalfp);
		    
/*		    for(int c=0;c<polymer->getNumMonomers()+1;c++){
			for(int m=0;m<polymer2->getNumMonomers()+1;m++){
			    r1[c].print();
			    r2[m].print();
			    printf("%i  %i  %.15le\n\n",c, m,(r1[c]-r2[m]).norm());
			}
		    }
*/
		    
		    //for Interaction
		    int site = 0;
		    int N12 = polymer->getNumMonomers() + polymer2->getNumMonomers() + 2 - site; // sum number of r-vectors in both chains + trash element
		    Vector* r12 = new Vector [N12];
	
		    //write the whole second chain:
		    for(int c=0; c<polymer2->getNumMonomers() + 1; c++)
			r12[c] = r2[c];

		    r12[polymer2->getNumMonomers()+1] = Vector::zero; // any number, never used;

		    //add this chain:
		    for(int c=polymer2->getNumMonomers()+2; c<N12; c++){
			r12[c] = r1[site+1+c-polymer2->getNumMonomers()-2];
			//r12[c].print();
		    }
		    double lj = interaction->energyIfSiteChanged(polymer2->getNumMonomers()+1, N12, r12);
		
		    fprintf(eFinalfp, "%.15le\t", lj);
		    fprintf(eFinalfp, "%.15le\t", e1+e2+lj);
		    
		    fflush(eFinalfp);
		    fprintf(eFinalfp, "\n");
		    
		    delete [] r12;
		
		}
		
	    }
	    
	}
	
	
	if(fakeCoreNumber==0)
	    polymer->printAcceptNumberR(logfp);
//printf("here 0\n"); fflush(stdout);	
        if(polymer2!=NULL)
    	    if(fakeCoreNumber==0)
		polymer2->printAcceptNumberR(logfp);
//printf("here 1\n"); fflush(stdout);
	if(myCoreNumber==0)
	    Timer::tock(k, "",logfp);
//printf("here 2\n");fflush(stdout);
	if(fakeCoreNumber==0){
	    fclose(tfp);
	    fclose(rgyrCurrentfp);
	    fclose(rgyrFinalfp);
	    fclose(eCurrentfp);
	    fclose(eFinalfp);
	}
//printf("here 3\n");fflush(stdout);
	delete polymer;
	fclose(ktfp);
	fclose(conffp);
//printf("here 4\n");fflush(stdout);	
	if(regime==Regime::twoChains){
	    delete polymer2;
	    fclose(ktfp2);
	    fclose(conffp2);
	}
    }
    
    if(myCoreNumber==0){
	fprintf(logfp,"END\n\n");
	fclose(logfp);
//	printf("done\n");
    }

    
}


void MonteCarlo::readFromParamFile(const char* fileName)
{
    std::string etalon;
    int number;
    
    ParamFileReader data(fileName);
    
    etalon = "REGIME";
    number = data.search(etalon);
    if(number>=0)
	regime = (MonteCarlo::Regime)data.value(number);
    else
	regime = MonteCarlo::Regime::normal;
    
    etalon = "LOOPS_PER_CORE";
    number = data.search(etalon);
    if(number>=0)
	loopsPerCore = (int)data.value(number);
    else
	loopsPerCore = 1;
    
    etalon = "MAX_LOG_T";
    number = data.search(etalon);
    if(number>=0)
	maxLogT = data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MIN_LOG_T";
    number = data.search(etalon);
    if(number>=0)
	minLogT = data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "MIN_LOG_T";
    number = data.search(etalon);
    if(number>=0)
	minLogT = data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "LOG_T_STEP";
    number = data.search(etalon);
    if(number>=0)
	logTstep = data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "SWEEPS_PER_STEP";
    number = data.search(etalon);
    if(number>=0)
	sweepsPerStep = (int)data.value(number);
	
    else{
	printf("Cannot find %s\n", etalon.c_str());
	exit(1);
    }
    
    etalon = "CORES";
    number = data.search(etalon);
    if(number>=0)
	cores = (int)data.value(number);
	
    else{
	cores = 0;
    }
    
}

}//end of namespace