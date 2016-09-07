/*                                                                            */
/*                           ----  SPARTA   ----                              */
/*     Shifts Prediction from Analogue of Residue type and Torsion Angle      */
/*                           Yang Shen and Ad Bax                             */
/*                    J. Biomol. NMR, 38, 289-302 (2007)                      */
/*                 NIH, NIDDK, Laboratory of Chemical Physics                 */
/*                     version, 1.00 (build 2010.0607.00)                     */
/*                                                                            */
/*                      for any problem, please contact                       */
/*                          shenyang@niddk.nih.gov                            */
/*                                                                            */
/******************************************************************************/

#include "SPARTA.h"
#include <ctime>
#include <sys/timeb.h>


int main( int argc, char ** argv )
{
	string argList;
	bool time = true;
	SPARTA sparta;

	if(argc == 1) {
		sparta.printSyntax();
	}

	for(int i = 1; i < argc; i++){
		argList+= (" ");
		argList+= (argv[i]);
	}

	vector<string> fields = GDB::split(" -", argList);


	map<string, string> argments;
	for(int i = 0; i < fields.size(); i++)
	{
		vector<string> temp;
		int pos = fields[i].find_first_of(' ');
		temp.push_back( fields[i].substr(0,pos) );
		temp.push_back( fields[i].substr(pos+1,fields[i].length()-pos-1 ));
		string arg = PDB::simplifyWhiteSpace(temp[1]);

		vector<string>::iterator it = find(sparta.argList.begin(), sparta.argList.end(), temp[0]);
		if( it == sparta.argList.end() ) //not an valid argument
		{
			cerr << "Invaid arg name -" << temp[0] << endl;
			exit(0);
		}

		if( temp[0] == "oldsparta" ) arg = "true";
		else if( temp[0] == "notime" ) time = false;

		argments[ temp[0] ] = arg;

	}

	clock_t start, finish;
	start = clock();

	sparta.setArgs(argments);

	if( argments["oldsparta"].length() > 0 )
		sparta.runPredict();
	else if( argments["ins"].length() > 0 )
		sparta.runANN_Predictions();
	else
		sparta.runANN_Prediction();

	finish = clock();
	if( time )
    		cout << "\tRunning time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;
}
