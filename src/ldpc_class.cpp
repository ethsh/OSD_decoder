// Copyright 2012 Xishuo Liu, Stark C. Draper, Benjamin Recht
//
// This program is distributed under the terms of the GNU General Public License.
//
// This file is part of ADMM Decoder.
//
//    ADMM Decoder is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ADMM Decoder is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ADMM Decoder.  If not, see <http://www.gnu.org/licenses/>.
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Name:				ADMM Decoder
// Files:				ldpc_simulator_example.cpp, ldpc_class.cpp,
//						ldpc_simulator_data_def.h, ldpc_class.h,
//						MersenneTwister.h
// Date:				8.30.2012
//
// Author:				Xishuo Liu, xliu94@wisc.edu
// Thanks to:			S. Barman, S. Draper and B. Recht.
//
// Papers:				1. S. Barman, X. Liu, S. Draper and B. Recht,  
//						"Decomposition Methods for Large Scale LP Decoding"
//						http://arxiv.org/abs/1204.0556
//						2. X. Liu, S. Draper and B. Recht,
//						"Suppressing Pseudocodewords by Penalizing the Objective of LP Decoding"
//						IEEE Information Theory Workshop (ITW), 2012. Lausanne: Switzerland, 2012
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
This file contains methods for the ADMM decoder classes. 
*/

#include "ldpc_class.h"
#include "timer.h"
#include <cstdint>
#include <algorithm>
using std::min;
using std::max;

// #include "Projection_algorithm.h"
//#include "configuration_parameters.h"
//#define __OVER_RELAX__
double G_rho = 1.0;

// Martzi: We don't use it for mRRD simulation
// Snippet 1 (from ldpc_class.cpp.bak)


// Decoder Class
Decoder::Decoder(string FileName,int nChecks, int BlockLength)
{
	// Learn and read parity check matrix
	int extraChecks = 30;
	mPCMatrixFileName = FileName;
	mNChecks = nChecks;
	mBlocklength = BlockLength;
	mPCheckMapSize = 0;
	mPA_Alex = 0;
	CheckDegree = new int[mNChecks + extraChecks];
	VariableDegree = new int[mBlocklength];
	OutputFromDecoder = new double[mBlocklength];
	_LogLikelihoodRatio = new double [mBlocklength];
	u0 = new double [mBlocklength];
	mF1V = new int[nChecks+1];

	mSetDefaultParameters();
	// learn the structure of parity check matrix: degree, # edges
	mLearnParityCheckMatrix(); 
	mPCheckMap = new TRIPLE[mPCheckMapSize + 500];
	mParityCheckMatrix = new TRIPLE[mPCheckMapSize + 500];
	u = new MESSAGE<double>[mPCheckMapSize + 500];
	v = new MESSAGE<double>[mPCheckMapSize + 500];
	mReadParityCheckMatrix();
	warmstartchecks = new double[mPCheckMapSize + 500];
	warmstartsave = new double[mPCheckMapSize + 500];

	RateOfCode = (double)(mBlocklength -mNChecks)/mBlocklength;
}
Decoder::~Decoder()
{
	delete [] mF1V;
	delete [] u;
	delete [] v;
	delete [] u0;
	delete [] CheckDegree;
	delete [] VariableDegree;
	delete [] warmstartsave;
	delete [] warmstartchecks;
	delete [] mPCheckMap;
	delete [] mParityCheckMatrix;
	delete [] OutputFromDecoder;
	delete [] _LogLikelihoodRatio;
	//
}

//Decodersettings
void Decoder::SetParameters(int mxIt, double feas,  double p_mu, double p_rho)
{
	maxIteration = mxIt;
	para_end_feas = feas;
	para_mu = p_mu;
	para_rho = p_rho;
}
void Decoder::SetParameters(int mxIt, double dmp_coef_alpha)
{
	maxIteration = mxIt;
	para_dmp_coef_alpha = dmp_coef_alpha;
}
void Decoder::SetBSCnormalized(bool normal)
{
	if(normal)
		normalizeBSC = 1;
	else
		normalizeBSC = 0;
}
void Decoder::SetBPSaturation(bool nonsat, double value)
{
	BP_non_sat = nonsat;
	BP_sat_value = value;
}
void Decoder::SetDecoder(string decoder)
{
	mRealDecoderName = decoder;
	if(decoder == "BP")
	{
		DecoderName = decoder;
		mInitializeMessages();
	}
	else if(decoder == "ADMMLP")
	{
		DecoderName = "ADMM";
		enumADMM = LinearProgram;
	}
	else if(decoder == "ADMML1")
	{
		DecoderName = "ADMM";
		enumADMM = L1Penalty;
	}
	else if(decoder == "ADMML2")
	{
		DecoderName = "ADMM";
		enumADMM = L2Penalty;
	}
	else if(decoder == "ADMMEntropy")
	{
		DecoderName = "ADMM";
		enumADMM = EntropyPenalty;
	}
	else if(decoder == "MSA")
	{
		DecoderName = "MSA";
	}
	else
	{
		cout<<"decoder not supported"<<endl;
		exit(0);
	}
}

// Decoder initializations
void Decoder::mSetDefaultParameters()
{
	// decoder settings
	maxIteration = 1000;
	para_end_feas = 1e-6;
	para_mu = 5;
	para_rho= 1.8;

	alpha = 1;
	for(int i = 0; i < mBlocklength; i++)
	{
		VariableDegree[i] = 0;
	}
	for(int i = 0; i < mNChecks; i++)
	{
		CheckDegree[i] = 0;
		mF1V[i] = 0;
	}
	normalizeBSC = 0.0;
}
void Decoder::mLearnParityCheckMatrix()// learn the degree of the parity check matrix
{
	ifstream myfile(mPCMatrixFileName.c_str());
	int tempvalue = 0;
	int i = 0;
	string line;
	if (myfile.is_open())
	{
		mF1V[0] = 0;
		//cout<<"here 2.1"<<endl;
		while(getline(myfile, line)) 
		{
		   istringstream is(line);
		   int curr_degree = 0;
		   while( is >> tempvalue ) 
		   {
			   
			   VariableDegree[tempvalue]++;
			   curr_degree++;
			   mPCheckMapSize++;
		   }
		   CheckDegree[i] = curr_degree;
		   //cout<<"here 2.11-"<<i<<endl;
		   i++;
		   mF1V[i] = mF1V[i-1] + curr_degree;
		}
		//cout<<"here 2.2"<<endl;
		if (myfile.is_open())
			myfile.close();
		//cout<<"here 2.3"<<endl;

	}
	else 
	{
		cout << "Unable to open file"; 
		exit(0);
	}
}
void Decoder::mReadParityCheckMatrix()// read the parity check matrix, initialize message passing
{
	ifstream myfile (mPCMatrixFileName.c_str());
	int tempvalue = 0;
	int i = 0;
	int count = 0;
	string line;
	//cout<<"here 3.1"<<endl;
	if (myfile.is_open())
	{
		while(getline(myfile, line)) 
		{
			istringstream is(line);
			while( is >> tempvalue ) 
			{
				mParityCheckMatrix[count].col = tempvalue;
				mParityCheckMatrix[count].row = i;

				mPCheckMap[count].col = tempvalue;
				mPCheckMap[count].row = count;
			
				u[count].SetDeg(CheckDegree[i]);

				u[count].Mcol = mParityCheckMatrix[count].col;
				u[count].Mrow = mParityCheckMatrix[count].row;
				u[count].Mprob = 0;

				v[count].SetDeg(VariableDegree[tempvalue]);
				v[count].Mcol = mParityCheckMatrix[count].col;
				v[count].Mrow = mParityCheckMatrix[count].row;
				v[count].Mprob = 0;
			
				count++;
			}
			i++;
		}
	}
	if (myfile.is_open())
		myfile.close();
}
void Decoder::mInitializeMessages() 
{
	for (int i = 0; i < mPCheckMapSize; i++)
	{
		int countu = 0, countv = 0;
		for (int j = 0; j < mPCheckMapSize; j++)
		{
			if (u[i].Mrow == u[j].Mrow && i!=j)
			{
				u[i].i[countu] = j;
				countu++;
			}
			if (v[i].Mcol == v[j].Mcol && i!=j)
			{
				v[i].i[countv] = j;
				countv++;
			}
			if(countu == u[i].degree - 1 && countv == v[i].degree - 1)
				break;
		}
		countu = 0; countv = 0;
			//cout<<(double)i/mPCheckMapSize *100<<endl;
	}
}

// Martzi: We don't use it for mRRD simulation
// Snippet 2 (from ldpc_class.cpp.bak)
double Decoder::NonSatUpd(double llr1, double llr2)	
{
	double pos_llr1, pos_llr2;
	int sign_llr1, sign_llr2;
	if(llr1 > 0)
	{
		pos_llr1 = llr1;
		sign_llr1 = 1;
	}
	else
	{
		pos_llr1 = -llr1;
		sign_llr1 = -1;
	}
	if(llr2 > 0)
	{
		pos_llr2 = llr2;
		sign_llr2 = 1;
	}
	else
	{
		pos_llr2 = -llr2;
		sign_llr2 = -1;
	}
	double llrout = sign_llr1 * sign_llr2 * min(pos_llr1, pos_llr2);
	double sumsum = llr1 + llr2;
	double minusminus = llr1 - llr2;
	llrout = llrout + log(1 + exp(- abs(sumsum))) - log(1 + exp(- abs(minusminus)));
	return llrout;
}


void Decoder::BPDecoder()
{
	double* uu = new double[mPCheckMapSize];
	//clock_t t1, t2, t3, t4;
	// uint64_t t1, t2, t3, t4, tCheck;
	double tVar = 0; double tCheck = 0;
	double ttVar = 0;
	double ttCheck = 0;

	bool error = false;
	//double feas_tol = para_end_feas;
	int maxIter = maxIteration;
	int numIter = 0;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ifstream myfile("_LLR.txt");
	double tempvalue = 0;
	int i = 0;
	string line;
	for (int i = 0; i < mBlocklength; i++)
	{
		_LogLikelihoodRatio[i] = 0;
	}
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			istringstream is(line);
			while (is >> tempvalue)
			{
				_LogLikelihoodRatio[i] = tempvalue;
			}
			i++;
		}
	}
	if (myfile.is_open())
		myfile.close();
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for (int i = 0; i < mPCheckMapSize; i++)
	{
		v[i].Mprob =  _LogLikelihoodRatio[v[i].Mcol];
	}
	for(int i = 0; i < mBlocklength; i++)
	{
		u0[i] = _LogLikelihoodRatio[i];
	}

	int itnum;
	for(itnum = 0; itnum < maxIter; itnum++)
	{
		mIteration = itnum + 1;
		//	cout<<"finished calculating u"<<endl;
		//t1 = clock();
		// t1 = __rdtsc();
		// Update Parity-Check Nodes
		for (int i = 0; i < mPCheckMapSize; i++)
		{
			double temp;
			if(BP_non_sat)
			{
				temp = 1.0;
				temp = v[u[i].i[0]].Mprob;
				for(int j = 1; j < u[i].degree - 1; j++)
				{
					temp = NonSatUpd(temp, v[u[i].i[j]].Mprob);
				}
				if(u[i].degree % 2 != 0)
					temp = temp;
				u[i].Mprob = temp;
			}
			else
			{
				temp = 1.0;
				for(int j = 0; j < u[i].degree - 1; j++)
				{
					if(abs(v[u[i].i[j]].Mprob) < BP_sat_value)
						temp = temp * tanh(0.5 * v[u[i].i[j]].Mprob);
					else
					{
						if(v[u[i].i[j]].Mprob > 0)
							temp = temp * tanh(0.5 * BP_sat_value);
						else
							temp = temp * tanh(- 0.5 * BP_sat_value);
					}
				}
				u[i].Mprob = log((1 + temp)/(1 - temp));
			}
			if(my_isinf(u[i].Mprob) || my_isnan(u[i].Mprob))
			{
				cout<<"is inf or is nan"<<endl;
			}
		}
//		t2 = clock();
//		tCheck += (double) (t2 - t1);
		// t2 = __rdtsc();
		// tCheck = (t2 - t1); //(double) (t2 - t1);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		for (int i = 0; i < mPCheckMapSize; i++)
		{
			uu[i] = u[i].Mprob;
		}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//cout<<"one iter"<<endl;
		// Update Variable Nodes
		for (int i=0;i<mPCheckMapSize;i++)
		{

				v[i].Mprob = u0[v[i].Mcol];
				for(int j = 0; j < v[i].degree - 1; j++)
					v[i].Mprob += u[v[i].i[j]].Mprob;			
		}
		// t3 = clock();
		// tVar += (double)(t3 - t2);
		// decision at each iteration
		for(int i = 0; i < mBlocklength; i++)
			OutputFromDecoder[i] = u0[i];
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			OutputFromDecoder[u[i].Mcol] += u[i].Mprob;
		}
		for(int i = 0; i < mBlocklength; i++)
		{
			if(my_isinf(OutputFromDecoder[i]))
			{
				cout<<"isinf"<<endl;
				break;
			}
			if(my_isnan(OutputFromDecoder[i]))
			{
				cout<<"isnan"<<endl;
				break;
			}
		}

		if (ValidateCodeword())
		{
			mValidCodeword = true;
			mAlgorithmConverge = true;
			break;
		}
	}
	ttVar = (double) (tVar / CLOCKS_PER_SEC);
	ttCheck = (double) (tCheck / CLOCKS_PER_SEC);	
	delete [] uu;
}

void Decoder::Eff_BPDecoder()
{
	//clock_t t1, t2, t3, t4;
	// uint64_t t1, t2, t3, t4, tCheck;
	double tVar = 0; double tCheck = 0;
	double ttVar = 0;
	double ttCheck = 0;
	int f1 = 0;
	bool error = false;
	//double feas_tol = para_end_feas;
	int maxIter = maxIteration;
	int numIter = 0;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	double *uu = new double[mPCheckMapSize];
	ifstream myfile("_LLR.txt");
	double tempvalue = 0;
	int i = 0;
	string line;
	for (int i = 0; i < mBlocklength; i++)
	{
		_LogLikelihoodRatio[i] = 0;
	}
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			istringstream is(line);
			while (is >> tempvalue)
			{
				_LogLikelihoodRatio[i] = tempvalue;
			}
			i++;
		}
	}
	if (myfile.is_open())
		myfile.close();

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	for (int i = 0; i < mPCheckMapSize; i++)
	{
		v[i].Mprob = _LogLikelihoodRatio[v[i].Mcol];
	}
	for (int i = 0; i < mBlocklength; i++)
	{
		u0[i] = _LogLikelihoodRatio[i];
	}

	int itnum;
	for (itnum = 0; itnum < maxIter; itnum++)
	{
		mIteration = itnum + 1;
		//	cout<<"finished calculating u"<<endl;
		// t1 = __rdtsc();
		///////////////////////////////
		// Update Parity-Check Nodes //
		///////////////////////////////
		
		// go over all Parity Check nodes (Rows of H)
		for (int j = 0; j < mNChecks; j++)
		{
			f1 = mF1V[j];
			double temp;
			// calc msg left on the edge coresponds to the pivot '1' in row j of H 
			if (BP_non_sat)
			{
				temp = 1.0;
				temp = tanh(0.5 * v[u[f1].i[0]].Mprob);
				for (int k = 1; k < (u[f1].degree - 1); k++)
				{
					temp = NonSatUpd(temp, v[u[f1].i[k]].Mprob);
				}
				u[f1].Mprob = temp;
			}
			else
			{
				temp = 1.0;  double raw;
				temp = tanh(0.5 * v[u[f1].i[0]].Mprob);
				for (int k = 1; k < (u[f1].degree - 1); k++)
				{
					raw = 
					temp = temp * tanh( clamp(0.5 * v[u[f1].i[k]].Mprob, -15.0 , 15.0) );
				}
				u[f1].Mprob = log((1 + temp) / (1 - temp)); // For real x<1 arctanh(x) = log((1 + x) / (1 - x))
			}
			if (my_isinf(u[f1].Mprob) || my_isnan(u[f1].Mprob))
			{
				cout << "is inf or is nan" << endl;
			}
			// calculate msgs left for the rest of the edges connected to Parity Check j
			if (BP_non_sat)
			{
				cout << "Implement NonSatUpd^(-1) !!!" << endl;
				system("pause");
			}
			else
			{
				
				temp = temp * tanh( clamp(0.5 * v[f1].Mprob, -15.0 , 15.0) );
				temp = temp / tanh( clamp(0.5 * v[u[f1].i[0]].Mprob, -15.0, 15.0));   // calc msg left for the "second" '1' of j
				u[u[f1].i[0]].Mprob = log((1 + temp) / (1 - temp));;											  // (the second '1' from left in the j-th row of H)
				if (my_isinf(u[u[f1].i[0]].Mprob) || my_isnan(u[u[f1].i[0]].Mprob))
				{
					cout << "is inf or is nan" << endl;
				}
				for (int k = 1; k < (u[f1].degree - 1); k++){
					temp = temp * tanh(clamp(0.5 * v[u[f1].i[k - 1]].Mprob, -15.0, 15.0));
					temp = temp / tanh(clamp(0.5 * v[u[f1].i[k]].Mprob, -15.0, 15.0));      // calc msg left for the "k-th" '1' of j 
					u[u[f1].i[k]].Mprob = log((1 + temp) / (1 - temp));;												// (the "k-th" '1' from left in the j-th row of H)
					if (my_isinf(u[u[f1].i[k]].Mprob) || my_isnan(u[u[f1].i[k]].Mprob))
					{
						cout << "is inf or is nan" << endl;
					}
				}

			}
		}// end of outer loop - running over the rows of H (Parity Check Nodes)
		// t2 = __rdtsc();
		// tCheck = (t2 - t1); //(double) (t2 - t1);
		
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		for (int i = 0; i < mPCheckMapSize; i++)
		{
			uu[i] = u[i].Mprob;
		}
		
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		

		



		//cout<<"one iter"<<endl;
		// Update Variable Nodes
		for (int i = 0; i<mPCheckMapSize; i++)
		{

			v[i].Mprob = u0[v[i].Mcol];
			for (int j = 0; j < v[i].degree - 1; j++)
				v[i].Mprob += u[v[i].i[j]].Mprob;
		}
		// t3 = clock();
		// tVar += (double) (t3 - t2);
		// decision at each iteration
		for (int i = 0; i < mBlocklength; i++)
			OutputFromDecoder[i] = u0[i];
		for (int i = 0; i < mPCheckMapSize; i++)
		{
			OutputFromDecoder[u[i].Mcol] += u[i].Mprob;
		}
		for (int i = 0; i < mBlocklength; i++)
		{
			if (my_isinf(OutputFromDecoder[i]))
			{
				cout << "isinf" << endl;
				break;
			}
			if (my_isnan(OutputFromDecoder[i]))
			{
				cout << "isnan" << endl;
				break;
			}
		}
		
		if (ValidateCodeword())
		{
			mValidCodeword = true;
			mAlgorithmConverge = true;
			break;
		}
	}
	

	ttVar = (double) (tVar / CLOCKS_PER_SEC);
	ttCheck = (double) (tCheck / CLOCKS_PER_SEC);
	delete[] uu;
}

void Decoder::Eff_TG_Decoder_a()
{
	//clock_t t1, t2, t3, t4;
	// uint64_t t1, t2, t3, t4, tCheck;
	double tVar = 0; double tCheck = 0;
	double ttVar = 0;
	double ttCheck = 0;
	int f1 = 0;
	bool error = false;
	//double feas_tol = para_end_feas;
	int maxIter = maxIteration;
	int numIter = 0;
	Timer t;
	double *v_prev_Mprob = new double[mPCheckMapSize];
/*
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double *uu = new double[mPCheckMapSize];
ifstream myfile("_LLR.txt");
double tempvalue = 0;
int i = 0;
string line;
for (int i = 0; i < mBlocklength; i++)
{
_LogLikelihoodRatio[i] = 0;
}
if (myfile.is_open())
{
while (getline(myfile, line))
{
istringstream is(line);
while (is >> tempvalue)
{
_LogLikelihoodRatio[i] = tempvalue;
}
i++;
}
}
if (myfile.is_open())
myfile.close();

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

	t.start();
	for (int i = 0; i < mPCheckMapSize; i++)
	{
		v[i].Mprob = _LogLikelihoodRatio[v[i].Mcol];
		u[i].Mprob = 0; 
		v_prev_Mprob[i] = 0;
	}
	for (int i = 0; i < mBlocklength; i++)
	{
		u0[i] = _LogLikelihoodRatio[i];
	}

	int itnum;
	for (itnum = 0; itnum < maxIter; itnum++)
	{
		mIteration = itnum + 1;
//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
		// UPDATE VARIABLE NODES VV											   VV
		//VVVVVVVVVVVVVVVVVVVVVVVVV                                            VV
		for (int i = 0; i<mPCheckMapSize; i++)                                 
		{

			v[i].Mprob = u0[v[i].Mcol];
			for (int j = 0; j < v[i].degree - 1; j++)
				v[i].Mprob += u[v[i].i[j]].Mprob;
			v[i].Mprob = v[i].Mprob*para_dmp_coef_alpha + (1.0 - para_dmp_coef_alpha)*v_prev_Mprob[i];
			v_prev_Mprob[i] = v[i].Mprob;
		}
//		t3 = clock();
//		tVar += (double) (t3 - t2);
//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV





		//	cout<<"finished calculating u"<<endl;
		
//PCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPC
		// Update Parity-Check Nodes //
		//CPCPCPCPCPCPCPCPCPCPCPCPCPCP
		// t1 = __rdtsc();
		// go over all Parity Check nodes (Rows of H)
		for (int j = 0; j < mNChecks; j++)
		{
			f1 = mF1V[j];
			double temp;
			// calc msg left on the edge coresponds to the pivot '1' in row j of H 
			if (BP_non_sat)
			{
				temp = 1.0;
				temp = tanh(0.5 * v[u[f1].i[0]].Mprob);
				for (int k = 1; k < (u[f1].degree - 1); k++)
				{
					temp = NonSatUpd(temp, v[u[f1].i[k]].Mprob);
				}
				u[f1].Mprob = temp;
			}
			else
			{
				temp = 1.0;  double raw;
				temp = tanh(0.5 * v[u[f1].i[0]].Mprob);
				for (int k = 1; k < (u[f1].degree - 1); k++)
				{
					raw =
						temp = temp * tanh(clamp(0.5 * v[u[f1].i[k]].Mprob, -15.0, 15.0));
				}
				u[f1].Mprob = log((1 + temp) / (1 - temp)); // For real x<1 arctanh(x) = log((1 + x) / (1 - x))
			}
			if (my_isinf(u[f1].Mprob) || my_isnan(u[f1].Mprob))
			{
				cout << "is inf or is nan" << endl;
			}
			// calculate msgs left for the rest of the edges connected to Parity Check j
			if (BP_non_sat)
			{
				cout << "Implement NonSatUpd^(-1) !!!" << endl;
				system("pause");
			}
			else
			{

				temp = temp * tanh(clamp(0.5 * v[f1].Mprob, -15.0, 15.0));
				temp = temp / tanh(clamp(0.5 * v[u[f1].i[0]].Mprob, -15.0, 15.0));   // calc msg left for the "second" '1' of j
				u[u[f1].i[0]].Mprob = log((1 + temp) / (1 - temp));;											  // (the second '1' from left in the j-th row of H)
				if (my_isinf(u[u[f1].i[0]].Mprob) || my_isnan(u[u[f1].i[0]].Mprob))
				{
					cout << "is inf or is nan" << endl;
				}
				for (int k = 1; k < (u[f1].degree - 1); k++){
					temp = temp * tanh(clamp(0.5 * v[u[f1].i[k - 1]].Mprob, -15.0, 15.0));
					temp = temp / tanh(clamp(0.5 * v[u[f1].i[k]].Mprob, -15.0, 15.0));      // calc msg left for the "k-th" '1' of j 
					u[u[f1].i[k]].Mprob = log((1 + temp) / (1 - temp));;												// (the "k-th" '1' from left in the j-th row of H)
					if (my_isinf(u[u[f1].i[k]].Mprob) || my_isnan(u[u[f1].i[k]].Mprob))
					{
						cout << "is inf or is nan" << endl;
					}
				}

			}
		}// end of outer loop - running over the rows of H (Parity Check Nodes)
		// t2 = __rdtsc();
		// tCheck = (t2 - t1); //(double) (t2 - t1);

/*
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for (int i = 0; i < mPCheckMapSize; i++)
{
uu[i] = u[i].Mprob;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
//PCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCPCP


//TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
		// TERMINATION //
		//TTTTTTTTTTTTTTT

		for (int i = 0; i < mBlocklength; i++)
			OutputFromDecoder[i] = u0[i];
		for (int i = 0; i < mPCheckMapSize; i++)
		{
			OutputFromDecoder[u[i].Mcol] += u[i].Mprob;
		}
		for (int i = 0; i < mBlocklength; i++)
		{
			if (my_isinf(OutputFromDecoder[i]))
			{
				cout << "isinf" << endl;
				break;
			}
			if (my_isnan(OutputFromDecoder[i]))
			{
				cout << "isnan" << endl;
				break;
			}
		}
	}

	t.stop();
	mTGDec_a_DecodeTime = t.getElapsedTimeInMicroSec();
	ttVar = (double) (tVar / CLOCKS_PER_SEC);
	ttCheck = (double) (tCheck / CLOCKS_PER_SEC);
	//delete[] uu;
	delete [] v_prev_Mprob;
}


// Martzi: We don't use it for mRRD simulation
// Snippet 3 (from ldpc_class.cpp.bak)

bool Decoder::ValidateCodeword()
{
	int *syndrome = new int[mNChecks];
	for(int i = 0; i < mNChecks; i++)
		syndrome[i] = 0;
	for(int i = 0; i < mPCheckMapSize; i++)
	{
		if(DecoderName == "ADMM")
		{
			syndrome[mParityCheckMatrix[i].row] += (int) floor(OutputFromDecoder[mParityCheckMatrix[i].col] + 0.5);
		}
		else // DecoderName == "BP" or "MSA"
		{
			if(OutputFromDecoder[mParityCheckMatrix[i].col] < 0)
				syndrome[mParityCheckMatrix[i].row]++;
		}
	}
	bool validcodeword = true;
	for(int i = 0; i < mNChecks; i++)
	{
		if (syndrome[i]%2 == 1)
		{
			validcodeword = false;
			break;
		}
	}
	delete [] syndrome;
	return validcodeword;
}
bool Decoder::ValidateCodeword(int input[])
{
	int *syndrome = new int[mNChecks];
	for(int i = 0; i < mNChecks; i++)
		syndrome[i] = 0;
	for(int i = 0; i < mPCheckMapSize; i++)
	{
		syndrome[mParityCheckMatrix[i].row] += input[mParityCheckMatrix[i].col];
	}
	bool validcodeword = true;
	for(int i = 0; i < mNChecks; i++)
	{
		if (syndrome[i]%2 == 1)
		{
			validcodeword = false;
			break;
		}
	}
	delete [] syndrome;
	return validcodeword;
}
void Decoder::UPDATE_L_CODE(int *row, char *CodeFileNameLiu, int m)
{
	/*
	row[0] = 0; 
	row[1] = 1; 
	row[2] = 1;
	row[3] = 0;
	row[4] = 1;
	row[5] = 1;
	row[6] = 1;
	
	const int new_m = (int) mNChecks + 1;
	int* new_CheckDegree = new int[new_m];
	for (int k = 0; k < new_m; k++)
	{
	new_CheckDegree[k] = 7;
	}
	CheckDegree = new_CheckDegree;
	*/
	int new_mPCheckMapSize = mPCheckMapSize;
	CheckDegree[m] = 0;
	for (int k = 0; k < mBlocklength; k++){
		if (row[k] == 1){
			
			mParityCheckMatrix[new_mPCheckMapSize].col = k;
			mParityCheckMatrix[new_mPCheckMapSize].row = m;

			mPCheckMap[new_mPCheckMapSize].col = k;
			mPCheckMap[new_mPCheckMapSize].row = new_mPCheckMapSize;

			VariableDegree[k]++;
			CheckDegree[m]++;

			new_mPCheckMapSize++;
		}
	}

	for (int count = 0; count < new_mPCheckMapSize; count++){
		u[count].SetDeg(CheckDegree[mParityCheckMatrix[count].row]);
		u[count].Mcol = mParityCheckMatrix[count].col;
		u[count].Mrow = mParityCheckMatrix[count].row;
		u[count].Mprob = 0;

		v[count].SetDeg(VariableDegree[mParityCheckMatrix[count].col]);
		v[count].Mcol = mParityCheckMatrix[count].col;
		v[count].Mrow = mParityCheckMatrix[count].row;
		v[count].Mprob = 0;
	}
	mPCheckMapSize = new_mPCheckMapSize;
	

	// try to update Liu's Code Data Stracture directly from row
	// or
	// 1. L_row <- convert row to a line in the format of CodeFileNameLiu
	// 2. write L_row in the end of CodeFileNameLiu
	// 3. call mLearnParityCheckMatrix and mReadParityCheckMatrix  
	return;
}

// Martzi: We don't use it for mRRD simulation
// Snippet 4 (from ldpc_class.cpp.bak)
