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
// Project:				ADMM Decoder
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
This file contains declarations for the ADMM decoder classes. 
*/
#include "ldpc_simulator_data_def.h"
//#include <algorithm>  
#define ABS(X) (((X) < (0)) ? (-X) : (X))
#define clamp(x, min_val, max_val) (min(max((x), (min_val)), (max_val)))

// #include <intrin.h>
// #pragma intrinsic(__rdtsc)

// Martzi: We don't use it for mRRD simulation
//!  Class for noisy sequence. 
/*!
  Could be replaced by just an array ... anyway.
*/
class NoisySeq{
public:
	int Blocklength; /*!< blocklength of code */
	double *NoisySequence; /*!< noisy sequence vector, can be R^n for AWGN or {0,1}^n for BSC */
	string ChannelName; /*!< channel type: AWGN or BSC */
	double ChannelParameter;
	//! Print function. Intended for debugging.
	void print()
	{
		for(int i = 0; i < Blocklength; i++)
		{
			cout<<NoisySequence[i]<<" ";
		}
		cout<<endl;
	}

	//! Constructor.
    /*!
      \param blocklength Blocklength of code
	  \param channel Channel name, can only be "AWGN" or "BSC"
    */
	NoisySeq(int blocklength);
	void SetProperties(string channel, double channel_para){ChannelName = channel; ChannelParameter = channel_para;}
	//! destructor.
	~NoisySeq();
};


//!  Class for channel.
/*!
  Takes input and generate output. There are two types of channel supported: AWGN and BSC.
  The random number generator uses the Mersenne Twister random number generator.
*/
class Channel{
private:

	int mBlocklength;  /*!< blocklength of code */
	string mChannelName; /*!< channel type: AWGN or BSC */

	MTRand channelrand; /*!< random number generator*/
	double mCrossoverProb;/*!< cross over probability, is used only when channel is BSC*/
	double mSTD;/*!< std, is used only when channel is AWGN */

	void mGenerateAWGN(int input[], NoisySeq &ChannelOutput); // AWGN by std
	void mGenerateBitFlip(int input[], NoisySeq &ChannelOutput); // bit flip channel by num of bit flips
public:
	//! Set Channel
	/*!
      \param blocklength Blocklength of code
	  \param channel channel name, can only be "AWGN" or "BSC"
	  \param parameter Parameter for the channel. std for AWGN or crossover prob for BSC
    */
	void SetChannel(int blocklength, string channel, double parameter);
	
	//! Generate output noisy sequence.
    /*!
      \param input Vector for input. Should contains only 0s and 1s.
	  \param ChannelOutput Object for output. 
	  \sa NoisySeq
    */
	void GenerateOutput(int input[], NoisySeq &ChannelOutput); 
	//void SetAWGNSTD(double std){mSTD = std;}
	//void SetBSCCrossoverProb(double p){mCrossoverProb = p;}
	//! Constructor.
	Channel(int blocklength){mBlocklength = blocklength;}
	//! destructor.
	~Channel()
	{
#ifdef DEBUG
	cout<<"~Channel()"<<endl;
#endif
	}

	//! Set fixed seed.
    /*!
      \param seed. Seed
    */
	void SetSeed(int seed);
};


//!  Class for decoded sequence. 
/*!
  Stores the sequence from decoder, also contains information about decoder status.
*/
class DecodedSeq{
public:
	int Blocklength; /*!< blocklength of code */

	double *SoftInfo; /*!< soft information */
	int *HardDecision; /*!< 0/1 decision */
	bool AlgorithmConverge; /*!< true if the decoder converged */
	bool GlobalMin; /*!< true if for ADMM with Penalty, decoder converges to a ML CW*/
	int Iteration; /*!< number of iteration used */
	bool ValidCodeword; /*!< true if the codeword is a valid codeword */
	bool LPIntegerSol; /*!< true if the codeword is an integral solution from ADMM-LP decoder */
	unsigned long int ExeTime; /*!< decoding time */
	bool IsNan, IsInf;
	//! Constructor.
    /*!
      \param blocklength Blocklength of code
	  \param channel Channel name, can only be "AWGN" or "BSC"
    */
	DecodedSeq(int blocklength)
	{
		Blocklength = blocklength;
		HardDecision = new int[blocklength];
		SoftInfo = new double[blocklength];
		IsNan = false;
		IsInf = false;
		LPIntegerSol = false;
	}

	//! destructor.
	~DecodedSeq()
	{
#ifdef DEBUG
	cout<<"~DecodedSeq()"<<endl;
#endif
		delete [] HardDecision;
		delete [] SoftInfo;
	}
};

//!  Class for decoder
/*!
  Main program for ADMM decoders. Also included a sum-product decoder and a min-sum decoder. These two are not well optimized.
*/
class Decoder{
public:
	
	string DecoderName;
	string ChannelName;
	string mRealDecoderName;

	//! Projection Algorithm. Use two-slice property. 
    /*!
      \param v Input vector 
	  \param length Length of the vector
	  \return Projection onto the parity polytope
    */
	double* mProjectionPPd(NODE v[], int length);

	// parameters for decoders
	string mPCMatrixFileName;  /*!< parity check matrix file */
	int mNChecks, mBlocklength, mPCheckMapSize; 
	int *CheckDegree; 
	int *VariableDegree;
	int *mF1V;
	int mPA_Alex;
	double RateOfCode; /*!< code rate, k/n */
	double mTGDec_a_DecodeTime;

	TRIPLE *mPCheckMap; /*!< parity check mapping */
	TRIPLE *mParityCheckMatrix; /*!< parity check matrix */
	
	MESSAGE<double> *u; /*!< messages in message passing */
	MESSAGE<double> *v;
	double *u0; 
	

	// algorithm parameter
	double para_end_feas;   /*!<  ADMM alg end tolerance, typically 1e-5 */
	double para_mu;   /*!< ADMM alg \mu, typically 5.5 */
	double para_rho; /*!< ADMM over relaxation para, typically 1.8, 1.9 */
	int maxIteration; /*!< max iteration */
	double para_dmp_coef_alpha; /*!< Damping coefficient - alpha, for Eff_BPDecoder  */
	bool BP_non_sat; /*!< true if nonsaturating BP is used */
	double BP_sat_value;  /*!< saturation value if saturating BP is used*/

	//! Initialize BP bipartite graph. Note this function is NOT efficient!!!
	/*!
	  Link u message with v messages. Simply search for all possible edges and set pointer.
	  The BP functions provided are not fully optimized. Please use with care.
	*/
	void mInitializeMessages(); 

	//! Learn the degree of the parity check matrix
	/*!
	  The parity check matrix should be in correct format. This function is useful for irregular LDPC codes.
	*/
	void mLearnParityCheckMatrix();
	//! Read and store the parity check matrix
	void mReadParityCheckMatrix();
	void mSetDefaultParameters();
	void UPDATE_L_CODE(int *row, char *CodeFileNameLiu,int m);
	//! Calculate log likelihood ratio using the output from the channel
	/*!
	  \param channeloutput Output from channel.
	  \sa NoisySeq
	*/
	void mGenerateLLR(NoisySeq &channeloutput);

	//Decoders: 
	void ADMMDecoder();
	void BPDecoder();
	void Eff_BPDecoder();
	void Eff_TG_Decoder_a();
	void IlanBPDecoder();
	void MSADecoder();
	// data
	double *_LogLikelihoodRatio; /*!< log-likelihood ratio from received vector */
	double *OutputFromDecoder; /*!<soft information from decoder. could be pseudocodeword or message from BP  */

	double alpha; /*!< the constant for penalty term */

	double normalizeBSC; 

	double *warmstartchecks;	/*!< warm start ADMM */
	double *warmstartsave;	/*!< save results for warm start*/
	int mNumOfGlobalMin;
	//decoder stats
	unsigned long int mExeTime; /*!< exevution time */
	bool mAlgorithmConverge; /*!< true if the decoder converges*/
	bool mGlobalMin;  /*!< true if for ADMM with Penalty, decoder converges to a ML CW*/
	int mIteration; /*!< number of iterations used for decoding */
	int mNumOfProj;
	bool mValidCodeword; /*!< true if the output is a valid codeword */
	// ADMM decoder select
	enum ADMMdecoding { 
                 LinearProgram, /*!< Enum value LP. */  
                 L1Penalty, /*!< Enum value L1. */  
                 L2Penalty,  /*!< Enum value L2. */  
				 EntropyPenalty, /*!< Enum value Entropy. */ 
				 GaussianPenalty /*!< Enum value Gaussian. */
               }
		enumADMM;
	// function for different x-updates
	double xUpdateADMMLP(double degree, double mu, double llr, double t);
	double xUpdateADMML1(double degree, double mu, double llr, double t, double alpha);
	double xUpdateADMML2(double degree, double mu, double llr, double t, double alpha);
	double xUpdateADMMEntropy(double prevx, double degree, double mu, double llr, double t, double alpha);
	double xUpdateADMMGauss(double prevx, double degree, double mu, double llr, double t, double alpha);
	// function for nonsaturating BP
	double NonSatUpd(double llr1, double llr2);


	// functions for entropy penalty
	double entropy_th; /*!< Threshold for entropy penalty */  
	double funcWithEntropy(double x, double t, double d, double mu, double alpha)
	{
		double upth = 1 - entropy_th, lowth = entropy_th;
		double y;
		if(x < upth && x > lowth)
		{
			y = t/d + alpha/d/mu*log(x) - alpha/d/mu*log(1-x) - x;
		}
		if(x >= upth)
		{
			y = t/d + alpha/d/mu*log(upth) - alpha/d/mu*log(1 - upth) - upth
				+ (alpha/d/mu/upth/(1 - upth) - 1)*(x - upth);
		}
		if(x <= lowth)
		{
			y = t/d + alpha/d/mu*log(lowth) - alpha/d/mu*log(1 - lowth) - lowth
				+ (alpha/d/mu/lowth/(1 - lowth) - 1)*(x - lowth);
		}
		return y;
	}
	double funcWithEntropyDrv(double x, double t, double d, double mu, double alpha)
	{
		double upth = 1 - entropy_th, lowth = entropy_th;
		double y;
		if(x < upth && x > lowth)
		{
			y = alpha/d/mu/x/(1 - x) - 1;
		}
		if(x >= upth)
		{
			y = alpha/d/mu/upth/(1 - upth) - 1;
		}
		if(x <= lowth)
		{
			y = alpha/d/mu/lowth/(1 - lowth) - 1;
		}
		return y;
	}
	double newtonsol(double y, double xinit, int iter, double t, double d, double mu, double alpha)
	{
		if(iter < 0)
			iter = 5;
		double x = xinit;
		for(int i = 0; i < iter; i++)
		{
			double a = funcWithEntropy(x,t,d,mu,alpha);
			double b = funcWithEntropyDrv(x,t,d,mu,alpha);
			double xnew = x - (a - y)/b;
			if(xnew > 1) xnew = 1;
			if(xnew < 0) xnew = 0;
			if(abs (x - xnew) < 0.0000001)
			{
				x = xnew;
				break;
			}
			x = xnew;
		}
		return x;
	}

//public:
	//! Set decoder parameters.
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000 
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	*/
	void SetParameters(int mxIt, double feas,  double p_mu, double p_rho);

	// for Eff_BPDecoder
	void SetParameters(int mxIt, double dmp_coef_alpha);
	//! Set BSC normalized log likelihood ratios.
	/*!
	 \param normal if true, then the LLRs are 1 or -1. if false, then it is log(p/(1-p)) or log((1-p)/p)
	*/
	void SetBSCnormalized(bool normal);
	//! Set non-saturating BP decoder.
	/*!
	 \param nonsat True if the decoder should be NON-Saturated decoder
	 \param value Saturation level for saturated decoder. nonsat should be false for this option to work.
	*/
	void SetBPSaturation(bool nonsat, double value);
	//! Set decoder type.
	/*!
	 \param decoder Decoder name. Currently support "BP", "MSA", "ADMMLP", "ADMML1", "ADMML2" and "ADMMEntropy"
	*/
	void SetDecoder(string decoder);
	//! Set channel type.
	/*!
	 \param channel Channel name. Currently support "BSC" and "AWGN"
	*/
	void SetChannel(string channel){ChannelName = channel;}
	//! Set entropy penalty constant
	/*!
	 \param alp Penalty constant
	*/
	void SetPenaltyConstant(double alp){alpha = alp;}
	//! Set entropy penalty threshold
	/*!
	 \param th Threshold
	*/
	void SetEntropyParameter(double th){entropy_th = th;}
	//! See if the output from decoder is valid codeword
	/*!
		Note that the output soft information from ADMM and message passing algorithms(BP and MSA) are different.
	*/
	bool ValidateCodeword();
	//! See if the input sequence is valid codeword
	/*!
		\param input Input sequence array.
	*/
	bool ValidateCodeword(int input[]);
	//! Invoke decode function
	/*!
		\param channeloutput Noisy sequence from channel.
		\param decoderesult Decode results. Contains decoded sequence and other information.
		\sa DecodedSeq
		\sa NoisySeq
	*/
	void Decode(NoisySeq &channeloutput, DecodedSeq &decoderesult);
	

	
	//! Constructor.
	/*!
		\param FileName File for parity check matrix.
		\param blocklength Blocklength of code (i.e. n)
		\param nChecks Number of checks (i.e. row number of parity check matrix or (n-k))
	*/
	Decoder(string FileName,int nChecks, int BlockLength);
	//! Destructor
	~Decoder();
};

//!  Class for simulation
/*!
	Run simulations. 
	Take input sequence, pass it through channel to get noisy sequence. Pass noisy sequence to decoder and then compare decoded sequence with original input.
*/
class Simulator{
public:
	typedef unsigned long uint32;
	int *mInput;
	NoisySeq mNoisySequence;
	Channel mChannel;
	DecodedSeq mDecodedSequence;
	Decoder mDecoder;

	string mChannelName;
	string mDecoderName;

	string mParityCheckFile;
	string mCodewordFile;

	int mBlocklength;
	int mNChecks;
	double mRateOfCode;
	double EbN0;
	double mSTD;
	double mCrossOverProb;
	double mChannelParameter;
	bool mErrorFlag;

	// Simulation control
	uint32 mOutputFileInterval;
	string mOutputFileName;
	ofstream mOutputFileStream;

	uint32 mOutputCommandInterval;
	uint32 mOutputCommandDetailLevel;

	// simulation stats
	uint32 mTargetErrors; 
	uint32 mTargetSims;
	uint32 mTargetExeTime;

	uint32 mTotalExeTime;

	uint32 mTotalErrors;
	uint32 mTotalSims;
	uint32 mTotalDecodingTime;
	uint32 mTotalIterations;

	uint32 mTotalCorrectIterations;
	uint32 mTotalWrongIterations;

	uint32 mTotalCorrectDecodingTime;
	uint32 mTotalWrongDecodingTime;

	uint32 mTotalMLErrors;
	uint32 mTotalUndetectedErrors;

	//! The input is Eb/N0 (dB). This should be converted to std in AWGN or crossover probability in BSC
	void mTranslateEbN0()
	{
		mSTD = (double)sqrt(1.0 / pow(10.0, EbN0/10) / 2.0 / mRateOfCode);
		double temp = sqrt(2 * mRateOfCode * pow(10.0 , EbN0/10));
		mCrossOverProb = 0.5 - 0.5 * myerf(temp/sqrt(2.0));
		if(mChannelName == "AWGN")
			mChannelParameter = mSTD;
		else
			mChannelParameter = mCrossOverProb;
	}
	void mUpdateStats();
	void mClear()
	{
		mTotalErrors = 0;
		mTotalSims = 0;
		mTotalExeTime = 0;
		mTotalDecodingTime = 0;
		mTotalIterations = 0;
		mTotalCorrectIterations = 0;
		mTotalWrongIterations = 0;
		mTotalCorrectDecodingTime = 0;
		mTotalWrongDecodingTime = 0;
		mTotalMLErrors = 0;
		mTotalUndetectedErrors = 0;
	}
	void mOutputCommand()
	{
		if (mTotalSims % mOutputCommandInterval == 0)
		{
			if(mOutputCommandDetailLevel == 1)
			{
				cout<<"TotalSims = "<<mTotalSims
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<<"*"<<(mTotalDecodingTime*100)/mTotalExeTime<<"%(s)  Total error = "<<mTotalErrors
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<endl;
			}
			else if(mOutputCommandDetailLevel == 2)
			{
				cout<<"TotalSims = "<<mTotalSims
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<<"*"<<(mTotalDecodingTime*100)/mTotalExeTime<<"%(s)"
					<<"  Total error = "<<mTotalErrors
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Avg # iteration correct = "<<mTotalCorrectIterations/(double)(mTotalSims - mTotalErrors)
					<<"  Avg # iteration wrong = "<<mTotalWrongIterations/(double)(mTotalErrors)
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<"  Avg decode time correct = "<<(double)mTotalCorrectDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims - mTotalErrors) <<"(s)"
					<<"  Avg decode time wrong = "<<(double)mTotalWrongDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalErrors) <<"(s)"
					<<endl;
			}
		}
	}
	void mWriteFile()
	{
		int RoundsCount = 1;
		if (mTotalSims == mOutputFileInterval)
		{

				mOutputFileStream<<"TotalSims = "<<mOutputFileInterval<<"*"<<RoundsCount
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<<"*"<<(mTotalDecodingTime*100)/mTotalExeTime<<"%(s)  Total error = "<<mTotalErrors
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<endl;
				RoundsCount++;
				mClear();
		}
	}
//public:

	//! Return total number of errors
	uint32 GetTotalErrors(){return mTotalErrors;}
	//! Return total number of simulations
	uint32 GetTotalSims(){return mTotalSims;}
	//! Return time used for decoding
	uint32 GetTotalDecodingTime(){return mTotalDecodingTime;}
	//! Return total number of iterations
	uint32 GetTotalIterations(){return mTotalIterations;}
	//! Return total number of iterations for correct decoding events
	uint32 GetTotalCorrectIterations(){return mTotalCorrectIterations;}
	//! Return total number of iterations for erroneous decoding events
	uint32 GetTotalWrongIterations(){return mTotalWrongIterations;}
	//! Return time used for correct decoding events
	uint32 GetTotalCorrectDecodingTime(){return mTotalCorrectDecodingTime;}
	//! Return time used for erroneous decoding events
	uint32 GetTotalWrongDecodingTime(){return mTotalWrongDecodingTime;}
	//! Return total number of errors that are sure to be ML errors.
	/*!
		The ML error is calculated whenever there is a decoded codeword but is not the ML error. This is a lower bound for number of ML errors.
	*/
	uint32 GetTotalMLErrors(){return mTotalMLErrors;}
	//! Return total number of undetected errors. These are valid codeword but not the transmitted codeword.
	uint32 GetTotalUndetectedErrors(){return mTotalUndetectedErrors;}
	//! Return time for simulation. This include decoding time and time for channels, etc.
	uint32 GetTotalExeTime(){return mTotalExeTime;}
	//! Set channel's random number seed
	/*!
		\param seed Seed for channel.
	*/
	void SetChannelSeed(int seed){mChannel.SetSeed(seed);}
	//! Set simulation targets
	/*!
		\param tError Target errors. Stop simulation if more than this number of errors are collected.
		\param tSim Target simulations. Stop simulation if more than this number of simulations are performed.
		\param tExeTime Target time. Stop simulation if more than this amount of time is consumed.
	*/
	void SetTargets(uint32 tError, uint32 tSim, uint32 tExeTime){mTargetErrors = tError; mTargetSims = tSim; mTargetExeTime = tExeTime;}
	//! Set codeword. 
	/*!
		Currently there is no encoding function. But you can choose a valid codeword to simulation if you have one other than the all zero codeword.
		\param filename Codeword file. The file should contain codeword. It should be 0,1 sequence, broken by space. e.g. 0 1 0 0 1 0
	*/
	void SetCodeword(string filename);
	//! Set command line output appearance. 
	/*!
		\param simInterval Interval for each command line output. 
		\param detaillevel Detail level for output. Only takes values 1 and 2.
	*/
	void SetCommandLineOutput(int simInterval, int detaillevel){mOutputCommandInterval = simInterval; mOutputCommandDetailLevel = detaillevel;}

	//! Set decoder name. 
	/*!
	\param decoder Decoder name.
	 \sa Decoder
	*/
	void SetDecoder(string decoder){mDecoderName = decoder; mDecoder.SetDecoder(mDecoderName);}

	//! Set decoder parameters. Sufficient for ADMM LP.
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000 
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho);
	}
	//! Set decoder parameters. Sufficient for ADMM L1 and ADMM L2
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000 
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param alp Penalty constant.
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, double alp)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetPenaltyConstant(alp);
	}
	//! Set decoder parameters. Sufficient for all ADMM algorithms
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000 
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param alp Penalty constant.
	 \param th Penalty thershold for entropy penalty
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, double alp, double th)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetPenaltyConstant(alp); mDecoder.SetEntropyParameter(th);
	}
	//! Set decoder parameters. Sufficient for BP decoding
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000 
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param nansat True if non-saturating BP is used
	 \param value Saturation value is original BP is used
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, bool nonsat, double value)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetBPSaturation(nonsat, value);
	}
	//! Set decoder parameters. Setting for normalized LLR for BSC
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000 
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param normal True if the LLRs from BSC are normalized to +1 or -1.
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, bool normal)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetBSCnormalized(normal);
	}
	void SetSNR(double SNR)
	{
		EbN0 = SNR;
		mTranslateEbN0();
		mNoisySequence.SetProperties(mChannelName, mChannelParameter);
		mChannel.SetChannel(mBlocklength, mChannelName, mChannelParameter);
	}
	//! Run simulations
	void RunSim(double *r, double *x, double *x_est, bool *b);

	//! Constructor
	/*!
		\param blocklength Blocklength of the code
		\param nchecks Number of checks
		\param pcfilename Parity check matrix file name. 
		\param channelname Channel name.
		\param eb_over_nzero Eb/N0 for the channel.
	*/
	Simulator(int blocklength, int nchecks, string pcfilename, string channelname, double eb_over_nzero) 
		:mNoisySequence(blocklength) 
		, mChannel(blocklength)
		, mDecodedSequence(blocklength)
		, mDecoder(pcfilename, nchecks, blocklength)
	{
		mChannelName = channelname;
		mBlocklength = blocklength;
		mNChecks = nchecks;
		mParityCheckFile = pcfilename;
		mRateOfCode = (double)(mBlocklength - mNChecks)/mBlocklength;
		EbN0 = eb_over_nzero;
		mTranslateEbN0();

		mNoisySequence.SetProperties(mChannelName, mChannelParameter);
		mChannel.SetChannel(mBlocklength, mChannelName, mChannelParameter);
		mDecoder.SetChannel(mChannelName);
		mInput = new int[mBlocklength];
		for(int i = 0; i < mBlocklength; i++)
			mInput[i] = 0;
	}
	~Simulator()
	{
		delete [] mInput;
	}
};
