
/*
  To run: g++ ex_bitwise_7.cpp -o ex_bitwise_7 && ./ex_bitwise_7
 */

// SF optimized to get minimum d(jetPtCorr)


/*
CHeck these:

https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/
https://www.delftstack.com/howto/cpp/read-csv-file-in-cpp/



 */


// bitset operators
#include <iostream>       // std::cout
#include <string>         // std::string
#include <bitset>         // std::bitset
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip> //#include <format>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <map>
#include <algorithm>

using namespace std;



const int printLevel = 0;
const int nLineToRead = -1;
const int nColsToRead = 7; // No. of columns to read from input lut_SFInDecimal file. File format: <col0: IEtaBin> <col1: QuantileMinPtInGeV> <col2: QuantileMaxPtInGeV> <col3: CalibrationAtPtInGeV> <col4: IEtaBin_forLUT> <col5: Pt_forLUT> <col6: SF_Decimal>
const unsigned int nCompBinMax = 2048;
const int IEtaBinOffsetForEtaCompressedLUT = 0; // 0: IEtaBin_forLUT = IEtaBin = [1, 41];  -1: IEtaBin_forLUT = IEtaBin - 1 = [0, 40]. # 0 give correct calibration.
const int Mode_calculateSFInBits = 1; // 0: calculate 'multiplier' such that 'addend' need not be zero ,
                                      // 1: calculate 'multiplier' such that 'addend=0'. Best performer with capLowPtSFAt2 = true;
                                      // 2: calculate 'multiplier' and 'addend' so as to have least (pTCorr - pTCorrTarget). Set 'capLowPtSFAt2 = false'. plot_checkJEC_v6p5. Not working properly.
const bool capLowPtSFAt2 = true; // if SF(pT<15 GeV) > 2 then SF(pT<15 GeV) = 2

std::string sInFile_SFs = "/home/siddhesh/Work/CMS/L1_Trigger_Work/L1T_ServiceTask/hcalPUsub/myAna/hcalPUsub_v5_20220311/run_1/makeLUTs/LUTs/Default_RawPUS_SF_RegressedTo_log_GenJetPt_v6/lut_calib_2022v1_ECALZS_decimal.txt";

std::string sOutFile = "/home/siddhesh/Work/CMS/L1_Trigger_Work/L1T_ServiceTask/hcalPUsub/myAna/hcalPUsub_v5_20220311/run_1/makeLUTs/LUTs/Default_RawPUS_SF_RegressedTo_log_GenJetPt_v6/lut_calib_2022v1_ECALZS.txt";






const unsigned int UnitySFInBits = 512;
// 1./512 = 0.0019531250 = least count of SF
const unsigned int MaxAdded      = pow(2, 8);  // 256
const unsigned int MaxMultiplier = pow(2, 10); // 2014
const int8_t       int8_t_Min = -128;
const int8_t       int8_t_Max =  127;




int8_t calculateAddend(unsigned int jet_hwPt, unsigned int multiplier, unsigned int jet_calibPt_target);



unsigned int calculateAddPlusMult(unsigned int multiplier, int8_t addend) {
  // addPlusMult_0 is calculated for 32-bit. The following procedure is done to calculate addPlusMult in 18-bits
  unsigned int    addPlusMult_0 = (addend << 10) + multiplier;
  std::bitset<18> addPlusMult_0_bitset(addPlusMult_0);
  unsigned int    addPlusMult   = (unsigned int)(addPlusMult_0_bitset.to_ulong());
  return addPlusMult;
}

unsigned int getMultiplier(unsigned int addPlusMult) {
  unsigned int multiplier = addPlusMult & 0x3ff; //  0x3ff = 1023 = 11 1111 1111
  return multiplier;
}

int8_t getAddend(unsigned int addPlusMult) {
  int8_t addend = (addPlusMult >> 10);
  return addend;
}

unsigned int calibrateJet(unsigned int jet_hwPt, unsigned int addPlusMult, int printLevel_i=0) {
  if (printLevel_i >= 2)
    printf("calibrateJet():: \n");
  
  std::bitset<18> addPlusMult_bitset(addPlusMult);
  unsigned int multiplier = addPlusMult & 0x3ff; //  0x3ff = 1023 = 11 1111 1111
  std::bitset<10> multiplier_bitset(multiplier);
  // handles -ve numbers correctly
  int8_t addend = (addPlusMult >> 10);
  std::bitset<8> addend_bitset(addend);

  unsigned int jetPtCorr = ((jet_hwPt * multiplier) >> 9) + addend;
  std::bitset<18> jet_hwPt_bitset(jet_hwPt);
  std::bitset<18> jetPtCorr_bitset(jetPtCorr);
  std::bitset<18> jetPtCorr1_bitset((jet_hwPt * multiplier));
  std::bitset<18> jetPtCorr2_bitset(((jet_hwPt * multiplier) >> 9));

  if (printLevel_i >= 2) {
    std::cout << "addPlusMult: " << addPlusMult << " " << addPlusMult_bitset
	      << ", multiplier:" << multiplier << " " << multiplier_bitset
	      << ", addend:" << static_cast<int16_t>(addend) << " -> " << (addPlusMult >> 10) << " " << addend_bitset
	      << "\n\t jet_hwPt: " << jet_hwPt << " " << jet_hwPt_bitset
	      << ", (jet_hwPt * multiplier): " << (jet_hwPt * multiplier) << " " << jetPtCorr1_bitset
	      << ", ((jet_hwPt * multiplier) >> 9): " << ((jet_hwPt * multiplier) >> 9) << " " << jetPtCorr2_bitset
	      << ", pTCorr: " << jetPtCorr << " " << jetPtCorr_bitset
      //<< ", :" <<
	      << "\n";
  }

  return jetPtCorr;
}


unsigned int calibrateJetPt(unsigned int jet_hwPt, unsigned int multiplier, int8_t addend, int printLevel_i=0) {
  unsigned int jetPtCorr = ((jet_hwPt * multiplier) >> 9) + addend;
  if (printLevel_i >= 2) {
    std::cout << " " << jet_hwPt << " (" << multiplier << ", " << static_cast<int16_t>(addend) << ") --> " << jetPtCorr << ", ";
  }
  return jetPtCorr;
}


unsigned int calculateTargetJetPtCorr(unsigned int jet_hwPt, double SFInDecimal) {
  unsigned int jet_calibPt_target = (unsigned int)( ( ((double)jet_hwPt) * SFInDecimal ) + 0.5); // 0.5 to round off fractional part of the number
  return jet_calibPt_target;
}


int calculateSumDeltaPt(unsigned int jet_hwPtMin, unsigned int jet_hwPtMax, double SFInDecimal, unsigned int SFInBits) {
  int SumDeltaPt = 0;
  for (unsigned int jet_hwPt_i = jet_hwPtMin; jet_hwPt_i <=jet_hwPtMax; jet_hwPt_i++) {
    unsigned int jetPtCorrTarget_i = calculateTargetJetPtCorr(jet_hwPt_i, SFInDecimal);
    unsigned int jetPtCorr_i       = calibrateJet(jet_hwPt_i, SFInBits);
    int absDPt_tmp = abs(((int)jetPtCorr_i) - ((int)jetPtCorrTarget_i));
    SumDeltaPt += absDPt_tmp;
  }
  return  SumDeltaPt;
}



unsigned int calculateSFInBits(unsigned int jet_hwPt, double SFInDecimal) {
  if (printLevel >= 2) printf("calculateSFInBits():: \n");
  
  //unsigned int jet_calibPt_target = (unsigned int)( ( ((double)jet_hwPt) * SFInDecimal ) + 0.5); // 0.5 to round off fractional part of the number
  unsigned int jet_calibPt_target = (unsigned int)( round( ((double)jet_hwPt) * SFInDecimal ) ); // 0.5 to round off fractional part of the number
  if (printLevel >= 2) {
    std::cout << ((double)jet_hwPt)
	      << ", " << SFInDecimal
	      << ", PtCorr: " <<  ( ((double)jet_hwPt) * SFInDecimal )
	      << " / " << jet_calibPt_target
	      << "\n";
  }
  
  unsigned int multiplier = (unsigned int)(SFInDecimal * ((double)UnitySFInBits));
  unsigned int multiplier_0 = multiplier;
  // multiplier is 10-bit, which is Power(2,10) = 1024 = MaxMultiplier
  if (multiplier >= MaxMultiplier ) { // truncate multiplier to 1023
    multiplier = MaxMultiplier - 1;
  }

  if (printLevel >= 2) {
  std::cout << "jet_hwPt: " << jet_hwPt
	    << ", SFInDecimal:" << SFInDecimal
	    << ", multiplier:" << multiplier_0 << " -> " << multiplier 
	    << ", ((double)UnitySFInBits):" << ((double)UnitySFInBits)
	    << ", (SFInDecimal * ((double)UnitySFInBits)):" << (SFInDecimal * ((double)UnitySFInBits))
	    << "\n";
  }

  int absDPt_Min = 99999;
  int8_t addend = 0;
  int dPt_iMinus1 = 0;
  //for (int8_t addend_i = 0; addend_i < (int8_t)UnitySFInBits; addend_i++) {
  //for (int8_t addend_i = 0; addend_i < 100; addend_i++) { // used for LUTs
  for (int8_t addend_i = 0; addend_i < MaxAdded; addend_i++) {
    unsigned int jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend_i;
    int dPt = (int)jetPtCorr_i - (int)jet_calibPt_target;
    int absDPt = std::abs(dPt);
    if (absDPt < absDPt_Min) {
      addend = addend_i;
      absDPt_Min = absDPt;
    }

    if (printLevel >= 2) {
      std::cout << "\tjet_hwPt: " << jet_hwPt
		<< ", jetPtCorr: " << jetPtCorr_i 
		<< " / " << jet_calibPt_target
		<< ", dPt: " << ((int)jetPtCorr_i - (int)jet_calibPt_target)
		<< ",  " << ((int)jetPtCorr_i - (int)jet_calibPt_target)/(double)jet_calibPt_target
		<< "%, multiplier: " << multiplier
		<< ", addend: " << static_cast<int16_t>(addend_i)
		<< "\n";
    }
    
    //if ((dPt * dPt_iMinus1) < 0) {
    if (dPt == 0) {
      break;
    }
    dPt_iMinus1 = dPt;    
  }

  

  unsigned int jetPtCorr = ((jet_hwPt * multiplier) >> 9) + addend;
  unsigned int addPlusMult = (addend << 10) + multiplier;

  if (printLevel >= 2) {
    std::cout << "jet_hwPt: " << jet_hwPt
	      << ", jetPtCorr: " << jetPtCorr 
	      << " / " << jet_calibPt_target
	      << ", dPt: " << ((int)jetPtCorr - (int)jet_calibPt_target)
	      << ",  " << (jetPtCorr - jet_calibPt_target)/(double)jet_calibPt_target
	      << "%, multiplier: " << multiplier
	      << ", addend: " << static_cast<int16_t>(addend)
	      << ", addPlusMult: " << addPlusMult
	      << "\n";
  }
  
  return addPlusMult;
}





unsigned int calculateSFInBits_addendZero(unsigned int jet_hwPt, double SFInDecimal) {
  if (printLevel >= 5)
    printf("calculateSFInBits_addendZero():: \n");
  
  unsigned int jet_calibPt_target = (unsigned int)( ( ((double)jet_hwPt) * SFInDecimal ) + 0.5); // 0.5 to round off fractional part of the number
  if (printLevel >= 5) {
  std::cout << ((double)jet_hwPt)
	    << ", " << SFInDecimal
	    << ", PtCorr: " <<  ( ((double)jet_hwPt) * SFInDecimal )
	    << " / " << jet_calibPt_target
	    << "\n";
  }
  
  unsigned int multiplier = (unsigned int)(SFInDecimal * ((double)UnitySFInBits));
  // multiplier is 10-bit, which is Power(2,10) = 1024 = MaxMultiplier
  if (multiplier >= MaxMultiplier ) { // truncate multiplier to 1023
    multiplier = MaxMultiplier - 1;
  }
  unsigned int multiplier_0 = multiplier;

  std::vector<unsigned int> multiplier_wAddendZero;
  for (unsigned int m_i = 0; m_i < MaxMultiplier; m_i++) {
    unsigned int m_tmp = multiplier + m_i;
    if (m_tmp >= MaxMultiplier) break;

    int8_t a_i = calculateAddend(jet_hwPt, m_tmp, jet_calibPt_target);
    if (a_i == 0) multiplier_wAddendZero.push_back( m_tmp );

    if (a_i < 0) break;
  }



    
  if ( multiplier_wAddendZero.size() > 0 ) {
    int idx_multiplier = int( multiplier_wAddendZero.size() / 2 );
    multiplier = multiplier_wAddendZero[idx_multiplier];
    int8_t addend = 0;
    
    unsigned int jetPtCorr = ((jet_hwPt * multiplier) >> 9) + addend;
    unsigned int addPlusMult = (addend << 10) + multiplier;
    return addPlusMult;
  }
  else {
  
    multiplier = multiplier_0;
  
    int absDPt_Min = 99999;
    int8_t addend = 0;
    for (int8_t addend_i = 0; addend_i < 127; addend_i++) { // int8_t range [-128, 127]
      unsigned int jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend_i;
      int dPt = (int)jetPtCorr_i - (int)jet_calibPt_target;
      int absDPt = std::abs(dPt);
      if (absDPt < absDPt_Min) {
	addend = addend_i;
	absDPt_Min = absDPt;
      }

      if (dPt == 0) {
	break;
      }
    }


    unsigned int addPlusMult = (addend << 10) + multiplier;
    return addPlusMult;
  }
    
}



int8_t calculateAddend(unsigned int jet_hwPt, unsigned int multiplier, unsigned int jet_calibPt_target) {
  int8_t addend = 0;
  unsigned int jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend;
  int dPt = (int)jetPtCorr_i - (int)jet_calibPt_target;
  if (dPt == 0) {
    return addend;
  }
  else if (dPt < 0) {
    int absDPt_Min = 99999;
    addend = 0;
    for (int8_t addend_i = 0; addend_i < 127; addend_i++) {
      jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend_i;
      dPt = (int)jetPtCorr_i - (int)jet_calibPt_target;
      int absDPt = std::abs(dPt);
      if (absDPt < absDPt_Min) {
	addend = addend_i;
	absDPt_Min = absDPt;
      }
    }
    return addend;
  }
  else if (dPt > 0) {
    int absDPt_Min = 99999;
    addend = 0;
    for (int8_t addend_i = -1; addend_i > -128; addend_i--) {
      jetPtCorr_i = ((jet_hwPt * multiplier) >> 9) + addend_i;
      dPt = (int)jetPtCorr_i - (int)jet_calibPt_target;
      int absDPt = std::abs(dPt);
      if (absDPt < absDPt_Min) {
	addend = addend_i;
	absDPt_Min = absDPt;
      }
    }
    return addend;
  }

  return addend;
}






unsigned int calculateSFInBits_stableSFsOverPtBin(unsigned int jet_hwPt, unsigned int jet_hwPtMin, unsigned int jet_hwPtMax, double SFInDecimal) {
  if (printLevel >= 2)
    printf("calculateSFInBits_stableSFsOverPtBin():: \n");


  unsigned int multiplier;
  int8_t       addend;
  double absDPtSum_Min = 99999999.0;
  for (unsigned int multiplier_i = 0;          multiplier_i  < MaxMultiplier; multiplier_i++) {
    if (printLevel >= 100) std::cout << " m_i: " <<  multiplier_i << "\n";
    
    for (int8_t     addend_i     = int8_t_Min; addend_i     < int8_t_Max;    addend_i++) {
      if (printLevel >= 100) std::cout << "\t a_i: " << static_cast<int16_t>(addend_i)  << "\n";
      
      unsigned int jetPtCorrTarget = calculateTargetJetPtCorr(jet_hwPt, SFInDecimal);
      unsigned int jetPtCorr       = calibrateJetPt(jet_hwPt, multiplier_i, addend_i);
	
      // Require perfect calibration at pT=PtCalibPoint
      if (jetPtCorr != jetPtCorrTarget) continue;
      
      // abs(pTCorr - pTCorrTarget) should be minimum over the full pT bin
      int absDPtSum = 0;
      for (unsigned int jet_hwPt_i = jet_hwPtMin; jet_hwPt_i <=jet_hwPtMax; jet_hwPt_i++) {
	if (jet_hwPt_i == 0) continue;
	
	unsigned int jetPtCorrTarget_i = calculateTargetJetPtCorr(jet_hwPt_i, SFInDecimal);
	unsigned int jetPtCorr_i       = calibrateJetPt(jet_hwPt_i, multiplier_i, addend_i);

	double absDPt_tmp = abs(((double)jetPtCorr_i) - ((double)jetPtCorrTarget_i));
	//double absDPt_tmp = abs(((double)jetPtCorr_i) - ((double)jetPtCorrTarget_i)) / ((double)jetPtCorrTarget_i);
	absDPtSum += absDPt_tmp;

	if (printLevel >= 30) {
	  std::cout << "\t\t\t m_i: " << multiplier_i << ", a_i: " << static_cast<int16_t>(addend_i)
		    << ",  jet_hwPt_i: " << jet_hwPt_i
		    << ", absDPt_tmp: " << absDPt_tmp << ", absDPtSum: " << absDPtSum
		    << "\n";
	}
      }

      if (printLevel >= 3) {
	std::cout << "\t\t m_i: " << multiplier_i << ", a_i: " << static_cast<int16_t>(addend_i)
		  << ", jet_hwPt: " << jet_hwPt
		  << ", jetPtCorr: " << jetPtCorr << " / " << jetPtCorrTarget
		  << ", absDPtSum: " << absDPtSum << ", absDPtSum_Min: " << absDPtSum_Min
		  << "\n";	  
      }     

      if (absDPtSum < absDPtSum_Min) {
	multiplier = multiplier_i;
	addend     = addend_i;
	absDPtSum_Min = absDPtSum;
      }
    }
  }

  //unsigned int addPlusMult = (addend << 10) + multiplier;
  unsigned int addPlusMult = calculateAddPlusMult(multiplier, addend); // correct way for 18bit addPlusMult

  
  if (printLevel >= 2) {
    unsigned int jet_calibPt_target = calculateTargetJetPtCorr(jet_hwPt, SFInDecimal);
    unsigned int jetPtCorr          = calibrateJetPt(jet_hwPt, multiplier, addend);

    
    std::cout << "jet_hwPt: " << jet_hwPt
	      << ", jetPtCorr: " << jetPtCorr
	      << " / " << jet_calibPt_target
	      << ", dPt: " << ((int)jetPtCorr - (int)jet_calibPt_target)
	      << ",  " << (jetPtCorr - jet_calibPt_target)/(double)jet_calibPt_target
	      << "%, multiplier: " << multiplier
	      << ", addend: " << static_cast<int16_t>(addend)
	      << ", addPlusMult: " << addPlusMult
	      << ", absDPtSum_Min: " << absDPtSum_Min
	      << "\n";
  }

  return addPlusMult;
}













int main ()
{


  /*
  std::string sIEtaColm    = "L1JetTowerIEtaAbs";
  std::string sPtColm      = "L1JetRawEnergy";
  std::string sBDTTypeColm = "ScaleFactor_Log(GenJet)";
  */


  std::cout << "sInFile_SFs:" << sInFile_SFs << "\n";

  std::ifstream infile(sInFile_SFs.c_str());


  std::vector<std::string> dataColmNames;
  std::vector<std::vector<double>> data;
  
  
  
  std::string line, word;
  std::vector<std::string> row;
  int count=0;


  // Read column names ----------------------------------
  if (! infile.good() ) {
    printf("%s could not open **** ERROR **** \n",sInFile_SFs.c_str());
    return 0;
  }

  // Read data from 2nd line ---------------------------
  printf("\n\nRead data: \n");
  while (std::getline(infile, line)) {
    // used for breaking words
    std::istringstream iss(line);
    std::stringstream ss(line);
    if (printLevel >= 5) {
      std::cout << "Line " << count << ": " << line << "\n";
    }
    
    
    std::vector<double> dataInRow;
    int iColm = 0;
    // read every column data of a row and
    // store it in a string variable, 'word'
    while (std::getline(ss, word, ' ')) {
      double n = std::stod( word );
      if (printLevel >= 5) {
	std::cout << ">" << word << " | " << n << "<" << std::endl;
      }
      //row.push_back(word);
      //data[iColm].push_back(n);
      dataInRow.push_back(n);
      iColm++;

      if (iColm >= nColsToRead) break;
    }
    
    data.push_back(dataInRow);
    count++;
    if (nLineToRead > 0 && count > nLineToRead) break;
  }


  std::cout << "data.size(): " << data.size()
	    << ", data[0].size(): " << data[0].size()
	    << std::endl;

  std::ofstream outfile(sOutFile.c_str());


  // LUT JEC header --------------------------------------------
  outfile << "# address to addend+multiplicative factor LUT\n"
	  << "# maps 11 bits to 18 bits\n"
	  << "# 18 bits = (addend<<10) + multiplier)\n"
	  << "# addend is signed 8 bits, multiplier is 10 bits\n"
	  << "# anything after # is ignored with the exception of the header\n"
	  << "# the header is first valid line starting with #<header> versionStr nrBitsAddress nrBitsData </header>\n"
	  << "#<header> v1 11 18 </header>\n";
    ;

  // -----------------------------------------------------------
  

  int IEtaBin_forLUT = -1;
  unsigned int compBin_last = 0;
  for (size_t iData=0; iData<data.size(); iData++) {
    unsigned int iEta          = static_cast<unsigned int>(data[iData][0]);
    unsigned int PtMin_iQuant  = static_cast<unsigned int>(data[iData][1] + 0.5);
    unsigned int PtMax_iQuant  = static_cast<unsigned int>(data[iData][2] + 0.5);
    unsigned int Pt            = static_cast<unsigned int>(data[iData][3] + 0.5);
    
    unsigned int etaBin        = static_cast<unsigned int>(data[iData][4]);
    unsigned int ptBin         = static_cast<unsigned int>(data[iData][5]);
    unsigned int compBin       = (etaBin << 4) | ptBin;
    
    double       SFInDecimal   = data[iData][6];
    
    unsigned int jetHwPt       = Pt * 2; // jetHwPt = 2*jetPtInGeV https://github.com/cms-sw/cmssw/blob/master/L1Trigger/L1TCalorimeter/src/firmware/Stage2Layer2JetAlgorithmFirmwareImp1.cc#L717 . Applying a factor of 2 does not affect SFInDecimal --> SFInBits calculation here.
    unsigned int jetHwPtMin    = PtMin_iQuant * 2;
    unsigned int jetHwPtMax    = PtMax_iQuant * 2;
    
    unsigned int SFInBits;
    if      (Mode_calculateSFInBits == 0)  SFInBits    = calculateSFInBits(jetHwPt, SFInDecimal);
    else if (Mode_calculateSFInBits == 1)  SFInBits    = calculateSFInBits_addendZero(jetHwPt, SFInDecimal);
    else if (Mode_calculateSFInBits == 2)  SFInBits    = calculateSFInBits_stableSFsOverPtBin(jetHwPt, jetHwPtMin, jetHwPtMax, SFInDecimal);
    unsigned int jetHwPtCorr     = calibrateJet(jetHwPt, SFInBits);
    unsigned int jetHwPtCorr_Exp = calculateTargetJetPtCorr(jetHwPt, SFInDecimal); //static_cast<unsigned int>( ((double)(jetHwPt) * SFInDecimal) + 0.5);
    //std::cout << ", SFInBits: " << SFInBits << " (m: " << (SFInBits & 0x3ff) << ", a: " << (SFInBits >> 10) << ") \n";
    
    std::string sComment = "";
    if ( ptBin == 0) {
      IEtaBin_forLUT++;
      //unsigned int etaBin_calo = etaBin;
      //if (IEtaBinOffsetForEtaCompressedLUT == 0) etaBin_calo = etaBin + 1;
      //sComment = Form("    # eta_bin %d, pt %d", IEtaBin_forLUT, Pt);
      sComment = "  # ieta " + std::to_string(iEta) + ", pt " + std::to_string(ptBin);  // + " "
    }

    unsigned int SFInBits_0 = SFInBits;
    std::string sComment1 = "";
    if ( capLowPtSFAt2 && (ptBin == 0) ) {
      unsigned int addPlusMult = SFInBits;
      unsigned int multiplier = addPlusMult & 0x3ff; //  0x3ff = 1023 = 11 1111 1111
      int8_t addend = (addPlusMult >> 10);
      if (multiplier == (MaxMultiplier - 1) && (addend > 0)) {
	//addend = 0;
	//SFInBits = multiplier; // Cap SFInBits at (MaxMultiplier - 1)=1023
	SFInBits = 2047; // Cap SFInBits at 2047, which is SF=2.008665 with multiplier=1023  and addend=1 
	sComment1 = "  >>> cap SF to 2: SFInBits: " + std::to_string(SFInBits_0) + " --> " + std::to_string(SFInBits) + ", ";
      }
    }

    int          AbsDeltaPt_SumOverPtQuant = calculateSumDeltaPt(jetHwPtMin, jetHwPtMax, SFInDecimal, SFInBits);
    
    if (printLevel >= 0) {
      std::cout << "etaBin: " << etaBin
		<< ", ptBin: " << ptBin
		<< ", compBin: " << compBin
		<< ", Pt: " << Pt << " (" << PtMin_iQuant << ", " <<  PtMax_iQuant << ")"
		<< ", SFInDecimal: " << SFInDecimal
		<< ", jetHwPt: " << jetHwPt
		<< ", SFInBits: " << SFInBits << " (m: " << getMultiplier(SFInBits) << ", a: " << static_cast<int16_t>(getAddend(SFInBits)) << ")"
		<< ", jetHwPtCorr: " << jetHwPtCorr
		<< " / " << jetHwPtCorr_Exp
		<< "\t\t error: " << (double)(jetHwPtCorr) / (double)(jetHwPtCorr_Exp)
		<< ",  sum(dPt) overPtBin: " << AbsDeltaPt_SumOverPtQuant
		<< sComment.c_str()
		<< sComment1.c_str()
		<< std::endl;
    }

    
    outfile << compBin 
	    << " " << SFInBits
	    << sComment.c_str()
	    << "\n";

    compBin_last = compBin;
  }

  std::cout  << "compBin_last: " << compBin_last
	      << ", nCompBinMax: " << nCompBinMax
	      << "\n";

  compBin_last++;
  
  for (compBin_last; compBin_last < nCompBinMax; compBin_last++) {
    outfile << compBin_last
	    << " 0 # dummy"
	    << "\n";
  }

  outfile.close();

  printf("\nlut_calib %s is being produced!!\n", sOutFile.c_str());

  











}
