

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
const int nColsToRead = 5; // No. of columns to read from input lut_SFInDecimal file. File format: <col0: IEtaBin> <col1: PtInGeV> <col2: IEtaBin_forLUT> <col3: Pt_forLUT> <col4: SF_Decimal>
const unsigned int nCompBinMax = 2048;
const int IEtaBinOffsetForEtaCompressedLUT = 0; // 0: IEtaBin_forLUT = IEtaBin = [1, 41];  -1: IEtaBin_forLUT = IEtaBin - 1 = [0, 40]. # 0 give correct calibration.
std::string sInFile_SFs = "/home/siddhesh/Work/CMS/L1_Trigger_Work/L1T_ServiceTask/hcalPUsub/myAna/hcalPUsub_v5_20220311/run_1/makeLUTs/LUTs/Default_RawPUS_SF_RegressedTo_log_GenJetPt_v6/lut_calib_2022v1_ECALZS_decimal.txt";

std::string sOutFile = "/home/siddhesh/Work/CMS/L1_Trigger_Work/L1T_ServiceTask/hcalPUsub/myAna/hcalPUsub_v5_20220311/run_1/makeLUTs/LUTs/Default_RawPUS_SF_RegressedTo_log_GenJetPt_v6/lut_calib_2022v1_ECALZS.txt";






const unsigned int UnitySFInBits = 512;
// 1./512 = 0.0019531250 = least count of SF
const unsigned int MaxAdded      = pow(2, 8);  // 256
const unsigned int MaxMultiplier = pow(2, 10); // 2014


unsigned int calibrateJet(unsigned int jet_hwPt, unsigned int addPlusMult) {
  if (printLevel >= 2)
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

  if (printLevel >= 2) {
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
  for (int8_t addend_i = 0; addend_i < 100; addend_i++) {
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
    
    // File format:
    // <col0: PtEtaCombineBinNumber> <col1: PtInGeV> <col2: SF_Decimal>
    
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
    unsigned int iEta        = static_cast<unsigned int>(data[iData][0]);
    unsigned int Pt          = static_cast<unsigned int>(data[iData][1] + 0.5);
    
    unsigned int etaBin      = static_cast<unsigned int>(data[iData][2]);
    unsigned int ptBin       = static_cast<unsigned int>(data[iData][3]);
    unsigned int compBin     = (etaBin << 4) | ptBin;
    
    double       SFInDecimal = data[iData][4];
    
    unsigned int jetHwPt     = Pt * 2; // jetHwPt = 2*jetPtInGeV https://github.com/cms-sw/cmssw/blob/master/L1Trigger/L1TCalorimeter/src/firmware/Stage2Layer2JetAlgorithmFirmwareImp1.cc#L717 . Applying a factor of 2 does not affect SFInDecimal --> SFInBits calculation here.
    unsigned int SFInBits    = calculateSFInBits(jetHwPt, SFInDecimal);
    unsigned int jetHwPtCorr = calibrateJet(jetHwPt, SFInBits);
    unsigned int jetHwPtCorr_Exp = static_cast<unsigned int>( ((double)(jetHwPt) * SFInDecimal) + 0.5);

    std::string sComment = "";
    if ( ptBin == 0) {
      IEtaBin_forLUT++;
      //unsigned int etaBin_calo = etaBin;
      //if (IEtaBinOffsetForEtaCompressedLUT == 0) etaBin_calo = etaBin + 1;
      //sComment = Form("    # eta_bin %d, pt %d", IEtaBin_forLUT, Pt);
      sComment = "  # ieta " + std::to_string(iEta) + ", pt " + std::to_string(ptBin);  // + " "
    }


    if (printLevel >= 0) {
      std::cout << "etaBin: " << etaBin
		<< ", ptBin: " << ptBin
		<< ", compBin: " << compBin
		<< ", Pt: " << Pt
		<< ", SFInDecimal: " << SFInDecimal
		<< ", jetHwPt: " << jetHwPt
		<< ", SFInBits: " << SFInBits
		<< ", jetHwPtCorr: " << jetHwPtCorr
		<< " / " << jetHwPtCorr_Exp
		<< "\t\t error: " << (double)(jetHwPtCorr) / (double)(jetHwPtCorr_Exp)
		<< sComment.c_str()
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
