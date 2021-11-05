#include <iostream>  //always, std, string
#include <vector>
#include <cmath>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
// *********** To adapt: *********** //
  const unsigned int n = 9;    // number of spins: only for function "print_basis" to use "bitset<n>()"
  
  //INPUT FILE:
  const string INPUT_folder = "INPUT/"; 
  //string INPUT_filename = "";

  //OUTPUT FILE:
  const string OUTPUT_directory = "OUTPUT/"; 
  //string filename = "";  // OUTPUT FILE NAME

  // Criteria at alpha*sigma  (ex.  1*sigma or 3*sigma)
  const unsigned int alpha = 3;

// *********** don't modify: *********** //
  //string OUTPUT_filename = OUTPUT_directory + filename;

  //string OUTPUT_filename_Data_VS_Model = OUTPUT_directory + filename + "_dataVSmodel.dat"; //"Birds_dataVSmodel_2ndOrder.dat";
  //const string OUTPUT_filename_ProbaSpace = OUTPUT_directory + filename + "_ProbaSpace.dat";

  //
  const uint32_t un = 1;
  const uint32_t NOp_tot = (un << n) - 1;                     // number of operators = 2^n - 1


/********************************************************************/
/**************************    STRUCTURES    ************************/
/********************************************************************/
struct Interaction
{
  uint32_t Op;      // binary operator associated to the interaction
  unsigned int k;   // order of the interaction
  double g;         // parameter of the interaction in current representation
  double g_Ising;   // parameter of the interaction in {-1,+1} representation  
  double g_ACE;     // parameter of the interaction in {0,1} representation
  double av_D;      // empirical average ("d" for data)
  double av_M;      // average in the Model
};

//Structure with the final information for the probability of appearance of each operator in the dataset
struct Operator
{
  uint32_t bin;     // binary representation of the operator
  mutable unsigned int layer;        // to which layer the operator belongs --> known only after the selection of the best Basis: by default, equal to n (=last layer)
  //unsigned int k1;  // nb of point where op = 1 --> it's a R.V.:  k1 = sum(op[s^i])
  mutable double p1_M;     // in the model: probability that op = 1 
  double p1_D;     // in the data: probability that op = 1 --> rem: it's a R.V. = sum(op[s^i]) / N

  double S;           // - [ p1*log(p1) + (1-p1)*log(1-p1) ]
  mutable double DKL;
  bool operator < (const Operator &other) const   // for ranking Operators from the most to the less likely
    { return (DKL > other.DKL || (DKL == other.DKL && bin < other.bin)); }
};





