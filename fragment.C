//////////////////////////////////////////////////////////////////////
// fragment.C  Copyright (c) 2014 Dario Ghersi                      //
// Version: 20140307                                                //
// Goal:    fragment small molecules with-user defined cleavage     //
//          rules (e.g., RECAP rules ), encoded as SMARTS patterns  //
// Usage:   fragmenter -r RULES -i SMALL_MOLECULES -o OUTFILE       //
//                     -n MIN_FRAGMENT_SIZE                         //
//                     -k MAX_COMBINED_FRAGMENTS                    //
//		       [-m MAX_FRAGMENT_SIZE] [-s SIZE_RATIO]	    		//
//		       [-w MAX_MOL_WEIGHT] [-x]			    				//
//								    								//
//          See User's Guide for more details                       //
//                                                                  //
// This file is part of the MOLBLOCKS suite.                        //
// MOLBLOCKS is free software: you can redistribute it and/or       //
// modify it under the terms of the GNU General Public License as   //
// published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.              //
//                                                                  //
// MOLBLOCKS is distributed in the hope that it will be useful,     //
// but WITHOUT ANY WARRANTY; without even the implied warranty of   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    //
// GNU General Public License for more details.                     //
//                                                                  //
// You should have received a copy of the GNU General Public        //
// License along with MOLBLOCKS.                                    //
// If not, see <http://www.gnu.org/licenses/>.                      //
//                                                                  //
// MODIFIED VERSION (Kathrin Heikamp)                               //
// Version: 20150408                                                //
// Modification: 4 additional input parameter & 1 removed (-e)      //
// -k maxCombinedFragments                                          //
//	Instead of deriving only terminal fragments, the smallest   	//
//	fragments can be combined up to the number of               	//
//	maxCombinedFragments.                                       	//
//	Example: -k 3 results in the generation of all fragments    	//
//	consisting of up to 3 smallest fragments.                  	 	//
// -m maxFrag		                                            	//
//	The parameter specifies the maximal size of fragments as    	//
//	the maximal number of atoms allowed in fragments.           	//
//	Example: -m 30 limits the number of atoms in the fragments  	//
//	to 30 heavy atoms.			                    				//
// -s sizeRatio                                                     //
//	The parameter specifies the maximal size of fragments as    	//
//	ratio of the size of the parent compound.                   	//
// 	Example: -s 0.5 filters out all compounds whose size        	//
//	(number of atoms) is larger than 50% of the atoms of the    	//
//	original molecule.                                          	//
// -w maxWeight		                                            	//
//	The parameter specifies the maximal molecular weight of     	//
//	fragments allowed.				            					//
//	Example: -w 300 limits the molecular weight of the	    		//
//	fragments to 300.			                    				//
// -x                                                               //
//	The parameter specifies whether dummy atoms should be       	//
//  attached to the atoms that formed cut bonds. The flag       	//
//	turns on writing dummy atoms to fragment smiles.            	//
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iostream>
#include <iterator>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/parsmart.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>

using namespace OpenBabel;
using namespace std;
#include "fragment.h"
#include "utilities.h"

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

void checkCommandLineArgs(char **argv, int argc)
{

  bool err = false;
  if (!cmdOptionExists(argv, argv+argc, "-i")) {
    cerr << "Input file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-r")) {
    cerr << "Using default RECAP rules\n";
  }
  if (!cmdOptionExists(argv, argv+argc, "-o")) {
    cerr << "Output file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-n")) {
    cerr << "Please specify the minimum fragment size\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-k")) {
    cerr << "Please specify the maximum number of combined (minimal) fragments\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-m")) {
    cerr << "No maximal fragment size specified\n";
  }
  if (!cmdOptionExists(argv, argv+argc, "-s")) {
    cerr << "No maximal fragment size relative to parent compound specified\n";
  }
  if (!cmdOptionExists(argv, argv+argc, "-w")) {
    cerr << "No maximal molecular weight specified\n";
  }

  if (err){
    cerr << USAGE;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////////////

void fragAllMols(const Parameters & p, vector<OBSmartsPattern> & rules)
{
  // main function -- fragment one molecule at a time, and print the
  // results to file
  // N.B. the functions expects the small molecules to be represented
  // as SMILES string, one molecule per line

  fstream smiFile, outFile;
  string line, smiles, name = "";
  OBConversion conv;
  unsigned long totNumMol = 0, count = 0;
  unsigned long oldPercentage = 0, percentage = 0;

  // set the input to SMILES, and the output to canonical SMILES
  conv.SetInAndOutFormats("SMI", "CAN");

  // open the small molecules file
  smiFile.open(p.smiFileName, fstream::in);

  // complain if the file doesn't exist
  if (! smiFile.good()) {
    cerr << "Can't open " << p.smiFileName << endl;
    exit(1);
  }

  // count how many molecules (lines) are in the file
  while (getline(smiFile, line)) {
    totNumMol++;
  }

  // rewind the small molecules file
  smiFile.clear();
  smiFile.seekg(0, ios::beg);

  // open the output file
  outFile.open(p.outFileName, fstream::out);
  cout << endl;

  // process each small molecule
  while (getline(smiFile, line)) {
    // reset the name of the molecule to the empty string
    name = "";

    // extract the SMILES string
    istringstream iss(line);
    iss >> smiles;

    // extract the name (if present)
    iss >> name;

    // build the molecule
    Molecule mol(conv, smiles, name);

    // fragment the molecule
    mol.fragment(rules, conv, p);

    // size of molecule
    int size = mol.moleculeSize();

    // print the results
    mol.printResults(outFile, conv, p.minFrag, p.maxFrag, size, p.sizeRatio, p.maxWeight);
    outFile.flush();

    count++;

    // call the progress bar function every 10 molecules
    if ((count % 10) == 0 || totNumMol < 10) {
      percentage = 100.0 * count / totNumMol;
      if (percentage > oldPercentage) {
	oldPercentage = percentage;
	printProgBar(percentage);
      }
    }
  }
  printProgBar(100.0 * count / totNumMol);
  cout << endl;

  // close the open streams
  smiFile.close();
  outFile.close();
}

//////////////////////////////////////////////////////////////////////

Molecule::Molecule(OBConversion conv, const string & smiles, const string & nn)
{
  // constructor for the Molecule class
  singletons.reserve(1000);
  conv.ReadString(&mol, smiles);

  // add the name, if not empty
  if (!nn.empty()) {
    name = nn;
  }
  else {
    name = "";
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::calculateDepMat(int minFrag, bool dummy)
{
  // calculate the dependence matrix between cleavable bonds
  // A[i, j] = 1 if bond i is independent from bond j
  // A[i, j] = 0 otherwise

  unsigned int n = matchPairs.size();

  // allocate memory for the dependence matrix
  depMat = (bool **) malloc(sizeof(bool *) * n);
  depMat[0] = NULL;
  for (unsigned int i = 0; i < n; i++) {
    depMat[i] = (bool *) malloc(sizeof(bool) * n);
  }

  // process each bonds
  bool isCleaved;
  int atom1, atom2;
  OBMol temp; // intermediate object to store the molecule after the
              // first cut

  // set the diagonal to "false"
  for (unsigned int i = 0; i < n; i++) {
    depMat[i][i] = false;
  }

  for (unsigned int i = 0; i < n; i++) {

    mol = original; // restore the original molecule

    // get the bond atoms
    atom1 = matchPairs[i][0];
    atom2 = matchPairs[i][1];

    // cleave the bond
    isCleaved = cleaveBond(atom1, atom2, minFrag, dummy);

    temp = mol;

    for (unsigned int j = i + 1; j < n; j++) {
      mol = temp; // restore the molecule after the first cut

      // get the bond atoms
      atom1 = matchPairs[j][0];
      atom2 = matchPairs[j][1];

      // attempt to cleave
      isCleaved = cleaveBond(atom1, atom2, minFrag, dummy);

      // fill in the dependence matrix
      if (isCleaved) {
	depMat[i][j] = depMat[j][i] = true;
      }
      else {
	depMat[i][j] = depMat[j][i] = false;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::go(int offset, int k)
{
  // get all possible combinations of cuttable bonds

  if (k == 0){
    independentSets.push_back(combinations);
    return;
  }
  for (int i = offset; i <= (int) matchPairs.size() - k; ++i) {
    combinations.push_back(i);
    go(i+1, k-1);
    combinations.pop_back();
  }
}

//////////////////////////////////////////////////////////////////////

inline bool Molecule::cleaveBond(int atom1, int atom2, int minFrag, bool dummy)
{
  // cleave the bond between two atoms

  // get the corresponding bond
  OBBond *bond = mol.GetBond(atom1, atom2);
  if (bond) {
    mol.DeleteBond(bond);

    // if dummy atoms should be written, create to new atoms and add 
    // to the atoms forming the deleted bond
    if (dummy) {
      mol.AddBond(atom1, mol.NewAtom()->GetIdx(), bond->GetBondOrder());
      mol.AddBond(atom2, mol.NewAtom()->GetIdx(), bond->GetBondOrder());
    }
    return true;
  }
  else {
    return false; // bond not cleaved
  }

  return false;
}

//////////////////////////////////////////////////////////////////////

inline void Molecule::cut(int minFrag, bool dummy)
{
  // cut the bonds between the atoms that match the patterns
  const unsigned int matchPairsSize = matchPairs.size();
  for (unsigned int i = 0; i < matchPairsSize; i++) {
    cleaveBond(matchPairs[i][0], matchPairs[i][1], minFrag, dummy);
  } 
}

//////////////////////////////////////////////////////////////////////

inline void Molecule::cut(int minFrag, int singleton, bool dummy)
{
  // cut the singleton bond

  cleaveBond(matchPairs[singleton][0], matchPairs[singleton][1],
	     minFrag, dummy);
}

//////////////////////////////////////////////////////////////////////

void Molecule::cut(int minFrag, const list<int> & perm, bool dummy)
{
  // cut the bonds between the atoms that match the patterns,
  // cleaving the bonds in the order specified in 'perm',
  for (list<int>::const_iterator iter = perm.begin(); iter != perm.end(); iter++) {
    cleaveBond(matchPairs[*iter][0], matchPairs[*iter][1],
	       minFrag, dummy);
  }

}

//////////////////////////////////////////////////////////////////////

void Molecule::findCleavableBonds( vector<OBSmartsPattern> & rules,
				  int minFrag) {
  // apply the cleavage rules and find the matching atom pairs

  // make sure that only uniqye matching atom pairs are found
  // i.e., hash of atom pair combinations
  std::map<std::string,int> uniquepairs;

  for (vector<OBSmartsPattern>::iterator rule = rules.begin();
       rule != rules.end(); rule++) {

    if ((*rule).Match(mol)) { // the pattern matches

      // find all pairs of matching atoms
      const vector<vector<int> > & maplist = (*rule).GetMapList();
      int atom1, atom2;

      for (unsigned int i = 0; i < maplist.size(); i++) {
        atom1 = maplist[i][0];
	atom2 = maplist[i][1];

        // get string representation of current atom pair
	std::string strpair;
	std::stringstream sstm;
	if (atom1 < atom2){
	  sstm << atom1 << "_" << atom2;
	} else {
	  sstm << atom2 << "_" << atom1;
	}
	strpair = sstm.str();

        // if current atom pair is not already in the list of atoms
	if (uniquepairs.find(strpair) == uniquepairs.end()){
	  
          // update hash of atom pair combinations
	  uniquepairs.insert(std::pair<std::string,int>(strpair,1));

	  // add the atoms to the list of atoms
	  vector<int> pair(2);
	  pair[0]=atom1; //pair.push_back(atom1);
	  pair[1]=atom2; //pair.push_back(atom2);
	  matchPairs.push_back(pair);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::fragment(  vector<OBSmartsPattern> & rules,
			const OBConversion & conv, const Parameters & p)
{
  // fragment the molecule according to 'rules'

  // find the cleavable bonds
  findCleavableBonds(rules, p.minFrag);

  // proceed if there is at least one cleavable bond
  if (matchPairs.size() > 0) {

    // first, make a copy of the molecule
    original = mol;
 
    // cluster the bonds to be cleaved into mutually exclusive sets
    getIndependentSets(p.minFrag, p.maxCombinedFragments, p.writeDummyAtoms);

    unsigned int indepSetSize = independentSets.size();
    // carry out the fragmentation for each independent set
    for (unsigned int i = 0; i < indepSetSize; i++) {
      mol = original;
      cut(p.minFrag, independentSets[i], p.writeDummyAtoms);
      storeFragments(conv);
    }

    unsigned int singeltonSize = singletons.size();
    // Take care of the singleton bonds
    for (unsigned int i = 0; i < singeltonSize; i++) {
      mol = original;
      cut(p.minFrag, singletons[i], p.writeDummyAtoms);
      storeFragments(conv);
    }

  }

  else { // no cleavable bonds, just print the results
    storeFragments(conv);
  }
}

//////////////////////////////////////////////////////////////////////

inline int Molecule::moleculeSize()
{
  OBMol temp = original;
  
  //int size = temp.NumAtoms();
  int size = 0;
  FOR_ATOMS_OF_MOL(a, temp){
    if ((! a->IsHydrogen()) && a->GetAtomicNum()>0) {
      size++;
    }
  }
  return size;
}

//////////////////////////////////////////////////////////////////////

inline int Molecule::fragmentSize(int posAtom1, int posAtom2)
{
  // calculate the size of the fragment comprised between atom1
  // and atom2  

  // get the atoms between atom1 and atom2
  vector<int> children;
  mol.FindChildren(children, posAtom1, posAtom2);

  // make sure we count only heavy atoms (i.e., ignore H)
  OBAtom *atom;
  int fragSize = 0;
  for (unsigned int i = 0; i < children.size(); i++) {
    atom = mol.GetAtom(children[i]);
    if (! atom->IsHydrogen()) {
      fragSize++;
    }
  }
  
  return fragSize + 1;
}

//////////////////////////////////////////////////////////////////////

inline int Molecule::fragmentSize( OBMol & fragment) const
{
  // calculate the size of the fragment 
  // make sure we count only heavy atoms (i.e., ignore H) and not
  // dummy atoms (i.e., ignore atoms with atomic number 0)
  int fragSize = 0;
  FOR_ATOMS_OF_MOL(a, fragment){
    if ((! a->IsHydrogen()) && a->GetAtomicNum()>0) {
      fragSize++;
    }
  }
  
  return fragSize;
}

//////////////////////////////////////////////////////////////////////

inline double Molecule::molecularWeight(   OBMol & fragment) const
{
  double weight = fragment.GetMolWt();
  return weight;
}

//////////////////////////////////////////////////////////////////////

inline void Molecule::freeDepMat()
{
  // free the memory allocated for the distance matrix

  unsigned int n = matchPairs.size();

  for (unsigned int i = 0; i < n; i++)
    free(depMat[i]);
  free(depMat);

}

//////////////////////////////////////////////////////////////////////

void Molecule::getIndependentSets(int minFrag, int maxCombinedFragments, bool dummy)
{
  // calculate the graph theoretic distance (shortest number of
  // bonds) between all atoms, and cluster the atoms pairs into either
  // dependent (if their distance < minFrag) or independent groups
  // using complete linkage

  // build the dependence matrix between the cleavable bonds
  calculateDepMat(minFrag, dummy);

  // NEW CODE
  // -> derive all possible combinations of bonds that should be cut up to a maximum
  // number of combined fragments

  // get start point of for loop
  int max = matchPairs.size() - maxCombinedFragments + 1;
  
  // if the number of combined fragments is larger than the number of bonds to be cut, 
  // generate all possible fragments
  if (max < 0){
      max = 0;
  }
  
  // if the number of combined fragments is negative or zero, generate all possible fragments
  if (maxCombinedFragments <= 0){
      max = 0;
  }
  
  // for loop
  // loop over all possible combinations of bonds that should be cut starting at max
  // idea: by doing this, all combinations of non-cuttable bonds (num <= maxCombinedFragments-1)
  // are generated 
  for (unsigned int i = max; i <= matchPairs.size(); ++i) {
    go(0,i);
  }

  // free memory for the distance matrix
  freeDepMat();

}

//////////////////////////////////////////////////////////////////////

void Molecule::printResults(fstream &outFile, OBConversion conv,
			    unsigned int minFrag, unsigned int maxFrag,
			    int size, float sizeRatio, double maxWeight)
{
  // print the molecule to file

  // further remove redundancy from fragments by converting them
  // to molecules and back to strings

  vector<string> nonRedundant;
  OBMol temp;
  string currFrag;

  std::map<string,bool>::const_iterator iter    = fragments.begin();
  std::map<string,bool>::const_iterator iterEnd = fragments.end();  
  for ( ; iter != iterEnd; iter++ ) {  
    conv.ReadString(&temp, iter->first);      

    unsigned int fragmentsize = fragmentSize(temp);
    double fragmentweight = molecularWeight(temp);

    // filter fragments by minimal fragment size
    if (fragmentsize >= minFrag) {
      // filter fragments by maximal fragment size (if defined)
      if ((maxFrag <= 0) || (maxFrag > 0 && fragmentsize <= maxFrag)){
	// filter fragments by relative size to parent ligand (if defined)
	if ((sizeRatio <= 0) || (sizeRatio > 0 && fragmentsize <= (float) size*sizeRatio)){
	  // filter fragments by molecular weight (if defined)
	  if ((maxWeight <= 0) || (maxWeight > 0 && fragmentweight <= maxWeight)){
            currFrag = conv.WriteString(&temp, true);
            if (find(nonRedundant.begin(), nonRedundant.end(), currFrag) ==
	        nonRedundant.end()) {
	      nonRedundant.push_back(currFrag);
	    }
	  }
        }
      }
    }
  }

  bool needDot = false;
  bool notEmpty = false;
  if (nonRedundant.size() > 0) {
    notEmpty = true;
  }
  for (vector<string>::iterator fragment = nonRedundant.begin();
       fragment != nonRedundant.end(); fragment++) {
    if (needDot) {
      outFile << "." << *fragment;
    }
    else {
      outFile << *fragment;
      needDot = true;
    }
  }

  // print the name (if not empty)
  if (notEmpty) {
    if (!name.empty()) {
      outFile << "\t" << name << endl;
    }
    else {
      outFile << endl;
    }
  }
  else {
    if (!name.empty()) {
      outFile << "<NA>\t" << name << endl;
    }
    else {
      outFile << "<NA>\t" << endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::storeFragments(OBConversion conv)
{
  // store the fragments (if not already stored)

  // convert the molecule into a string
  const string & molString = conv.WriteString(&mol, true);

  // the '.' character separate the fragments...
  istringstream iss(molString);
  string fragment;
  OBMol temp;
  while (getline(iss, fragment, '.')) {
    // include only "new" fragments
         
    if(fragments.find(fragment) == fragments.end()) {
      fragments.insert(std::pair<string,bool>(fragment,true));
    }
  }
}

//////////////////////////////////////////////////////////////////////

Parameters::Parameters(char **argv, int argc)
{
  // parse the command-line arguments

  // get file names
  rulesFileName = getCmdOption(argv, argv + argc, "-r");
  smiFileName = getCmdOption(argv, argv + argc, "-i");
  outFileName = getCmdOption(argv, argv + argc, "-o");

  // get minimum size of fragmetns
  stringstream temp(getCmdOption(argv, argv + argc, "-n"));
  temp >> minFrag;

  // sanity check
  if (minFrag < 2) {
    cerr << "The minimum fragment size (-n) should be > 1\n";
    exit(1);
  }

  // get maximum number of combined minimal fragments
  stringstream temp1(getCmdOption(argv, argv + argc, "-k"));
  temp1 >> maxCombinedFragments;

  // get maximum size of fragments
  if (cmdOptionExists(argv, argv + argc, "-m")){
    stringstream temp2(getCmdOption(argv, argv + argc, "-m"));
    temp2 >> maxFrag;

    // sanity check
    if (maxFrag <= 0) {
      cerr << "The maximum fragment size is <= 0, i.e., no filtering by maximum size is applied\n";
    }
    else if (maxFrag < minFrag) {
      cerr << "The maximum fragment size (-m) should be > the minimum fragment size (-n)\n";
      exit(1);
    }
  }
  else {
    maxFrag = -1;
  }

  // get maximal size ratio of fragments
  if (cmdOptionExists(argv, argv + argc, "-s")){
    stringstream temp3(getCmdOption(argv, argv + argc, "-s"));
    temp3 >> sizeRatio;

    // sanity check
    if (sizeRatio <= 0) {
      cerr << "The relative size (-s) is <= 0, i.e., no filtering by relative size is applied\n";
    }
  }
  else {
    sizeRatio = -1;
  }

  // get maximal molecular weight of fragments
  if (cmdOptionExists(argv, argv + argc, "-w")){
    stringstream temp4(getCmdOption(argv, argv + argc, "-w"));
    temp4 >> maxWeight;

    // sanity check
    if (maxWeight <= 0) {
      cerr << "The maximum molecular weight (-w) is <= 0, i.e., no filtering by molecular weight is applied\n";
    }
  }
  else {
    maxWeight = -1;
  }

  // get optional paramter for writing cutting points
  writeDummyAtoms = cmdOptionExists(argv, argv + argc, "-x");
}

//////////////////////////////////////////////////////////////////////

vector<OBSmartsPattern> storeCleavageRules(const char *rulesFileName)
{
  // store the cleavage rules as a vector of SMARTS PATTERN objects

  // declare the pattern vector
  vector<OBSmartsPattern> rules;

  // open the rules file
  fstream rulesFile;
  rulesFile.open(rulesFileName, fstream::in);

  // complain if the file doesn't exist
  if (! rulesFile.good()) {
    cerr << "Can't open " << rulesFileName << endl;
    exit(1);
  }

  // store the rules
  string line, rule;
  unsigned int count = 0;
  while (getline(rulesFile, line)) {

    // add an element to the pattern vectors
    rules.push_back(OBSmartsPattern());

    // process the pattern
    istringstream iss(line);
    iss >> rule;
    rules[count++].Init(rule);
  }

  rulesFile.close();

  return rules;
}

//////////////////////////////////////////////////////////////////////

vector<OBSmartsPattern> storeCleavageRules()
{
  // store the cleavage rules as a vector of SMARTS PATTERN objects

  // declare the pattern vector
  vector<OBSmartsPattern> rules;


  // store the rules
  string line, rule;
  for (unsigned int i = 0; i < NUM_RECAP; i++) {

    // add an element to the pattern vectors
    rules.push_back(OBSmartsPattern());

    // process the pattern
    rules[i].Init(RECAP[i]);
  }

  return rules;
}

//////////////////////////////////////////////////////////////////////
// MAIN PROGRAM                                                     //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // check the command-line arguments
  checkCommandLineArgs(argv, argc);

  // get the parameters
  Parameters p(argv, argc);
  vector<OBSmartsPattern> rules;
  if (p.rulesFileName) {  // load the SMARTS cleavage rulse

    rules = storeCleavageRules(p.rulesFileName);
  }
  else { // use default RECAP rules
    rules = storeCleavageRules();
  }

  // process the molecules
  fragAllMols(p, rules);
  cout << endl;

  return 0;
}
