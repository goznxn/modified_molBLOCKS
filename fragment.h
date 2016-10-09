//////////////////////////////////////////////////////////////////////
// fragment.h  Copyright (c) 2014 Dario Ghersi                      //
// Version: 20140307                                                //
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
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// CONSTANTS                                                        //
//////////////////////////////////////////////////////////////////////


const string USAGE = "fragmenter -r RULES -i SMALL_MOLECULES\n\
                    -o OUTFILE -n MIN_FRAGMENT_SIZE \n\
		    -k MAX_COMBINED_FRAGMENTS \n\
		    [-m MAX_FRAGMENT_SIZE] [-s SIZE_RATIO] \n\
		    [-w MAX_MOL_WEIGHT] [-x]\n";

// default RECAP rules
const unsigned int NUM_RECAP = 11;
const char *RECAP[NUM_RECAP] = {
  "[$([C!$(C([#7])(=O)[!#1!#6])](=[O]))]!@[#7!$([#7][!#1!#6])]",
  "[$(C=!@O)]!@[$([O;+0])]",
  "[#6]!@[N;!$(N=*);!$(N[#6]=[!#6]);!$(N~[!#1!#6])!X4]",
  "[$(C(=!@O)([#7;+0;D2,D3])!@[#7;+0;D2,D3])]!@[$([#7;+0;D2,D3])]",
  "[O!$(O[#6]~[!#1!#6])]([#6])!@[#6]",
  "[#6]=!@[#6]",
  "[N;+1;D4]!@[#6]",
  "[$([n;+0])]-!@[#6!$([#6]=[!#6])]",
  "[$(N(@-C(=O)))]!@-[#6!$([#6]=[!#6])]",
  "[c]!@[c]",
  "[$([#7;+0;D2,D3])]-!@[$([S](=[O])=[O])]"
};

//////////////////////////////////////////////////////////////////////
// CLASSES                                                          //
//////////////////////////////////////////////////////////////////////

class Parameters {

 public:
  char *rulesFileName;
  char *smiFileName;
  char *outFileName;
  int minFrag;
  int maxFrag;
  float sizeRatio;
  double maxWeight;
  int maxCombinedFragments;
  bool writeDummyAtoms;

  Parameters(char **, int);
};

//////////////////////////////////////////////////////////////////////

class Molecule {

 private:
  OBMol mol;
  OBMol original; // copy of the original (for extensive fragment.)
  string name;
  int size;
  int fragmentsize;
  double fragmentweight;
  vector<vector<int> > matchPairs;
  std::map<string,bool> fragments;
  bool **depMat;
  vector<int> singletons;

  std::map<std::string,int> uniquepairs;
  std::string strpair;
  
  vector< list<int>   > independentSets;
  list<int>             combinations;
  
  void cut(int, const list<int> &, bool);
  void calculateDepMat(int, bool);
  void getIndependentSets(int, int, bool);
  bool cleaveBond(int, int, int, bool);
  void cut(int, bool);
  void cut(int, int, bool);
  

  void dealWithAmbiguousCuts();
  void findCleavableBonds( vector<OBSmartsPattern>&, int);
  int fragmentSize(int, int);
  int fragmentSize( OBMol&) const;
  double molecularWeight( OBMol&) const;
  void freeDepMat();
  void storeFragments(OBConversion);
  void go(int, int);

 public:
  Molecule (OBConversion, const string&, const string&);
  void fragment(vector<OBSmartsPattern> &, const OBConversion &,
		const Parameters &);
  void getMatchAtoms();
  void printResults(fstream &, OBConversion, unsigned int, unsigned int, int, float, double);
  int moleculeSize();
};

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

void checkCommandLineArgs(char **, int);
void fragAllMols(const Parameters&, vector<OBSmartsPattern>&);
vector<OBSmartsPattern> storeCleavageRules(const char *);
vector<OBSmartsPattern> storeCleavageRules();
