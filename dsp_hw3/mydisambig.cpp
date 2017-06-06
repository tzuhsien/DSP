#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <array>
#include <vector>
#include <unordered_map>

#include "Ngram.h"

using namespace std;

// Get P(W2 | W1) -- bigram
double getBigramProb(Vocab& voc, Ngram& lm, const char *w1, const char *w2)
{
	VocabIndex wid1 = voc.getIndex(w1);
	VocabIndex wid2 = voc.getIndex(w2);

	if(wid1 == Vocab_None)  //OOV
		wid1 = voc.getIndex(Vocab_Unknown);
	if(wid2 == Vocab_None)  //OOV
		wid2 = voc.getIndex(Vocab_Unknown);

	VocabIndex context[] = { wid1, Vocab_None };
	return lm.wordProb( wid2, context);
}

void doViterbi(const string& line, unordered_map < string , vector < string > >& map, Vocab& voc, Ngram& lm)
{
	//fill seq
	stringstream str(line);
	vector< string > seq;
	seq.push_back("<s>");
	string buffer;
	while (str >> buffer) {
		seq.push_back(buffer);
	}
	seq.push_back("</s>");
	
//	cerr << "finish fill seq" << endl;
	//Viterbi
	vector < vector <double> > delta(seq.size());
	vector < vector <unsigned int> > backtrack(seq.size());

	delta[0].push_back(0.0);
	for (unsigned int i = 1; i < seq.size(); i++) {
		for (unsigned int j = 0; j < map[seq[i]].size(); j++) {
//			cerr <<  getBigramProb(voc, lm, map[seq[i - 1]][0].c_str(), map[seq[i]][j].c_str()) << endl;
//			cerr << "i = " << i << " j = " << j << endl;
//			cerr << delta[i-1].size() << endl;
//			cerr << delta[i- 1][0] << endl;
			double maxprobability = getBigramProb(voc, lm, map[seq[i - 1]][0].c_str(), map[seq[i]][j].c_str()) + delta[i - 1][0];
			unsigned int maxIndex = 0;
			for (unsigned int k = 1; k < map[seq[i - 1]].size(); k++) {
				double probability = getBigramProb(voc, lm, map[seq[i - 1]][k].c_str(), map[seq[i]][j].c_str())+ delta[i - 1][k];
				if (probability > maxprobability) {
					maxprobability = probability;
					maxIndex = k;
				}
			}
			delta[i].push_back(maxprobability);
			backtrack[i].push_back(maxIndex);

		}
	}

//	cerr << "finish Viterbi" << endl;
	//fill result
	vector < string > result;
	unsigned int Index = 0;
	for (unsigned int i = seq.size() - 1; i > 0; i--) {
		result.push_back(map[seq[i]][Index]);
		Index = backtrack[i][Index];
	}
	result.push_back("<s>");
	for(auto rit = result.rbegin();rit != result.rend();rit++){
		cout << *rit;
		if(rit+1 != result.rend()){
			cout << " ";
		}
		else{
			cout << endl;
		}
	}
//	cerr << "finish fill result" << endl;
}

void ReadMap(unordered_map < string, vector < string > >& map, fstream& mapfile)
{
	string line;
	while (getline(mapfile, line)) {
		stringstream str(line);
		string key;
		str >> key;
		vector < string > chinesechars;
		string chinesechar;
		while(str >> chinesechar){
			chinesechars.push_back(chinesechar);
		}
		map.insert(make_pair(key, chinesechars));
	}

	//insert lineBegin to ZhtoBig5
	{
		string key("<s>");
		vector < string > linebegin;
		linebegin.push_back(key);
		map.insert(make_pair(key, linebegin));
	}
	//insert lineEnd to ZhtoBig5
	{
		string key("</s>");
		vector < string > lineend;
		lineend.push_back(key);
		map.insert(make_pair(key, lineend));
	}
	return;
}

int main(int argc, char** argv)
{
	fstream inputfile, mapfile;
	Vocab voc;
	inputfile.open(argv[2], ios::in);
	mapfile.open(argv[4], ios::in);
	int ngram_order = atoi(argv[8]);
	Ngram lm( voc, ngram_order );

	{
		File lmFile( argv[6] , "r" );
		lm.read(lmFile);
		lmFile.close();
	}
//	cerr << "before read map" << endl;
	unordered_map < string , vector < string > > map;
	ReadMap(map,mapfile);
	
//	cerr << "finish read map" << endl; 
	string line;
	while(getline(inputfile,line)){
		doViterbi(line,map,voc,lm);
	}
}
