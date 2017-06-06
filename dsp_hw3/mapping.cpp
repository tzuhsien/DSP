#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <unordered_map>
#include <set>
#include <string>
#include <fstream>
#define BUFFER_SIZE 300

using namespace std;

int main(int argc, char **argv)
{
	fstream input;
	fstream file;
	input.open(argv[1],ios::in);
	file.open(argv[2], ios::out | fstream::trunc);
	char buffer[BUFFER_SIZE] = {0};
	unordered_map<string, string> map;
	while (input.getline(buffer, BUFFER_SIZE)) {
		string textspeak;
		string chinese;
		set<string> setoftextspeak;
		
		char* pch = strchr(buffer, ' ');
		while (pch != NULL) {
			textspeak.clear();
			textspeak.push_back(pch[1]);
			textspeak.push_back(pch[2]);
			setoftextspeak.insert(textspeak);
			pch = strchr(pch + 1, '/');
		}
		
		chinese.append(buffer, 2);
		setoftextspeak.insert(chinese);
		for(set<string>::iterator it = setoftextspeak.begin(); it!=setoftextspeak.end(); it++) {
//			file << "textspeak : " << *it << endl;
			unordered_map<string, string>::iterator got = map.find(*it);
			if (got == map.end()) {
//				file << "in\n";
				map.insert(make_pair(*it, chinese));
				got = map.find(*it);
//				file << got->first << " " << got->second << endl;
			}
			else {
				got->second.push_back(' ');
				got->second.append(chinese);
//				file << got->first << " " << got->second << endl;
			}
		}
		memset(buffer, 0, BUFFER_SIZE);
	}
//	file << "--------------------------------------" << endl;
	for (auto& it : map) 
		file << it.first << " " << it.second << endl;
}
