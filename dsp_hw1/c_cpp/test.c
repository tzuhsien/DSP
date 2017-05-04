#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmm.h"

#define MAX_NUM 10

double DoViterbi(HMM* hmm, char* seq)
{
	int len = strlen(seq);
	double delta[hmm->state_num][len];
	for (int i = 0; i < hmm->state_num; i++) {
		delta[i][0] = hmm->initial[i] * hmm->observation[seq[0] - 'A'][i];
	}
	for (int t = 0; t < len - 1; t++) {
		for (int j = 0; j < hmm->state_num; j++) {
			double max_probability = 0;
			for (int i = 0; i < hmm->state_num; i++) {
				double probability = delta[i][t] * hmm->transition[i][j];
				if (probability > max_probability)
					max_probability = probability;
			}
			delta[j][t + 1] = max_probability * hmm->observation[seq[t + 1] - 'A'][j];
		}
	}
	double max_probability = 0;
	for (int i = 0; i < hmm->state_num; i++) {
//		printf("%g\n", delta[i][len - 1]);
		if (delta[i][len - 1] > max_probability)
			max_probability = delta[i][len - 1];
	}
	return max_probability;
}

int main(int argc, char* argv[])
{
	HMM hmm[MAX_NUM];
	int count = load_models(argv[1], hmm, MAX_NUM);
	FILE* read_fp = open_or_die(argv[2], "r" );
	FILE* write_fp = open_or_die(argv[3], "w");
	char seq[MAX_SEQ] = {0};
	while (fscanf(read_fp, "%s", seq) > 0) {	
		int max_index;
		double max_probability = 0;
		for (int i = 0; i < count; i++) {
			double probability = DoViterbi(&hmm[i], seq);
//			printf("%g\n", probability);
			if (probability > max_probability) {
				max_probability = probability;
				max_index = i;
			}
		}
		memset(seq, 0, MAX_SEQ);
		fprintf(write_fp, "model_%02d.txt %e\n", max_index + 1, max_probability);
	}
	return 0;
}
