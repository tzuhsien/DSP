#include <stdio.h>
#include <string.h>
#include "hmm.h"

		
double** CreateDoubleArray(int N, int T)
{
	double** array = (double**)malloc(N * sizeof(int*));
	for (int i = 0; i < N; i++) {
		array[i] = (double*)malloc(T * sizeof(double));
	}
	return array;
}

void FreeDoubleArray(double** array, int N)
{
	for (int i = 0; i < N; i++) {
		free(array[i]);
	}
	free(array);
	return;
}

double*** CreateThreeArray(int N, int T)
{
	double*** array = (double***)malloc(T * sizeof(double**));
	for (int i = 0; i < T; i++) {
		array[i] = (double**)malloc(N * sizeof(double*));
		for (int j = 0; j < N; j++) {
			array[i][j] = (double*)malloc(N * sizeof(double));
		}
	}
	return array;
}
void FreeThreeArray(double*** array, int N, int T)
{
	for (int i = 0; i < T; i++) {
		for (int j = 0; j < N; j++) {
			free(array[i][j]);
		}
		free(array[i]);
	}
	free(array);
}
void CalculateAlpha(HMM* hmm, double** alpha, char* seq, int N, int T)
{
	for (int i = 0; i < N; i++) {
		alpha[i][0] = hmm->initial[i] * hmm->observation[seq[0] - 'A'][i];
	}
	for (int t = 0; t < T - 1; t++) {
		for (int j = 0; j < N; j++) {
			alpha[j][t + 1] = 0;
			for (int i = 0; i < N; i++) {
				alpha[j][t + 1] += alpha[i][t] * hmm->transition[i][j];
			}
			alpha[j][t + 1] *= hmm->observation[seq[t + 1] - 'A'][j];
		}
	}
}
void CalculateBeta(HMM* hmm, double** beta, char* seq, int N, int T)
{
	for (int i = 0; i < N; i++) {
		beta[i][T - 1] = 1;
	}
	for (int t = T - 2; t >= 0; t--) {
		for (int i = 0; i < N; i++) {
			beta[i][t] = 0;
			for (int j = 0; j < N; j++) {
				beta[i][t] += hmm->transition[i][j] * hmm->observation[seq[t + 1] - 'A'][j] * beta[j][t + 1]
			}
		}
	}
}
void CalculateGamma(HMM* hmm, double** gamma, double** alpha, double** beta, char* seq, int N, int T)
{
	for (int t = 0; t < T; t++) {
		double denominator = 0;
		for (int i = 0; i < N; i++) {
			denominator += alpha[i][t] * beta[i][t];
		}

		for (int i = 0; i < N; i++) {
			gamma[i][t] = alpha[i][t] * beta[i][t] / denominator;
		}
	}
	return;
}
void CalculateEpsilon(Hmm* hmm, double*** epsilon, double** alpha, double** beta, char* seq, int N, int T)
{
	for (int t = 0; t < T; t++) {
		double denominator = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				denominator += alpha[i][t] * hmm->transition[i][j] * hmm->observation[seq[t + 1] - 'A'][j] * beta[j][t + 1];
			}
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				epsilon[t][i][j] = alpha[i][t] * hmm->transition[i][j] * hmm->observation[seq[t + 1] - 'A'][j] * beta[j][t + 1] / denominator;
			}
		}
	}
	return;
}
	
void Training(HMM* hmm, double** alpha, double** beta, double** gamma, double*** epsilon, int N, int T)
{
	for (int i = 0; i < N; i++) {
		hmm->initial[i] = gamma[i][0];
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; n++) {
			double numerator = 0;
			double denominator = 0;
			for (int t = 0; t < T - 1; t++) {
				numerator += epsilon[t][i][j];
				denominator += gamma[i][t];
			}
			hmm->transition[i][j] = numerator / denominator;
		}
	}

	
}
int main(int argc,char* argv[])
{
	HMM hmm;
	FILE* fp;
	int times = atoi(argv[1]);
	loadHMM(&hmm, argv[2]);
	fp = open_or_die(argv[3], "r");
	char seq[MAX_SEQ] = {0};
	while (fscanf(fp, "%s", seq) > 0) {
		int N = hmm.state_num;
		int T = strlen(seq);
		double **alpha = CreateDoubleArray(N, T);
		double **beta = CreateDoubleArray(N, T);
		double **gamma = CreateDoubleArray(N, T);
		double ***epsilon = CreateThreeArray(N, T);
		CalculateAlpha(&hmm, alpha, seq, N, T);
		CalculateBeta(&hmm, beta, seq, N, T);
		CalculateGamma(&hmm, gamma, alpha, beta, seq, N, T);
		CalculateEpsilon(&hmm, epsilon, alpha, beta, seq, N, T);
		Training(&hmm, alpha, beta, gamma, epsilon, N, T);
	}
	return 0;
}
