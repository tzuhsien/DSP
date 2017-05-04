#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmm.h"

		
double** CreateDoubleArray(int N, int T)
{
	double** array = (double**)malloc(N * sizeof(double*));
	for (int i = 0; i < N; i++) {
		array[i] = (double*)malloc(T * sizeof(double));
		memset(array[i], 0, T * sizeof(double));
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
			memset(array[i][j], 0, N * sizeof(double));
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
void CalculateAlpha(HMM* hmm, double** alpha, char* seq, int len)
{
	for (int i = 0; i < hmm->state_num; i++) {
		alpha[i][0] = hmm->initial[i] * hmm->observation[seq[0] - 'A'][i];
	}
	for (int t = 0; t < len - 1; t++) {
		for (int j = 0; j < hmm->state_num; j++) {
			alpha[j][t + 1] = 0;
			for (int i = 0; i < hmm->state_num; i++) {
				alpha[j][t + 1] += alpha[i][t] * hmm->transition[i][j];
			}
			alpha[j][t + 1] *= hmm->observation[seq[t + 1] - 'A'][j];
		}
	}
}
void CalculateBeta(HMM* hmm, double** beta, char* seq, int len)
{
	for (int i = 0; i < hmm->state_num; i++) {
		beta[i][len - 1] = 1;
	}
	for (int t = len - 2; t >= 0; t--) {
		for (int i = 0; i < hmm->state_num; i++) {
			beta[i][t] = 0;
			for (int j = 0; j < hmm->state_num; j++) {
				beta[i][t] += hmm->transition[i][j] * hmm->observation[seq[t + 1] - 'A'][j] * beta[j][t + 1];
			}
		}
	}
}
void CalculateGamma(HMM* hmm, double** gamma, double** alpha, double** beta, char* seq, int len)
{
	for (int t = 0; t < len; t++) {
		double denominator = 0;
		for (int i = 0; i < hmm->state_num; i++) {
			denominator += alpha[i][t] * beta[i][t];
		}

		for (int i = 0; i < hmm->state_num; i++) {
			gamma[i][t] = alpha[i][t] * beta[i][t] / denominator;
		}
	}
	return;
}
void CalculateEpsilon(HMM* hmm, double*** epsilon, double** alpha, double** beta, char* seq, int len)
{
	for (int t = 0; t < len - 1; t++) {
		double denominator = 0;
		for (int i = 0; i < hmm->state_num; i++) {
			for (int j = 0; j < hmm->state_num; j++) {
				denominator += alpha[i][t] * hmm->transition[i][j] * hmm->observation[seq[t + 1] - 'A'][j] * beta[j][t + 1];
			}
		}
		for (int i = 0; i < hmm->state_num; i++) {
			for (int j = 0; j < hmm->state_num; j++) {
				epsilon[t][i][j] = alpha[i][t] * hmm->transition[i][j] * hmm->observation[seq[t + 1] - 'A'][j] * beta[j][t + 1] / denominator;
			}
		}
	}
	return;
}
	
void CalculatePi(HMM* hmm, double* pi, double** gamma)
{
	for (int i = 0; i < hmm->state_num; i++) {
		pi[i] += gamma[i][0];
	}
}
void CalculateSigmaGamma(HMM* hmm, double** sigma_gamma, double** gamma, char* seq, int len, double* sigma_gamma_T)
{
	for (int i = 0; i < hmm->state_num; i++) {
		for (int t = 0; t < len; t++) {
			sigma_gamma[i][seq[t] - 'A'] += gamma[i][t];
		}
		sigma_gamma_T[i] += gamma[i][len - 1];
	}
}
void CalculateSigmaEpsilon(HMM* hmm, double** sigma_epsilon, double*** epsilon, int len)
{
	for (int i = 0; i < hmm->state_num; i++) {
		for (int j = 0; j < hmm->state_num; j++) {
			for (int t = 0; t < len - 1; t++) {
				sigma_epsilon[i][j] += epsilon[t][i][j];
			}
		}
	}
}
void Training(HMM* hmm, double* pi, double** sigma_gamma, double** sigma_epsilon, int N, double* sigma_gamma_T) 
{
	for (int i = 0; i < hmm->state_num; i++) {
		hmm->initial[i] = pi[i] / (double)N;
	}

	double gamma[hmm->state_num];
	for (int i = 0; i < hmm->state_num; i++) {
		gamma[i] = 0;
		for (int j = 0; j < hmm->observ_num; j++) {
			gamma[i] += sigma_gamma[i][j];
		}
	}

	for (int i = 0; i < hmm->state_num; i++) {
		for (int j = 0; j < hmm->state_num; j++) {
			hmm->transition[i][j] = sigma_epsilon[i][j] / (gamma[i] - sigma_gamma_T[i]);
		}
	}

	for (int i = 0; i < hmm->state_num; i++) {
		for (int j = 0; j < hmm->observ_num; j++) {
			hmm->observation[j][i] = sigma_gamma[i][j] / gamma[i];
		}
	}
}
int main(int argc,char* argv[])
{
	HMM hmm;
	FILE* fp;
	int times = atoi(argv[1]);
	loadHMM(&hmm, argv[2]);
	for (int i = 0; i < times; i++) {
		fp = open_or_die(argv[3], "r");
		char seq[MAX_SEQ] = {0};
		double pi[hmm.state_num];
		memset(pi, 0, hmm.state_num * sizeof(double));
		double** sigma_gamma = CreateDoubleArray(hmm.state_num, hmm.observ_num);
		double** sigma_epsilon = CreateDoubleArray(hmm.state_num, hmm.state_num);
		double sigma_gamma_T[hmm.state_num];
		memset(sigma_gamma_T, 0, hmm.state_num * sizeof(double));
		int N = 0;
		while (fscanf(fp, "%s", seq) > 0) {
			int len = strlen(seq);
			double **alpha = CreateDoubleArray(hmm.state_num, len);
			double **beta = CreateDoubleArray(hmm.state_num, len);
			double **gamma = CreateDoubleArray(hmm.state_num, len);
			double ***epsilon = CreateThreeArray(hmm.state_num, len);
			CalculateAlpha(&hmm, alpha, seq, len);
			CalculateBeta(&hmm, beta, seq, len);
			CalculateGamma(&hmm, gamma, alpha, beta, seq, len);
			CalculateEpsilon(&hmm, epsilon, alpha, beta, seq, len);
			CalculatePi(&hmm, pi, gamma);
			CalculateSigmaGamma(&hmm, sigma_gamma, gamma, seq, len, sigma_gamma_T);
			CalculateSigmaEpsilon(&hmm, sigma_epsilon, epsilon, len);
			FreeDoubleArray(alpha, hmm.state_num);
			FreeDoubleArray(beta, hmm.state_num);
			FreeDoubleArray(gamma, hmm.state_num);
			FreeThreeArray(epsilon, hmm.state_num, len);
			N++;
		}
		Training(&hmm, pi, sigma_gamma, sigma_epsilon, N, sigma_gamma_T);
		fclose(fp);	
		FreeDoubleArray(sigma_epsilon, hmm.state_num);
		FreeDoubleArray(sigma_gamma, hmm.state_num);
	}
	fp = open_or_die(argv[4], "w");
	dumpHMM(fp, &hmm);
	fclose(fp);
		return 0;
}
