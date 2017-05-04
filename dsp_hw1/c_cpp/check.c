#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	FILE* fp1 = fopen(argv[1], "r");
	FILE* fp2 = fopen(argv[2], "r");
	char seq1[256] = {0};
	char seq2[256] = {0};
	char buffer[256];
	int right = 0;
	int num = 0;
	while (fscanf(fp1, "%s", seq1) > 0 && fscanf(fp2, "%s", seq2) > 0) {
		fscanf(fp1, "%s", buffer);
		if (!strcmp(seq1, seq2))
			right++;
		num++;
		memset(seq1, 0, 256);
		memset(seq2, 0, 256);
	}
	printf("%lf\n", (double)right / (double)num);
}
