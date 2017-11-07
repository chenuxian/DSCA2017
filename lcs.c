#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

char s1[101], s2[101], output[200000][101], temp[101]; // last char store \0
int LCS[101][101], len_s1, len_s2, same_len, count = 0;

int compare(const void*a, const void *b) {
    return(strcmp((char *)a, (char *)b));
}

void findAll(int index_i, int index_j, int curr_len) {
    char k;
    int i, j;
    if(curr_len == same_len) {
		for (i = 0; i < same_len; ++i) {
			output[count][i] = temp[i];
		}
		count++;
        return;
    }
    if(index_i == len_s1 || index_j == len_s2) {
        return;
	}

    for(i = index_i; i < len_s1; ++i) {
        k = s1[i];
        for(j = index_j; j < len_s2; ++j) {
            if (k == s2[j] && LCS[i+1][j+1] == curr_len + 1) {
                temp[curr_len] = k;
                findAll(i+1, j+1, curr_len+1);
            }
    	}
	}
}

int main() {
    int i,j;
    scanf("%s", s1);
    scanf("%s", s2);
    len_s1 = strlen(s1);
    len_s2 = strlen(s2);

    for (i = 0; i <= len_s1; ++i) {
        LCS[i][0] = 0;
	}
    for (i = 0; i <= len_s2; ++i) {
        LCS[0][i] = 0;
    }
    for (i = 1; i <= len_s1; ++i) {
        for (j = 1; j <= len_s2; ++j) {
            if (s1[i-1] == s2[j-1]) {
                LCS[i][j] = LCS[i-1][j-1] + 1;
			} else {
                LCS[i][j] = MAX(LCS[i-1][j], LCS[i][j-1]);
            }
		}
	}
    
	same_len = LCS[len_s1][len_s2];
    findAll(0, 0, 0);
	printf("%d %d\n", same_len, count);
	
	// sort
    qsort(output, count, sizeof output[0], compare);
    
	for(i = 0; i < count; ++i) {
        printf("%s\n", output[i]);
	}
    return 0;
}
