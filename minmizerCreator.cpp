#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>


using namespace std;

#define A 0
#define C 1
#define G 2
#define T 3
#define K 12
#define W 4
#define N 7

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#define READ_LENGTH                150
#define ERROR_THRESHOLD            3
#define REF_SUB_SEQUENCE_LENGTH    READ_LENGTH + 2 * ERROR_THRESHOLD

int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score,  bool backtraching, int wop, int wex, int wsub);

int computeMinimum3(int a, int b, int c) {
    int min_val = a;
    if (b < min_val) {
        min_val = b;
    }
    if (c < min_val) {
        min_val = c;
    }
    return min_val;
}

int main() {
    int score = 0;
    string read = "002011023132211001032232000111103133301300000301000300330213222323223223222101132300331102130131222002320021022020031010320011132202210202203210232023";
    string sub  = "002011023132211001032232000111103133301300000301000300330213222323223223222101132300331102130131222002320021022020031010320011132202210202203210232023020220";
    for (int i = 0; i < 1; i++) {
        wagnerFischerAffineGap2(read, sub, &score, false, 1, 1, 1);
        cout << "score is: " << score;
    }
    return 0;
}

int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score,  bool backtraching, int wop, int wex, int wsub) {
    int n = S1.size();
    int m = S2.size();
    vector<int> contenders;
    int minCon = 0;
    string seq1_align = "";
    string seq2_align = "";

    int D[REF_SUB_SEQUENCE_LENGTH][READ_LENGTH] = {0};
    int M1[REF_SUB_SEQUENCE_LENGTH][READ_LENGTH] = {0};
    int M2[REF_SUB_SEQUENCE_LENGTH][READ_LENGTH] = {0};


    // Initialize matrices with appropriate values
    for (int i = 1; i <= n; ++i) {
        D[i][0] = i * wex;
        M1[i][0] = i * wex;
    }
    for (int j = 1; j <= m; ++j) {
        D[0][j] = j * wex;
        M2[0][j] = j * wex;
    }
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (abs(i - j) > 5) {
                D[i][j] = 10;
                M1[i][j] = 10;
                M2[i][j] = 10;
            }
        }
    }

    // Fill in the matrices using dynamic programming
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (abs(i-j) <= 6) {
                M1[i][j] = min(M1[i - 1][j] + wex, D[i - 1][j] + wop + wex);
                M2[i][j] = min(M2[i][j - 1] + wex, D[i][j - 1] + wop + wex);
                if (S1[i - 1] == S2[j - 1]) {
                    D[i][j] = D[i - 1][j - 1];
                } else {
                    D[i][j] = computeMinimum3(M1[i][j], M2[i][j], D[i - 1][j - 1] + wsub);
                    //D[i][j] = min({M1[i][j], M2[i][j], D[i - 1][j - 1] + wsub});
                }
            }
        }
    }
    // Return the optimal alignment score
    *score = D[n][m];

    // find the alignment of the read according to sub reference sequence
    if(backtraching) {
        int i = n, j = m;
        while (i > 0 && j > 0 && D[i][j] > 0) {
            contenders = {D[i - 1][j - 1], D[i - 1][j], D[i][j - 1]};
            minCon = *(std::min_element(contenders.begin(), contenders.end()));

            if (minCon == contenders[0]) {
                if (S1[i - 1] != S2[j - 1]) {
                    seq1_align += S2[i - 1];
                    seq2_align += S2[j - 1];
                } else {
                    seq1_align += S1[i - 1];
                    seq2_align += S2[j - 1];
                }
                i -= 1;
                j -= 1;
            } else if (minCon == contenders[1]) {
                seq1_align += S1[i - 1];
                seq2_align += "-";
                i -= 1;
            } else {
                seq1_align += "-";
                seq2_align += S2[j - 1];
                j -= 1;
            }

        }

        std::reverse(seq1_align.begin(), seq1_align.end());
        std::reverse(seq2_align.begin(), seq2_align.end());
    }

    return *score;
}

