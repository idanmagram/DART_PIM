#include <iostream>
#include <vector>
#include <algorithm>


using namespace std;

#define MAX 10

#define READ_LENGTH                150
#define ERROR_THRESHOLD            3
#define REF_SUB_SEQUENCE_LENGTH    READ_LENGTH + 2 * ERROR_THRESHOLD
#define MAX_GAP                    2 * ERROR_THRESHOLD

int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score,  bool backtraching, int wop, int wex, int wsub);

int main() {
    int score = 0;
    string read = "002011023132211001032232000111103133301300000301000300330213222323223223222101132300331102130131222002320021022020031010320011132202210202203210232023";
    string sub  = "002011023132211001032232000111103133301300000301000300330213222323223223222101132300331102130131222002320021022020031010320011132202210202203210232023020220";
    int i = 0;
    int res = 0;
    for (i = 0; i < 10; i++) {
        res = wagnerFischerAffineGap2(read, sub, &score, false, 1, 1, 1);
    }
    return res;
}

int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score,  bool backtraching, int wop, int wex, int wsub) {
    vector<int> contenders;
    int minCon = 0;
    string seq1_align = "";
    string seq2_align = "";
    const int n = S1.size();
    const int m = S2.size();

    // Initialize matrices with appropriate values
    int* D = new int[(n + 1) * (m + 1)]{};
    int* M1 = new int[(n + 1) * (m + 1)]{};
    int* M2 = new int[(n + 1) * (m + 1)]{};

    // Fill the DP tables using dynamic programming
    for (int i = 1; i <= n; ++i) {
        D[i * (m + 1)] = i;
        M1[i * (m + 1)] = i;
    }
    for (int j = 1; j <= m; ++j) {
        D[j] = j;
        M2[j] = j;
    }

    const int max_gap = MAX_GAP;
    const int max_gap_penalty = max_gap + wop;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            const int idx = i * (m + 1) + j;
            if (abs(i - j) > max_gap) {
                D[idx] = max_gap_penalty;
                M1[idx] = max_gap_penalty;
                M2[idx] = max_gap_penalty;
                continue;
            }

            const int idx_diag = (i - 1) * (m + 1) + (j - 1);
            const int idx_up = (i - 1) * (m + 1) + j;
            const int idx_left = i * (m + 1) + (j - 1);

            M1[idx] = min(M1[idx_up] + wex, D[idx_up] + wop + wex);
            M2[idx] = min(M2[idx_left] + wex, D[idx_left] + wop + wex);
            D[idx] = (S1[i - 1] == S2[j - 1]) ? D[idx_diag] : min({M1[idx], M2[idx], D[idx_diag] + wsub});
        }
    }

    *score = D[n * (m + 1) + m];

    delete[] D;
    delete[] M1;
    delete[] M2;

    // find the alignment of the read according to sub reference sequence
    /*
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
    */
    return *score;
}

