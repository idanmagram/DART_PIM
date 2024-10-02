#include <iostream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <fstream>


using namespace std;


#define MAX 10

#define READ_LENGTH                150
#define ERROR_THRESHOLD            3
#define REF_SUB_SEQUENCE_LENGTH    READ_LENGTH + 2 * ERROR_THRESHOLD
#define MAX_GAP                    2 * ERROR_THRESHOLD

int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score,  bool backtraching, int wop, int wex, int wsub);
void getReadsFromFile(ifstream& readsFile);

int main(int argc, char* argv[]) {

    int score = 0;
    string read = "300302032130100103230320011332000010300300032000130231223103000030313103033030320330103310303220032311020010222002010223001003010302033313300221322022";
    string sub = "002301300302032130100103230320011332000010300300032000130231223103000030313103033030320330103310303220032311020010222002010223001003010302033313300221322022";
    //string read = "mat";
    //string sub = "mat" ;
    int i = 0;
    int res = 0;
    for (i = 0; i < 1; i++) {
        res = wagnerFischerAffineGap2(read, sub, &score, false, 1, 1, 1);
        cout << res;
    }
    return 0;
}


int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score, bool backtracing, int wop, int wex, int wsub) {
    const int n = S1.size();
    const int m = S2.size();
    vector<int> contenders;
    int minCon = 0;
    string seq1_align = "";
    string seq2_align = "";


    //int D[REF_SUB_SEQUENCE_LENGTH+1][READ_LENGTH+1] = {0};

    int M1[REF_SUB_SEQUENCE_LENGTH+1][READ_LENGTH+1] = {0};
    int M2[REF_SUB_SEQUENCE_LENGTH+1][READ_LENGTH+1] = {0};
    int D[READ_LENGTH][REF_SUB_SEQUENCE_LENGTH] = {
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155},  // First row
            {1},  // First column for row 1
            {2},  // First column for row 2
            {3},  // First column for row 3
            {4},  // First column for row 4
            {5},  // First column for row 5
            {6},  // First column for row 6
            {7},  // First column for row 7
            {8},  // First column for row 8
            {9},  // First column for row 9
            {10}, // First column for row 10
            {11}, // First column for row 11
            {12}, // First column for row 12
            {13}, // First column for row 13
            {14}, // First column for row 14
            {15}, // First column for row 15
            {16}, // First column for row 16
            {17}, // First column for row 17
            {18}, // First column for row 18
            {19}, // First column for row 19
            {20}, // First column for row 20
            {21}, // First column for row 21
            {22}, // First column for row 22
            {23}, // First column for row 23
            {24}, // First column for row 24
            {25}, // First column for row 25
            {26}, // First column for row 26
            {27}, // First column for row 27
            {28}, // First column for row 28
            {29}, // First column for row 29
            {30}, // First column for row 30
            {31}, // First column for row 31
            {32}, // First column for row 32
            {33}, // First column for row 33
            {34}, // First column for row 34
            {35}, // First column for row 35
            {36}, // First column for row 36
            {37}, // First column for row 37
            {38}, // First column for row 38
            {39}, // First column for row 39
            {40}, // First column for row 40
            {41}, // First column for row 41
            {42}, // First column for row 42
            {43}, // First column for row 43
            {44}, // First column for row 44
            {45}, // First column for row 45
            {46}, // First column for row 46
            {47}, // First column for row 47
            {48}, // First column for row 48
            {49}, // First column for row 49
            {50}, // First column for row 50
            {51}, // First column for row 51
            {52}, // First column for row 52
            {53}, // First column for row 53
            {54}, // First column for row 54
            {55}, // First column for row 55
            {56}, // First column for row 56
            {57}, // First column for row 57
            {58}, // First column for row 58
            {59}, // First column for row 59
            {60}, // First column for row 60
            {61}, // First column for row 61
            {62}, // First column for row 62
            {63}, // First column for row 63
            {64}, // First column for row 64
            {65}, // First column for row 65
            {66}, // First column for row 66
            {67}, // First column for row 67
            {68}, // First column for row 68
            {69}, // First column for row 69
            {70}, // First column for row 70
            {71}, // First column for row 71
            {72}, // First column for row 72
            {73}, // First column for row 73
            {74}, // First column for row 74
            {75}, // First column for row 75
            {76}, // First column for row 76
            {77}, // First column for row 77
            {78}, // First column for row 78
            {79}, // First column for row 79
            {80}, // First column for row 80
            {81}, // First column for row 81
            {82}, // First column for row 82
            {83}, // First column for row 83
            {84}, // First column for row 84
            {85}, // First column for row 85
            {86}, // First column for row 86
            {87}, // First column for row 87
            {88}, // First column for row 88
            {89}, // First column for row 89
            {90}, // First column for row 90
            {91}, // First column for row 91
            {92}, // First column for row 92
            {93}, // First column for row 93
            {94}, // First column for row 94
            {95}, // First column for row 95
            {96}, // First column for row 96
            {97}, // First column for row 97
            {98}, // First column for row 98
            {99}, // First column for row 99
            {100}, // First column for row 100
            {101}, // First column for row 101
            {102}, // First column for row 102
            {103}, // First column for row 103
            {104}, // First column for row 104
            {105}, // First column for row 105
            {106}, // First column for row 106
            {107}, // First column for row 107
            {108}, // First column for row 108
            {109}, // First column for row 109
            {110}, // First column for row 110
            {111}, // First column for row 111
            {112}, // First column for row 112
            {113}, // First column for row 113
            {114}, // First column for row 114
            {115}, // First column for row 115
            {116}, // First column for row 116
            {117}, // First column for row 117
            {118}, // First column for row 118
            {119}, // First column for row 119
            {120}, // First column for row 120
            {121}, // First column for row 121
            {122}, // First column for row 122
            {123}, // First column for row 123
            {124}, // First column for row 124
            {125}, // First column for row 125
            {126}, // First column for row 126
            {127}, // First column for row 127
            {128}, // First column for row 128
            {129}, // First column for row 129
            {130}, // First column for row 130
            {131}, // First column for row 131
            {132}, // First column for row 132
            {133}, // First column for row 133
            {134}, // First column for row 134
            {135}, // First column for row 135
            {136}, // First column for row 136
            {137}, // First column for row 137
            {138}, // First column for row 138
            {139}, // First column for row 139
            {140}, // First column for row 140
            {141}, // First column for row 141
            {142}, // First column for row 142
            {143}, // First column for row 143
            {144}, // First column for row 144
            {145}, // First column for row 145
            {146}, // First column for row 146
            {147}, // First column for row 147
            {148}, // First column for row 148
            {149}  // First column for row 149
    };




    const int max_gap = MAX_GAP;
    const int max_gap_penalty = max_gap + wop;

    // Optimized initialization for both first row and first column
    /*
    for (int i = 1; i <= max(n, m); ++i) {

        if (i <= n) { // Initialize first column (gaps in S1 aligned to empty S2)
            D[i][0] = wop + (i - 1) * wex; // Gap penalty for aligning S1 to empty S2
            M1[i][0] = D[i][0];            // Same gap penalty for M1
            M2[i][0] = max_gap_penalty;     // M2 irrelevant here (no gaps in S1)
        }
        if (i <= m) { // Initialize first row (gaps in S2 aligned to empty S1)
            D[0][i] = wop + (i - 1) * wex; // Gap penalty for aligning S2 to empty S1
            M1[0][i] = max_gap_penalty;     // M1 irrelevant here (no gaps in S2)
            M2[0][i] = D[0][i];            // Same gap penalty for M2
        }
    }
     */

    // Fill the DP tables using dynamic programming
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (abs(i - j) <= max_gap) {
                if (abs(i - 1 - j) > max_gap) {
                    M1[i - 1][j] = max_gap_penalty;
                    D[i - 1][j] = max_gap_penalty;
                }
                if (abs(i - (j - 1)) > max_gap) {
                    M2[i][j - 1] = max_gap_penalty;
                    D[i][j - 1] = max_gap_penalty;
                }

                M1[i][j] = min(M1[i - 1][j] + wex, D[i - 1][j] + wop + wex);
                M2[i][j] = min(M2[i][j - 1] + wex, D[i][j - 1] + wop + wex);

                if (S1[i - 1] == S2[j - 1])
                    D[i][j] = D[i - 1][j - 1];
                else
                    D[i][j] = min({M1[i][j], M2[i][j], D[i - 1][j - 1] + wsub});
            }
        }
    }

    *score = D[n][m]; // Final alignment score
    return *score;
}