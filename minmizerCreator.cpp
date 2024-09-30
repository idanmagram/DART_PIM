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
    string sub  = "002301300302032130100103230320011332000010300300032000130231223103000030313103033030320330103310303220032311020010222002010223001003010302033313300221322022";
    int i = 0;
    int res = 0;
    for (i = 0; i < 100; i++) {
        res = wagnerFischerAffineGap2(read, sub, &score, false, 1, 1, 1);
        cout << score;
    }
    return 0;

/*
    ifstream readsFile;
    bool readsFileOpen = false;
    readsFile = ifstream(argv[2]);
    readsFileOpen = readsFile.is_open();
    if(!readsFileOpen){
        std::cout << "ERROR: Can't open file " << string(argv[2]) << endl;
        return 1;
    }

    getReadsFromFile(readsFile);
*/
}
/*
void getReadsFromFile(ifstream& readsFile){
    cout << "start getReadsFromFile" << endl;
    string line;
    //skip first line
    if(!getline(readsFile, line)){
        std::cout << "MSG: Reads file is empty." << line << endl;
        return;
    }
    while(getline(readsFile, line)){
        //convertSeq2Nums(line), The conversion is after find_minimizers because the function gets read of letters
        cout << "idan magram ";
        //skip two line
    }
}
*/

int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score, bool backtracing, int wop, int wex, int wsub) {
    const int n = S1.size();
    const int m = S2.size();
    vector<int> contenders;
    int minCon = 0;
    string seq1_align = "";
    string seq2_align = "";

    int D[REF_SUB_SEQUENCE_LENGTH+1][READ_LENGTH+1] = {0};
    int M1[REF_SUB_SEQUENCE_LENGTH+1][READ_LENGTH+1] = {0};
    int M2[REF_SUB_SEQUENCE_LENGTH+1][READ_LENGTH+1] = {0};

    const int max_gap = MAX_GAP;
    const int max_gap_penalty = max_gap + wop;

    // Optimized initialization for both first row and first column
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
