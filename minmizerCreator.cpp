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

class Kmer {
public:
    int position;
    string kmerSeq;

    Kmer(int pos, const string &kmerSeq) : position(pos), kmerSeq(kmerSeq) {}
};

class RowMinimizer {
public:
    Kmer min;
    string refSeq;

    RowMinimizer(Kmer min, string refSeq) : min(min), refSeq(refSeq) {}
};


bool compareMinimizers(const Kmer* a, const Kmer* b);
int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score,  bool backtraching, int wop, int wex, int wsub);


class Heap {
public:
    std::vector<Kmer*> elements;

    void heapify(size_t index) {
        size_t size = elements.size();
        size_t largest = index;
        size_t left = 2 * index + 1;
        size_t right = 2 * index + 2;

        if (left < size && (elements[left]->kmerSeq < elements[largest]->kmerSeq)) {
            largest = left;
        }

        if (right < size && (elements[right]->kmerSeq < elements[largest]->kmerSeq)) {
            largest = right;
        }

        if (largest != index) {
            std::swap(elements[index], elements[largest]);
            heapify(largest);
        }
    }

    Heap(const std::vector<Kmer*>& values) : elements(values) {
        std::make_heap(elements.begin(), elements.end(), compareMinimizers);
    }


    void insert(Kmer* value) {
        elements.push_back(value);
        size_t index = elements.size() - 1;

        while (index > 0) {
            size_t parent = (index - 1) / 2;
            if (elements[index]->kmerSeq > elements[parent]->kmerSeq) {
                break;
            }

            std::swap(elements[index], elements[parent]);
            index = parent;
        }
        heapify(0);
    }

    void eraseLastMin() {
        elements.pop_back();
        heapify(0);
    }

    Kmer* minElement() {
        return elements[0];
    }
};

std::vector<Kmer*> createMinimizers(const string &seq) {
    int i;
    std::vector<Kmer*> heapKmers;
    std::vector<Kmer*> kmersInOrder;
    Kmer* minMinimizer;
    std::vector<Kmer*> outMinimizers;

    for(i = 0; (i < seq.size()) && (i < W); i++) {
        std::string subString = seq.substr(i, K);

        minMinimizer = new Kmer(i, subString);
        heapKmers.push_back(minMinimizer);
        kmersInOrder.push_back(minMinimizer);
    }
    Heap heapMin(heapKmers);

    for(; (i <= seq.size() - K + 1); i++) {
        auto currentMin = heapMin.minElement();
        if (outMinimizers.empty()) {
            outMinimizers.push_back(new Kmer(currentMin->position, currentMin->kmerSeq));
        }
        else if(currentMin->position != outMinimizers.back()->position) {
            outMinimizers.push_back(new Kmer(currentMin->position, currentMin->kmerSeq));
        }
        std::string subString = seq.substr(i, K);
        auto newMinimizer = new Kmer(i, subString);
        kmersInOrder.push_back(newMinimizer);
        auto elem1 = kmersInOrder[kmersInOrder.size() - W - 1];
        auto iteratorToDelete = std::find(heapMin.elements.begin(), heapMin.elements.end(), elem1);
        heapMin.elements.erase(iteratorToDelete);

        heapMin.insert(newMinimizer);
    }
    return outMinimizers;
}

bool compareMinimizers(const Kmer* a, const Kmer* b) {
    return b->kmerSeq < a->kmerSeq;
}

int wagnerFischerAffineGap(const string& S1, const string& S2, int wop, int wex, int wsub) {
    int n = S1.size();
    int m = S2.size();

    // Initialize matrices D, M1, and M2 with appropriate dimensions
    vector<vector<int>> D(n + 1, vector<int>(m + 1, 0));
    vector<vector<int>> M1(n + 1, vector<int>(m + 1, 0));
    vector<vector<int>> M2(n + 1, vector<int>(m + 1, 0));

    // Initialize matrices with appropriate values
    for (int i = 1; i <= n; ++i) {
        D[i][0] = i * wex;
        M1[i][0] = i * wex;
    }
    for (int j = 1; j <= m; ++j) {
        D[0][j] = j * wex;
        M2[0][j] = j * wex;
    }

    // Fill in the matrices using dynamic programming
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            M1[i][j] = min(M1[i - 1][j] + wex, D[i - 1][j] + wop + wex);
            M2[i][j] = min(M2[i][j - 1] + wex, D[i][j - 1] + wop + wex);

            if (S1[i - 1] == S2[j - 1]) {
                D[i][j] = D[i - 1][j - 1];
            } else {
                D[i][j] = min({M1[i][j], M2[i][j], D[i - 1][j - 1] + wsub});
            }

        }
    }
    // Return the optimal alignment score
    return D[n][m];
}
string backtrackingWF(const string& S1, const string& S2) {
    string seq1_align = "";
    string seq2_align = "";
    int n = S1.size();
    int m = S2.size();
    int wex = 1;
    int wop = 1;
    int wsub = 1;
    vector<int> contenders;
    int minCon = 0;

    vector<vector<int>> D(n + 1, vector<int>(m + 1, 0));
    vector<vector<int>> M1(n + 1, vector<int>(m + 1, 0));
    vector<vector<int>> M2(n + 1, vector<int>(m + 1, 0));

    // Initialize matrices with appropriate values
    for (int i = 1; i <= n; ++i) {
        D[i][0] = i * wex;
        M1[i][0] = i * wex;
    }
    for (int j = 1; j <= m; ++j) {
        D[0][j] = j * wex;
        M2[0][j] = j * wex;
    }

    // Fill in the matrices using dynamic programming
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            M1[i][j] = min((M1[i - 1][j] + wex), (D[i - 1][j] + wop + wex));
            M2[i][j] = min(M2[i][j - 1] + wex, D[i][j - 1] + wop + wex);

            if (S1[i - 1] == S2[j - 1]) {
                D[i][j] = D[i - 1][j - 1];
            } else {
                D[i][j] = min({M1[i][j], M2[i][j], D[i - 1][j - 1] + wsub});
            }
        }
    }

    int i = m, j = n;
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

    return seq1_align;
}

int calculateReadLocation(string read,  vector<RowMinimizer*> &MinimizersListForCPU, int eth) {
    vector<Kmer*> minimizers = createMinimizers(read);
    int EdScore = -1;
    int minEdScore = MinimizersListForCPU[0]->refSeq.size();
    int address = 0;
    map<Kmer*,int> mapKmerToAdressInReq;
    string minSubRef;

    for (Kmer* minimizer : minimizers) {
        for (const auto& row : MinimizersListForCPU) {
            if (row->min.kmerSeq == minimizer->kmerSeq) {
                string subRef = row->refSeq.substr(N - K + minimizer->position - eth, N + 2*eth);
                EdScore = wagnerFischerAffineGap(read, subRef, 1, 1, 1);
                if(EdScore < minEdScore) {
                    address = row->min.position; // location of minimizer in the genome
                    minEdScore = EdScore;
                    minSubRef = subRef;
                }
            }
        }
    }
    return address;
}

int main() {
    //createMinimizers("10213230312320233");

    int score = 0;
    string read = "002011023132211001032232000111103133301300000301000300330213222323223223222101132300331102130131222002320021022020031010320011132202210202203210232023";
    string sub  = "002011023132211001032232000111103133301300000301000300330213222323223223222101132300331102130131222002320021022020031010320011132202210202203210232023020220";
    //string sub  = "211133202121211002111102210333230221022333113131112102222120323013122120321310212231322211102111323021131301103222131021321021321101132001032210123221020020";

    //string read = "111123";
    //string sub  = "1111230000";

    //string sub= "0010132110210202021323222330101000222230210222330003323313212010313033321132023001112232112221201233112123221101300311330222111030030102001000232";
    //vector<vector<int>> D(n + 1, vector<int>(m + 1, 0));
    for (int i = 0; i < 1; i++) {
        wagnerFischerAffineGap2(read, sub, &score, true, 1, 1, 1);
        cout << "score is: " << score;
    }
    //backtrackingWF("AGGCCTA", "TAGCTTA");
    //vector<RowMinimizer*> MinimizersListForCPU;
    //Kmer min1 = Kmer(0,"123");
    //Kmer min2 = Kmer(6,"231");
    //Kmer min3 = Kmer(8,"023");
    //MinimizersListForCPU.push_back(new RowMinimizer(min1, "0231032123233101"));
    //MinimizersListForCPU.push_back(new RowMinimizer(min2, "222032123133101"));
    //MinimizersListForCPU.push_back(new RowMinimizer(min3, "231032100323221"));
    //calculateReadLocation("231032101233101", MinimizersListForCPU, 2);

    return 0;
}

int wagnerFischerAffineGap2(const string& S1, const string& S2, int* score,  bool backtraching, int wop, int wex, int wsub) {
    int n = S1.size();
    int m = S2.size();
    vector<int> contenders;
    int minCon = 0;
    string seq1_align = "";
    string seq2_align = "";

    // Initialize matrices D, M1, and M2 with appropriate dimensions
    vector<vector<int>> D(n + 1, vector<int>(m + 1, 0));
    vector<vector<int>> M1(n + 1, vector<int>(m + 1, 0));
    vector<vector<int>> M2(n + 1, vector<int>(m + 1, 0));

    // Initialize matrices with appropriate values
    for (int i = 1; i <= n; ++i) {
        D[i][0] = i * wex;
        M1[i][0] = i * wex;
    }
    for (int j = 1; j <= m; ++j) {
        D[0][j] = j * wex;
        M2[0][j] = j * wex;
    }

    // Fill in the matrices using dynamic programming
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (abs(i-j) < 8) {
                M1[i][j] = min(M1[i - 1][j] + wex, D[i - 1][j] + wop + wex);
                M2[i][j] = min(M2[i][j - 1] + wex, D[i][j - 1] + wop + wex);
                if (S1[i - 1] == S2[j - 1]) {
                    D[i][j] = D[i - 1][j - 1];
                } else {
                    D[i][j] = min({M1[i][j], M2[i][j], D[i - 1][j - 1] + wsub});
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
