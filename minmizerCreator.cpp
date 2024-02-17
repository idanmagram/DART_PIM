#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>


using namespace std;

#define A 0
#define C 1
#define G 2
#define T 3
#define K 3
#define W 4

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

void calculateReadLocation(string read,  vector<RowMinimizer*> &MinimizersListForCPU) {
    vector<Kmer*> minimizers = createMinimizers(read);
    int EdScore = -1;
    int minEdScore = MinimizersListForCPU[0]->refSeq.size();
    int address = 0;
    map<Kmer*,int> mapKmerToAdressInReq;

    for (Kmer* minimizer : minimizers) {
        for (const auto& row : MinimizersListForCPU) {
            if (row->min.kmerSeq == minimizer->kmerSeq) {
                EdScore = wagnerFischerAffineGap(read, row->refSeq, 1, 1, 1);
                if(EdScore < minEdScore) {
                    address = row->min.position;
                    minEdScore = EdScore;
                }
            }
        }
        if (EdScore != -1) {
            mapKmerToAdressInReq[minimizer] = address;
        }
        EdScore = -1;
    }
}

int main() {
    createMinimizers("231032101233101");
    wagnerFischerAffineGap("231", "231", 1, 1, 1);
    vector<RowMinimizer*> MinimizersListForCPU;
    Kmer min1 = Kmer(0,"123");
    Kmer min2 = Kmer(6,"231");
    Kmer min3 = Kmer(8,"023");
    MinimizersListForCPU.push_back(new RowMinimizer(min1, "0231032221233101"));
    MinimizersListForCPU.push_back(new RowMinimizer(min2, "222032101233101"));
    MinimizersListForCPU.push_back(new RowMinimizer(min3, "231032101233221"));

    calculateReadLocation("231032101233101", MinimizersListForCPU);

    return 0;
}
