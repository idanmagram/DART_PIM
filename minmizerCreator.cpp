#include <iostream>
#include <vector>
#include <algorithm>

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

class Kmer {
public:
    int position;
    vector<int> kmerSeq;

    Kmer(int pos, const vector<int> &kmerSeq) : position(pos), kmerSeq(kmerSeq) {}

};


bool isLexBefore(const vector<int> &seq1, const vector<int> &seq2);
bool compareMinimizers(const Kmer* a, const Kmer* b);


class Heap {
public:
    std::vector<Kmer*> elements;

    void heapify(size_t index) {
        size_t size = elements.size();
        size_t largest = index;
        size_t left = 2 * index + 1;
        size_t right = 2 * index + 2;

        if (left < size && isLexBefore(elements[left]->kmerSeq, elements[largest]->kmerSeq)) {
            largest = left;
        }

        if (right < size && isLexBefore(elements[right]->kmerSeq, elements[largest]->kmerSeq)) {
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
            if (!isLexBefore(elements[index]->kmerSeq, elements[parent]->kmerSeq)) {
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
        return (elements[0]);
    }

};


vector<Kmer*> createMinimizers(const vector<int> &seq) {
    int i;
    vector<Kmer*> heapKmers;
    vector<Kmer*> kmersInOrder;
    Kmer* minMinimizer;
    vector<Kmer*> outMinimizers;

    for(i = 0; (i < seq.size()) && (i < W); i++) {
        std::vector<int> subVector(seq.begin() + i, seq.begin() + i + K);

        minMinimizer = new Kmer(i, subVector);
        heapKmers.push_back(minMinimizer);
        kmersInOrder.push_back(minMinimizer);
    }
    Heap heapMin(heapKmers);

    for(; (i <= seq.size() - K + 1); i++) {
        auto currentMin = heapMin.minElement();
        if (outMinimizers.empty()) {
            auto minElement = reinterpret_cast<Kmer *const>(heapMin.minElement());
            outMinimizers.push_back(new Kmer(minElement->position, minElement->kmerSeq));
        }
        else if(currentMin->position != outMinimizers.back()->position) {
                outMinimizers.push_back(reinterpret_cast<Kmer *const>(heapMin.minElement()));
        }
        vector<int> subVector(seq.begin() + i, seq.begin() + i + K);
        auto newMinimizer = new Kmer(i, subVector);
        kmersInOrder.push_back(newMinimizer);
        auto elem1 = kmersInOrder[kmersInOrder.size() - W - 1];
        auto iteratorToDelete = std::find(heapMin.elements.begin(), heapMin.elements.end(), elem1);
        heapMin.elements.erase(iteratorToDelete);

        minMinimizer = reinterpret_cast<Kmer*>(heapMin.minElement());
        heapMin.insert(newMinimizer);

        if(isLexBefore(newMinimizer->kmerSeq, minMinimizer->kmerSeq)) {
            minMinimizer = newMinimizer;
        }
    }
    return outMinimizers;
}

bool isLexBefore(const vector<int> &seq1, const vector<int> &seq2) {
    for (int i = 0; i < seq1.size(); i++) {
        if(seq1[i] < seq2[i]) {
            return true;
        }
        else if (seq1[i] > seq2[i]) {
            return false;
        }
    }
    return false;
}

bool compareMinimizers(const Kmer* a, const Kmer* b) {
    return !isLexBefore(a->kmerSeq, b->kmerSeq);
}

int main() {
    createMinimizers({2,3,1,0,3,2,1,0,1,2,3,3,1,0,1});
    return 0;
}

