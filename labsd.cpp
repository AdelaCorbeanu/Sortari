#pragma GCC optimize("O3", "unroll-loops")

#include <vector>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <chrono>
using namespace std;

void count_sort10 (vector<long  long> &v, int p, int baza = 10) {
    vector<int> ap(baza, 0);
    vector<long long> aux(v.size());

    for (auto it : v) {
        ap[(it / p) % baza]++;
    }

    for (int i = 1; i < baza; ++i)
        ap[i] += ap[i - 1];        //retine pozitia pe care apare prima data numarul cu cifra i

    int n = (int)v.size() - 1;
    for (int i = n; i >= 0; --i) {
        long long index = (v[i] / p) % baza;  //ultima cifra
        aux[ap[index] - 1] = v[i];
        ap[index]--;
    }
    v = move(aux);
}

void radix_sort10 (vector<long long> &v, int baza = 10) {
    long long x = v[0];
    for (auto it : v) if (it > x) x = it;     //gasim maximul

    int cif = 0;
    while (x) {
        cif++;              //nr de cifre al maximului
        x /= baza;
    }

    int p = 1;
    for (int i = 0; i < cif; ++i) {
        count_sort10 (v, p, baza);
        p *= baza;
    }
}

void count_sort (vector<long long> &v, int p, int baza2) {
    vector<long long> aux(v.size());
    int baza = 1 << baza2;         //2^(baza2)
    vector<int> ap(baza, 0);

    for (auto it : v) {
        ap[(it >> p) & (baza - 1)]++;    //x % 2^n = x & (2^n - 1)
    }

    for (int i = 1; i < baza; ++i)
        ap[i] += ap[i - 1];        //retine pozitia pe care apare prima data numarul cu cifra i

    int n = (int)v.size() - 1;
    for (int i = n; i >= 0; --i) {
        int index = (v[i] >> p) & (baza - 1);
        aux[ap[index] - 1] = v[i];
        ap[index]--;
    }
    v = move(aux);
}

void radix_sort (vector<long long> &v, int baza2) {     //baza este 2^(baza2)
    long long x = v[0];
    for (auto it : v) if (it > x) x = it;     //gasim maximul

    int cif = 0;
    while (x) {
        cif++;
        x >>= baza2;     //shiftam la dreapta cu exponentul lui 2 (impartire)
    }

    int p = 0;
    for (int i = 0; i < cif; ++i) {
        count_sort (v, p, baza2);
        p += baza2;
    }
}



void merge(vector<long long>&a, vector<long long>&b, vector <long long> &v) {
    static vector<long long> c;
    c.resize(0);

    int i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) c.push_back(a[i++]);
        else c.push_back(b[j++]);
    }

    while (i < a.size())
        c.push_back(a[i++]);

    while (j < b.size())
        c.push_back(b[j++]);
    swap(v, c);

}

void merge_sort (vector<long long> &v) {
    if (v.size() > 1) {
        int size = (int)v.size() / 2;

        vector<long long> b(v.begin() + size, v.end());
        v.resize(size);
        merge_sort(v);
        merge_sort(b);

        merge(v, b, v);
    }
}

void merge2(vector<long long>&a, vector<long long>&b, vector <long long> &v) {
    static vector<long long> c;
    c.resize(0);

    auto i = a.begin();
    auto j = b.begin();

    while (i != a.end() && j != b.end()) {
        if (*i < *j) c.push_back(*i++);
        else c.push_back(*j++);
    }

    while (i != a.end()) c.push_back(*i++);
    while (j != b.end()) c.push_back(*j++);

    swap(v, c);
}

void merge_sort2 (vector<long long> &v) {
    if (v.size() > 1) {
        int size = v.size() / 2;

        vector<long long> b(v.begin() + size, v.end());
        v.resize(size);
        merge_sort2(v);
        merge_sort2(b);

        merge2(v, b, v);
    }
}



void shell_sort (vector<long long> &v) {
    int n = (int)v.size();
    for (int gap = n / 2; gap >= 1; gap /= 2) {
        for (int i = gap; i < n; ++i) {
            long long aux = v[i], j;
            for (j = i; j >= gap && v[j - gap] > aux; j -= gap)
                v[j] = v[j - gap];
            v[j] = aux;
        }
    }
}

void shell_sortC (vector<long long> &v) {
    vector <int> gaps = {1,4,10,23,57,132,301,701};
    int n = (int)v.size();

    for (int h = 7; h >= 0; h--) {
        int gap = gaps[h];
        for (int i = gap; i < n; ++i) {
            long long aux = v[i], j;

            for (j = i; j >= gap && v[j - gap] > aux; j -= gap)
                v[j] = v[j - gap];

            v[j] = aux;
        }
    }
}



void insertion_sort (vector<long long> &v) {
    int n = (int)v.size();

    for (int i = 1; i < n; ++i) {
        long long aux = v[i];
        int j = i - 1;

        while (aux < v[j] && j >= 0) {
            v[j + 1] = v[j];
            j--;
        }
        v[j + 1] = aux;
    }
}

void bucket_sort (vector<long long> &v) {
    int x = min((int)v.size() / 100, 10000);
    vector<long long> buckets[x];

    long long Max = 0;
    for (auto it : v)
        if (it > Max) Max = it;
    Max++;

    for (auto it : v)
        buckets[(it * x) / Max].push_back(it);

    for (auto & bucket : buckets)
        insertion_sort(bucket);

    int i = 0;
    for (const auto& bucket : buckets) {
        for (auto it : bucket) v[i++] = it;
    }
}



void heapify (vector<long long> &v, const int &n, const int &tata) {
    int Max = tata;
    int copil1 = 2 * tata + 1;
    int copil2 = 2 * tata + 2;

    if (copil1 < n && v[copil1] > v[Max]) Max = copil1;
    if (copil2 < n && v[copil2] > v[Max]) Max = copil2;

    if (Max != tata) {
        swap(v[tata], v[Max]);
        heapify(v, n, Max);
    }
}

void heap_sort (vector<long long> &v) {
    int n = (int)v.size();
    for (int i = n / 2 - 1; i >= 0; i--) heapify(v, n, i);

    for (int i = n - 1; i >= 0; i--) {
        swap(v[0], v[i]);
        heapify(v, i, 0);
    }
}



bool sortat_bine (vector<long long> &v, vector<long long> &a) {
    for (int i = 0; i < v.size(); ++i)
        if (v[i] != a[i]) return false;
    return true;
}


void testeaza_cu_nr_random () {
    int p[] = {1000, 1000000, 100000000};
    long long p2[] = {1000, 1000000000000, 1000000000000000000};
    for (long long Max : p2) {
        for (int n : p) {
            cout << "N = " << n << "  Max = " << Max << '\n';
            srand(unsigned(time(nullptr)));
            vector<long long> v(n);
            generate(v.begin(), v.end(), rand);
            for (auto &it: v)
                it = (it * it * it * it) % Max;
            vector<long long> a = v;
            vector<long long> b = v;

            auto t1 = std::chrono::high_resolution_clock::now();
            sort(a.begin(), a.end());
            auto t2 = std::chrono::high_resolution_clock::now();

            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "c++ sort:  " << ms_int.count() << "ms" << endl;

            b = v;
            t1 = std::chrono::high_resolution_clock::now();
            radix_sort10(b);
            t2 = std::chrono::high_resolution_clock::now();

            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "radix sort cu baza 10: " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';

            b = v;
            t1 = std::chrono::high_resolution_clock::now();
            radix_sort10(b, 256);
            t2 = std::chrono::high_resolution_clock::now();

            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "radix sort cu baza 2^8 fara operatii pe biti:  " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';


            b = v;
            t1 = std::chrono::high_resolution_clock::now();
            radix_sort(b, 8);
            t2 = std::chrono::high_resolution_clock::now();

            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "radix sort cu baza 2^8: " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';

            b = v;
            t1 = std::chrono::high_resolution_clock::now();
            radix_sort(b, 16);
            t2 = std::chrono::high_resolution_clock::now();

            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "radix sort cu baza 2^16: " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';


            b = v;
            t1 = std::chrono::high_resolution_clock::now();
            merge_sort(b);
            t2 = std::chrono::high_resolution_clock::now();

            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "merge sort:  " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';


            b = v;
            t1 = std::chrono::high_resolution_clock::now();
            shell_sort(b);
            t2 = std::chrono::high_resolution_clock::now();

            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "shell sort basic:  " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';


            b = v;
            t1 = std::chrono::high_resolution_clock::now();
            shell_sortC(b);
            t2 = std::chrono::high_resolution_clock::now();

            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "shell sort ciura:  " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';

            if (Max != p2[2]) {
                b = v;
                t1 = std::chrono::high_resolution_clock::now();
                bucket_sort(b);
                t2 = std::chrono::high_resolution_clock::now();

                ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
                cout << "bucket sort:  " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';
            }


            b = v;
            t1 = std::chrono::high_resolution_clock::now();
            heap_sort(b);
            t2 = std::chrono::high_resolution_clock::now();

            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            cout << "heap sort:  " << ms_int.count() << "ms" << "  " << (sortat_bine(b, a) ? "sortat bine" : "nesortat") << '\n';


            cout << '\n' << '\n';
        }
    }
}

int main()
{
    testeaza_cu_nr_random();
    return 0;
}