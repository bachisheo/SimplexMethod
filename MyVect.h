//
// Created by kesa on 23.05.2021.
//

#ifndef SIMPLEXMETHOD_MYVECT_H
#define SIMPLEXMETHOD_MYVECT_H

#include <utility>
#include <vector>
#include <iostream>
#include <cassert>

#define db double
#define ll long long int


class MyVect {
    std::vector<db> dat;

public:
    MyVect(std::vector<double> dat) : dat(std::move(dat)) {}

    MyVect slice(int first, int last) {
        assert(first >= 0 && first < last && last <= dat.size());
        return MyVect(std::vector<db>(dat.begin() + first, dat.begin() + last));
    }

    static std::vector<db> input(std::istream &in, int n) {
        std::vector<db> vec = std::vector<db>(n);
        for (int i = 0; i < n; ++i)
            in >> vec[i];
        return vec;
    }

    MyVect() {
        dat = std::vector<db>();
    }

    db &operator[](const int index) {
        return dat[index];
    }

    int getMaxId(int l, int r) {
        int mxmId = l;
        for (int i = l + 1; i < r; i++)
            if (dat[i] > dat[mxmId])
                mxmId = i;
        return mxmId;
    }

    int getMaxId() {
        return getMaxId(0, dat.size());
    }

    db getMax() {
        db mxm = dat[0];
        for (int i = 1; i < dat.size(); i++)
            if (dat[i] > mxm) mxm = dat[i];
        return mxm;
    }
};


#endif //SIMPLEXMETHOD_MYVECT_H
