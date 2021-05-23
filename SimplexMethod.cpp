//
// Created by kesa on 23.05.2021.
//

#include "MyMatrix.h"

using namespace std;
/// <summary>
/// индекс наибольшего элемента вектора
/// </summary>
/// <param name="b">исходный вектор</param>
/// <returns>индекс наибольшего элемента в векторе</returns>
struct GameModel {
    int rowCount{}, colCount{};
    vector<db> b, zFuncCoeff;
    Matrix A;
    bool toInf = false;
};

struct SimplexSolution {
    vector<db> zFuncCoeff, x;
};

void printSystemOfEquation(ostream &out, GameModel gm) {
    for (int i = 1; i < gm.A.rowNumb; i++) {
        cout << endl;
        for (int j = 1; j < gm.A.colNumb; j++) {
            if (gm.A[i][j] == 0) cout << "\t";
            else cout << gm.A[i][j] << "*X" << j << " + ";
        }
        cout << "1*X" << gm.A.colNumb - 1 + i << " = " << gm.A[i][0] << ";";
    }
}

ostream &operator<<(ostream &out, const SimplexSolution &ss) {
    vector<db> c = ss.zFuncCoeff;
    out << "\nTarget function formula expressed in terms of free variables:\n";
    out << "Z(X) = " << c[c.size() - 1];
    for (int i = 0; i < c.size() - 1; ++i)
        out << " + " << c[i] << "*X" << i + 1;

    out << "\nX = {" << ss.x[0];
    for (int i = 1; i < ss.x.size(); ++i)
        out << ", " << ss.x[i];
    out << " }";
    return out;
}

class SimplexMethod {
private:
    vector<int> baseId;
    vector<int> freeId;
    Matrix exA;

    void createExtendedMatrix(GameModel gm) {
        exA = Matrix(gm.rowCount + 1, gm.colCount + 1);
        exA[0] = gm.zFuncCoeff;

        for (int i = 1; i < exA.rowNumb; i++)
            exA[i] = gm.A[i - 1];
        for (int i = 0; i < exA.rowNumb; i++)
            exA[i].push_back(gm.b[i]);

    }

    /// look for positive coefficients in resolved column - resolved string
    /// \param resCol index of resolved column
    /// \return index of resolved string
    int getResolvedStrInd(int resCol) {
        bool hasPositive = false;
        int resStr = 0, bInd = exA.colNumb - 1;
        db mnm = 0;
        for (int i = 1; i < exA.rowNumb; i++)
            if (exA[i][resCol] > 0)
                if (!hasPositive || mnm > exA[i][bInd] / exA[i][resCol]) {
                    hasPositive = true;
                    resStr = i;
                    mnm = exA[i][bInd] / exA[i][resCol];
                }
        if (!hasPositive)
            throw invalid_argument("\nThe function is unlimited. Optimal solution doesn't exist");
        return resStr;
    }

    void findOptimalResult() {
        //look for positive coefficients of the target function
        //resolved column is the largest of them
        int n = exA.rowNumb, m = exA.colNumb;
        MyVect c = MyVect(exA[0]);
        for (int resCol = c.getMaxId(0, m - 1); c[resCol] > eps; resCol = c.getMaxId(0, m - 1)) {
            int resStr = getResolvedStrInd(resCol);
            //swap the resolved row and column (Jordan's transformation)
            swap(freeId[resCol], baseId[resStr - 1]);
            Matrix A_ = exA;
            db resEl = A_[resStr][resCol] = 1 / exA[resStr][resCol];

            //transform a resolved string
            for (int j = 0; j < A_.colNumb; j++)
                if (j != resCol)
                    A_[resStr][j] = exA[resStr][j] * resEl;

            //transform a resolved column
            for (int i = 0; i < A_.rowNumb; i++)
                if (i != resStr)
                    A_[i][resCol] = exA[i][resCol] * resEl * -1.;

            //transform all other elements according to the rectangle rule
            for (int i = 0; i < A_.rowNumb; i++)
                if (i != resStr)
                    for (int j = 0; j < A_.colNumb; j++)
                        if (j != resCol)
                            A_[i][j] = exA[i][j] + exA[resStr][j] * A_[i][resCol];
            exA = A_;
            c = MyVect(exA[0]);
        }
    }

public:
    SimplexMethod(const GameModel &gm) {
        createExtendedMatrix(gm);
        baseId = vector<int>(gm.rowCount);
        freeId = vector<int>(gm.colCount);
        //set correspondence between free variable and table column
        for (int j = 0; j < gm.colCount; j++)
            freeId[j] = j + 1;
        //set correspondence between base variable and table row
        for (int i = 0; i < gm.rowCount; i++)
            baseId[i] = i + 1 + gm.colCount;
    }

    SimplexSolution calculateOptSolution() {
        int n = exA.rowNumb, m = exA.colNumb;
        findOptimalResult();
        SimplexSolution sp;
        sp.zFuncCoeff = vector<db>(n + m - 1, 0);
        sp.x = vector<db>(n + m - 2, 0);
        for (int i = 0; i < baseId.size(); ++i)
            sp.x[baseId[i] - 1] = exA[i + 1][m - 1];
        for (int i = 0; i < freeId.size(); ++i)
            sp.zFuncCoeff[freeId[i] - 1] = exA[0][i];
        sp.zFuncCoeff[n + m - 2] = exA[0][m - 1];
        return sp;
    }
};

istream &operator>>(istream &in, GameModel &gm) {
    int n, m;
    in >> n >> m;
    gm.zFuncCoeff = MyVect::input(in, m);
    gm.b = MyVect::input(in, n + 1);
    gm.A = Matrix(n, m);
    gm.A.input(in);
    gm.colCount = m;
    gm.rowCount = n;
    return in;
}

int main() {
    setlocale(LC_ALL, "rus");
    freopen("..\\input.dat", "r", stdin);
    GameModel gm;
    cin >> gm;
    SimplexMethod sm = SimplexMethod(gm);
    SimplexSolution ss = sm.calculateOptSolution();
    cout << ss;
}

