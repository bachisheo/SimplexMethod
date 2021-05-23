
#include <fstream>
#include "MyMatrix.h"
#pragma warning (disable:4996)
using namespace std;
    void Matrix::output(ostream& fout) {
        fout << endl;
        for (int i = 0; i < rowNumb; i++)
        {
            for (int j = 0; j < colNumb; j++)
                fout << dat[i][j] << "\t";
            fout << endl;
        }
    }
    //ввод матрицы
    void Matrix::input(istream& fin) {
        for (int i = 0; i < rowNumb; i++)
            for (int j = 0; j < colNumb; j++)
                fin >> dat[i][j];
    }
    ///
    /// \param cstr - количество строк
    /// \param ccol количество столбцов
    Matrix::Matrix(int cstr, int ccol) : rowNumb(cstr), colNumb(ccol) {
        size = min(rowNumb, colNumb);
        dat = vector<vector<double>>(rowNumb);
        for (int i = 0; i < rowNumb; i++)
            dat[i] = vector<double>(colNumb);
    };

    /// конструктор квадратной матрицы
    /// \param size
    Matrix::Matrix(int size) : Matrix(size, size) {}

    /// конструктор матрицы из двумерного вектора
    /// \param M
    Matrix::Matrix(vector<vector<double>> M) : Matrix(M.size(), M[0].size())
    {
        for (int i = 0; i < M.size(); i++)
            for (int j = 0; j < M[i].size(); j++)
                dat[i][j] = M[i][j];
    }

    /// Конструктор одномерной матрицы из вектора
    /// \param M
    Matrix::Matrix(vector<double> M) : Matrix(M.size(), 1)
    {
        for (int i = 0; i < rowNumb; i++)
            dat[i][0] = M[i];
    }

    /// сумма двух матриц
    /// \param otherMatrix
    /// \return
    Matrix Matrix::operator+(const Matrix &otherMatrix) const {
                Matrix c = Matrix(otherMatrix);
        for (int i = 0; i < this->rowNumb; i++)
            for (int j = 0; j < colNumb; j++)
                c.dat[i][j] = otherMatrix.dat[i][j] + dat[i][j];
        return c;
    }
    Matrix Matrix::operator +(db d) const {
        Matrix c = Matrix(rowNumb, colNumb);
        for (int i = 0; i < rowNumb; i++)
            for (int j = 0; j < colNumb; j++)
                c.dat[i][j] = dat[i][j] + d;
        return c;
    }
    Matrix Matrix::operator * (const Matrix& b) const {
        Matrix c = Matrix(0);
        if (colNumb == b.rowNumb) {
            c = Matrix(rowNumb, b.colNumb);
            for (int i = 0; i < rowNumb; i++)
                for (int j = 0; j < b.colNumb; j++)
                {
                    c.dat[i][j] = 0;
                    for (int k = 0; k < colNumb; k++)
                        c.dat[i][j] += dat[i][k] * b.dat[k][j];
                }
        }
        return c;
    }




Matrix Matrix::operator *(db d) const {
        Matrix c = Matrix(rowNumb, colNumb);
        for (int i = 0; i < rowNumb; i++)
            for (int j = 0; j < colNumb; j++)
                c.dat[i][j] = dat[i][j] *d;
        return c;
    }

    ///умножение матрицы на вектор
    ///
    const vector<db> Matrix::operator *(vector<db> v) const {
        Matrix b = Matrix(v);
        return ((*this) * b).trans().dat[0];
    }
    vector<db>& Matrix::operator [](const int rowIndex)  {
        return dat[rowIndex];
    }


    Matrix Matrix::trans() const{
        Matrix B = Matrix(colNumb, rowNumb);
        for (int i = 0; i < rowNumb; i++)
            for (int j = 0; j < colNumb; j++)
                B.dat[j][i] = dat[i][j];
        return B;
}

void Matrix::push_back(vector<double> vec) {
        dat.push_back(vec);
        rowNumb++;
}


void TriangleMatrix::swapString(int a, int b) {
        //зафиксировать изменение в порядке строк
        swap(orderOfString[a], orderOfString[b]);
        //поменять местами содержимое строк
        vector<double> tmp = dat[a];
        dat[a] = dat[b];
        dat[b] = tmp;
        numbOfTurns++;
    }
TriangleMatrix::TriangleMatrix(vector<vector<double>> M) : Matrix(M) {
        //приветси матрицу к треугольному виду
        orderOfString = vector<int>(rowNumb);
        for (int i = 0; i < rowNumb; i++)
        {
            orderOfString[i] = i;
        }

        for (int k = 0; k < rowNumb; k++)
        {
            //найти главный (ведущий элемент) в столбце
            int indOfMxm = k;
            for (int p = k; p < rowNumb; p++)
            {
                if (abs(dat[p][k]) > abs(dat[indOfMxm][k]))
                    indOfMxm = p;
            }
            //если все эелемнты 0 - пропустить итерацию
            if (abs(dat[indOfMxm][k]) < eps) {
                det = 0;
                //exit(12);
            }
            //поместить ведущий элемент на главную диагональ
            if (indOfMxm != k) swapString(indOfMxm, k);
            //в каждой строке, начиная с k+1
            for (int p = k + 1; p < rowNumb; p++)
            {
                //найти коэффициент  С[p][k]
                double C = dat[p][k] / dat[k][k];
                //вычесть из текущей строки, k-тую строку
                //умноженную на С
                for (int l = k; l < colNumb; l++)
                {
                    dat[p][l] = dat[p][l] - dat[k][l] * C;
                }
                //сохранить коэффициент
                dat[p][k] = C;
            }
        };
    }
    //вычислить определитель матрицы
    double TriangleMatrix::getDet() {
        det = 1;
        //вычислить определитель треугольной матрицы
        //как произведение элементов, стоящих на главной диагонали
        for (int k = 0; k < size; k++)
        {
            if (abs(dat[k][k]) < eps) return 0;
            det *= dat[k][k];
        }
        //изменить знак определителя, если строки
        //были переставлены четное число раз
        if (numbOfTurns % 2) det *= (-1);
        return det;
    }
    //преобразовать вектор b к трекгольному виду матрицы
    vector<double> TriangleMatrix::transformVector(vector <double> b) {
        vector<double> bk;
        if (b.size() != size) exit(40);
        bk = vector<double>(size);
        //переставить строки в том же порядке, что и у треугольной матрицы А
        for (int k = 0; k < size; k++)
        {
            bk[k] = b[orderOfString[k]];
        }
        //провети над вектором преобразовния с коэффициентами треугольной матрицы
        for (int k = 0; k < size - 1; k++)
        {
            for (int l = k + 1; l < size; l++)
            {
                bk[l] = bk[l] - bk[k] * dat[l][k];
            }
        }
        return bk;
    }

// норма разности двух векторов
db norm(vector<db> a, vector<db> b) {
    db mxm = abs(a[0] - b[0]);
    for (int i = 1; i < a.size(); i++)
    {
        if (mxm < abs(b[i] - a[i])) mxm = abs(b[i] - a[i]);
    }
    return mxm;
}

//обратный ход
vector <double> revStep(vector<double> b, vector<vector<double>> a) {
    int m = b.size();
    vector <double> x = vector<double>(m);
    //найти каждый х, начиная с i-того
    for (int i = m - 1; i >= 0; i--)
    {
        //посчитать сумму произведений известных х
        //на соответствующие элементы матрицы А
        double sum = 0;
        for (int l = i + 1; l < m; l++)
        {
            sum += x[l] * a[i][l];
        }
        //вычислить текущий х
        x[i] = b[i] - sum;
        x[i] = x[i] / a[i][i];
    }
    return x;
}
//поиск решения СЛАУ методом Гаусса
vector<double> GaussMethod(Matrix B, vector <double> b) {
    TriangleMatrix A = TriangleMatrix(B.dat);
    int m = b.size();
    vector<double> x;
    //проверить определитель матрицы
    if (abs(A.getDet()) < eps) {
        return x;
    }
    //преобразовать правую часть уравнения
    vector<double> bk = A.transformVector(b);
    //обратный ход метода Гаусса
    return revStep(bk, A.dat);
}
vector <db> prodVectOnDigit(vector<db> b, db a) {
    for (int i = 0; i < b.size(); i++)
    {
        b[i] = b[i] * a;
    }
    return b;
}
//вычисление невязки r = Ae1 - L1e1
vector<db> residual(Matrix A, vector <db> e1, db l1) {
    int m = e1.size();
    vector<db> res = vector<db>(m);
    vector<db> Ae1 = A*e1;
    vector<db> L1e1 = prodVectOnDigit(e1, l1);
    //для каждой строки
    for (int i = 0; i < m; i++)
    {
        res[i] = Ae1[i] - L1e1[i];
    }
    return res;
}

//скалярное произведение векторов
db scalarProd(vector<db> a, vector<db> b) {
    db res = 0;
    for (int i = 0; i < a.size(); i++)
    {
        res += a[i] * b[i];
    }
    return res;
}
db getNorm(vector<db> a) {
    db mxm = abs(a[0]);
    for (int i = 1; i < a.size(); i++)
    {
        mxm = max(abs(a[i]), mxm);
    }
    return mxm;
}

//нормирование вектора
vector<db> doNorm(vector<db> x) {
    //найти норму вектора
    db vNorm = getNorm(x);
    //разделить вектор на его норму
    for (int i = 0; i < x.size(); i++)
    {
        x[i] = x[i] / vNorm;
    }
    return x;
}

// максимальное по модулю собственное числа матрицы А
db searchMaxAbsL(Matrix A, vector<db>& x0, int& k) {
    int delt = 15;
    vector<db> x = x0;
    x = A*x0;
    db l1 = scalarProd(x, x0) / scalarProd(x0, x0);
    db l2 = l1 * 1000;
    //пока разница между последовательно найденными числами больше
    //заданной точности или число итераций не достигло макисмальное
    //- вычислять новое приближение вектора
    //и собственного числа
    for (k = 0; abs(l1 - l2) > eps && k <= maxIterationNumber; k++) {
        x0 = x, l2 = l1;
        //x(rowNumb+1) = A*x(rowNumb)
        x = A* x0;
        //нормируем вектора через заданное число итераций
        l1 = scalarProd(x, x0) / scalarProd(x0, x0);
        if (k % delt == 0) {
            x = doNorm(x);
        }
    }
    //нормировать последний вычесленный вектор, получить собственный
    //вектор для вычисленного собственного числа
    x0 = doNorm(x);
    return l1;
}


db getL(Matrix A, vector <db>& x0, int& k, long long c) {
    int sz = A.size;
    Matrix Ec = Matrix(sz);
    for (int i = 0; i < sz; i++)
    {
        for (int j = 0; j < sz; j++)
        {
            Ec.dat[i][j] = 0;
        }
        Ec.dat[i][i] = c * 1.0;
    }
    Matrix B = A + Ec;
    db l = searchMaxAbsL(B, x0, k);
    return l - c;
}

//наибольшее собственное число матрицы А
db getMaxL(Matrix A, vector<db>& x0, int& k) {
    return getL(A, x0, k, C);
}

//наименьшее собственное число матрицы А
db getMinL(Matrix A, vector<db>& x0, int& k) {
    return getL(A, x0, k, C * (-1));
}

//собственное число матрицы А, ближайшее к заданному l0
db nearbyLambda(Matrix A, db l0, vector<db>& x0, int& k) {
    int sz = A.size;
    Matrix El = Matrix(sz);
    for (int i = 0; i < sz; i++)
    {
        for (int j = 0; j < sz; j++)
        {
            El.dat[i][j] = 0;
        }
        El.dat[i][i] = l0 * (-1);
    }
    El = A + El;
    int delt = 15;
    vector<db> x = x0;
    //x = itherMethod(El, x0, x0, eps,0);
    x = GaussMethod(El, x0);
    if (x.size() == 0)
        return nearbyLambda(A, l0 + eps, x0, k);
    db l1 = l0 + scalarProd(x0, x0) / scalarProd(x, x0);
    db l2 = l1 * 100;
    //пока разница между последовательно найденными числами больше
    //заданной точности или число итераций не достигло макисмальное
    //- вычислять новое приближение вектора
    //и собственного числа
    for (k = 0; abs(l1 - l2) > eps && k <= maxIterationNumber; k++)
    {
        l2 = l1;
        x0 = doNorm(x);
        //x = itherMethod(El, x0, x0, eps, 0);
        x = GaussMethod(El, x0);

        l1 = l0 + scalarProd(x0, x0) / scalarProd(x, x0);

    }
    x0 = x;
    return l1;
}


//печать вектора
void printVector(vector<double> vect, ostream& out) {
    out << endl;
    if (vect.empty()) out << "\nThis vector is empty";
    for (int i = 0; i < vect.size(); i++)
    {
        out << vect[i] << endl;
    }
}
//ввод вектора
vector <double> inputVector(int n, istream& in) {
    vector<double> b = vector<double>(n);
    for (int i = 0; i < n; i++) in >> b[i];
    return b;
}
//норма вектора
Matrix getHilbert(int n) {
    vector<vector<double>> h = vector<vector<double>>(n);

    for (int i = 0; i < n; i++)
    {
        h[i] = vector<double>(n);
        for (int j = 0; j < n; j++)
        {
            h[i][j] = 1.0 / (i + j + 1);
        }
    }
    return Matrix(h);
}