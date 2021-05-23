#ifndef MYMATR
#define MYMATR

#include "MyVect.h"
const ll C = 10000;
using namespace std;
const double eps = 1e-7;
const int maxIterationNumber = 50000;
class Matrix {
public:
    //rowNumb - количество строк, colNumb - столбцов
    int size, colNumb, rowNumb;
    vector<vector<double>> dat;
    //печать матрицы
    void output(ostream& fout) ;
    //ввод матрицы
    void input(istream& fin) ;
    Matrix(){}
    //TODO add def for matrix() constructor

    Matrix(int cstr);
   ///
   /// \param cstr - количество строк
   /// \param ccol количество столбцов
    Matrix(int cstr, int ccol) ;
    /// конструктор матрицы из двумерного вектора
    /// \param M
    Matrix(vector<vector<double>> M) ;
   /// Конструктор одномерной матрицы из вектора
   /// \param M
    Matrix(vector<double> M) ;
    void push_back(vector<db> vec);
    /// сумма двух матриц
    /// \param otherMatrix
    /// \return
    //const после параметров - предупреждает о том, что метод не изменяет
    //поля класса
    Matrix operator +(const Matrix& otherMatrix) const;
    Matrix operator +(db d) const ;
    Matrix operator * (const Matrix& b) const ;
    Matrix operator *(db d) const ;
    ///
    ///умножение матрицы на вектор
    ///
    const vector<db> operator *(vector<db> v) const ;
    vector<db>& operator [](const int rowIndex)  ;

    Matrix trans() const;
};


class TriangleMatrix : public Matrix {
    //количество перестановок строк
    int numbOfTurns = 0;
    double det = 0;
    //измененный порядок строк
    vector<int> orderOfString;
    //поменять строки треугольной матрицы местами
    void swapString(int a, int b) ;
public:
    TriangleMatrix(vector<vector<double>> M) ;
    //вычислить определитель матрицы
    double getDet() ;
    //преобразовать вектор b к трекгольному виду матрицы
    vector<double> transformVector(vector <double> b);
};
// норма разности двух векторов
db norm(vector<db> a, vector<db> b) ;
//обратный ход
vector <double> revStep(vector<double> b, vector<vector<double>> a) ;
//поиск решения СЛАУ методом Гаусса
vector<double> GaussMethod(Matrix B, vector <double> b) ;
vector <db> prodVectOnDigit(vector<db> b, db a);
vector<db> residual(Matrix A, vector <db> e1, db l1);
//скалярное произведение векторов
db scalarProd(vector<db> a, vector<db> b) ;
db getNorm(vector<db> a) ;

//нормирование вектора
vector<db> doNorm(vector<db> x) ;

// максимальное по модулю собственное числа матрицы А
db searchMaxAbsL(Matrix A, vector<db>& x0, int& k);


db getL(Matrix A, vector <db>& x0, int& k, long long c) ;

//наибольшее собственное число матрицы А
db getMaxL(Matrix A, vector<db>& x0, int& k) ;


//наименьшее собственное число матрицы А
db getMinL(Matrix A, vector<db>& x0, int& k) ;
//собственное число матрицы А, ближайшее к заданному l0
db nearbyLambda(Matrix A, db l0, vector<db>& x0, int& k) ;
//печать вектора
void printVector(vector<double> vect, ostream& out) ;
//ввод вектора
vector <double> inputVector(int n, istream& in);
//норма вектора
Matrix getHilbert(int n) ;

#endif