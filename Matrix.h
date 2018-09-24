#include <iostream>
#include <cassert>
#include <cmath>

using std::cout;
using std::endl;
using std::cin;

class Matrix {

private:
    size_t m, n;
    double *data;

public:
    Matrix(size_t rows, size_t columns, double value = 0.0); //default constructor
    Matrix(const Matrix &other); //copy constructor
    ~Matrix(); // destructor
    Matrix &operator=(Matrix other); //copy assignment and move assignment: copy and swap idiom

/*---------------------------------------------------------------------------------------------------------------*/

    //Trying move constructors (It worked)

    void swap(Matrix &other); //for efficient swapping
    Matrix(Matrix &&other); //move constructor

    /*Important reference video regarding copy and move constructors basics: https://www.youtube.com/watch?v=cO1lb2MiDr8
    The video clearly explains what is rvalue reference, move constructor, std::swap and std::move functions, how copy
    assignment and move assignment operators can be combined without using call by value method without any overhead
    losses of call by value method.*/

/*---------------------------------------------------------------------------------------------------------------*/

    //Data access operator overloading
    double &operator()(size_t i, size_t j); // ()operator overloading for returning value
    const double &
    operator()(size_t i, size_t j) const; // () operator overloading for assigning value to an element (i,j)

/*---------------------------------------------------------------------------------------------------------------*/

    //Arithmetic operators overloading
    Matrix operator+(const Matrix &other) const; //+ operator overloading
    Matrix operator-(const Matrix &other) const; // - operator overloading
    void operator+=(const Matrix &other); // += operator overloading
    void operator-=(const Matrix &other); // -= operator overloading
    Matrix operator*(const Matrix &other) const; // * operator overloading
    void operator*=(const Matrix &other); // *= operator overloading

    bool operator==(const Matrix &other) const; // == operator overloading
    bool operator!=(const Matrix &other) const; // != operator overloading

/*---------------------------------------------------------------------------------------------------------------*/

    //iostream operators overloading
    friend std::ostream &operator<<(std::ostream &os, const Matrix &a); // ostream << operator overloading
    friend std::istream &operator>>(std::istream &, Matrix &input); // >> istream >> operator overloading

/*---------------------------------------------------------------------------------------------------------------*/

    //Implemented as per question
    size_t rows() const; //returns the number of rows
    size_t cols() const; //return the number of columns
};

Matrix::Matrix(size_t rows, size_t columns, double value) : m(rows), n(columns), data(new double[m * n]) {
    for (size_t i = 0; i < m * n; i++)
        data[i] = value;
}

Matrix::Matrix(const Matrix &other) : m(other.m), n(other.n), data(new double[m * n]) {
    for (size_t i = 0; i < m * n; i++)
        data[i] = other.data[i];
}

void Matrix::swap(Matrix &other) {
    std::swap(m, other.m);
    std::swap(n, other.n);
    std::swap(data, other.data);
}

Matrix &Matrix::operator=(Matrix other) {
    swap(other);
    return *this;
}

Matrix::Matrix(Matrix &&other) : m(0), n(0), data(nullptr) {
    swap(other);
}

double &Matrix::operator()(size_t i, size_t j) {
    assert((i < m) && (j < n));
    return data[i * n + j];
}

const double &Matrix::operator()(size_t i, size_t j) const {
    assert((i <= m) && (j <= n));
    return data[i * n + j];
}

Matrix Matrix::operator+(const Matrix &other) const {
    assert(m == other.m && n == other.n);
    Matrix temp = *this;
    for (size_t i = 0; i < m * n; i++)
        temp.data[i] = data[i] + other.data[i];
    return temp;
}

void Matrix::operator+=(const Matrix &other) {
    assert(m == other.m && n == other.n);
    *this = *this + other;
}

Matrix Matrix::operator-(const Matrix &other) const {
    assert(m == other.m && n == other.n);
    Matrix temp = *this;
    for (size_t i = 0; i < m * n; i++)
        temp.data[i] = data[i] - other.data[i];
    return temp;
}

void Matrix::operator-=(const Matrix &other) {
    assert(m == other.m && n == other.n);
    *this = *this - other;
}

Matrix Matrix::operator*(const Matrix &other) const {
    assert(n == other.m);
    Matrix result(m, other.n, 0.0);
    for (size_t i = 0; i < m; i++) {
        for (size_t k = 0; k < other.n; k++) {
            double sum = 0.0;
            for (size_t j = 0; j < other.m; j++) {
                sum += ((*this)(i, j) * other(j, k));
            }
            result(i, k) = sum;
        }
    }
    return result;
}

void Matrix::operator*=(const Matrix &other) {
    assert(n == other.m);
    Matrix prod = *this * other;
    *this = prod;
}

Matrix::~Matrix() {
    delete[] data;
}

std::ostream &operator<<(std::ostream &os, const Matrix &output) {
    os << std::endl; //formatting
    for (size_t i = 0; i < output.m; i++) {
        for (size_t j = 0; j < output.n; j++)
            os << output(i, j) << " ";
        os << std::endl;
    }
    os << std::endl; //formatting
    return os;
}

std::istream &operator>>(std::istream &is, Matrix &input) {
    for (size_t i = 0; i < input.m; i++)
        for (size_t j = 0; j < input.n; j++) {
            if (is >> input(i, j))
                continue;
            else {
                cout << "The elements of the matrix can only be double" << endl;
                abort();
            }
        }
    return is;
}


bool Matrix::operator==(const Matrix &other) const {
    if (m == other.m && n == other.n) {
        for (size_t i = 0; i < m; i++) {
            for (size_t j = 0; j < n; j++) {
                if ((*this)(i, j) == other(i, j))
                    continue;
                else
                    return 0;
            }
        }
        return 1;
    } else
        return 0;
}

bool Matrix::operator!=(const Matrix &other) const {
    return !(*this == other);
}

size_t Matrix::rows() const {
    return m;
}

size_t Matrix::cols() const {
    return n;
}
