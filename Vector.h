#include <iostream>
#include <cassert>
#include <cmath>

using std::cout;
using std::endl;
using std::cin;

class Vector {

private:
    size_t n;
    double *data;

public:
    Vector(size_t rows, double value = 0.0); //default constructor
    Vector(const Vector &other); //copy constructor
    ~Vector(); // destructor
    Vector &operator=(Vector other); //copy assignment and move assignment: copy and swap idiom

/*---------------------------------------------------------------------------------------------------------------*/

    //Trying move constructors (It worked)

    void swap(Vector &other); //for efficient swapping
    Vector(Vector &&other); //move constructor

    /*Important reference video regarding copy and move constructors basics: https://www.youtube.com/watch?v=cO1lb2MiDr8
    The video clearly explains what is rvalue reference, move constructor, std::swap and std::move functions, how copy
    assignment and move assignment operators can be combined without using call by value method without any overhead
    losses of call by value method.*/

/*---------------------------------------------------------------------------------------------------------------*/

    //Data access operator overloading
    double &operator()(size_t i); // ()operator overloading for returning value
    const double &
    operator()(size_t i) const; // () operator overloading for assigning value to an element (i,j)

/*---------------------------------------------------------------------------------------------------------------*/

    //Arithmetic operators overloading
    Vector operator+(const Vector &other) const; //+ operator overloading
    Vector operator-(const Vector &other) const; // - operator overloading
    void operator+=(const Vector &other); // += operator overloading
    void operator-=(const Vector &other); // -= operator overloading


    bool operator==(const Vector &other) const; // == operator overloading
    bool operator!=(const Vector &other) const; // != operator overloading

/*---------------------------------------------------------------------------------------------------------------*/

    //iostream operators overloading
    friend std::ostream &operator<<(std::ostream &os, const Vector &a); // ostream << operator overloading
    friend std::istream &operator>>(std::istream &, Vector &input); // >> istream >> operator overloading

/*---------------------------------------------------------------------------------------------------------------*/

    //Implemented as per question
    size_t rows() const; //returns the number of rows
    
};

Vector::Vector( size_t columns, double value) :  n(columns), data(new double[ n]) {
    for (size_t i = 0; i <  n; i++)
        data[i] = value;
}

Vector::Vector(const Vector &other) :  n(other.n), data(new double[ n]) {
    for (size_t i = 0; i <  n; i++)
        data[i] = other.data[i];
}

void Vector::swap(Vector &other) {
    std::swap(n, other.n);
    std::swap(data, other.data);
}

Vector &Vector::operator=(Vector other) {
    swap(other);
    return *this;
}

Vector::Vector(Vector &&other) : n(0), data(nullptr) {
    swap(other);
}

double &Vector::operator()(size_t i) {
    assert( i < n);
    return data[i ];
}

const double &Vector::operator()(size_t i) const {
    assert(i <= n) ;
    return data[i];
}

Vector Vector::operator+(const Vector &other) const {
    assert( n == other.n);
    Vector temp = *this;
    for (size_t i = 0; i <  n; i++)
        temp.data[i] = data[i] + other.data[i];
    return temp;
}

void Vector::operator+=(const Vector &other) {
    assert( n == other.n);
    *this = *this + other;
}

Vector Vector::operator-(const Vector &other) const {
    assert( n == other.n);
    Vector temp = *this;
    for (size_t i = 0; i < n; i++)
        temp.data[i] = data[i] - other.data[i];
    return temp;
}

void Vector::operator-=(const Vector &other) {
    assert( n == other.n);
    *this = *this - other;
}


Vector::~Vector() {
    delete[] data;
}

std::ostream &operator<<(std::ostream &os, const Vector &output) {
    os << std::endl; //formatting
    for (size_t i = 0; i < output.n; i++) 
            os << output(i) << " ";
    os << std::endl; //formatting
    return os;
}

std::istream &operator>>(std::istream &is, Vector &input) {
        for (size_t j = 0; j < input.n; j++) {
            if (is >> input( j))
                continue;
            else {
                cout << "The elements of the Vector can only be double" << endl;
                abort();
            }
        }
    return is;
}


bool Vector::operator==(const Vector &other) const {
    if (n == other.n) {
        for (size_t j = 0; j < n; j++) {
            if ((*this)( j) == other(j))
                continue;
            else
                return 0;
        }
        return 1;
    } else
        return 0;
}

bool Vector::operator!=(const Vector &other) const {
    return !(*this == other);
}

size_t Vector::rows() const {
    return n;
}

