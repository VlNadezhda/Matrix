
#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <cassert>
 #include <cstring>

template<typename T>
class Matrix;

template<typename T>
std::ostream & operator<< (std::ostream &, const Matrix<T> &);

template<typename T>
Matrix<T> operator+(const Matrix<T> M1, const Matrix<T> M2);

template<typename T>
Matrix<T> operator-(const Matrix<T> M1, const Matrix<T> M2);

template<typename T>
Matrix<T> operator*(const Matrix<T> M1, const Matrix<T> M2);

template<typename T>
class Matrix
{
private:
    int row;
    int col;
    std::vector<std::vector<T>> matrix;

public:

    Matrix();
    Matrix(const std::initializer_list<T> &list);
    Matrix(const std::initializer_list<std::initializer_list<T>> &list);

    Matrix(const std::initializer_list<char> &list)                        = delete;
    Matrix(const std::initializer_list<std::initializer_list<char>> &list) = delete;

    Matrix(size_t r = 0, size_t c = 0):row(r),col(c){
        matrix.reserve(row);
    };

    void fillMatrix(){
        for(uint i = 0; i < row; ++i){
            matrix.emplace_back(std::vector(col,0));
        }
    }
    void fillMatrixRand(int k);
    Matrix operator=(const Matrix&);

    std::vector<T> getRow(int i){
        return  matrix[i];
    };
    std::vector<T> getColumn(int j){
        std::vector<T> col;
        for(auto a: matrix){
            col.push_back(a.at(j));
        }
        return col;
    };

    void setRow(const std::vector<T>&,int);
    void setColumn(const std::vector<T>&,int);
    void setElement(uint,uint, const int);

    int getRowsNum() const noexcept{
        return this->row;
    }
    int getColumnsNum() const noexcept{
        return this->col;
    }

    bool operator==(const Matrix&);
    T operator()(uint i, uint j) const{
        assert(i <= this->row && j <= this->col);
        return matrix[i][j];
    }
    T &operator()(uint i, uint j){
        assert(i <= this->row && j <= this->col);
        return matrix[i][j];
    }
    Matrix& operator+=(const Matrix&);
    Matrix& operator+=(double);

    Matrix& operator-=(const Matrix&);
    Matrix& operator-=(double);

    Matrix& operator*=(const Matrix&);
    Matrix& operator*=(double);

    Matrix& operator/=(double);
    Matrix&  operator^(int);

    Matrix inverse();
    double determinant();

    friend std::ostream& operator<< <>(std::ostream&, const Matrix&);
    ~Matrix(){}
};

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<T> &list):
        Matrix(list.size(), 0){
    matrix.push_back(list);
}

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>> &list):
        Matrix(list.size(), list.size() ? list.begin()->size() : 0){

    for (const auto &l : list){
       matrix.push_back(l);
    }
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &mat) {
    for (uint i = 0; i < mat.getRowsNum(); ++i) {
        for (uint j = 0; j < mat.getColumnsNum(); ++j){
            out<<std::setw(4)<<mat.matrix[i][j]<<" ";
        }
        std::cout<<"\n";
    }
    return out;
}

template<typename T>
Matrix<T> Matrix<T>::operator=(const Matrix<T> &mat) {
//    assert(this->row == mat.getRowsNum());
//    assert(this->col == mat.getColumnsNum());
    if(matrix.size() == mat.matrix.size() &&
            std::equal(matrix.begin(), matrix.end(), mat.matrix.begin())){
        return *this;
    }
    this->row = mat.getRowsNum();
    this->col = mat.getColumnsNum();

    for (uint i = 0; i < this->row; ++i) {
            for (uint j = 0; j < this->col; ++j){
                matrix[i].clear();
                this->matrix[i][j] = mat.matrix[i][j];
            }
        }

    return *this;
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T> &mat){
    assert(this->row == mat.getRowsNum());
    assert(this->col == mat.getColumnsNum());
     for(uint i = 0; i < this->row; ++i){
            for(uint j = 0; j < this->col; ++j){
                if (mat.matrix[i][j] != this->matrix[i][j])
                    return false;
            }
        }
        return true;
}

template<typename T>
Matrix<T> operator+(const Matrix<T> m1,const Matrix<T>  m2){
   int row = m2.getRowsNum();
   int col = m2.getColumnsNum();
   assert(m1.getColumnsNum() == col && m1.getRowsNum() == row);

   Matrix<T> m(row,col);
   m.fillMatrix();
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            m.setElement(i,j,m1(i,j) + m2(i,j));
        }
    }
    return m;
}

template<typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &mat){
    Matrix<T> result = matrix + mat;
    (*this) = result;
}


template<typename T>
Matrix<T> operator-(const Matrix<T> m1,const Matrix<T>  m2){
   int row = m2.getRowsNum();
   int col = m2.getColumnsNum();
   assert(m1.getColumnsNum() == col && m1.getRowsNum() == row);

   Matrix<T> m = m1;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            m.setElement(i,j, m1(i,j) - m2(i,j));
        }
    }
    return m;
}

template<typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &mat) {
    Matrix<T> result = matrix - mat.matrix;
    (*this) = result;
    return *this;
}


template<typename T>
Matrix<T> operator*(const Matrix<T> m1,const Matrix<T>  m2){
   int row = m2.getRowsNum();
   int col = m2.getColumnsNum();
   assert(m1.getRowsNum() == col && m1.getColumnsNum() == row);

   Matrix<T> m(col,row);
   m.fillMatrix();
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            for(int k = 0; k < col;  ++k){
                 m.setElement(i,j,m1(i,k)*m2(k,j));
            }
        }
    }
    return m;
}

template<typename T>
Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &mat) {
    Matrix<T> result = matrix * mat;
    (*this) = result;
    return this;
}


template<typename T>
void Matrix<T>::setRow(const std::vector<T>& row, int i){
    assert(row.size() == col);
        matrix[i] = row;
}

template<typename T>
void Matrix<T>::setColumn(const std::vector<T>& col, int j){
    assert(col.size() == row);
    for(uint i = 0; i < this->getRowsNum(); ++i){
        matrix.at(i).at(j) = col.at(i);
    }
}

template<typename T>
void Matrix<T>::fillMatrixRand(int k){
    std::vector<T> tmp;
    for(uint i = 0; i < this->getRowsNum(); ++i){
        for(uint j = 0; j < this->getColumnsNum(); ++j){
            tmp.push_back(rand()%k);
        }
        matrix.push_back(std::move(tmp));
    }

}

template<typename T>
void Matrix<T>::setElement(uint i,uint j, const int num){
    assert(i <= row && j <= col);
    matrix[i][j] = num;
//    std::cout<<*this;
};

template<typename T>
Matrix<T> &Matrix<T>::operator^(const int n){
    for(uint i = 1; i <= n; ++i){
        matrix*= matrix;
    }
};

template<typename T>
Matrix<T> Matrix<T>::inverse(){};

template<typename T>
double Matrix<T>::determinant(){};

#endif //MATRIX_MATRIX_H
