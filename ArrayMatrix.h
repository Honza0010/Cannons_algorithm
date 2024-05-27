#ifndef __ARRAYMATRIX_
#define __ARRAYMATRIX_



#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>


template< typename Real >
class ArrayMatrix
{
public:
    ArrayMatrix(int size)
    {
        this->matrix = new Real[size * size];
        this->size = size;
    }

    void setValue(const Real& v)
    {
        for (int i = 0; i < size * size; i++)
            this->matrix[i] = v;
    }

    void setElement(int i, int j, const Real& v)
    {
        this->matrix[i * size + j] = v;
    }

    ArrayMatrix& operator = (const ArrayMatrix& m)
    {
        memcpy(this->matrix, m.matrix, size * size * sizeof(Real));
        return (*this);
    }

    bool operator == (const ArrayMatrix& m) const
    {
        if (memcmp(this->matrix, m.matrix, size * size * sizeof(Real)) == 0)
            return true;
        return false;
    }

    bool operator != (const ArrayMatrix& m) const
    {
        return !this->operator == (m);
    }

    void multiplyMatrices(const ArrayMatrix< Real >& m1,
        const ArrayMatrix< Real >& m2)
    {
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
            {
                Real aux(0.0);
                for (int k = 0; k < size; k++)
                    aux += m1.matrix[i * size + k] * m2.matrix[k * size + j];
                this->matrix[i * size + j] = aux;
            }
    }

    void multiplyMatricesWithTransposition(const ArrayMatrix< Real >& m1,
        const ArrayMatrix< Real >& m2)
    {
        Real* transposedMatrix = new Real[size * size];
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                transposedMatrix[i * size + j] =
                m2.matrix[j * size + i];

        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
            {
                Real aux(0.0);
                for (int k = 0; k < size; k++)
                    aux += m1.matrix[i * size + k] * transposedMatrix[j * size + k];
                this->matrix[i * size + j] = aux;
            }

        delete[] transposedMatrix;
    }

    template< int blockSize = 16 >
    void multiplyMatricesBlockwise1(const ArrayMatrix< Real >& m1,
        const ArrayMatrix< Real >& m2)
    {
        memset(this->matrix, 0, size * size * sizeof(Real));
        for (int i = 0; i < size; i += blockSize)
            for (int j = 0; j < size; j += blockSize)
                for (int k = 0; k < size; k += blockSize)
                {
                    for (int i1 = i; i1 < i + blockSize; i1++)
                        for (int j1 = j; j1 < j + blockSize; j1++)
                            for (int k1 = k; k1 < k + blockSize; k1++) {
                                this->matrix[i1 * size + j1] +=
                                    m1.matrix[i1 * size + k1] *
                                    m2.matrix[k1 * size + j1];
                            }

                }
    }
    template< int blockSize = 16 >
    void multiplyMatricesBlockwise2(const ArrayMatrix< Real >& m1,
        const ArrayMatrix< Real >& m2)
    {
        memset(this->matrix, 0, size * size * sizeof(Real));
        for (int i = 0; i < size; i += blockSize)
            for (int j = 0; j < size; j += blockSize)
                for (int k = 0; k < size; k += blockSize)
                {
                    Real* res = &this->matrix[i * size + j];
                    Real* _m1 = &m1.matrix[i * size + k];
                    for (int i1 = 0; i1 < blockSize; i1++, res += size, _m1 += size)
                        for (int j1 = 0; j1 < blockSize; j1++)
                        {
                            Real* _m2 = &m2.matrix[k * size + j];
                            for (int k1 = 0; k1 < blockSize; k1++, _m2 += size)
                                res[j1] += _m1[k1] * _m2[j1];
                        }
                }
    }

    template< int blockSize = 16 >
    void multiplyMatricesBlockwise3(const ArrayMatrix< Real >& m1,
        const ArrayMatrix< Real >& m2)
    {
        memset(this->matrix, 0, size * size * sizeof(Real));
        for (int i = 0; i < size; i += blockSize)
            for (int j = 0; j < size; j += blockSize)
                for (int k = 0; k < size; k += blockSize)
                {
                    Real* res = &this->matrix[i * size + j];
                    Real* _m1 = &m1.matrix[i * size + k];
                    for (int i1 = 0; i1 < blockSize; i1++, res += size, _m1 += size)
                    {
                        Real* _m2 = &m2.matrix[k * size + j];
                        for (int k1 = 0; k1 < blockSize; k1++, _m2 += size)
                            for (int j1 = 0; j1 < blockSize; j1++)
                                res[j1] += _m1[k1] * _m2[j1];
                    }
                }
    }

    bool check(const ArrayMatrix< Real >& m, Real tolerance = 1.0e-6)
    {
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                if (this->matrix[i * size + j] != m.matrix[i * size + j])
                {
                    std::cerr << "Error of position: " << i << " " << j << " -> "
                        << this->matrix[i * size + j] << " x " << m.matrix[i * size + j] << std::endl;
                    std::cerr << "this = " << std::endl << *this << std::endl << "m = " << std::endl << m << std::endl;
                    return false;
                }
        return true;
    }

    void print(std::ostream& str) const
    {
        str << "Ptr = " << matrix << std::endl;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
                str << matrix[i * size + j] << "\t";
            str << std::endl;
        }
    }

    ~ArrayMatrix()
    {
        delete[] this->matrix;
    }
protected:
    Real* matrix = nullptr;

    int size;
};

template< typename Real >
std::ostream& operator<<(std::ostream& str, const ArrayMatrix< Real >& m)
{
    m.print(str);
    return str;
}


#endif // !__ARRAYMATRIX_