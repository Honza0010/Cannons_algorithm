#ifndef __ARRAYMATRIX_
#define __ARRAYMATRIX_



#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>


template< typename Real >
class ArrayMatrix
{
public:
    ArrayMatrix(int size)
    {
        this->matrix = new Real[size * size];
        this->size = size;
    }

    ArrayMatrix(const std::string& filename)
    {
        std::ifstream file(filename, std::ios::in);

        if (!file.is_open())
        {
            throw std::ios_base::failure("The could not be opened");
        }

        file >> this->size;
        this->matrix = new Real[size * size];

        Real val;
        for (int i = 0; i < this->size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                file >> val;
                this->setElement(i, j, val);
                if ((i < size - 1 || j < size - 1) && file.eof())
                {
                    throw std::out_of_range("Too few elements in the file");
                    file.close();
                }
                if (file.fail())
                {
                    throw std::invalid_argument("Invalid data type");
                    file.close();
                }
            }
        }

        file.close();
    }

    void randFillMatrix()
    {
        std::random_device rd;
        std::mt19937 gen(rd());

        if constexpr (std::is_integral_v<Real>) {
            // Use uniform_int_distribution for integer types
            std::uniform_int_distribution<Real> dis(0, 100); // Range [0, 100]
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    this->setElement(i, j, dis(gen));
                }
            }
        }
        else if constexpr (std::is_floating_point_v<Real>) {
            // Use uniform_real_distribution for floating-point types
            std::uniform_real_distribution<Real> dis(0.0, 100.0); // Range [0.0, 1.0]
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    this->setElement(i, j, roundToDecimalPlaces(dis(gen), 4));
                }
            }
        }
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