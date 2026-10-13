#include <csignal>
#include <iostream>
#include <set>
#include <vector>
#include <limits>

using namespace std;

#define EPS 1e-9

template<typename E>
class Matrix { // not robust to irregular matrix
private:
    int rows, cols;
    vector<E> v;
public:
    Matrix(const vector<vector<E>>& v1) {
        rows = v1.size();
        if(rows <= 0) throw 100;
        cols = v1[0].size();
        if(cols <= 0) throw 101;
        v = vector<E>(rows * cols, 0.0f);
        for(int i = 0; i < v1.size(); i++) {
            for(int j = 0; j < v1[0].size(); j++) {
                v[i * cols + j] = v1[i][j];
            }
        }
    }

    Matrix(const vector<E>& v1) {
        rows = v1.size();
        cols = 1;
        v = v1;
    }

    Matrix(int r, int c, E def = 0.0f) {
        if(r <= 0 || c <= 0) throw 102;
        rows = r;
        cols = c;
        v = vector<E>(rows * cols, def);
    }

    E& operator()(int i, int j) {
        return v[i * cols + j];
    }

    const E& operator()(int i, int j) const {
        return v[i * cols + j];
    }

    E& operator()(int i) {
        return v[i];
    }

    const E& operator()(int i) const {
        return v[i];
    }

    E& operator[](int i) {
        return v[i];
    }

    const E& operator[](int i) const {
        return v[i];
    }

    Matrix<E> T() const {
        Matrix res(cols, rows);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                res(j, i) = (*this)(i, j);
            }
        }
        return res;
    }

    Matrix<E>& operator+=(const Matrix<E>& mat) {
        if(mat.rows != rows || mat.cols != cols) throw 103;
        for(int i = 0; i < rows * cols; i++) {
            v[i] += mat.v[i];
        }
        return *this;
    }

    Matrix<E> operator+(const Matrix<E>& mat) const {
        Matrix<E> res = *this;
        res += mat;
        return res;
    }

    Matrix<E>& operator*=(E scalar) {
        for(int i = 0; i < rows * cols; i++) {
            v[i] *= scalar;
        }
        return *this;
    }

    Matrix<E> operator*(E scalar) const {
        Matrix<E> res = *this;
        res *= scalar;
        return res;
    }

    vector<E> operator*(const vector<E>& v1) const {
        if(v1.size() != cols) throw 104;
        vector<E> res(rows, 0.0f);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                res[i] += (*this)(i, j) * v1[j];
            }
        }
        return res;
    }

    Matrix<E> operator*(const Matrix<E>& mat) const {
        if(mat.rows != cols) throw 105;
        Matrix<E> res(rows, mat.cols, 0.0f);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < mat.cols; j++) {
                for(int k = 0; k < cols; k++) {
                    res(i, j) += (*this)(i, k) * mat(k, j);
                }
            }
        }
        return res;
    }

    int size() const {
        return rows;
    }

    int size(int dim) const {
        return (dim == 0) ? rows : cols;
    }

    Matrix<E> row(int j) const {
        if(j < 0 || j >= rows) throw 107;
        Matrix<E> res(1, cols);
        for(int i = 0; i < cols; i++) {
            res(i) = (*this)(j, i);
        }
        return res;
    }

    Matrix<E> col(int j) const {
        if(j < 0 || j >= cols) throw 108;
        Matrix<E> res(rows, 1);
        for(int i = 0; i < rows; i++) {
            res(i) = (*this)(i, j);
        }
        return res;
    }

    void resize(int new_rows, int new_cols, E def = 0.0f) {
        if(new_rows <= 0 || new_cols <= 0) throw 106;
        vector<E> new_v(new_rows * new_cols, def);
        for(int i = 0; i < min(rows, new_rows); i++) {
            for(int j = 0; j < min(cols, new_cols); j++) {
                new_v[i * new_cols + j] = (*this)(i, j);
            }
        }
        v = std::move(new_v);
        rows = new_rows;
        cols = new_cols;
    }

    friend ostream& operator<<(ostream& os, const Matrix<E>& mat) {
        for(int i = 0; i < mat.rows; i++) {
            for(int j = 0; j < mat.cols; j++) {
                os << mat(i, j) << " ";
            }
            os << '\n';
        }
        return os;
    }
};

template<typename T>
ostream& operator<<(ostream& os, const vector<T>& vec) {
    for(const T& el: vec) {
        os << el << ' ';
    }
    os << '\n';
    return os;
}

Matrix<double> simplex(
    const Matrix<double>& A,
    Matrix<double>& b_,
    Matrix<int>& basis_idx,
    const Matrix<double>& c,
    Matrix<double>& B1
) {
    /*
     * 1) only the first time: B (= I)
     * 2) B1, u = c[basis_idx]B1
     * 3) c_ = c - uF
     * 4) if c_[i] >= 0 for every i: STOP -> return x* = ( B1b | 0)
     * 5) find h s.t. c_[h] > 0
     * 6) b_ = B1b, A_h = B1Ah
     * 7) t = argmin{b_[i]/A_[i][h]: A_[i][h] > 0}, theta = argmin{b_[i]/A_[i][h]: A_[i][h] > 0}
     * 8) pivoting to find B1, z1, b_
     */
    int m = A.size(0), n = A.size(1);
    set<int> basis_set;
    for(int i = 0; i  < basis_idx.size(1); i++) basis_set.insert(basis_idx(i));
    Matrix<double> cB(m, 1, 0.0f);
    for(int i = 0; i < m; i++) {
        cB(i) = c(basis_idx(i));
    }
    Matrix<double> u = cB.T() * B1;
    Matrix<double> c_ = c.T() + (u * A) * -1;
    int h = -1;
    for(int i = 0; i < n; i++) {
        if(basis_set.find(i) == basis_set.end()) {
            if(c_(i) < -EPS) {
                h = i;
                break;
            }
        }
    }
    if(h == -1) {
        cout << "STOP! c_(i) <= 0 for all idxs." << endl;
        Matrix<double> x(A.size(1), 1, 0.0f);
        for(int i = 0; i < basis_idx.size(1); i++) {
            x(basis_idx(i)) = b_(i);
        }
        return x;
    }
    Matrix<double> A_h = B1 * A.col(h);
    int t = -1;
    double theta = numeric_limits<double>::infinity();;
    for(int i = 0; i < b_.size(); i++) {
        if(A_h(i) <= EPS) continue;
        double val = b_(i) / A_h(i);
        if(val < theta) {
            theta = val;
            t = i;
        }
    }
    if(t == -1) {
        cout << "STOP! Unbounded problem." << endl;
        throw 200;
    }
    // PIVOTING
    double p = A_h(t);
    for(int i = 0; i < m; i++) {
        B1(t, i) /= p;
    }
    for(int i = 0; i < A_h.size(); i++) {
        if(i != t) {
            for(int j = 0; j < m; j++) {
                B1(i, j) -= A_h(i) * B1(t, j);
            }
        }
    }
    b_(t) *= 1/p;
    for(int i = 0; i < b_.size(); i++) {
        if(i != t) b_(i) -= A_h(i) * b_(t);
    }
    basis_idx(t) = h;
    Matrix<double> x(A.size(1), 1, 0.0f);
    for(int i = 0; i < basis_idx.size(1); i++) {
        x(basis_idx(i)) = b_(i);
    }
    double z0 = (c.T() * x)(0);
    cout << "z0: " << z0 << endl;
    return simplex(A, b_, basis_idx, c, B1);
}


int main() {
    int m, n;
    cout << "Num. of equations (excluding xi >= 0): ";
    cin >> m;
    cout << "Num. of variables: ";
    cin >> n;
    Matrix<double> A(m, n, 0.0f);
    Matrix<double> b(m, 1, 0);
    cout << "Upload A below" << endl;
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            cin >> A(i, j);
        }
    }
    cout << "Upload b below" << endl;
    for(int i = 0; i < m; i++) cin >> b(i);
    // scan columns to find k*e_i (all elements == 0 except one)
    Matrix<int> basis_idx(1, m, -1);
    for(int i = 0; i < n; i++) {
        int zero_count = 0;
        int nonzero_row;
        int j;
        for(j = 0; j < m; j++) {
            if(A(j, i) == 0) zero_count++;
            else {
                if(basis_idx[j] >= 0) break;
                if(j - zero_count >= 1) break;
                nonzero_row = j;
            }
        }
        if(j == m && zero_count != m) {
            basis_idx(nonzero_row) = i;
        }
    }
    vector<int> artificial_vars;
    for(int i = 0; i < m; i++) {
        if(basis_idx(i) < 0) {
            basis_idx(i) = (n + artificial_vars.size());
            artificial_vars.push_back(i);
        }
    }
    A.resize(m, n + artificial_vars.size(), 0.0f);
    for(int i = 0; i < artificial_vars.size(); i++) {
        A(artificial_vars[i], n + i) = 1.0;
    }
    cout << A;
    cout << "Basis idx: " << basis_idx;
    // min z = ∑ x_i to remove artificial variables
    Matrix<double> c(n + artificial_vars.size(), 1, 0.0f);
    for(int i = n; i < n + artificial_vars.size(); i++) {
        c(i, 0) = 1;
    }
    Matrix<double> B1(m, m, 0.0f);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) B1(j, i) = 0.0f;
        B1(i, i) = 1 / A(i, basis_idx(i));
    }
    cout << B1;
    Matrix<double> b_ = B1 * b;
    cout << b_;
    Matrix<double> x(n + artificial_vars.size(), 1, 0.0f);
    for(int i = 0; i < basis_idx.size(1); i++) {
        x(basis_idx(i)) = b_(i);
    }
    cout << x;
    double z0 = (c.T() * x)(0);
    cout << "z0: " << z0 << endl;
    x = simplex(A, b_, basis_idx, c, B1);
    cout << "Solution removing artificial vars: " << x;
    A.resize(m, n);
    cout << "Upload c below" << endl;
    c.resize(n, 1);
    for(int i = 0; i < n; i++) {
        cin >> c(i);
    }
    x = simplex(A, b_, basis_idx, c, B1);
    cout << "Final solution: " << x;
}
