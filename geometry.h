#pragma once
#include <cmath>
#include <cassert>
#include <iostream>

// ================================================================
//  Generic vec<n>
// ================================================================
template<int n>
struct vec {
    double data[n] = { 0 };

    double& operator[](const int i) {
        assert(i >= 0 && i < n);
        return data[i];
    }

    double operator[](const int i) const {
        assert(i >= 0 && i < n);
        return data[i];
    }
};

// vec2 specialization
template<>
struct vec<2> {
    double x = 0, y = 0;

    double& operator[](const int i) {
        assert(i >= 0 && i < 2);
        return (i == 0 ? x : y);
    }
    double operator[](const int i) const {
        assert(i >= 0 && i < 2);
        return (i == 0 ? x : y);
    }
};

// vec3 specialization
template<>
struct vec<3> {
    double x = 0, y = 0, z = 0;

    double& operator[](const int i) {
        assert(i >= 0 && i < 3);
        return (i == 0 ? x : (i == 1 ? y : z));
    }
    double operator[](const int i) const {
        assert(i >= 0 && i < 3);
        return (i == 0 ? x : (i == 1 ? y : z));
    }
};

// vec4 specialization (compatible with generic vec storage)
template<>
struct vec<4> {
    double data[4] = { 0,0,0,0 };

    double& operator[](int i) {
        assert(i >= 0 && i < 4);
        return data[i];
    }
    double operator[](int i) const {
        assert(i >= 0 && i < 4);
        return data[i];
    }

    // convenience accessors (note: not references to avoid aliasing issues)
    double x() const { return data[0]; }
    double y() const { return data[1]; }
    double z() const { return data[2]; }
    double w() const { return data[3]; }

    vec<2> xy() const { return vec<2>{ data[0], data[1] }; }
    vec<3> xyz() const { return vec<3>{ data[0], data[1], data[2] }; }
};

typedef vec<2> vec2;
typedef vec<3> vec3;
typedef vec<4> vec4;

// ================================================================
//  Basic vec operations
// ================================================================
template<int n>
vec<n> operator+(const vec<n>& a, const vec<n>& b) {
    vec<n> r;
    for (int i = 0; i < n; ++i) r[i] = a[i] + b[i];
    return r;
}

template<int n>
vec<n> operator-(const vec<n>& a, const vec<n>& b) {
    vec<n> r;
    for (int i = 0; i < n; ++i) r[i] = a[i] - b[i];
    return r;
}

template<int n>
vec<n> operator*(const vec<n>& a, double v) {
    vec<n> r;
    for (int i = 0; i < n; ++i) r[i] = a[i] * v;
    return r;
}

template<int n>
vec<n> operator*(double v, const vec<n>& a) { return a * v; }

template<int n>
vec<n> operator/(const vec<n>& a, double v) {
    vec<n> r;
    for (int i = 0; i < n; ++i) r[i] = a[i] / v;
    return r;
}

// dot (generic)
template<int n>
double dot(const vec<n>& a, const vec<n>& b) {
    double s = 0;
    for (int i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

// length
template<int n>
double norm(const vec<n>& a) {
    return std::sqrt(dot(a, a));
}

// normalized for vec3
inline vec3 normalized(const vec3& v) {
    double l = norm<3>(v);
    if (l == 0) return v;
    return v / l;
}

// cross product for vec3
inline vec3 cross(const vec3& a, const vec3& b) {
    return vec3{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

// ================================================================
//  Matrix
// ================================================================
template<int nrows, int ncols>
struct mat {
    vec<ncols> rows[nrows];

    vec<ncols>& operator[](int r) { assert(r >= 0 && r < nrows); return rows[r]; }
    const vec<ncols>& operator[](int r) const { assert(r >= 0 && r < nrows); return rows[r]; }

    static mat<nrows, ncols> identity() {
        mat<nrows, ncols> M;
        for (int r = 0; r < nrows; ++r)
            for (int c = 0; c < ncols; ++c)
                M[r][c] = (r == c ? 1.0 : 0.0);
        return M;
    }

    mat<ncols, nrows> transpose() const {
        mat<ncols, nrows> R;
        for (int r = 0; r < nrows; ++r)
            for (int c = 0; c < ncols; ++c)
                R[c][r] = rows[r][c];
        return R;
    }
};

// mat * vec
template<int nrows, int ncols>
vec<nrows> operator*(const mat<nrows, ncols>& M, const vec<ncols>& v) {
    vec<nrows> r;
    for (int i = 0; i < nrows; ++i)
        r[i] = dot<ncols>(M[i], v);
    return r;
}

// mat * mat
template<int R1, int C1, int C2>
mat<R1, C2> operator*(const mat<R1, C1>& A, const mat<C1, C2>& B) {
    mat<R1, C2> R;
    for (int i = 0; i < R1; ++i)
        for (int j = 0; j < C2; ++j) {
            R[i][j] = 0;
            for (int k = 0; k < C1; ++k)
                R[i][j] += A[i][k] * B[k][j];
        }
    return R;
}

// stream output
template<int r, int c>
std::ostream& operator<<(std::ostream& out, const mat<r, c>& M) {
    for (int i = 0; i < r; ++i) out << M[i] << "\n";
    return out;
}
