// ================================================================
// geometry.h - Геометрические типы и операции (векторы, матрицы)
// ================================================================
#pragma once
#include <cmath>
#include <cassert>
#include <iostream>
#include <array>  // Добавлено для std::array

// ================================================================
//  Generic vec<n> - Обобщённый вектор размерности n
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

// Специализация для vec2
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

// Специализация для vec3
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

// Специализация для vec4 (совместима с generic хранилищем)
template<>
struct vec<4> {
    double data[4] = { 0, 0, 0, 0 };

    double& operator[](int i) {
        assert(i >= 0 && i < 4);
        return data[i];
    }
    double operator[](int i) const {
        assert(i >= 0 && i < 4);
        return data[i];
    }

    // Удобные методы доступа (возвращают копии, чтобы избежать проблем с алиасингом)
    double x() const { return data[0]; }
    double y() const { return data[1]; }
    double z() const { return data[2]; }
    double w() const { return data[3]; }

    vec<2> xy() const { return vec<2>{ data[0], data[1] }; }
    vec<3> xyz() const { return vec<3>{ data[0], data[1], data[2] }; }
};

// Псевдонимы для удобства
typedef vec<2> vec2;
typedef vec<3> vec3;
typedef vec<4> vec4;

// ================================================================
//  Базовые операции с векторами
// ================================================================
template<int n>
vec<n> operator+(const vec<n>& a, const vec<n>& b) {
    vec<n> result;
    for (int i = 0; i < n; ++i) result[i] = a[i] + b[i];
    return result;
}

template<int n>
vec<n> operator-(const vec<n>& a, const vec<n>& b) {
    vec<n> result;
    for (int i = 0; i < n; ++i) result[i] = a[i] - b[i];
    return result;
}

template<int n>
vec<n> operator*(const vec<n>& a, double scalar) {
    vec<n> result;
    for (int i = 0; i < n; ++i) result[i] = a[i] * scalar;
    return result;
}

template<int n>
vec<n> operator*(double scalar, const vec<n>& a) { return a * scalar; }

template<int n>
vec<n> operator/(const vec<n>& a, double scalar) {
    vec<n> result;
    for (int i = 0; i < n; ++i) result[i] = a[i] / scalar;
    return result;
}

// Скалярное произведение (обобщённое)
template<int n>
double dot(const vec<n>& a, const vec<n>& b) {
    double sum = 0;
    for (int i = 0; i < n; ++i) sum += a[i] * b[i];
    return sum;
}

// Длина вектора
template<int n>
double norm(const vec<n>& a) {
    return std::sqrt(dot(a, a));
}

// Нормализация для vec3
inline vec3 normalized(const vec3& v) {
    double length = norm<3>(v);
    if (length == 0) return v;
    return v / length;
}

// Векторное произведение для vec3
inline vec3 cross(const vec3& a, const vec3& b) {
    return vec3{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

// ================================================================
//  Матрицы
// ================================================================
template<int n_rows, int n_cols>
struct mat {
    vec<n_cols> rows[n_rows];

    vec<n_cols>& operator[](int row) {
        assert(row >= 0 && row < n_rows);
        return rows[row];
    }

    const vec<n_cols>& operator[](int row) const {
        assert(row >= 0 && row < n_rows);
        return rows[row];
    }

    static mat<n_rows, n_cols> identity() {
        mat<n_rows, n_cols> result;
        for (int r = 0; r < n_rows; ++r)
            for (int c = 0; c < n_cols; ++c)
                result[r][c] = (r == c ? 1.0 : 0.0);
        return result;
    }

    mat<n_cols, n_rows> transpose() const {
        mat<n_cols, n_rows> result;
        for (int r = 0; r < n_rows; ++r)
            for (int c = 0; c < n_cols; ++c)
                result[c][r] = rows[r][c];
        return result;
    }
};

// Умножение матрицы на вектор
template<int n_rows, int n_cols>
vec<n_rows> operator*(const mat<n_rows, n_cols>& M, const vec<n_cols>& v) {
    vec<n_rows> result;
    for (int i = 0; i < n_rows; ++i)
        result[i] = dot<n_cols>(M[i], v);
    return result;
}

// Умножение матрицы на матрицу
template<int R1, int C1, int C2>
mat<R1, C2> operator*(const mat<R1, C1>& A, const mat<C1, C2>& B) {
    mat<R1, C2> result;
    for (int i = 0; i < R1; ++i)
        for (int j = 0; j < C2; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < C1; ++k)
                result[i][j] += A[i][k] * B[k][j];
        }
    return result;
}

// Вывод матрицы в поток
template<int r, int c>
std::ostream& operator<<(std::ostream& out, const mat<r, c>& M) {
    for (int i = 0; i < r; ++i) out << M[i] << "\n";
    return out;
}

// ----------------------------------------------------------------
// Вспомогательные конструкторы и унарный минус для векторов
// ----------------------------------------------------------------

inline vec2 make_vec2(double x, double y) {
    vec2 result;
    result.x = x;
    result.y = y;
    return result;
}

inline vec3 make_vec3(double x, double y, double z) {
    vec3 result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
}

inline vec4 make_vec4(double x, double y, double z, double w) {
    vec4 result;
    result[0] = x;
    result[1] = y;
    result[2] = z;
    result[3] = w;
    return result;
}

// Унарный минус для любого vec<n>
template<int n>
vec<n> operator-(const vec<n>& a) {
    return a * -1.0;
}

// ================================================================
// Frustum Culling - Плоскости и AABB
// ================================================================

// Плоскость в виде Ax + By + Cz + D = 0
struct Plane {
    vec3 normal;
    double d;

    Plane() : normal(vec3{ 0,0,1 }), d(0) {}
    Plane(const vec3& n, const vec3& point) {
        normal = normalized(n);
        d = -dot(normal, point);
    }

    // Расстояние от точки до плоскости
    double distance(const vec3& point) const {
        return dot(normal, point) + d;
    }
};

// Axis-Aligned Bounding Box
struct AABB {
    vec3 min;
    vec3 max;

    AABB() : min(vec3{ 0,0,0 }), max(vec3{ 0,0,0 }) {}
    AABB(const vec3& min_val, const vec3& max_val) : min(min_val), max(max_val) {}

    vec3 getCenter() const {
        return (min + max) * 0.5;
    }

    vec3 getSize() const {
        return max - min;
    }

    vec3 getHalfSize() const {
        return getSize() * 0.5;
    }

    // Проверка пересечения с другим AABB
    bool intersects(const AABB& other) const {
        return (min.x <= other.max.x && max.x >= other.min.x) &&
            (min.y <= other.max.y && max.y >= other.min.y) &&
            (min.z <= other.max.z && max.z >= other.min.z);
    }

    // Трансформация AABB матрицей (приблизительно)
    AABB transform(const mat<4, 4>& matrix) const {
        // Массив углов AABB
        std::array<vec3, 8> corners;
        corners[0] = vec3{min.x, min.y, min.z};
        corners[1] = vec3{max.x, min.y, min.z};
        corners[2] = vec3{min.x, max.y, min.z};
        corners[3] = vec3{max.x, max.y, min.z};
        corners[4] = vec3{min.x, min.y, max.z};
        corners[5] = vec3{max.x, min.y, max.z};
        corners[6] = vec3{min.x, max.y, max.z};
        corners[7] = vec3{max.x, max.y, max.z};

        vec3 new_min = vec3{ 1e9, 1e9, 1e9 };
        vec3 new_max = vec3{ -1e9, -1e9, -1e9 };

        for (int i = 0; i < 8; i++) {
            const vec3& corner = corners[i];
            vec4 transformed = matrix * make_vec4(corner.x, corner.y, corner.z, 1.0);
            vec3 pos = transformed.xyz() / transformed.w();

            new_min.x = std::min(new_min.x, pos.x);
            new_min.y = std::min(new_min.y, pos.y);
            new_min.z = std::min(new_min.z, pos.z);

            new_max.x = std::max(new_max.x, pos.x);
            new_max.y = std::max(new_max.y, pos.y);
            new_max.z = std::max(new_max.z, pos.z);
        }

        return AABB(new_min, new_max);
    }
};