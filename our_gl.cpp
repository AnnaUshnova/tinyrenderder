// ================================================================
// our_gl.cpp - Реализация OpenGL-подобного конвейера
// ================================================================
#include "our_gl.h"
#include <limits>
#include <algorithm>
#include <iostream>
#include <climits>
#include <cmath>

// Глобальные переменные
mat<4, 4> ModelView = mat<4, 4>::identity();
mat<4, 4> Perspective = mat<4, 4>::identity();
mat<4, 4> Viewport = mat<4, 4>::identity();
std::vector<double> zbuffer;

// Диагностические счётчики
static std::size_t fragments_drawn = 0;
static std::size_t triangles_rasterized = 0;
static int min_x = INT_MAX, min_y = INT_MAX, max_x = INT_MIN, max_y = INT_MIN;
static double min_z = std::numeric_limits<double>::infinity();
static double max_z = -std::numeric_limits<double>::infinity();

// Инициализация камеры (look-at матрица)
void lookat(const vec3 eye, const vec3 center, const vec3 up) {
    vec3 z_axis = normalized(eye - center);  // Направление взгляда (к камере)
    vec3 x_axis = normalized(cross(up, z_axis));
    vec3 y_axis = cross(z_axis, x_axis);

    ModelView = mat<4, 4>::identity();

    // Поворот
    ModelView[0][0] = x_axis[0]; ModelView[0][1] = x_axis[1]; ModelView[0][2] = x_axis[2];
    ModelView[1][0] = y_axis[0]; ModelView[1][1] = y_axis[1]; ModelView[1][2] = y_axis[2];
    ModelView[2][0] = z_axis[0]; ModelView[2][1] = z_axis[1]; ModelView[2][2] = z_axis[2];

    // Трансляция
    ModelView[0][3] = -dot(x_axis, eye);
    ModelView[1][3] = -dot(y_axis, eye);
    ModelView[2][3] = -dot(z_axis, eye);
}

// Инициализация перспективной проекции
void init_perspective(double fov_deg, double aspect, double znear, double zfar) {
    double fov_rad = fov_deg * M_PI / 180.0;
    double tan_half_fov = std::tan(fov_rad / 2.0);

    Perspective = mat<4, 4>::identity();
    Perspective[0][0] = 1.0 / (aspect * tan_half_fov);
    Perspective[1][1] = 1.0 / tan_half_fov;
    // Стандартная перспективная проекция (стиль OpenGL, NDC z в [-1,1])
    Perspective[2][2] = (zfar + znear) / (znear - zfar);
    Perspective[2][3] = (2.0 * zfar * znear) / (znear - zfar);
    Perspective[3][2] = -1.0;
    Perspective[3][3] = 0.0;
}

// Инициализация видового экрана
void init_viewport(int x, int y, int w, int h) {
    Viewport = mat<4, 4>::identity();
    Viewport[0][0] = w / 2.0;
    Viewport[1][1] = h / 2.0;
    Viewport[0][3] = x + w / 2.0;
    Viewport[1][3] = y + h / 2.0;

    // Z-преобразование можно настроить здесь при необходимости
    Viewport[2][2] = 1.0;
    Viewport[2][3] = 0.0;
}

// Инициализация Z-буфера
void init_zbuffer(int width, int height) {
    zbuffer.assign(width * height, std::numeric_limits<double>::infinity());
}

// Барицентрические координаты (метод векторного произведения)
static vec3 barycentric(const vec2& A, const vec2& B, const vec2& C, const vec2& P) {
    vec3 s0 = { C[0] - A[0], B[0] - A[0], A[0] - P[0] };
    vec3 s1 = { C[1] - A[1], B[1] - A[1], A[1] - P[1] };
    vec3 u = cross(s0, s1);

    if (std::abs(u[2]) < 1e-12)
        return vec3{ -1, 1, 1 };  // Вырожденный треугольник

    return vec3{ 1.0 - (u[0] + u[1]) / u[2], u[1] / u[2], u[0] / u[2] };
}

// Основная функция растризации треугольника
void rasterize(const Triangle& clip, const IShader& shader, TGAImage& framebuffer) {
    ++triangles_rasterized;

    // Вершины в clip-space (из шейдера)
    vec4 v0 = clip[0], v1 = clip[1], v2 = clip[2];

    // Проверка: если какая-либо w == 0, нельзя делить
    if (std::abs(v0[3]) < 1e-12 || std::abs(v1[3]) < 1e-12 || std::abs(v2[3]) < 1e-12)
        return;

    // Перспективное деление -> NDC
    vec4 ndc4[3] = { v0 / v0[3], v1 / v1[3], v2 / v2[3] };

    // Проверка на NaN/Inf после деления
    for (int i = 0; i < 3; ++i) {
        for (int c = 0; c < 4; ++c) {
            double val = ndc4[i][c];
            if (!std::isfinite(val)) return;
        }
    }

    // Экранные координаты (только x, y через Viewport)
    vec2 screen[3] = {
        (Viewport * ndc4[0]).xy(),
        (Viewport * ndc4[1]).xy(),
        (Viewport * ndc4[2]).xy()
    };

    // Отсечение задних граней в экранном пространстве
    vec2 edge1 = screen[1] - screen[0];
    vec2 edge2 = screen[2] - screen[0];
    double cross_product = edge1.x * edge2.y - edge1.y * edge2.x;
    if (cross_product <= 0) return;

    // Ограничивающий прямоугольник
    int min_x_px = std::max(0, (int)std::floor(std::min({ screen[0][0], screen[1][0], screen[2][0] })));
    int max_x_px = std::min(framebuffer.width() - 1, (int)std::ceil(std::max({ screen[0][0], screen[1][0], screen[2][0] })));
    int min_y_px = std::max(0, (int)std::floor(std::min({ screen[0][1], screen[1][1], screen[2][1] })));
    int max_y_px = std::min(framebuffer.height() - 1, (int)std::ceil(std::max({ screen[0][1], screen[1][1], screen[2][1] })));

    if (min_x_px > max_x_px || min_y_px > max_y_px) return;

    // Обновление диагностики ограничивающего прямоугольника
    min_x = std::min(min_x, min_x_px);
    min_y = std::min(min_y, min_y_px);
    max_x = std::max(max_x, max_x_px);
    max_y = std::max(max_y, max_y_px);

    // Предвычисление clip.w для перспективно-корректной интерполяции
    double w0 = v0[3], w1 = v1[3], w2 = v2[3];

    // Растеризация по ограничивающему прямоугольнику
    for (int x = min_x_px; x <= max_x_px; ++x) {
        for (int y = min_y_px; y <= max_y_px; ++y) {
            vec2 pixel_center{ (double)x + 0.5, (double)y + 0.5 };
            vec3 barycentric_coords = barycentric(screen[0], screen[1], screen[2], pixel_center);

            if (barycentric_coords[0] < 0 || barycentric_coords[1] < 0 || barycentric_coords[2] < 0)
                continue;

            // Интерполяция NDC z
            double z_ndc = barycentric_coords[0] * ndc4[0].xyz()[2]
                + barycentric_coords[1] * ndc4[1].xyz()[2]
                + barycentric_coords[2] * ndc4[2].xyz()[2];

            if (!std::isfinite(z_ndc)) continue;

            int idx = x + y * framebuffer.width();

            // Z-тест: меньший z означает ближе (в NDC near=-1 < far=1)
            if (!(z_ndc < zbuffer[idx])) continue;

            // Перспективно-корректные барицентрические координаты
            double inv_w0 = (std::abs(w0) > 1e-12) ? (1.0 / w0) : 0.0;
            double inv_w1 = (std::abs(w1) > 1e-12) ? (1.0 / w1) : 0.0;
            double inv_w2 = (std::abs(w2) > 1e-12) ? (1.0 / w2) : 0.0;

            double denom = barycentric_coords[0] * inv_w0
                + barycentric_coords[1] * inv_w1
                + barycentric_coords[2] * inv_w2;

            vec3 perspective_correct_bary;
            if (std::abs(denom) < 1e-15) {
                // Вырожденный случай — используем обычные барицентрические координаты
                perspective_correct_bary = barycentric_coords;
            }
            else {
                perspective_correct_bary[0] = (barycentric_coords[0] * inv_w0) / denom;
                perspective_correct_bary[1] = (barycentric_coords[1] * inv_w1) / denom;
                perspective_correct_bary[2] = (barycentric_coords[2] * inv_w2) / denom;
            }

            auto [discard, color] = shader.fragment(perspective_correct_bary);
            if (discard) continue;

            // Запись глубины и цвета
            zbuffer[idx] = z_ndc;
            framebuffer.set(x, y, color);

            ++fragments_drawn;

            // Обновление диагностики диапазона z
            min_z = std::min(min_z, z_ndc);
            max_z = std::max(max_z, z_ndc);
        }
    }
}

// Вывод статистики рендеринга
void print_render_stats() {
    std::cerr << "DEBUG: triangles=" << triangles_rasterized
        << " fragments_drawn=" << fragments_drawn
        << " bbox=[" << min_x << "," << min_y << "] - [" << max_x << "," << max_y << "]"
        << " z-range=[" << (std::isfinite(min_z) ? std::to_string(min_z) : "inf") << ","
        << (std::isfinite(max_z) ? std::to_string(max_z) : "-inf") << "]\n";
}

Frustum Frustum::createFromMatrix(const mat<4, 4>& matrix) {
    Frustum frustum;

    // Извлечение плоскостей из матрицы
    // Left plane
    frustum.planes[LEFT].normal.x = matrix[0][3] + matrix[0][0];
    frustum.planes[LEFT].normal.y = matrix[1][3] + matrix[1][0];
    frustum.planes[LEFT].normal.z = matrix[2][3] + matrix[2][0];
    frustum.planes[LEFT].d = matrix[3][3] + matrix[3][0];

    // Right plane
    frustum.planes[RIGHT].normal.x = matrix[0][3] - matrix[0][0];
    frustum.planes[RIGHT].normal.y = matrix[1][3] - matrix[1][0];
    frustum.planes[RIGHT].normal.z = matrix[2][3] - matrix[2][0];
    frustum.planes[RIGHT].d = matrix[3][3] - matrix[3][0];

    // Bottom plane
    frustum.planes[BOTTOM].normal.x = matrix[0][3] + matrix[0][1];
    frustum.planes[BOTTOM].normal.y = matrix[1][3] + matrix[1][1];
    frustum.planes[BOTTOM].normal.z = matrix[2][3] + matrix[2][1];
    frustum.planes[BOTTOM].d = matrix[3][3] + matrix[3][1];

    // Top plane
    frustum.planes[TOP].normal.x = matrix[0][3] - matrix[0][1];
    frustum.planes[TOP].normal.y = matrix[1][3] - matrix[1][1];
    frustum.planes[TOP].normal.z = matrix[2][3] - matrix[2][1];
    frustum.planes[TOP].d = matrix[3][3] - matrix[3][1];

    // Near plane
    frustum.planes[NEAR].normal.x = matrix[0][3] + matrix[0][2];
    frustum.planes[NEAR].normal.y = matrix[1][3] + matrix[1][2];
    frustum.planes[NEAR].normal.z = matrix[2][3] + matrix[2][2];
    frustum.planes[NEAR].d = matrix[3][3] + matrix[3][2];

    // Far plane
    frustum.planes[FAR].normal.x = matrix[0][3] - matrix[0][2];
    frustum.planes[FAR].normal.y = matrix[1][3] - matrix[1][2];
    frustum.planes[FAR].normal.z = matrix[2][3] - matrix[2][2];
    frustum.planes[FAR].d = matrix[3][3] - matrix[3][2];

    // Нормализация плоскостей
    for (int i = 0; i < 6; i++) {
        double length = norm(frustum.planes[i].normal);
        if (length > 0.0) {
            frustum.planes[i].normal = frustum.planes[i].normal / length;
            frustum.planes[i].d /= length;
        }
    }

    return frustum;
}

bool Frustum::intersects(const AABB& aabb) const {
    for (int i = 0; i < 6; i++) {
        const Plane& plane = planes[i];

        // Положительная вершина (по нормали)
        vec3 positive = aabb.min;
        if (plane.normal.x >= 0) positive.x = aabb.max.x;
        if (plane.normal.y >= 0) positive.y = aabb.max.y;
        if (plane.normal.z >= 0) positive.z = aabb.max.z;

        // Если положительная вершина снаружи плоскости, AABB снаружи
        if (plane.distance(positive) < 0) {
            return false;
        }
    }
    return true;
}
