// ================================================================
// our_gl.h - Основной интерфейс OpenGL-подобного конвейера
// ================================================================
#pragma once
#include "geometry.h"
#include "tgaimage.h"
#include <vector>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ================================================================
//  Глобальные матрицы и буфер глубины
// ================================================================
extern mat<4, 4> ModelView;     // Модельно-видовая матрица
extern mat<4, 4> Perspective;   // Матрица проекции
extern mat<4, 4> Viewport;      // Матрица видового экрана
extern std::vector<double> zbuffer;  // Z-буфер

// ================================================================
//  Функции инициализации / утилиты
// ================================================================
void lookat(const vec3 eye, const vec3 center, const vec3 up);
void init_perspective(double fov_deg, double aspect, double znear, double zfar);
void init_viewport(int x, int y, int w, int h);

// Инициализация Z-буфера размером width x height
// Буфер заполняется std::numeric_limits<double>::infinity() (значение "пусто")
void init_zbuffer(int width, int height);

// ================================================================
//  Интерфейс шейдера / растризации
// ================================================================
struct IShader {
    // Вспомогательная функция: выборка из 2D текстуры (uv в [0,1])
    static TGAColor sample2D(const TGAImage& img, const vec2& uv) {
        int x = std::min<int>(img.width() - 1,
            std::max<int>(0, int(uv.x * img.width())));
        int y = std::min<int>(img.height() - 1,
            std::max<int>(0, int(uv.y * img.height())));
        return img.get(x, y);
    }

    // Вершинный шейдер
    virtual vec4 vertex(int face_index, int vertex_index) { return vec4(); }

    // Фрагментный шейдер должен вернуть pair<discard, color>
    // Если discard == true, пиксель не записывается
    virtual std::pair<bool, TGAColor> fragment(const vec3 barycentric) const = 0;
};

// Треугольник в clip-space (4D координаты)
typedef vec<4> Triangle[3];

// Основная функция растризации треугольника
void rasterize(const Triangle& clip, const IShader& shader, TGAImage& framebuffer);

// Отладочная информация
void print_render_stats();