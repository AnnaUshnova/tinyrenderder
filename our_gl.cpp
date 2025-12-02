#include "our_gl.h"
#include <limits>
#include <algorithm>
#include <iostream>
#include <climits>

// globals
mat<4, 4> ModelView = mat<4, 4>::identity();
mat<4, 4> Perspective = mat<4, 4>::identity();
mat<4, 4> Viewport = mat<4, 4>::identity();
std::vector<double> zbuffer;

// diagnostics
static std::size_t g_fragments_drawn = 0;
static int g_min_x = INT_MAX, g_min_y = INT_MAX, g_max_x = INT_MIN, g_max_y = INT_MIN;
static double g_min_z = 1e9, g_max_z = -1e9;

// helpers
static mat<4, 4> identity4() { return mat<4, 4>::identity(); }

// setup
void lookat(const vec3 eye, const vec3 center, const vec3 up) {
    vec3 z = normalized(eye - center);  // Направление ВИДЕНИЯ (от центра к глазу)
    vec3 x = normalized(cross(up, z));
    vec3 y = cross(z, x);

    ModelView = mat<4, 4>::identity();

    // Поворот
    ModelView[0][0] = x[0]; ModelView[0][1] = x[1]; ModelView[0][2] = x[2];
    ModelView[1][0] = y[0]; ModelView[1][1] = y[1]; ModelView[1][2] = y[2];
    ModelView[2][0] = z[0]; ModelView[2][1] = z[1]; ModelView[2][2] = z[2];

    // Смещение (камера в начале координат)
    ModelView[0][3] = -dot(x, eye);
    ModelView[1][3] = -dot(y, eye);
    ModelView[2][3] = -dot(z, eye);
}

void init_perspective(double fov_deg, double aspect, double znear, double zfar) {
    double fov = fov_deg * M_PI / 180.0;
    double t = std::tan(fov / 2.0);

    Perspective = mat<4, 4>::identity();
    Perspective[0][0] = 1.0 / (aspect * t);
    Perspective[1][1] = 1.0 / t;
    Perspective[2][2] = (zfar + znear) / (znear - zfar);  // ИЗМЕНЕНО!
    Perspective[2][3] = (2.0 * zfar * znear) / (znear - zfar);  // ИЗМЕНЕНО!
    Perspective[3][2] = -1.0;
    Perspective[3][3] = 0.0;

}

void init_viewport(int x, int y, int w, int h) {
    Viewport = mat<4, 4>::identity();
    Viewport[0][0] = w / 2.0;
    Viewport[1][1] = h / 2.0;
    Viewport[0][3] = x + w / 2.0;
    Viewport[1][3] = y + h / 2.0;

    // Z преобразование: из [-1, 1] в [0, depth_range]
    // Используем 0.5 чтобы было в середине диапазона [0,1]
    Viewport[2][2] = 0.5;    // масштабирование
    Viewport[2][3] = 0.5;    // смещение
}

void init_zbuffer(int width, int height) {
    zbuffer.assign(width * height, 1.0);  // ИЗМЕНЕНО!
}

// barycentric as cross-product method
static vec3 barycentric(const vec2& A, const vec2& B, const vec2& C, const vec2& P) {
    vec3 s0 = { C[0] - A[0], B[0] - A[0], A[0] - P[0] };
    vec3 s1 = { C[1] - A[1], B[1] - A[1], A[1] - P[1] };
    vec3 u = cross(s0, s1);
    if (std::abs(u[2]) < 1e-9) return vec3{ -1,1,1 };
    return vec3{ 1.0 - (u[0] + u[1]) / u[2], u[1] / u[2], u[0] / u[2] };
}

void rasterize(const Triangle& clip, const IShader& shader, TGAImage& framebuffer) {

    static int triangle_count = 0;
    triangle_count++;


    // clip are clip-space vec4 (from shader.vertex)
    vec4 v0 = clip[0], v1 = clip[1], v2 = clip[2];

    // perspective divide -> NDC
    vec4 ndc4[3] = { v0 / v0[3], v1 / v1[3], v2 / v2[3] };
    vec2 screen[3] = { (Viewport * ndc4[0]).xy(), (Viewport * ndc4[1]).xy(), (Viewport * ndc4[2]).xy() };

    vec2 edge1 = screen[1] - screen[0];
    vec2 edge2 = screen[2] - screen[0];
    double cross_product = edge1.x * edge2.y - edge1.y * edge2.x;
    if (cross_product <= 0) return;


    // bounding box
    int minx = std::max(0, (int)std::floor(std::min({ screen[0][0], screen[1][0], screen[2][0] })));
    int maxx = std::min(framebuffer.width() - 1, (int)std::ceil(std::max({ screen[0][0], screen[1][0], screen[2][0] })));
    int miny = std::max(0, (int)std::floor(std::min({ screen[0][1], screen[1][1], screen[2][1] })));
    int maxy = std::min(framebuffer.height() - 1, (int)std::ceil(std::max({ screen[0][1], screen[1][1], screen[2][1] })));

    if (triangle_count < 10) {
        std::cerr << "Triangle " << triangle_count << ": ";
        std::cerr << "screen[0] = (" << screen[0][0] << ", " << screen[0][1] << "), ";
        std::cerr << "bbox = [" << minx << "," << miny << "] to [" << maxx << "," << maxy << "]" << std::endl;
    }

    if (minx > maxx || miny > maxy) return;

    for (int x = minx; x <= maxx; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            vec2 P{ (double)x + 0.5, (double)y + 0.5 };
            vec3 bc = barycentric(screen[0], screen[1], screen[2], P);
            if (bc[0] < 0 || bc[1] < 0 || bc[2] < 0) continue;

            // interpolate depth (use ndc z)
            double z = bc[0] * ndc4[0].xyz()[2] + bc[1] * ndc4[1].xyz()[2] + bc[2] * ndc4[2].xyz()[2];
            int idx = x + y * framebuffer.width();

            if (triangle_count < 10 && x == minx && y == miny) {
                std::cerr << "First pixel: z = " << z
                    << ", zbuffer = " << zbuffer[idx]
                    << ", ndc z values: "
                        << ndc4[0].xyz()[2] << ", "
                        << ndc4[1].xyz()[2] << ", "
                        << ndc4[2].xyz()[2] << std::endl;
            }

            // ИЗМЕНЕНИЕ ЗДЕСЬ: в OpenGL NDC z [-1, 1], где -1 = ближе
            // Мы хотим, чтобы БОЛЬШЕЕ значение z (ближе к 1) было ДАЛЬШЕ
            if (z >= zbuffer[idx]) continue;  // МЕНЯЕМ знак сравнения!

            auto [discard, color] = shader.fragment(bc);
            if (discard) continue;

            zbuffer[idx] = z;
            framebuffer.set(x, y, color);

            // diagnostics
            if (triangle_count < 10 && x == minx && y == miny) {
                std::cerr << "First pixel: z = " << z << ", zbuffer = " << zbuffer[idx] << std::endl;
            }
        }
    }
}

void print_render_stats() {
    std::cerr << "DEBUG: fragments=" << g_fragments_drawn
        << " bbox=[" << g_min_x << "," << g_min_y << "] - ["
        << g_max_x << "," << g_max_y << "]"
        << " z-range=[" << g_min_z << "," << g_max_z << "]\n";
}
