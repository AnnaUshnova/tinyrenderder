#include "our_gl.h"
#include <limits>
#include <algorithm>
#include <iostream>
#include <climits>
#include <cmath>

// globals
mat<4, 4> ModelView = mat<4, 4>::identity();
mat<4, 4> Perspective = mat<4, 4>::identity();
mat<4, 4> Viewport = mat<4, 4>::identity();
std::vector<double> zbuffer;

// diagnostics
static std::size_t g_fragments_drawn = 0;
static std::size_t g_triangles_rasterized = 0;
static int g_min_x = INT_MAX, g_min_y = INT_MAX, g_max_x = INT_MIN, g_max_y = INT_MIN;
static double g_min_z = std::numeric_limits<double>::infinity(), g_max_z = -std::numeric_limits<double>::infinity();

// helpers
static mat<4, 4> identity4() { return mat<4, 4>::identity(); }

// setup
void lookat(const vec3 eye, const vec3 center, const vec3 up) {
    vec3 z = normalized(eye - center);  // view direction (from center to eye)
    vec3 x = normalized(cross(up, z));
    vec3 y = cross(z, x);

    ModelView = mat<4, 4>::identity();

    // rotation
    ModelView[0][0] = x[0]; ModelView[0][1] = x[1]; ModelView[0][2] = x[2];
    ModelView[1][0] = y[0]; ModelView[1][1] = y[1]; ModelView[1][2] = y[2];
    ModelView[2][0] = z[0]; ModelView[2][1] = z[1]; ModelView[2][2] = z[2];

    // translation
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
    // Standard perspective projection (OpenGL style, NDC z in [-1,1])
    Perspective[2][2] = (zfar + znear) / (znear - zfar);
    Perspective[2][3] = (2.0 * zfar * znear) / (znear - zfar);
    Perspective[3][2] = -1.0;
    Perspective[3][3] = 0.0;
}

void init_viewport(int x, int y, int w, int h) {
    Viewport = mat<4, 4>::identity();
    Viewport[0][0] = w / 2.0;
    Viewport[1][1] = h / 2.0;
    Viewport[0][3] = x + w / 2.0;
    Viewport[1][3] = y + h / 2.0;

    // Z transform from NDC [-1,1] to [0,1] if needed:
    // Here мы чаще работаем с NDC z при тесте, поэтому Viewport z-слой оставляем нейтральным.
    // Если хотите, можно записывать post-viewport depth: Viewport[2][2] = 0.5; Viewport[2][3] = 0.5;
    Viewport[2][2] = 1.0;
    Viewport[2][3] = 0.0;
}

void init_zbuffer(int width, int height) {
    // инициализируем sentinel = +inf (означает "пусто / ничего отрисовано")
    zbuffer.assign(width * height, std::numeric_limits<double>::infinity());
}

// barycentric as cross-product method
static vec3 barycentric(const vec2& A, const vec2& B, const vec2& C, const vec2& P) {
    vec3 s0 = { C[0] - A[0], B[0] - A[0], A[0] - P[0] };
    vec3 s1 = { C[1] - A[1], B[1] - A[1], A[1] - P[1] };
    vec3 u = cross(s0, s1);
    if (std::abs(u[2]) < 1e-12) return vec3{ -1,1,1 };
    return vec3{ 1.0 - (u[0] + u[1]) / u[2], u[1] / u[2], u[0] / u[2] };
}

void rasterize(const Triangle& clip, const IShader& shader, TGAImage& framebuffer) {
    ++g_triangles_rasterized;

    // clip are clip-space vec4 (from shader.vertex)
    vec4 v0 = clip[0], v1 = clip[1], v2 = clip[2];

    // guard: if any w is zero -> cannot perspective divide reliably
    if (std::abs(v0[3]) < 1e-12 || std::abs(v1[3]) < 1e-12 || std::abs(v2[3]) < 1e-12)
        return;

    // perspective divide -> NDC
    vec4 ndc4[3] = { v0 / v0[3], v1 / v1[3], v2 / v2[3] };

    // guard against NaN/Inf after divide
    for (int i = 0; i < 3; ++i) {
        for (int c = 0; c < 4; ++c) {
            double val = ndc4[i][c];
            if (!std::isfinite(val)) return;
        }
    }

    // screen coordinates (using Viewport for x,y only)
    vec2 screen[3] = { (Viewport * ndc4[0]).xy(), (Viewport * ndc4[1]).xy(), (Viewport * ndc4[2]).xy() };

    // backface culling in screen space (consistent winding)
    vec2 edge1 = screen[1] - screen[0];
    vec2 edge2 = screen[2] - screen[0];
    double cross_product = edge1.x * edge2.y - edge1.y * edge2.x;
    if (cross_product <= 0) return;

    // bounding box
    int minx = std::max(0, (int)std::floor(std::min({ screen[0][0], screen[1][0], screen[2][0] })));
    int maxx = std::min(framebuffer.width() - 1, (int)std::ceil(std::max({ screen[0][0], screen[1][0], screen[2][0] })));
    int miny = std::max(0, (int)std::floor(std::min({ screen[0][1], screen[1][1], screen[2][1] })));
    int maxy = std::min(framebuffer.height() - 1, (int)std::ceil(std::max({ screen[0][1], screen[1][1], screen[2][1] })));

    if (minx > maxx || miny > maxy) return;

    // update bbox diagnostics
    g_min_x = std::min(g_min_x, minx);
    g_min_y = std::min(g_min_y, miny);
    g_max_x = std::max(g_max_x, maxx);
    g_max_y = std::max(g_max_y, maxy);

    // precompute clip.w for perspective-correct interpolation
    double w0 = v0[3], w1 = v1[3], w2 = v2[3];

    for (int x = minx; x <= maxx; ++x) {
        for (int y = miny; y <= maxy; ++y) {
            vec2 P{ (double)x + 0.5, (double)y + 0.5 };
            vec3 bc = barycentric(screen[0], screen[1], screen[2], P);
            if (bc[0] < 0 || bc[1] < 0 || bc[2] < 0) continue;

            // interpolate NDC z
            double z_ndc = bc[0] * ndc4[0].xyz()[2] + bc[1] * ndc4[1].xyz()[2] + bc[2] * ndc4[2].xyz()[2];

            if (!std::isfinite(z_ndc)) continue;

            int idx = x + y * framebuffer.width();

            // DEPTH TEST:
            // we store currently best z (smaller z => ближе, since NDC near=-1 < far=1)
            // sentinel is +inf so first fragment always passes
            if (!(z_ndc < zbuffer[idx])) continue;

            // Perspective-correct barycentric: divide by clip.w (v.w)
            double invw0 = (std::abs(w0) > 1e-12) ? (1.0 / w0) : 0.0;
            double invw1 = (std::abs(w1) > 1e-12) ? (1.0 / w1) : 0.0;
            double invw2 = (std::abs(w2) > 1e-12) ? (1.0 / w2) : 0.0;

            double denom = bc[0] * invw0 + bc[1] * invw1 + bc[2] * invw2;

            vec3 pcbar;
            if (std::abs(denom) < 1e-15) {
                // вырожденный случай — упрощаем до обычных барицентриков
                pcbar = bc;
            }
            else {
                pcbar[0] = (bc[0] * invw0) / denom;
                pcbar[1] = (bc[1] * invw1) / denom;
                pcbar[2] = (bc[2] * invw2) / denom;
            }

            auto [discard, color] = shader.fragment(pcbar);
            if (discard) continue;

            // write depth & color
            zbuffer[idx] = z_ndc;
            framebuffer.set(x, y, color);

            ++g_fragments_drawn;

            // update z-range diagnostics
            g_min_z = std::min(g_min_z, z_ndc);
            g_max_z = std::max(g_max_z, z_ndc);

            // optional: early debug print for first few triangles/pixels
            if (g_triangles_rasterized < 10 && g_fragments_drawn < 20) {
                std::cerr << "Debug: tri#" << g_triangles_rasterized << " px(" << x << "," << y << ") z=" << z_ndc
                    << " prev_z=" << zbuffer[idx] << " fragments=" << g_fragments_drawn << "\n";
            }
        }
    }
}

void print_render_stats() {
    std::cerr << "DEBUG: triangles=" << g_triangles_rasterized
        << " fragments_drawn=" << g_fragments_drawn
        << " bbox=[" << g_min_x << "," << g_min_y << "] - [" << g_max_x << "," << g_max_y << "]"
        << " z-range=[" << (std::isfinite(g_min_z) ? std::to_string(g_min_z) : "inf") << ","
        << (std::isfinite(g_max_z) ? std::to_string(g_max_z) : "-inf") << "]\n";
}
