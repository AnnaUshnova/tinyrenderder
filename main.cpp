// main.cpp — Phong render + SSAO (AO) post-process
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

#include "our_gl.h"
#include "model.h"
#include "tgaimage.h"
#include "camera.h"

// externs from your GL helpers
extern mat<4, 4> ModelView, Perspective;
extern std::vector<double> zbuffer;

// ---------------------------------------------------------------------------
// Phong Shader
struct PhongShader : IShader {
    const Model& model;

    vec4 light_dir_eye;
    vec2 varying_uv[3];
    vec3 varying_pos[3];
    vec3 varying_nrm[3];

    PhongShader(const vec3 light, const Model& m) : model(m) {
        vec4 L = ModelView * vec4{ light[0], light[1], light[2], 0. };

        vec3 l3 = normalized(L.xyz());
        light_dir_eye = vec4{ l3[0], l3[1], l3[2], 0.0 };
    }


    virtual vec4 vertex(const int face, const int vert) {
        varying_uv[vert] = model.uv(face, vert);
        vec3 v = model.vert(face, vert);
        varying_nrm[vert] = model.normal(face, vert);

        vec4 gl_Position = ModelView * vec4{ v[0], v[1], v[2], 1. };
        varying_pos[vert] = gl_Position.xyz();

        vec4 clip = Perspective * gl_Position;

        // ---- Debug print (PATCH F) ----
        static int dbg_counter = 0;
        if (dbg_counter < 5 * 3) {
            std::cerr << "[vertex] clip = "
                << clip[0] << " "
                << clip[1] << " "
                << clip[2] << " "
                << clip[3] << std::endl;
        }
        dbg_counter++;
        // ---------------------------------

        return clip;
    }

    virtual std::pair<bool, TGAColor> fragment(const vec3 bar) const {
        // Получаем позиции в пространстве камеры
        vec3 p = varying_pos[0] * bar[0] + varying_pos[1] * bar[1] + varying_pos[2] * bar[2];

        // Вычисляем геометрическую нормаль из интерполированных позиций
        // Но нужно иметь доступ к трем вершинам треугольника...
        // Вместо этого используем интерполированные нормали
        vec3 n = normalized(varying_nrm[0] * bar[0] +
            varying_nrm[1] * bar[1] +
            varying_nrm[2] * bar[2]);

        vec3 l = normalized(light_dir_eye.xyz());
        double diff = std::max(0., dot(n, l));

        vec2 uv = varying_uv[0] * bar[0] + varying_uv[1] * bar[1] + varying_uv[2] * bar[2];
        TGAColor color = model.diffuse(uv);

        // Простое освещение как в старом коде
        double ambient = 0.3;
        double spec = 0.0; // Отключим блики для глаз

        for (int ch = 0; ch < 3; ch++) {
            color[ch] = std::min(255.0, color[ch] * (ambient + 0.7 * diff));
        }

        return { false, color };
    }
};

struct EyeShader : IShader {
    const Model& model;
    vec2 varying_uv[3];

    virtual vec4 vertex(const int face, const int vert) {
        varying_uv[vert] = model.uv(face, vert);
        vec3 v = model.vert(face, vert);
        return Perspective * ModelView * vec4{ v[0], v[1], v[2], 1. };
    }

    virtual std::pair<bool, TGAColor> fragment(const vec3 bar) const {
        vec2 uv = varying_uv[0] * bar[0] + varying_uv[1] * bar[1] + varying_uv[2] * bar[2];
        TGAColor color = model.diffuse(uv);
        // Никакого освещения, только текстура
        return { false, color };
    }
};

// ---------------------------------------------------------------------------
// AO utilities
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr int WIDTH = 800;
constexpr int HEIGHT = 800;
constexpr double AO_MAX_DISTANCE = 1000.0;

double max_elevation_angle_from_zbuffer(const std::vector<double>& zb, int px, int py, double dirx, double diry) {
    const int w = WIDTH;
    const int h = HEIGHT;

    double maxangle = 0.0;
    double z0 = zb[px + py * w];

    for (double t = 1.0; t < AO_MAX_DISTANCE; t += 1.0) {
        double cx = px + dirx * t;
        double cy = py + diry * t;

        int ix = int(cx);
        int iy = int(cy);
        if (ix < 0 || ix >= w || iy < 0 || iy >= h) return maxangle;

        double zc = zb[ix + iy * w];

        double dist = std::hypot(cx - px, cy - py);
        if (dist < 1.0) continue;

        double elev = zc - z0;
        double ang = std::atan2(elev, dist);

        if (ang > maxangle) maxangle = ang;
    }
    return maxangle;
}

// ---------------------------------------------------------------------------
// MAIN
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    constexpr vec3 light{ 1, 2, 3 };
    constexpr vec3 eye{ -1, 0.5, 3 };    // Камера справа-сверху-спереди
    constexpr vec3 center{ 0, 0, 0 }; // Смотрим на центр
    constexpr vec3 up{ 0, 1, 0 };

    Camera cam(eye, center, up);
    cam.setPerspective(60.0, (double)WIDTH / HEIGHT, 0.1, 100.0);

    ModelView = cam.view;
    Perspective = cam.proj;

    init_viewport(WIDTH / 16, HEIGHT / 16, WIDTH * 7 / 8, HEIGHT * 7 / 8);
    init_zbuffer(WIDTH, HEIGHT);

    TGAImage framebuffer(WIDTH, HEIGHT, TGAImage::RGB);

    std::cout << "Start rendering (Phong pass)..." << std::endl;

    for (int m = 1; m < argc; m++) {
        Model model(argv[m]);
        std::cerr << "Loaded model: verts=" << model.nverts()
            << " faces=" << model.nfaces() << std::endl;

        PhongShader shader(light, model);

        for (int f = 0; f < model.nfaces(); f++) {
            Triangle tri = {
                shader.vertex(f, 0),
                shader.vertex(f, 1),
                shader.vertex(f, 2)
            };
            rasterize(tri, shader, framebuffer);
        }
    }
    print_render_stats();

    std::cout << "Phong pass done. Computing AO..." << std::endl;

    TGAImage ao(WIDTH, HEIGHT, TGAImage::GRAYSCALE);

    const int samples = 8;
    const double step = M_PI / 4.0;

    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {

            double z = zbuffer[x + y * WIDTH];

            if (z < -1e5) {
                ao.set(x, y, TGAColor(255));
                continue;
            }

            double total = 0.0;

            for (double a = 0; a < 2 * M_PI - 1e-6; a += step) {
                double ang = max_elevation_angle_from_zbuffer(
                    zbuffer, x, y, std::cos(a), std::sin(a)
                );
                total += (M_PI / 2.0 - ang);
            }

            total /= (M_PI / 2.0) * samples;
            total = std::pow(std::max(0.0, std::min(1.0, total)), 100.0);

            unsigned char v = static_cast<unsigned char>(std::min(1.0, total) * 255.0);
            ao.set(x, y, TGAColor(v));
        }
    }

    std::cout << "AO computed. Compositing final image..." << std::endl;

    TGAImage final_img(WIDTH, HEIGHT, TGAImage::RGB);

    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            TGAColor base = framebuffer.get(x, y);
            TGAColor aoc = ao.get(x, y);

            double ao_factor = aoc[0] / 255.0;

            TGAColor out;
            out[2] = std::min(255.0, base[2] * ao_factor);
            out[1] = std::min(255.0, base[1] * ao_factor);
            out[0] = std::min(255.0, base[0] * ao_factor);
            out[3] = 255;

            final_img.set(x, y, out);
        }
    }

    framebuffer.write_tga_file("phong_framebuffer.tga");
    ao.write_tga_file("ao_map.tga");
    final_img.write_tga_file("final.tga");

    std::cout << "Done. Outputs: phong_framebuffer.tga, ao_map.tga, final.tga" << std::endl;
    return 0;
}
