// main.cpp — Phong render + SSAO (AO) post-process
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

#include "our_gl.h"
#include "model.h"
#include "tgaimage.h"
#include "camera.h"

// externs from GL helpers
extern mat<4, 4> ModelView, Perspective;
extern std::vector<double> zbuffer;

// ============================================================================
// LIGHTING PARAMETERS
// ============================================================================
struct PhongLight {
    vec3 light_dir_world;

    double ambient;
    double diffuse_k;
    double specular_k;
    double shininess;

    double roughness;   // 0 = гладко, 1 = матово
};

// ============================================================================
// Phong Shader
// ============================================================================
struct PhongShader : IShader {
    const Model& model;
    const PhongLight& light;

    vec4 light_dir_eye;
    vec2 varying_uv[3];
    vec3 varying_pos[3];
    vec3 varying_nrm[3];

    PhongShader(const PhongLight& lp, const Model& m)
        : model(m), light(lp)
    {
        // convert light direction into eye-space using *current* ModelView
        vec4 L = ModelView * vec4{ lp.light_dir_world[0], lp.light_dir_world[1], lp.light_dir_world[2], 0.0 };
        vec3 l3 = normalized(L.xyz());
        light_dir_eye = vec4{ l3[0], l3[1], l3[2], 0.0 };
    }

    virtual vec4 vertex(const int face, const int vert) {
        varying_uv[vert] = model.uv(face, vert);
        vec3 v = model.vert(face, vert);
        varying_nrm[vert] = model.normal(face, vert);

        vec4 gl_Position = ModelView * vec4{ v[0], v[1], v[2], 1.0 };
        varying_pos[vert] = gl_Position.xyz();

        return Perspective * gl_Position;
    }

    virtual std::pair<bool, TGAColor> fragment(const vec3 bar) const {
        vec3 p =
            varying_pos[0] * bar[0] +
            varying_pos[1] * bar[1] +
            varying_pos[2] * bar[2];

        vec3 n = normalized(
            varying_nrm[0] * bar[0] +
            varying_nrm[1] * bar[1] +
            varying_nrm[2] * bar[2]
        );

        vec2 uv =
            varying_uv[0] * bar[0] +
            varying_uv[1] * bar[1] +
            varying_uv[2] * bar[2];

        TGAColor color = model.diffuse(uv);

        // Light vectors
        vec3 l = normalized(light_dir_eye.xyz());
        vec3 v = normalized(p * -1.0);
        vec3 r = normalized(2 * dot(n, l) * n - l);

        double diff = std::max(0.0, dot(n, l));

        // Roughness-modified shininess
        double smooth_shininess = light.shininess * (1.0 - light.roughness);
        smooth_shininess = std::max(smooth_shininess, 1.0);

        double spec = 0.0;
        double rv = std::max(dot(r, v), 0.0);
        if (rv > 0.0)
            spec = std::pow(rv, smooth_shininess);

        // Final shading
        for (int c = 0; c < 3; c++) {
            double base = color[c];
            double shaded =
                base * (light.ambient + diff * light.diffuse_k) +
                255.0 * spec * light.specular_k;

            color[c] = std::min(255.0, shaded);
        }

        return { false, color };
    }
};

// ============================================================================
// AO utilities
// ============================================================================
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr int WIDTH = 800;
constexpr int HEIGHT = 800;
constexpr double AO_MAX_DISTANCE = 1000.0;

double max_elevation_angle_from_zbuffer(const std::vector<double>& zb, int px, int py, double dirx, double diry) {
    double maxangle = 0.0;
    double z0 = zb[px + py * WIDTH];

    for (double t = 1.0; t < AO_MAX_DISTANCE; t += 1.0) {
        double cx = px + dirx * t;
        double cy = py + diry * t;

        int ix = int(cx);
        int iy = int(cy);
        if (ix < 0 || ix >= WIDTH || iy < 0 || iy >= HEIGHT) return maxangle;

        double zc = zb[ix + iy * WIDTH];

        double dist = std::hypot(cx - px, cy - py);
        if (dist < 1.0) continue;

        double elev = zc - z0;
        double ang = std::atan2(elev, dist);

        if (ang > maxangle) maxangle = ang;
    }
    return maxangle;
}

// ============================================================================
// MAIN
// ============================================================================
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj\n";
        return 1;
    }

    // ------------------------------------------------------------------------
    // 1) LIGHT PARAMETERS (YOUR EDITABLE BLOCK!)
    // ------------------------------------------------------------------------
    PhongLight light_cfg;
    light_cfg.light_dir_world = { 1, 2, 3 };
    light_cfg.ambient = 0.25;
    light_cfg.diffuse_k = 1.0;
    light_cfg.specular_k = 0.2;
    light_cfg.shininess = 64.0;
    light_cfg.roughness = 0.8;    // ← мотируем поверхность

    // ------------------------------------------------------------------------
    // 2) CAMERA SETUP
    // ------------------------------------------------------------------------
    constexpr vec3 eye{ -1, 0.5, 3 };
    constexpr vec3 center{ 0, 0, 0 };
    constexpr vec3 up{ 0, 1, 0 };

    Camera cam(eye, center, up);
    cam.setPerspective(60.0, (double)WIDTH / HEIGHT, 0.1, 100.0);

    ModelView = cam.view;
    Perspective = cam.proj;

    // ------------------------------------------------------------------------
    // 3) RENDER SETUP
    // ------------------------------------------------------------------------
    init_viewport(WIDTH / 16, HEIGHT / 16, WIDTH * 7 / 8, HEIGHT * 7 / 8);
    init_zbuffer(WIDTH, HEIGHT);

    TGAImage framebuffer(WIDTH, HEIGHT, TGAImage::RGB);

    std::cout << "Start rendering (Phong pass)...\n";

    // ------------------------------------------------------------------------
    // 4) PHONG PASS
    // ------------------------------------------------------------------------
    for (int m = 1; m < argc; m++) {
        Model model(argv[m]);
        std::cerr << "Loaded model: verts=" << model.nverts()
            << " faces=" << model.nfaces() << "\n";

        // NOTE: shader created *after* ModelView was set → OK now
        PhongShader shader(light_cfg, model);

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

    // ------------------------------------------------------------------------
    // 5) SSAO PASS (unchanged)
    // ------------------------------------------------------------------------
    std::cout << "Phong pass done. Computing AO...\n";

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
                total += (M_PI / 2.0 -
                    max_elevation_angle_from_zbuffer(
                        zbuffer, x, y, std::cos(a), std::sin(a)
                    )
                    );
            }

            total /= (M_PI / 2.0) * samples;
            total = std::pow(std::clamp(total, 0.0, 1.0), 100.0);

            ao.set(x, y, TGAColor((unsigned char)(total * 255.0)));
        }
    }

    // ------------------------------------------------------------------------
    // 6) COMPOSITE AO + PHONG
    // ------------------------------------------------------------------------
    std::cout << "AO computed. Compositing final image...\n";

    TGAImage final_img(WIDTH, HEIGHT, TGAImage::RGB);

    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            TGAColor base = framebuffer.get(x, y);
            TGAColor aoc = ao.get(x, y);

            double ao_factor = aoc[0] / 255.0;

            TGAColor out;
            out[2] = base[2] * ao_factor;
            out[1] = base[1] * ao_factor;
            out[0] = base[0] * ao_factor;
            out[3] = 255;

            final_img.set(x, y, out);
        }
    }

    framebuffer.write_tga_file("phong_framebuffer.tga");
    ao.write_tga_file("ao_map.tga");
    final_img.write_tga_file("final.tga");

    std::cout << "Done.\n";
    return 0;
}
