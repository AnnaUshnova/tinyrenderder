// main.cpp
// Phong rendering + normal map (with eye-protection) + SSAO-like AO
// Outputs: phong.tga, zbuffer.tga, ao.tga, final.tga

#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>

#include "our_gl.h"
#include "model.h"
#include "tgaimage.h"

// externals from our_gl
extern mat<4, 4> ModelView;
extern mat<4, 4> Perspective;
extern mat<4, 4> Viewport;
extern std::vector<double> zbuffer;

// ------------------------------------------------------------------
// Settings
// ------------------------------------------------------------------
const int WIDTH = 800;
const int HEIGHT = 800;
const char* OBJ_PATH = "obj/african_head/african_head.obj";
const char* EYES_OBJ_PATH = "obj/african_head/african_head_eye_inner.obj";

// AO params
constexpr int AO_NUM_DIRECTIONS = 8;
constexpr int AO_STEPS_PER_DIR = 8;
constexpr double AO_STEP_RADIUS = 16.0; // pixels
constexpr double AO_OCCLUSION_THRESHOLD = 1e-3;

// Light (world space)
const vec3 LIGHT_WORLD = normalized(vec3{ 1.0, 1.0, 1.0 });

// Eye detection thresholds (to avoid normalmap on eyes)
constexpr double EYE_DIFFUSE_BRIGHTNESS = 0.85; // normalized [0..1]
constexpr double EYE_SPECULAR_THRESHOLD = 5.0;  // if specular < this -> likely eye (no shiny highlight needed)

// ------------------------------------------------------------------
// Phong shader with normalmap blending + eye protection
// (unchanged)
// ------------------------------------------------------------------
struct PhongShader : public IShader {
    const Model* model;
    vec4 light_dir_eye;              // light direction in eye-space
    vec2 varying_uv[3];
    vec3 varying_pos_eye[3];         // position in eye-space before perspective
    vec3 varying_nrm_eye[3];         // geometric normal transformed to eye-space

    // normalmap strength (0..1)
    double normalmap_strength = 1.0; // use full normalmap by default, blend below for eyes

    PhongShader(const Model* m, const vec3& light_world) : model(m) {
        // transform light dir into eye-space (ModelView must be set)
        vec4 Lw = vec4{ light_world[0], light_world[1], light_world[2], 0.0 };
        vec4 Le = ModelView * Lw;
        vec3 le3 = normalized(Le.xyz());
        light_dir_eye = vec4{ le3[0], le3[1], le3[2], 0.0 };
    }

    virtual vec4 vertex(int iface, int nth) {
        vec3 v = model->vert(iface, nth);
        vec3 n = model->normal(iface, nth);   // geometric normal (from OBJ)

        vec2 uv = model->uv(iface, nth);
        varying_uv[nth] = uv;

        // position in eye-space
        vec4 pos_eye4 = ModelView * vec4{ v[0], v[1], v[2], 1.0 };
        varying_pos_eye[nth] = pos_eye4.xyz();

        // transform geometric normal to eye-space (w=0)
        vec4 n_eye4 = ModelView * vec4{ n[0], n[1], n[2], 0.0 };
        varying_nrm_eye[nth] = n_eye4.xyz();

        // return clip-space position
        vec4 clip = Perspective * pos_eye4;
        return clip;
    }

    // pcbar = perspective-correct barycentric
    virtual std::pair<bool, TGAColor> fragment(const vec3 pcbar) const {
        // interpolate attributes
        vec3 p_eye = varying_pos_eye[0] * pcbar[0] + varying_pos_eye[1] * pcbar[1] + varying_pos_eye[2] * pcbar[2];
        vec3 n_eye_geom = varying_nrm_eye[0] * pcbar[0] + varying_nrm_eye[1] * pcbar[1] + varying_nrm_eye[2] * pcbar[2];
        vec2 uv = varying_uv[0] * pcbar[0] + varying_uv[1] * pcbar[1] + varying_uv[2] * pcbar[2];

        // base diffuse color & specular power from model
        TGAColor base = model->diffuse(uv);
        double spec_power = std::max(1.0, (double)model->specular(uv));

        // decide if this pixel is an eye region:
        double brightness = ((double)base[0] + (double)base[1] + (double)base[2]) / (3.0 * 255.0);
        bool likely_eye = (brightness >= EYE_DIFFUSE_BRIGHTNESS) && (spec_power <= EYE_SPECULAR_THRESHOLD);

        // sample normal map (if present)
        vec3 nm = model->normal(uv); // note: model.normal(uv) returns (x,y,z) mapped from map
        // Transform sampled normal to eye-space too (treating as vector in object space)
        vec4 nm_eye4 = ModelView * vec4{ nm[0], nm[1], nm[2], 0.0 };
        vec3 nm_eye = nm_eye4.xyz();

        // Blend geometric normal and normalmap — but protect eyes
        vec3 n_eye;
        if (likely_eye) {
            // use geometric normal only on eyes
            n_eye = n_eye_geom;
        }
        else {
            // blend: more weight to normalmap but keep some geometric base to avoid bad tangents
            double nm_strength = normalmap_strength;
            // If normalmap is absent, model.normal returns (0,0,1) — blending will be harmless
            n_eye = normalized(n_eye_geom * (1.0 - nm_strength) + nm_eye * nm_strength);
        }

        // normalize final normal
        n_eye = normalized(n_eye);

        // lighting (eye-space)
        vec3 l = normalized(light_dir_eye.xyz());
        vec3 v = normalized(vec3{ -p_eye.x, -p_eye.y, -p_eye.z }); // view vector (camera at origin)
        vec3 r = normalized((2.0 * dot(n_eye, l)) * n_eye - l);

        double diff = std::max(0.0, dot(n_eye, l));
        double spec = 0.1;
        double rv = std::max(0.0, dot(r, v));
        if (rv > 0.0) spec = std::pow(rv, spec_power);

        // Compose final color
        const double ambient = 0.1;
        const double spec_strength = 0.3;

        TGAColor out = base;
        for (int c = 0; c < 3; ++c) {
            double col = base[c];
            double shaded = col * (ambient + diff * (1.0 - ambient)) + 255.0 * spec_strength * spec;
            if (shaded > 255.0) shaded = 255.0;
            out[c] = (unsigned char)(shaded);
        }

        return { false, out };
    }
};

// ------------------------------------------------------------------
// Simple Eye shader (separate, glossy but simple)
// - uses same light direction so skin lighting remains consistent
// - simpler fragment (no normalmap blending), stronger specular
// ------------------------------------------------------------------
struct EyeShader : public IShader {
    const Model* model;
    vec4 light_dir_eye;
    vec2 varying_uv[3];
    vec3 varying_pos_eye[3];
    vec3 varying_nrm_eye[3];

    EyeShader(const Model* m, const vec3& light_world) : model(m) {
        vec4 Lw = vec4{ light_world[0], light_world[1], light_world[2], 0.0 };
        vec4 Le = ModelView * Lw;
        vec3 le3 = normalized(Le.xyz());
        light_dir_eye = vec4{ le3[0], le3[1], le3[2], 0.0 };
    }

    virtual vec4 vertex(int iface, int nth) {
        vec3 v = model->vert(iface, nth);
        vec3 n = model->normal(iface, nth);

        varying_uv[nth] = model->uv(iface, nth);

        vec4 pos_eye4 = ModelView * vec4{ v[0], v[1], v[2], 1.0 };
        varying_pos_eye[nth] = pos_eye4.xyz();

        vec4 n_eye4 = ModelView * vec4{ n[0], n[1], n[2], 0.0 };
        varying_nrm_eye[nth] = n_eye4.xyz();

        return Perspective * pos_eye4;
    }

    virtual std::pair<bool, TGAColor> fragment(const vec3 pcbar) const {
        vec3 p_eye = varying_pos_eye[0] * pcbar[0] + varying_pos_eye[1] * pcbar[1] + varying_pos_eye[2] * pcbar[2];
        vec3 n_eye = varying_nrm_eye[0] * pcbar[0] + varying_nrm_eye[1] * pcbar[1] + varying_nrm_eye[2] * pcbar[2];
        vec2 uv = varying_uv[0] * pcbar[0] + varying_uv[1] * pcbar[1] + varying_uv[2] * pcbar[2];

        // base color from eye texture (if exists)
        TGAColor base = model->diffuse(uv);

        // make eyes glossier: amplify specular power
        double spec_map = std::max(1.0, (double)model->specular(uv));
        double spec_power = spec_map * 8.0; // stronger highlight for eyes

        vec3 l = normalized(light_dir_eye.xyz());
        vec3 v = normalized(vec3{ -p_eye.x, -p_eye.y, -p_eye.z });
        vec3 rn = normalized(n_eye);

        double diff = std::max(0.0, dot(rn, l));

        vec3 r = normalized((2.0 * dot(rn, l)) * rn - l);
        double rv = std::max(0.0, dot(r, v));
        double spec = 0.0;
        if (rv > 0.0) spec = std::pow(rv, spec_power);

        const double ambient = 0.1;
        const double spec_strength = 1.5;

        TGAColor out = base;
        for (int c = 0; c < 3; ++c) {
            double col = base[c];
            double shaded = col * (ambient + diff * (1.0 - ambient)) + 255.0 * spec_strength * spec;
            if (shaded > 255.0) shaded = 255.0;
            out[c] = (unsigned char)(shaded);
        }

        return { false, out };
    }
};

// ------------------------------------------------------------------
// Save zbuffer to grayscale image (no flip)
// ------------------------------------------------------------------
static void save_zbuffer_image(const std::vector<double>& zb, int w, int h, const char* filename) {
    TGAImage img(w, h, TGAImage::RGB);

    // find finite min/max
    double zmin = std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < w * h; ++i) {
        double z = zb[i];
        if (!std::isfinite(z)) continue;
        zmin = std::min(zmin, z);
        zmax = std::max(zmax, z);
    }
    if (zmin == std::numeric_limits<double>::infinity()) {
        // empty -> white
        for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) img.set(x, y, TGAColor(255, 255, 255));
        img.write_tga_file(filename);
        return;
    }
    if (zmax - zmin < 1e-7) zmax = zmin + 1e-7;

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double z = zb[x + y * w];
            if (!std::isfinite(z)) {
                img.set(x, y, TGAColor(255, 255, 255));
            }
            else {
                double t = (z - zmin) / (zmax - zmin);
                unsigned char v = (unsigned char)(std::round(255.0 * (1.0 - t)));
                img.set(x, y, TGAColor(v, v, v));
            }
        }
    }
    img.write_tga_file(filename);
}

// ------------------------------------------------------------------
// Simple SSAO-ish screen-space occlusion (no flip)
// ------------------------------------------------------------------
static double compute_ssao_at(const std::vector<double>& zb, int w, int h, int px, int py) {
    int idx0 = px + py * w;
    double z0 = zb[idx0];
    if (!std::isfinite(z0)) return 1.0;

    int occluded = 0;
    int total = 0;

    for (int d = 0; d < AO_NUM_DIRECTIONS; ++d) {
        double ang = (2.0 * M_PI * d) / AO_NUM_DIRECTIONS;
        double dx = std::cos(ang);
        double dy = std::sin(ang);
        for (int s = 1; s <= AO_STEPS_PER_DIR; ++s) {
            double t = (double)s / (double)AO_STEPS_PER_DIR;
            double radius = t * AO_STEP_RADIUS;
            int sx = (int)std::round(px + dx * radius);
            int sy = (int)std::round(py + dy * radius);
            if (sx < 0 || sx >= w || sy < 0 || sy >= h) continue;
            int sidx = sx + sy * w;
            double zs = zb[sidx];
            if (!std::isfinite(zs)) {
                total++;
                continue;
            }
            if (zs < z0 - AO_OCCLUSION_THRESHOLD) occluded++;
            total++;
        }
    }

    if (total == 0) return 1.0;
    double occ = (double)occluded / (double)total;
    return 1.0 - occ;
}

// ------------------------------------------------------------------
// Main
// ------------------------------------------------------------------
int main(int argc, char** argv) {
    const char* model_path = (argc > 1) ? argv[1] : OBJ_PATH;
    const char* eyes_path = EYES_OBJ_PATH;
    Model model(model_path);

    // load eyes model located in the same folder as the head (as requested)
    Model eye_model(eyes_path);

    TGAImage framebuffer(WIDTH, HEIGHT, TGAImage::RGB);

    // init transforms + zbuffer (our_gl uses +inf sentinel)
    init_zbuffer(WIDTH, HEIGHT);
    vec3 eye = vec3{ 0.7, 1.0, 3.0 };
    vec3 center = vec3{ 0.0, 0.0, 0.0 };
    vec3 up = vec3{ 0.0, 1.0, 0.0 };
    lookat(eye, center, up);
    init_perspective(60.0, (double)WIDTH / (double)HEIGHT, 0.1, 100.0);
    init_viewport(0, 0, WIDTH, HEIGHT);

    PhongShader shader(&model, LIGHT_WORLD);

    // render main model (updates zbuffer)
    for (int f = 0; f < model.nfaces(); ++f) {
        vec4 clip[3];
        for (int v = 0; v < 3; ++v) clip[v] = shader.vertex(f, v);
        rasterize(clip, shader, framebuffer);
    }

    // Backup zbuffer before drawing eyes so AO/zbuffer remain unchanged by eyes.
    std::vector<double> zbuffer_backup = zbuffer;

    // Render eyes into the same framebuffer using a separate, simple shader.
    EyeShader eye_shader(&eye_model, LIGHT_WORLD);
    for (int f = 0; f < eye_model.nfaces(); ++f) {
        vec4 clip[3];
        for (int v = 0; v < 3; ++v) clip[v] = eye_shader.vertex(f, v);
        rasterize(clip, eye_shader, framebuffer);
    }
    std::cerr << "Info: eyes rendered into framebuffer (phong) using EyeShader.\n";

    // Save PHONG image (with eyes)
    TGAImage phong_img = framebuffer;
    phong_img.write_tga_file("phong.tga");

    // restore zbuffer so AO and zbuffer image do not include the eyes
    zbuffer = zbuffer_backup;

    // ZBUFFER (no flip) — now without eyes
    save_zbuffer_image(zbuffer, WIDTH, HEIGHT, "zbuffer.tga");

    // AO map (no flip)
    TGAImage ao_img(WIDTH, HEIGHT, TGAImage::RGB);
    for (int y = 0; y < HEIGHT; ++y) {
        for (int x = 0; x < WIDTH; ++x) {
            double ao = compute_ssao_at(zbuffer, WIDTH, HEIGHT, x, y);
            unsigned char v = (unsigned char)std::round(255.0 * ao);
            ao_img.set(x, y, TGAColor(v, v, v));
        }
    }
    ao_img.write_tga_file("ao.tga");

    // Final compose: phong * ao (final IS flipped to match previous output style)
    TGAImage final_img(WIDTH, HEIGHT, TGAImage::RGB);
    for (int y = 0; y < HEIGHT; ++y) {
        for (int x = 0; x < WIDTH; ++x) {
            TGAColor p = phong_img.get(x, y);
            TGAColor a = ao_img.get(x, y);
            double af = a[0] / 255.0;
            TGAColor out;
            out[0] = (unsigned char)std::min(255.0, std::round(p[0] * af));
            out[1] = (unsigned char)std::min(255.0, std::round(p[1] * af));
            out[2] = (unsigned char)std::min(255.0, std::round(p[2] * af));
            final_img.set(x, y, out);
        }
    }
    final_img.write_tga_file("final.tga");

    print_render_stats();
    std::cerr << "Saved: phong.tga (with eyes), zbuffer.tga (no eyes), ao.tga, final.tga (with eyes)\n";
    return 0;
}
