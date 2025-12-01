#pragma once
#include "geometry.h"
#include "tgaimage.h"
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// Globals
extern mat<4, 4> ModelView;
extern mat<4, 4> Perspective;
extern mat<4, 4> Viewport;
extern std::vector<double> zbuffer;

// matrix helpers used by main
void lookat(const vec3 eye, const vec3 center, const vec3 up);
void init_perspective(double fov_deg, double aspect, double znear, double zfar);
void init_viewport(int x, int y, int w, int h);
void init_zbuffer(int width, int height);

// rasterizer / shader interface (tiny renderer style)
struct IShader {
    // helper: sample 2D texture (uv in [0,1])
    static TGAColor sample2D(const TGAImage& img, const vec2& uv) {
        int x = std::min<int>(img.width() - 1, std::max<int>(0, int(uv.x * img.width())));
        int y = std::min<int>(img.height() - 1, std::max<int>(0, int(uv.y * img.height())));
        return img.get(x, y);
    }

    virtual vec4 vertex(int iface, int nth) { return vec4(); }
    virtual std::pair<bool, TGAColor> fragment(const vec3 bar) const = 0;
};

typedef vec<4> Triangle[3];

void rasterize(const Triangle& clip, const IShader& shader, TGAImage& framebuffer);

// debug
void print_render_stats();
