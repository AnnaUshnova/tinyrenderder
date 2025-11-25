#pragma once
#include "tgaimage.h"
#include "geometry.h"
#include <vector>

// ------------------------------------------------------------
// Global transformation matrices
// ------------------------------------------------------------
extern mat<4, 4> ModelView;
extern mat<4, 4> Projection;
extern mat<4, 4> Viewport;
extern mat<4, 4> Perspective;

extern std::vector<double> zbuffer;

// ------------------------------------------------------------
// Matrix helpers
// ------------------------------------------------------------
void lookat(const vec3 eye, const vec3 center, const vec3 up);
void init_perspective(double f);
void init_viewport(int x, int y, int w, int h);
void init_zbuffer(int width, int height);

// ------------------------------------------------------------
// Shader interface
// ------------------------------------------------------------
struct IShader {
    // Used by texture-based shaders
    static TGAColor sample2D(const TGAImage& img, const vec2& uv) {
        return img.get(
            uv.x * img.width(),
            uv.y * img.height()
        );
    }

    // Vertex shader
    virtual vec4 vertex(int face, int nth) { return vec4{ 0,0,0,1 }; }

    // Fragment shader
    virtual std::pair<bool, TGAColor> fragment(const vec3 bar) const = 0;
};

// ------------------------------------------------------------
// Rasterizer
// ------------------------------------------------------------

typedef vec4 Triangle[3];

void rasterize(const Triangle& clip, const IShader& shader, TGAImage& framebuffer);
