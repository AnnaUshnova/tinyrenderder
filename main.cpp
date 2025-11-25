#include "our_gl.h"
#include "model.h"

extern mat<4,4> ModelView, Perspective; // "OpenGL" state matrices and
extern std::vector<double> zbuffer;     // the depth buffer

struct PhongShader : IShader {
    const Model& model;

    vec4 light_dir_eye;     // направление света в eye-space
    vec2 varying_uv[3];     // UV на вершинах
    vec3 varying_pos[3];    // позиция в eye-space (для вычисления view direction)
    vec3 varying_nrm[3];    // геометрическая нормаль (для построения TBN)

    PhongShader(const vec3 light, const Model& m) : model(m) {
        light_dir_eye = normalized(ModelView * vec4{ light.x, light.y, light.z, 0. });
    }

    // --- VERTEX SHADER -------------------------------------------------------
    virtual vec4 vertex(const int face, const int vert) {

        varying_uv[vert] = model.uv(face, vert);

        // исходная позиция в мире
        vec3 v = model.vert(face, vert);

        // нормаль вершины (из .obj)
        varying_nrm[vert] = model.normal(face, vert);

        // перевод позиции в eye-space
        vec4 gl_Position = ModelView * vec4{ v.x, v.y, v.z, 1. };
        varying_pos[vert] = gl_Position.xyz();

        return Perspective * gl_Position;
    }

    // --- FRAGMENT SHADER -----------------------------------------------------
    virtual std::pair<bool, TGAColor> fragment(const vec3 bar) const {

        vec2 uv = varying_uv[0] * bar[0] +
            varying_uv[1] * bar[1] +
            varying_uv[2] * bar[2];

        // diffuse
        TGAColor color = model.diffuse(uv);

        // normal from normal map (tinyrenderer format = already in view space)
        vec3 n = normalized(model.normal(uv));

        vec3 l = normalized(light_dir_eye.xyz());

        // phong components
        double ambient = 0.3;
        double diff = std::max(0., n * l);

        // reflection vector
        vec3 r = normalized(n * (2.0 * (n * l)) - l);

        // SPECULAR FROM SPEC MAP (important!!)
        float spec_power = model.specular(uv);
        if (spec_power < 1) spec_power = 1;

        // *** key trick: mix old r.z hack with correct r*l ***
        double spec_old = std::pow(std::max(r.z, 0.), 35.0);      // bright punchy highlight
        double spec_real = std::pow(std::max(0., r * l), spec_power);  // correct phong

        // blend them: keeps realism + gives nice bright highlight
        double spec = spec_real * 0.5 + spec_old * 0.7;

        // final color
        for (int ch : {0, 1, 2})
            color[ch] = std::min(255.,
                color[ch] * (ambient + diff * 0.4 + spec * 0.9));

        return { false, color };
    }


};



int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    constexpr int width  = 800;      // output image size
    constexpr int height = 800;
    constexpr vec3  light{ -1, 1, 3}; // light source
    constexpr vec3    eye{-1, 0, 2}; // camera position
    constexpr vec3 center{ 0, 0, 0}; // camera direction
    constexpr vec3     up{ 0, 1, 0}; // camera up vector

    lookat(eye, center, up);                                   // build the ModelView   matrix
    init_perspective(norm(eye-center));                        // build the Perspective matrix
    init_viewport(width/16, height/16, width*7/8, height*7/8); // build the Viewport    matrix
    init_zbuffer(width, height);
    TGAImage framebuffer(width, height, TGAImage::RGB);

    std::cout << "Start rendering..." << std::endl;


    for (int m=1; m<argc; m++) {                    // iterate through all input objects
        Model model(argv[m]);                       // load the data
        PhongShader shader(light, model);
        for (int f=0; f<model.nfaces(); f++) {      // iterate through all facets
            Triangle clip = { shader.vertex(f, 0),  // assemble the primitive
                              shader.vertex(f, 1),
                              shader.vertex(f, 2) };
            rasterize(clip, shader, framebuffer);   // rasterize the primitive
        }
    }

    std::cout << "Done rendering!" << std::endl;


    framebuffer.write_tga_file("framebuffer.tga");
    return 0;
}

