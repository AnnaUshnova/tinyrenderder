#include <fstream>
#include <sstream>
#include "model.h"

Model::Model(const std::string filename) {
    std::ifstream in(filename);
    if (in.fail()) {
        std::cerr << "Failed to open OBJ: " << filename << std::endl;
        return;
    }

    // ensure tangents vector size
    tangents.assign(verts.size(), vec3{ 0,0,0 });

    // accumulate face tangents
    for (int f = 0; f < nfaces(); ++f) {
        int i0 = facet_vrt[f * 3 + 0];
        int i1 = facet_vrt[f * 3 + 1];
        int i2 = facet_vrt[f * 3 + 2];

        vec3 v0 = verts[i0];
        vec3 v1 = verts[i1];
        vec3 v2 = verts[i2];

        vec2 uv0 = tex[facet_tex[f * 3 + 0]];
        vec2 uv1 = tex[facet_tex[f * 3 + 1]];
        vec2 uv2 = tex[facet_tex[f * 3 + 2]];

        vec3 deltaPos1 = v1 - v0;
        vec3 deltaPos2 = v2 - v0;

        vec2 deltaUV1 = uv1 - uv0;
        vec2 deltaUV2 = uv2 - uv0;

        double r = (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);
        if (std::abs(r) < 1e-8) continue;
        double invr = 1.0 / r;

        vec3 tangent = (deltaPos1 * deltaUV2.y - deltaPos2 * deltaUV1.y) * invr;

        tangents[i0] = tangents[i0] + tangent;
        tangents[i1] = tangents[i1] + tangent;
        tangents[i2] = tangents[i2] + tangent;
    }

    // normalize tangents
    for (size_t i = 0; i < tangents.size(); ++i) {
        tangents[i] = normalized(tangents[i]);
    }

    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        char trash;

        // ----------------------------------------------------------------------
        // vertices
        // ----------------------------------------------------------------------
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            vec3 v;
            iss >> v.x >> v.y >> v.z;
            verts.push_back(v);
        }

        // ----------------------------------------------------------------------
        // vertex normals
        // ----------------------------------------------------------------------
        else if (!line.compare(0, 3, "vn ")) {
            iss >> trash >> trash;
            vec3 n;
            iss >> n.x >> n.y >> n.z;
            norms.push_back(normalized(n));
        }

        // ----------------------------------------------------------------------
        // texture coordinates
        // ----------------------------------------------------------------------
        else if (!line.compare(0, 3, "vt ")) {
            iss >> trash >> trash;
            vec2 uv;
            iss >> uv.x >> uv.y;
            tex.push_back({ uv.x, 1 - uv.y });   // flip V
        }

        // ----------------------------------------------------------------------
        // faces: f v/t/n v/t/n v/t/n
        // ----------------------------------------------------------------------
        else if (!line.compare(0, 2, "f ")) {
            iss >> trash;
            int v, t, n;
            int cnt = 0;
            while (iss >> v >> trash >> t >> trash >> n) {
                facet_vrt.push_back(v - 1);
                facet_tex.push_back(t - 1);
                facet_nrm.push_back(n - 1);
                cnt++;
            }
            if (cnt != 3) {
                std::cerr << "Error: OBJ file must be triangulated\n";
                return;
            }
        }
    }

    std::cerr << "# v# " << nverts() << " f# " << nfaces() << std::endl;

    // compute per-vertex tangents (after parsing!)
    tangents.assign(verts.size(), vec3{ 0,0,0 });

    for (int f = 0; f < nfaces(); f++) {
        int i0 = facet_vrt[3 * f + 0];
        int i1 = facet_vrt[3 * f + 1];
        int i2 = facet_vrt[3 * f + 2];

        vec3 v0 = verts[i0];
        vec3 v1 = verts[i1];
        vec3 v2 = verts[i2];

        vec2 uv0 = tex[facet_tex[3 * f + 0]];
        vec2 uv1 = tex[facet_tex[3 * f + 1]];
        vec2 uv2 = tex[facet_tex[3 * f + 2]];

        vec3 e1 = v1 - v0;
        vec3 e2 = v2 - v0;
        vec2 duv1 = uv1 - uv0;
        vec2 duv2 = uv2 - uv0;

        double r = (duv1.x * duv2.y - duv1.y * duv2.x);
        if (std::abs(r) < 1e-9) continue;
        r = 1.0 / r;

        vec3 t = (e1 * duv2.y - e2 * duv1.y) * r;

        tangents[i0] = tangents[i0] + t;
        tangents[i1] = tangents[i1] + t;
        tangents[i2] = tangents[i2] + t;
    }

    for (size_t i = 0; i < tangents.size(); i++)
        tangents[i] = normalized(tangents[i]);


    // load textures ------------------------------------------------------------
    auto load_texture = [&](const std::string& suffix, TGAImage& img) {
        size_t dot = filename.find_last_of(".");
        if (dot == std::string::npos) return;
        std::string texfile = filename.substr(0, dot) + suffix;
        bool ok = img.read_tga_file(texfile.c_str());
        std::cerr << "texture file " << texfile << " loading "
            << (ok ? "ok" : "failed") << std::endl;
       // img.flip_vertically();
        };

    load_texture("_diffuse.tga", diffusemap);
    load_texture("_nm.tga", normalmap);
    load_texture("_spec.tga", specmap);
}

// -----------------------------------------------------------------------------
// getters
// -----------------------------------------------------------------------------

vec3 Model::vert(int i) const {
    return verts[i];
}

vec3 Model::vert(int iface, int nthvert) const {
    return verts[facet_vrt[iface * 3 + nthvert]];
}

vec3 Model::normal(int iface, int nthvert) const {
    return norms[facet_nrm[iface * 3 + nthvert]];
}

vec3 Model::normal(const vec2& uv) const {
    TGAColor c = normalmap.get(
        uv.x * normalmap.width(),
        uv.y * normalmap.height()
    );

    vec3 n;
    n.x = (double)c[2] / 255.0 * 2.0 - 1.0;
    n.y = (double)c[1] / 255.0 * 2.0 - 1.0;
    n.z = (double)c[0] / 255.0 * 2.0 - 1.0;

    return normalized(n);
}

vec2 Model::uv(int iface, int nthvert) const {
    return tex[facet_tex[iface * 3 + nthvert]];
}

// diffuse lookup --------------------------------------------------------------

TGAColor Model::diffuse(const vec2& uv) const {
    return diffusemap.get(
        uv.x * diffusemap.width(),
        uv.y * diffusemap.height()
    );
}

// specular coefficient --------------------------------------------------------

float Model::specular(const vec2& uv) const {
    TGAColor c = specmap.get(
        uv.x * specmap.width(),
        uv.y * specmap.height()
    );
    return c[0] / 1.0f;   // takes R-channel as shininess
}


