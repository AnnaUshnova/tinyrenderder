#pragma once
#include <vector>
#include <string>
#include "geometry.h"
#include "tgaimage.h"

class Model {
public:
    Model(const std::string filename);

    int nverts() const { return verts.size(); }
    int nfaces() const { return facet_vrt.size() / 3; }

    // geometry
    vec3 vert(int i) const;
    vec3 vert(int iface, int nthvert) const;

    // normals
    vec3 normal(int iface, int nthvert) const;
    vec3 normal(const vec2& uv) const;

    // uv coordinates
    vec2 uv(int iface, int nthvert) const;

    // textures
    TGAColor diffuse(const vec2& uv) const;
    float specular(const vec2& uv) const;

private:
    // vertex attributes
    std::vector<vec3> verts;
    std::vector<vec3> norms;
    std::vector<vec2> tex;

    // tangents per-vertex (object space)
    std::vector<vec3> tangents;

    // face indices
    std::vector<int> facet_vrt;
    std::vector<int> facet_tex;
    std::vector<int> facet_nrm;

    // texture maps (must exist!)
    TGAImage diffusemap;
    TGAImage normalmap;
    TGAImage specmap;
};
