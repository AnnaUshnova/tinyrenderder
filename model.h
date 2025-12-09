#pragma once
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "tgaimage.h"

// Включаем заголовки Assimp
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <assimp/Importer.hpp>

class Model {
public:
    Model(const std::string filename);
    ~Model() = default;

    int nverts() const { return verts.size(); }
    int nfaces() const { return facet_vrt.size() / 3; }

    // geometry
    vec3 vert(int i) const;
    vec3 vert(int iface, int nthvert) const;

    // normals
    vec3 normal(int iface, int nthvert) const;
    vec3 normal(const vec2& uv) const;
    bool hasNormalMap() const { return normalmap.width() > 0; }

    // uv coordinates
    vec2 uv(int iface, int nthvert) const;

    // textures
    TGAColor diffuse(const vec2& uv) const;
    float specular(const vec2& uv) const;

private:
    // Загрузка через Assimp
    void loadModel(const std::string& path);
    void processNode(aiNode* node);
    void processMesh(aiMesh* mesh, const aiScene* scene);

    // Векторные данные
    std::vector<vec3> verts;
    std::vector<vec3> norms;
    std::vector<vec2> tex;
    std::vector<vec3> tangents;

    // Индексы
    std::vector<int> facet_vrt;
    std::vector<int> facet_tex;
    std::vector<int> facet_nrm;

    // Текстуры
    TGAImage diffusemap;
    TGAImage normalmap;
    TGAImage specmap;

    // Импортер Assimp
    Assimp::Importer importer;
    const aiScene* scene = nullptr;
    std::string baseDirectory; // базовая директория для текстур
};