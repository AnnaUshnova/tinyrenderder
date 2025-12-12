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

// Структура вершины (унифицированная)
struct Vertex {
    vec3 position;
    vec3 normal;
    vec2 texcoord;
    vec3 tangent;
    vec3 bitangent;  // Добавляем битангент для полноты
};

// Подмеш с материалом
struct SubMesh {
    std::string name;
    unsigned int startIndex;
    unsigned int indexCount;
    int materialIndex;
    bool hasNormals = false;
    bool hasTexCoords = false;
    bool hasTangents = false;
};

// Текстурный набор
struct MaterialTextures {
    TGAImage diffuse;
    TGAImage normal;
    TGAImage specular;
    TGAImage emission;

    bool hasDiffuse() const { return diffuse.width() > 0; }
    bool hasNormal() const { return normal.width() > 0; }
    bool hasSpecular() const { return specular.width() > 0; }
    bool hasEmission() const { return emission.width() > 0; }
};

class Model {
public:
    Model(const std::string& filename);
    ~Model() = default;

    // Основные методы
    bool load();
    void unload();

    // Доступ к данным
    int getVertexCount() const { return vertices.size(); }
    int getIndexCount() const { return indices.size(); }
    int getSubMeshCount() const { return subMeshes.size(); }
    int getMaterialCount() const { return materials.size(); }

    // Доступ к вершинам (обратная совместимость)
    vec3 vert(int i) const;
    vec3 vert(int iface, int nthvert) const;
    vec3 normal(int iface, int nthvert) const;
    vec2 uv(int iface, int nthvert) const;
    vec3 getCenter() const {
        return localAABB.getCenter();
    }
    vec3 getSize() const {
        return localAABB.max - localAABB.min;
    }

    // Для шейдеров
    TGAColor diffuse(const vec2& uv) const;
    vec3 normal(const vec2& uv) const;
    float specular(const vec2& uv) const;
    TGAColor emission(const vec2& uv) const;

    // Новые методы
    const SubMesh& getSubMesh(int index) const { return subMeshes[index]; }
    const MaterialTextures& getMaterial(int index) const { return materials[index]; }
    const std::vector<Vertex>& getVertices() const { return vertices; }
    const std::vector<unsigned int>& getIndices() const { return indices; }

    // Старая совместимость (для существующего кода)
    int nverts() const { return getVertexCount(); }
    int nfaces() const { return getIndexCount() / 3; }
    bool hasNormalMap() const { return materials.size() > 0 && materials[0].hasNormal(); }

    // AABB модели в локальных координатах
    const AABB& getLocalAABB() const { return localAABB; }

    // Вычисление AABB в мировых координатах
    AABB getWorldAABB(const mat<4, 4>& modelMatrix) const {
        return localAABB.transform(modelMatrix);
    }

private:
    // Загрузка через Assimp
    bool loadModel(const std::string& path);
    void processNode(aiNode* node);
    void processMesh(aiMesh* mesh);

    // Загрузка текстур
    bool loadTexturesFromMaterial(aiMaterial* mat, MaterialTextures& textures);
    bool loadTexture(aiMaterial* mat, aiTextureType type, TGAImage& image,
        const std::string& defaultName = "");

    // Утилиты
    void computeTangentsIfNeeded();
    void generateNormalsIfNeeded();

    // Основные данные
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    std::vector<SubMesh> subMeshes;
    std::vector<MaterialTextures> materials;

    // Assimp данные
    Assimp::Importer importer;
    const aiScene* scene = nullptr;
    std::string directory;      // Директория модели
    std::string filename;       // Имя файла модели
    bool isLoaded = false;

    // Дружественный класс ModelManager
    friend class ModelManager;

    AABB localAABB;  // AABB в локальных координатах модели
    void computeAABB();  // Вычисление AABB при загрузке
};