#include "model.h"
#include <iostream>
#include <filesystem>
#include <algorithm>
#include "model_manager.h"  // Для дружественного доступа

// Конструктор - только инициализация, реальная загрузка в load()
Model::Model(const std::string& filename)
    : filename(filename), isLoaded(false) {
    std::filesystem::path path(filename);
    directory = path.parent_path().string();
}


void Model::computeAABB() {
    if (vertices.empty()) {
        localAABB = AABB(vec3{ 0,0,0 }, vec3{ 0,0,0 });
        return;
    }

    vec3 min = vec3{ 1e9, 1e9, 1e9 };
    vec3 max = vec3{ -1e9, -1e9, -1e9 };

    for (const auto& vertex : vertices) {
        min.x = std::min(min.x, vertex.position.x);
        min.y = std::min(min.y, vertex.position.y);
        min.z = std::min(min.z, vertex.position.z);

        max.x = std::max(max.x, vertex.position.x);
        max.y = std::max(max.y, vertex.position.y);
        max.z = std::max(max.z, vertex.position.z);
    }

    // Добавляем небольшой запас
    vec3 margin = (max - min) * 0.01;
    localAABB = AABB(min - margin, max + margin);

    std::cout << "Model AABB: min(" << min.x << ", " << min.y << ", " << min.z
        << "), max(" << max.x << ", " << max.y << ", " << max.z << ")" << std::endl;
}


// Основной метод загрузки
bool Model::load() {
    if (isLoaded) {
        std::cerr << "Model already loaded: " << filename << std::endl;
        return true;
    }

    std::cout << "Loading model: " << filename << std::endl;

    // Загружаем модель через Assimp
    if (!loadModel(filename)) {
        std::cerr << "Failed to load model: " << filename << std::endl;
        return false;
    }

    // Генерируем нормали если их нет
    generateNormalsIfNeeded();

    // Вычисляем тангенсы если их нет
    computeTangentsIfNeeded();

    computeAABB();  // Добавить эту строку

    isLoaded = true;
    std::cout << "Model loaded successfully. Vertices: " << vertices.size()
        << ", Indices: " << indices.size()
        << ", SubMeshes: " << subMeshes.size()
        << ", Materials: " << materials.size() << std::endl;

    return true;
}

void Model::unload() {
    vertices.clear();
    indices.clear();
    subMeshes.clear();
    materials.clear();

    // Освобождаем Assimp сцену
    importer.FreeScene();
    scene = nullptr;

    isLoaded = false;
    std::cout << "Model unloaded: " << filename << std::endl;
}

bool Model::loadModel(const std::string& path) {
    // Флаги для импорта
    unsigned int flags =
        aiProcess_Triangulate |           // Конвертируем в треугольники
        aiProcess_FlipUVs |               // Переворачиваем UV по Y
        aiProcess_GenNormals |            // Генерируем нормали если их нет
        aiProcess_CalcTangentSpace |      // Вычисляем тангентное пространство
        aiProcess_JoinIdenticalVertices | // Объединяем одинаковые вершины
        aiProcess_ImproveCacheLocality |  // Оптимизируем для кэша
        aiProcess_OptimizeMeshes |        // Оптимизируем меши
        aiProcess_ValidateDataStructure;  // Проверяем структуру данных

    scene = importer.ReadFile(path, flags);

    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "ASSIMP ERROR: " << importer.GetErrorString() << std::endl;
        return false;
    }

    // Обрабатываем все меши
    processNode(scene->mRootNode);

    // Загружаем материалы
    if (scene->HasMaterials()) {
        materials.resize(scene->mNumMaterials);
        for (unsigned int i = 0; i < scene->mNumMaterials; ++i) {
            aiMaterial* mat = scene->mMaterials[i];
            if (!loadTexturesFromMaterial(mat, materials[i])) {
                std::cerr << "Failed to load textures for material " << i << std::endl;
            }
        }
    }
    else {
        // Создаем дефолтный материал
        materials.emplace_back();
        std::cout << "No materials found, using default" << std::endl;
    }

    return true;
}

void Model::processNode(aiNode* node) {
    // Обрабатываем все меши текущего узла
    for (unsigned int i = 0; i < node->mNumMeshes; i++) {
        aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
        processMesh(mesh);
    }

    // Рекурсивно обрабатываем дочерние узлы
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        processNode(node->mChildren[i]);
    }
}

void Model::processMesh(aiMesh* mesh) {
    SubMesh subMesh;
    subMesh.name = mesh->mName.C_Str();
    subMesh.startIndex = indices.size();
    subMesh.materialIndex = mesh->mMaterialIndex;
    subMesh.hasNormals = mesh->HasNormals();
    subMesh.hasTexCoords = mesh->HasTextureCoords(0);
    subMesh.hasTangents = mesh->HasTangentsAndBitangents();

    // Начальный индекс вершин для этого меша
    unsigned int vertexStart = vertices.size();

    // Обрабатываем вершины
    for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
        Vertex vertex;

        // Позиция
        vertex.position.x = mesh->mVertices[i].x;
        vertex.position.y = mesh->mVertices[i].y;
        vertex.position.z = mesh->mVertices[i].z;

        // Нормаль
        if (mesh->HasNormals()) {
            vertex.normal.x = mesh->mNormals[i].x;
            vertex.normal.y = mesh->mNormals[i].y;
            vertex.normal.z = mesh->mNormals[i].z;
        }

        // Текстурные координаты
        if (mesh->HasTextureCoords(0)) {
            vertex.texcoord.x = mesh->mTextureCoords[0][i].x;
            vertex.texcoord.y = mesh->mTextureCoords[0][i].y;
        }

        // Тангенты и битангенты
        if (mesh->HasTangentsAndBitangents()) {
            vertex.tangent.x = mesh->mTangents[i].x;
            vertex.tangent.y = mesh->mTangents[i].y;
            vertex.tangent.z = mesh->mTangents[i].z;

            vertex.bitangent.x = mesh->mBitangents[i].x;
            vertex.bitangent.y = mesh->mBitangents[i].y;
            vertex.bitangent.z = mesh->mBitangents[i].z;
        }

        vertices.push_back(vertex);
    }

    // Обрабатываем индексы (грани)
    for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
        aiFace face = mesh->mFaces[i];
        for (unsigned int j = 0; j < face.mNumIndices; j++) {
            indices.push_back(vertexStart + face.mIndices[j]);
        }
    }

    subMesh.indexCount = indices.size() - subMesh.startIndex;
    subMeshes.push_back(subMesh);

    std::cout << "Processed mesh: " << subMesh.name
        << " (vertices: " << mesh->mNumVertices
        << ", faces: " << mesh->mNumFaces << ")" << std::endl;
}

bool Model::loadTexturesFromMaterial(aiMaterial* mat, MaterialTextures& textures) {
    bool loadedAny = false;

    // Загружаем диффузную текстуру
    if (loadTexture(mat, aiTextureType_DIFFUSE, textures.diffuse, "_diffuse.tga")) {
        loadedAny = true;
        std::cout << "  Loaded diffuse texture" << std::endl;
    }

    // Загружаем нормал-мап (может быть в aiTextureType_NORMALS или aiTextureType_HEIGHT)
    if (!loadTexture(mat, aiTextureType_NORMALS, textures.normal, "_nm.tga")) {
        loadTexture(mat, aiTextureType_HEIGHT, textures.normal, "_nm.tga");
    }

    // Загружаем specular карту
    loadTexture(mat, aiTextureType_SPECULAR, textures.specular, "_spec.tga");

    // Загружаем emission карту
    loadTexture(mat, aiTextureType_EMISSIVE, textures.emission, "_emission.tga");

    return loadedAny;
}

bool Model::loadTexture(aiMaterial* mat, aiTextureType type, TGAImage& image,
    const std::string& defaultName) {
    aiString texturePath;
    if (mat->GetTextureCount(type) > 0 &&
        mat->GetTexture(type, 0, &texturePath) == AI_SUCCESS) {

        std::filesystem::path fullPath = std::filesystem::path(directory) / texturePath.C_Str();

        // Пробуем несколько расширений
        std::string extensions[] = { ".tga", ".png", ".jpg", ".bmp" };
        for (const auto& ext : extensions) {
            std::string testPath = fullPath.string();
            if (!testPath.empty() && testPath.find('.') == std::string::npos) {
                testPath += ext;
            }

            if (image.read_tga_file(testPath.c_str())) {
                return true;
            }
        }

        std::cerr << "Failed to load texture: " << fullPath.string() << std::endl;
        return false;
    }

    // Fallback: пробуем загрузить по имени модели
    if (!defaultName.empty()) {
        std::filesystem::path modelPath(filename);
        std::string stem = modelPath.stem().string();
        std::string fallbackPath = directory + "/" + stem + defaultName;

        if (image.read_tga_file(fallbackPath.c_str())) {
            return true;
        }
    }

    return false;
}

void Model::generateNormalsIfNeeded() {
    bool needsNormals = false;
    for (auto& vertex : vertices) {
        if (norm(vertex.normal) < 0.001) {
            needsNormals = true;
            break;
        }
    }

    if (!needsNormals) return;

    std::cout << "Generating normals..." << std::endl;

    // Инициализируем все нормали нулями
    for (auto& vertex : vertices) {
        vertex.normal = vec3{ 0, 0, 0 };
    }

    // Вычисляем нормали для каждой грани и добавляем к вершинам
    for (size_t i = 0; i < indices.size(); i += 3) {
        unsigned int i0 = indices[i];
        unsigned int i1 = indices[i + 1];
        unsigned int i2 = indices[i + 2];

        vec3 v0 = vertices[i0].position;
        vec3 v1 = vertices[i1].position;
        vec3 v2 = vertices[i2].position;

        vec3 edge1 = v1 - v0;
        vec3 edge2 = v2 - v0;
        vec3 faceNormal = cross(edge1, edge2);

        // Добавляем к каждой вершине грани
        vertices[i0].normal = vertices[i0].normal + faceNormal;
        vertices[i1].normal = vertices[i1].normal + faceNormal;
        vertices[i2].normal = vertices[i2].normal + faceNormal;
    }

    // Нормализуем
    for (auto& vertex : vertices) {
        if (norm(vertex.normal) > 0.001) {
            vertex.normal = normalized(vertex.normal);
        }
        else {
            vertex.normal = vec3{ 0, 0, 1 };
        }
    }
}

void Model::computeTangentsIfNeeded() {
    bool needsTangents = false;
    for (auto& vertex : vertices) {
        if (norm(vertex.tangent) < 0.001) {
            needsTangents = true;
            break;
        }
    }

    if (!needsTangents) return;

    std::cout << "Computing tangents..." << std::endl;

    // Инициализируем все тангенты нулями
    for (auto& vertex : vertices) {
        vertex.tangent = vec3{ 0, 0, 0 };
        vertex.bitangent = vec3{ 0, 0, 0 };
    }

    // Вычисляем тангенты для каждой грани
    for (size_t i = 0; i < indices.size(); i += 3) {
        unsigned int i0 = indices[i];
        unsigned int i1 = indices[i + 1];
        unsigned int i2 = indices[i + 2];

        Vertex& v0 = vertices[i0];
        Vertex& v1 = vertices[i1];
        Vertex& v2 = vertices[i2];

        vec3 deltaPos1 = v1.position - v0.position;
        vec3 deltaPos2 = v2.position - v0.position;

        vec2 deltaUV1 = v1.texcoord - v0.texcoord;
        vec2 deltaUV2 = v2.texcoord - v0.texcoord;

        double r = deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y;
        if (std::abs(r) < 1e-8) continue;

        double invr = 1.0 / r;

        vec3 tangent = (deltaPos1 * deltaUV2.y - deltaPos2 * deltaUV1.y) * invr;
        vec3 bitangent = (deltaPos2 * deltaUV1.x - deltaPos1 * deltaUV2.x) * invr;

        v0.tangent = v0.tangent + tangent;
        v1.tangent = v1.tangent + tangent;
        v2.tangent = v2.tangent + tangent;

        v0.bitangent = v0.bitangent + bitangent;
        v1.bitangent = v1.bitangent + bitangent;
        v2.bitangent = v2.bitangent + bitangent;
    }

    // Ортогонализация по Граму-Шмидту
    for (auto& vertex : vertices) {
        if (norm(vertex.tangent) > 0.001 && norm(vertex.normal) > 0.001) {
            // Ортогонализируем тангент относительно нормали
            vec3 n = normalized(vertex.normal);
            vec3 t = normalized(vertex.tangent);

            // t = normalize(t - n * dot(n, t))
            vertex.tangent = normalized(t - n * dot(n, t));

            // Пересчитываем битангент
            vertex.bitangent = cross(vertex.normal, vertex.tangent);
        }
        else {
            vertex.tangent = vec3{ 1, 0, 0 };
            vertex.bitangent = vec3{ 0, 1, 0 };
        }
    }
}

// Методы доступа для обратной совместимости
vec3 Model::vert(int i) const {
    if (i < 0 || i >= (int)vertices.size()) return vec3{ 0, 0, 0 };
    return vertices[i].position;
}

vec3 Model::vert(int iface, int nthvert) const {
    int idx = iface * 3 + nthvert;
    if (idx < 0 || idx >= (int)indices.size()) return vec3{ 0, 0, 0 };
    return vertices[indices[idx]].position;
}

vec3 Model::normal(int iface, int nthvert) const {
    int idx = iface * 3 + nthvert;
    if (idx < 0 || idx >= (int)indices.size()) return vec3{ 0, 0, 1 };
    return vertices[indices[idx]].normal;
}

vec2 Model::uv(int iface, int nthvert) const {
    int idx = iface * 3 + nthvert;
    if (idx < 0 || idx >= (int)indices.size()) return vec2{ 0, 0 };
    return vertices[indices[idx]].texcoord;
}

// Методы для текстур (используют первый материал)
TGAColor Model::diffuse(const vec2& uv) const {
    if (materials.empty() || !materials[0].hasDiffuse()) {
        return TGAColor(255, 255, 255, 255);
    }

    int x = std::clamp(int(uv.x * materials[0].diffuse.width()),
        0, materials[0].diffuse.width() - 1);
    int y = std::clamp(int(uv.y * materials[0].diffuse.height()),
        0, materials[0].diffuse.height() - 1);

    return materials[0].diffuse.get(x, y);
}

vec3 Model::normal(const vec2& uv) const {
    if (materials.empty() || !materials[0].hasNormal()) {
        return vec3{ 0, 0, 1 };
    }

    int x = std::clamp(int(uv.x * materials[0].normal.width()),
        0, materials[0].normal.width() - 1);
    int y = std::clamp(int(uv.y * materials[0].normal.height()),
        0, materials[0].normal.height() - 1);

    TGAColor c = materials[0].normal.get(x, y);
    vec3 n;
    n.x = (double)c[2] / 255.0 * 2.0 - 1.0;  // R
    n.y = (double)c[1] / 255.0 * 2.0 - 1.0;  // G  
    n.z = (double)c[0] / 255.0 * 2.0 - 1.0;  // B

    return normalized(n);
}

float Model::specular(const vec2& uv) const {
    if (materials.empty() || !materials[0].hasSpecular()) {
        return 1.0f;
    }

    int x = std::clamp(int(uv.x * materials[0].specular.width()),
        0, materials[0].specular.width() - 1);
    int y = std::clamp(int(uv.y * materials[0].specular.height()),
        0, materials[0].specular.height() - 1);

    TGAColor c = materials[0].specular.get(x, y);
    return c[0] / 255.0f;
}

TGAColor Model::emission(const vec2& uv) const {
    if (materials.empty() || !materials[0].hasEmission()) {
        return TGAColor(0, 0, 0, 255);
    }

    int x = std::clamp(int(uv.x * materials[0].emission.width()),
        0, materials[0].emission.width() - 1);
    int y = std::clamp(int(uv.y * materials[0].emission.height()),
        0, materials[0].emission.height() - 1);

    return materials[0].emission.get(x, y);
}