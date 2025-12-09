#include "model.h"
#include <iostream>
#include <filesystem>

Model::Model(const std::string filename) {
    // Сохраняем базовую директорию для текстур
    std::filesystem::path path(filename);
    baseDirectory = path.parent_path().string();

    // Загружаем модель через Assimp
    loadModel(filename);

    // Вычисляем тангенсы (остается как было, но после загрузки данных)
    tangents.assign(verts.size(), vec3{ 0,0,0 });

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

    // Нормализуем тангенсы
    for (size_t i = 0; i < tangents.size(); i++) {
        tangents[i] = normalized(tangents[i]);
    }
}

void Model::loadModel(const std::string& path) {
    // Загружаем модель через Assimp
    scene = importer.ReadFile(path,
        aiProcess_Triangulate |           // Треугольники
        aiProcess_FlipUVs |               // Переворачиваем UV
        aiProcess_GenNormals |            // Генерируем нормали
        aiProcess_CalcTangentSpace |      // Вычисляем тангенты
        aiProcess_JoinIdenticalVertices); // Объединяем одинаковые вершины

    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "ASSIMP ERROR: " << importer.GetErrorString() << std::endl;
        return;
    }

    // Рекурсивно обрабатываем все меши
    processNode(scene->mRootNode);

    std::cerr << "# v# " << nverts() << " f# " << nfaces() << std::endl;

    // Загружаем текстуры
    // Пробуем стандартные имена текстур (как в старом коде)
    std::filesystem::path modelPath(path);
    std::string baseName = modelPath.stem().string();
    std::string dir = modelPath.parent_path().string();

    if (!dir.empty()) dir += "/";

    // Диффузная текстура
    if (!diffusemap.read_tga_file((dir + baseName + "_diffuse.tga").c_str())) {
        std::cerr << "Diffuse texture not found: " << (dir + baseName + "_diffuse.tga") << std::endl;
    }

    // Нормал-мап
    if (!normalmap.read_tga_file((dir + baseName + "_nm.tga").c_str())) {
        std::cerr << "Normal map not found: " << (dir + baseName + "_nm.tga") << std::endl;
    }

    // Specular карта
    if (!specmap.read_tga_file((dir + baseName + "_spec.tga").c_str())) {
        std::cerr << "Specular map not found: " << (dir + baseName + "_spec.tga") << std::endl;
    }
}

void Model::processNode(aiNode* node) {
    // Обрабатываем все меши текущего узла
    for (unsigned int i = 0; i < node->mNumMeshes; i++) {
        aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
        processMesh(mesh, scene);
    }

    // Рекурсивно обрабатываем дочерние узлы
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        processNode(node->mChildren[i]);
    }
}

void Model::processMesh(aiMesh* mesh, const aiScene* scene) {
    // Начальный индекс для текущего меша
    size_t vertexStart = verts.size();

    // Вершины
    for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
        vec3 vertex;
        vertex.x = mesh->mVertices[i].x;
        vertex.y = mesh->mVertices[i].y;
        vertex.z = mesh->mVertices[i].z;
        verts.push_back(vertex);

        // Нормали
        if (mesh->mNormals) {
            vec3 normal;
            normal.x = mesh->mNormals[i].x;
            normal.y = mesh->mNormals[i].y;
            normal.z = mesh->mNormals[i].z;
            norms.push_back(normalized(normal));
        }
        else {
            norms.push_back(vec3{ 0,0,1 });
        }

        // Текстурные координаты (берем первый набор UV)
        if (mesh->mTextureCoords[0]) {
            vec2 uv;
            uv.x = mesh->mTextureCoords[0][i].x;
            uv.y = mesh->mTextureCoords[0][i].y;
            tex.push_back(uv);
        }
        else {
            tex.push_back(vec2{ 0,0 });
        }
    }

    // Индексы (грани) - преобразуем в плоский массив
    for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
        aiFace face = mesh->mFaces[i];
        // Меш уже триангулирован (благодаря aiProcess_Triangulate)
        for (unsigned int j = 0; j < face.mNumIndices; j++) {
            // Смещаем индексы на начало текущего меша
            int index = vertexStart + face.mIndices[j];
            facet_vrt.push_back(index);
            facet_tex.push_back(index);
            facet_nrm.push_back(index);
        }
    }
}

// Методы доступа остаются БЕЗ ИЗМЕНЕНИЙ (как в вашем старом коде)

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
    if (normalmap.width() == 0) {
        return vec3{ 0, 0, 1 };
    }

    TGAColor c = normalmap.get(
        uv.x * normalmap.width(),
        uv.y * normalmap.height()
    );

    vec3 n;
    n.x = (double)c[2] / 255.0 * 2.0 - 1.0;  // R
    n.y = (double)c[1] / 255.0 * 2.0 - 1.0;  // G  
    n.z = (double)c[0] / 255.0 * 2.0 - 1.0;  // B

    return normalized(n);
}

vec2 Model::uv(int iface, int nthvert) const {
    return tex[facet_tex[iface * 3 + nthvert]];
}

TGAColor Model::diffuse(const vec2& uv) const {
    return diffusemap.get(
        uv.x * diffusemap.width(),
        uv.y * diffusemap.height()
    );
}

float Model::specular(const vec2& uv) const {
    TGAColor c = specmap.get(
        uv.x * specmap.width(),
        uv.y * specmap.height()
    );
    return c[0] / 1.0f;
}