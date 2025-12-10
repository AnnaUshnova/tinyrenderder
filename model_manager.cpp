#include "model_manager.h"
#include <iostream>
#include <filesystem>
#include <iomanip>

std::shared_ptr<Model> ModelManager::loadModel(const std::string& path) {
    std::lock_guard<std::mutex> lock(cacheMutex);

    // Нормализуем путь
    std::filesystem::path fsPath(path);
    std::string canonicalPath = std::filesystem::canonical(fsPath).string();

    // Проверяем кэш
    auto it = modelCache.find(canonicalPath);
    if (it != modelCache.end()) {
        if (auto model = it->second.lock()) {
            std::cout << "Model cache hit: " << canonicalPath << std::endl;
            return model;
        }
        // weak_ptr истек, удаляем из кэша
        modelCache.erase(it);
    }

    // Создаем новую модель
    auto model = std::make_shared<Model>(canonicalPath);
    if (!model->load()) {
        std::cerr << "Failed to load model: " << canonicalPath << std::endl;
        return nullptr;
    }

    // Сохраняем в кэш
    modelCache[canonicalPath] = model;
    std::cout << "Model loaded and cached: " << canonicalPath << std::endl;

    return model;
}

std::shared_ptr<Model> ModelManager::getModel(const std::string& path) {
    return loadModel(path);
}

bool ModelManager::unloadModel(const std::string& path) {
    std::lock_guard<std::mutex> lock(cacheMutex);

    std::filesystem::path fsPath(path);
    std::string canonicalPath = std::filesystem::canonical(fsPath).string();

    auto it = modelCache.find(canonicalPath);
    if (it != modelCache.end()) {
        if (auto model = it->second.lock()) {
            model->unload();
        }
        modelCache.erase(it);
        std::cout << "Model unloaded from cache: " << canonicalPath << std::endl;
        return true;
    }

    return false;
}

void ModelManager::unloadAll() {
    std::lock_guard<std::mutex> lock(cacheMutex);

    for (auto& pair : modelCache) {
        if (auto model = pair.second.lock()) {
            model->unload();
        }
    }

    modelCache.clear();
    std::cout << "All models unloaded from cache" << std::endl;
}

void ModelManager::printStats() const {
    std::lock_guard<std::mutex> lock(cacheMutex);

    std::cout << "\n=== Model Manager Statistics ===" << std::endl;
    std::cout << "Cached models: " << modelCache.size() << std::endl;

    size_t loadedCount = 0;
    for (const auto& pair : modelCache) {
        if (auto model = pair.second.lock()) {
            loadedCount++;
            std::cout << "  - " << std::filesystem::path(pair.first).filename().string()
                << " (refs: " << model.use_count() - 1 << ")" << std::endl;
        }
    }

    std::cout << "Loaded models: " << loadedCount << std::endl;
    std::cout << "=============================\n" << std::endl;
}