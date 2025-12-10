#pragma once
#include "model.h"
#include <memory>
#include <unordered_map>
#include <mutex>
#include <string>

class ModelManager {
public:
    // Singleton
    static ModelManager& getInstance() {
        static ModelManager instance;
        return instance;
    }

    // Запрещаем копирование
    ModelManager(const ModelManager&) = delete;
    ModelManager& operator=(const ModelManager&) = delete;

    // Основные методы
    std::shared_ptr<Model> loadModel(const std::string& path);
    std::shared_ptr<Model> getModel(const std::string& path);
    bool unloadModel(const std::string& path);
    void unloadAll();

    // Утилиты
    size_t getModelCount() const {
        std::lock_guard<std::mutex> lock(cacheMutex);
        return modelCache.size();
    }

    void printStats() const;

private:
    ModelManager() = default;
    ~ModelManager() { unloadAll(); }

    std::unordered_map<std::string, std::weak_ptr<Model>> modelCache;
    mutable std::mutex cacheMutex;
};