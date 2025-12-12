// ================================================================
// camera.h - Камера для учебного рендерера (фиксированный рендеринг)
// ================================================================
#pragma once
#include "geometry.h"
#include <cmath>
#include <algorithm>

class Camera {
public:
    // Стандартные пресеты для разных типов сцен
    enum Preset {
        SPONZA_SCENE,     // Большая сцена (архитектура)
        CHARACTER_CLOSEUP, // Крупный план персонажа
        OVERVIEW,         // Обзор сверху
        DEFAULT           // Стандартная настройка
    };

    // Основные параметры камеры
    struct Params {
        vec3 eye = { 0, 0, 10 };
        vec3 target = { 0, 0, 0 };
        vec3 up = { 0, 1, 0 };

        double fov = 60.0;           // в градусах
        double aspect = 16.0 / 9.0;
        double nearPlane = 0.1;
        double farPlane = 1000.0;
    };

    Camera() = default;

    // Создание камеры с пресетом
    Camera(Preset preset, double aspectRatio = 16.0 / 9.0) {
        setPreset(preset, aspectRatio);
    }

    // Установка пресета
    void setPreset(Preset preset, double aspectRatio = 16.0 / 9.0) {
        params.aspect = aspectRatio;

        switch (preset) {
        case SPONZA_SCENE:
            // Для больших сцен типа Sponza
            params.eye = { 0, 15, 40 };
            params.target = { 0, 10, 0 };
            params.fov = 55.0;
            params.nearPlane = 0.5;
            params.farPlane = 500.0;
            break;

        case CHARACTER_CLOSEUP:
            // Для крупных планов персонажей
            params.eye = { 0, 5, 12 };
            params.target = { 0, 4, 0 };
            params.fov = 45.0;
            params.nearPlane = 0.1;
            params.farPlane = 100.0;
            break;

        case OVERVIEW:
            // Вид сверху
            params.eye = { 0, 50, 0 };
            params.target = { 0, 0, 0 };
            params.up = { 0, 0, -1 };
            params.fov = 60.0;
            params.nearPlane = 1.0;
            params.farPlane = 200.0;
            break;

        case DEFAULT:
        default:
            params.eye = { 0, 0, 10 };
            params.target = { 0, 0, 0 };
            params.fov = 60.0;
            params.nearPlane = 0.1;
            params.farPlane = 200.0;
            break;
        }

        updateMatrices();
    }

    // Автоматическая настройка под сцену (самое важное!)
    void autoSetupForScene(const AABB& sceneBounds, double aspectRatio = 16.0 / 9.0) {
        params.aspect = aspectRatio;

        // Центр сцены
        vec3 center = (sceneBounds.min + sceneBounds.max) * 0.5;
        vec3 size = sceneBounds.max - sceneBounds.min;

        // Наибольший размер сцены
        double maxDim = std::max({ size.x, size.y, size.z });

        // Рассчитываем расстояние для обзора всей сцены
        double fovRad = params.fov * M_PI / 180.0;
        double requiredDistance = (maxDim * 1.5) / (2.0 * tan(fovRad / 2.0));

        // Корректируем с учетом соотношения сторон
        if (params.aspect > 1.0) {
            requiredDistance *= params.aspect;
        }

        // Ограничиваем расстояние
        requiredDistance = std::max(5.0, std::min(requiredDistance, 200.0));

        // Позиция камеры - смотрим на сцену под углом 30 градусов сверху
        params.eye = center + vec3{ 0, requiredDistance * 0.5, requiredDistance };
        params.target = center;

        // Автоматически настраиваем far plane
        double sceneRadius = maxDim * 0.5;
        params.farPlane = std::max(100.0, requiredDistance + sceneRadius * 3.0);

        updateMatrices();
    }

    // Настройка для нескольких моделей
    void setupForMultipleModels(const std::vector<AABB>& modelBounds, double aspectRatio = 16.0 / 9.0) {
        if (modelBounds.empty()) {
            setPreset(DEFAULT, aspectRatio);
            return;
        }

        // Объединяем все AABB
        vec3 overallMin = modelBounds[0].min;
        vec3 overallMax = modelBounds[0].max;

        for (size_t i = 1; i < modelBounds.size(); ++i) {
            overallMin.x = std::min(overallMin.x, modelBounds[i].min.x);
            overallMin.y = std::min(overallMin.y, modelBounds[i].min.y);
            overallMin.z = std::min(overallMin.z, modelBounds[i].min.z);

            overallMax.x = std::max(overallMax.x, modelBounds[i].max.x);
            overallMax.y = std::max(overallMax.y, modelBounds[i].max.y);
            overallMax.z = std::max(overallMax.z, modelBounds[i].max.z);
        }

        AABB sceneBounds(overallMin, overallMax);
        autoSetupForScene(sceneBounds, aspectRatio);
    }

    // Обновление матриц
    void updateMatrices() {
        updateViewMatrix();
        updateProjectionMatrix();
    }

    // Получение матриц (совместимость со старым кодом)
    mat<4, 4> getViewMatrix() const { return viewMatrix; }
    mat<4, 4> getProjectionMatrix() const { return projectionMatrix; }
    mat<4, 4> getViewProjectionMatrix() const { return projectionMatrix * viewMatrix; }

    // Геттеры для параметров
    const Params& getParams() const { return params; }
    vec3 getEye() const { return params.eye; }
    vec3 getTarget() const { return params.target; }
    vec3 getUp() const { return params.up; }
    double getFOV() const { return params.fov; }
    double getAspect() const { return params.aspect; }
    double getNear() const { return params.nearPlane; }
    double getFar() const { return params.farPlane; }

    // Ручная настройка (если нужно)
    void setEye(const vec3& eye) { params.eye = eye; updateViewMatrix(); }
    void setTarget(const vec3& target) { params.target = target; updateViewMatrix(); }
    void setUp(const vec3& up) { params.up = up; updateViewMatrix(); }
    void setFOV(double fov) { params.fov = fov; updateProjectionMatrix(); }
    void setAspect(double aspect) { params.aspect = aspect; updateProjectionMatrix(); }
    void setClipping(double near, double far) {
        params.nearPlane = near;
        params.farPlane = far;
        updateProjectionMatrix();
    }

    // Для отладки
    void printInfo() const {
        std::cout << "Camera Info:" << std::endl;
        std::cout << "  Eye: (" << params.eye.x << ", " << params.eye.y << ", " << params.eye.z << ")" << std::endl;
        std::cout << "  Target: (" << params.target.x << ", " << params.target.y << ", " << params.target.z << ")" << std::endl;
        std::cout << "  FOV: " << params.fov << " degrees" << std::endl;
        std::cout << "  Aspect: " << params.aspect << std::endl;
        std::cout << "  Clipping: " << params.nearPlane << " - " << params.farPlane << std::endl;
        std::cout << "  Distance to target: " << norm(params.eye - params.target) << std::endl;
    }

private:
    Params params;
    mat<4, 4> viewMatrix;
    mat<4, 4> projectionMatrix;

    void updateViewMatrix() {
        vec3 zAxis = normalized(params.eye - params.target);
        vec3 xAxis = normalized(cross(params.up, zAxis));
        vec3 yAxis = cross(zAxis, xAxis);

        viewMatrix = mat<4, 4>::identity();
        viewMatrix[0][0] = xAxis.x; viewMatrix[0][1] = xAxis.y; viewMatrix[0][2] = xAxis.z;
        viewMatrix[1][0] = yAxis.x; viewMatrix[1][1] = yAxis.y; viewMatrix[1][2] = yAxis.z;
        viewMatrix[2][0] = zAxis.x; viewMatrix[2][1] = zAxis.y; viewMatrix[2][2] = zAxis.z;

        viewMatrix[0][3] = -dot(xAxis, params.eye);
        viewMatrix[1][3] = -dot(yAxis, params.eye);
        viewMatrix[2][3] = -dot(zAxis, params.eye);
    }

    void updateProjectionMatrix() {
        double fovRad = params.fov * M_PI / 180.0;
        double tanHalfFov = tan(fovRad / 2.0);

        projectionMatrix = mat<4, 4>::identity();
        projectionMatrix[0][0] = 1.0 / (params.aspect * tanHalfFov);
        projectionMatrix[1][1] = 1.0 / tanHalfFov;
        projectionMatrix[2][2] = (params.farPlane + params.nearPlane) / (params.nearPlane - params.farPlane);
        projectionMatrix[2][3] = (2.0 * params.farPlane * params.nearPlane) / (params.nearPlane - params.farPlane);
        projectionMatrix[3][2] = -1.0;
        projectionMatrix[3][3] = 0.0;
    }
};

// ================================================================
// Вспомогательные функции для совместимости со старым кодом
// ================================================================

// Заменяет старый lookat (для обратной совместимости)
inline void lookat(const vec3& eye, const vec3& center, const vec3& up) {
    // Эта функция теперь не рекомендуется, используйте класс Camera
    // Оставлена для совместимости со старым кодом
}

// Функция для быстрой настройки камеры в main.cpp
inline void setupCameraForRendering(Camera& camera, const std::vector<AABB>& modelBounds,
    int width, int height, bool autoAdjust = true) {
    if (autoAdjust && !modelBounds.empty()) {
        camera.setupForMultipleModels(modelBounds, (double)width / height);
    }
    else {
        camera.setPreset(Camera::SPONZA_SCENE, (double)width / height);
    }

    camera.printInfo();
}