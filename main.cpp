// ================================================================
// main.cpp - Основной файл рендеринга с Phong освещением и SSAO
// ================================================================
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>
#include <memory>

#include "tgaimage.h"
#include "geometry.h"
#include "our_gl.h"
#include "model.h"
#include "model_manager.h"
#include "camera.h"

// Явное объявление глобальных переменных
extern mat<4, 4> ModelView;
extern mat<4, 4> Perspective;
extern mat<4, 4> Viewport;
extern std::vector<double> zbuffer;

// Константы рендеринга
const int WIDTH = 1200;  // Увеличим размер окна для Sponza
const int HEIGHT = 800;
const char* DEFAULT_MODEL_PATH = "obj/african_head/african_head.obj";
const char* EYES_MODEL_PATH = "obj/african_head/african_head_eye_inner.obj";
const char* SPONZA_MODEL_PATH = "obj/sponza/sponza.obj";

// Настройки защиты глаз от нормал-мапа
constexpr double EYE_DIFFUSE_BRIGHTNESS_THRESHOLD = 0.85;
constexpr double EYE_SPECULAR_POWER_THRESHOLD = 5.0;

// ================================================================
//  Шейдер Phong освещения (Key, Fill, Rim + нормал-мап)
// ================================================================
struct PhongShader : public IShader {
    const Model* model;

    vec3 key_light_dir_eye;   // Направление ключевого света в пространстве камеры
    vec3 fill_light_dir_eye;  // Направление заполняющего света
    vec3 rim_light_dir_eye;   // Направление контурного света

    // Интерполируемые атрибуты
    vec2 varying_uv[3];
    vec3 varying_position_eye[3];
    vec3 varying_normal_eye[3];

    double normal_map_strength = 1.0;

    PhongShader(const Model* model_ptr) : model(model_ptr) {}

    void initLightDirections(const vec3& key_dir_world,
        const vec3& fill_dir_world,
        const vec3& rim_dir_world) {
        // Преобразование мировых направлений света в пространство камеры
        mat<3, 3> normal_matrix;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                normal_matrix[i][j] = ModelView[i][j];
            }
        }

        key_light_dir_eye = normalized(normal_matrix * key_dir_world);
        fill_light_dir_eye = normalized(normal_matrix * fill_dir_world);
        rim_light_dir_eye = normalized(normal_matrix * rim_dir_world);
    }

    vec4 vertex(int face_index, int vertex_index) override {
        vec3 vertex_position = model->vert(face_index, vertex_index);
        vec3 vertex_normal = model->normal(face_index, vertex_index);

        varying_uv[vertex_index] = model->uv(face_index, vertex_index);

        vec4 position_eye = ModelView * make_vec4(vertex_position[0],
            vertex_position[1],
            vertex_position[2],
            1.0);
        varying_position_eye[vertex_index] = position_eye.xyz();

        vec4 normal_eye = ModelView * make_vec4(vertex_normal[0],
            vertex_normal[1],
            vertex_normal[2],
            0.0);
        varying_normal_eye[vertex_index] = normal_eye.xyz();

        return Perspective * position_eye;
    }

    std::pair<bool, TGAColor> fragment(const vec3 barycentric) const override {
        // Интерполяция атрибутов
        vec3 position_eye = varying_position_eye[0] * barycentric[0]
            + varying_position_eye[1] * barycentric[1]
            + varying_position_eye[2] * barycentric[2];

        vec3 geometry_normal = varying_normal_eye[0] * barycentric[0]
            + varying_normal_eye[1] * barycentric[1]
            + varying_normal_eye[2] * barycentric[2];

        vec2 uv = varying_uv[0] * barycentric[0]
            + varying_uv[1] * barycentric[1]
            + varying_uv[2] * barycentric[2];

        TGAColor base_color = model->diffuse(uv);
        double specular_power = std::max(1.0, (double)model->specular(uv));

        // Определение, является ли пиксель частью глаза
        double brightness = (base_color[0] + base_color[1] + base_color[2]) / (3.0 * 255.0);
        bool is_eye_pixel = (brightness >= EYE_DIFFUSE_BRIGHTNESS_THRESHOLD)
            && (specular_power <= EYE_SPECULAR_POWER_THRESHOLD);

        // Нормаль из нормал-мапа
        vec3 normal_map_value = model->normal(uv);
        vec3 normal_map_eye = (ModelView * make_vec4(normal_map_value[0],
            normal_map_value[1],
            normal_map_value[2],
            0.0)).xyz();

        // Финальная нормаль: смешивание геометрической и нормал-мап нормали
        vec3 final_normal = is_eye_pixel
            ? geometry_normal
            : normalized(geometry_normal * (1.0 - normal_map_strength)
                + normal_map_eye * normal_map_strength);

        vec3 view_direction = normalized(-position_eye);  // Камера в начале координат

        // Интенсивности источников света
        const double KEY_DIFFUSE_INTENSITY = 1.0;
        const double KEY_SPECULAR_INTENSITY = 1.0;
        const double FILL_DIFFUSE_INTENSITY = 0.35;
        const double RIM_DIFFUSE_INTENSITY = 0.6;

        // Ключевой свет (диффузный + зеркальный)
        double key_diffuse = std::max(0.0, dot(final_normal, key_light_dir_eye))
            * KEY_DIFFUSE_INTENSITY;
        double key_specular = 0.0;

        if (KEY_SPECULAR_INTENSITY > 0.0) {
            vec3 reflect_dir = normalized(final_normal * (2.0 * dot(final_normal, key_light_dir_eye))
                - key_light_dir_eye);
            double reflect_view_dot = std::max(0.0, dot(reflect_dir, view_direction));
            key_specular = (reflect_view_dot > 0.0 ? std::pow(reflect_view_dot, specular_power) : 0.0)
                * KEY_SPECULAR_INTENSITY;
        }

        // Заполняющий свет (только диффузный)
        double fill_diffuse = std::max(0.0, dot(final_normal, fill_light_dir_eye))
            * FILL_DIFFUSE_INTENSITY;

        // Контурный свет (для акцента силуэта)
        double rim_diffuse = std::max(0.0, dot(final_normal, rim_light_dir_eye))
            * RIM_DIFFUSE_INTENSITY;

        double total_diffuse = key_diffuse + fill_diffuse + rim_diffuse;
        double total_specular = key_specular;
        double ambient = 0.10;

        // Финальный цвет
        TGAColor result = base_color;
        for (int channel = 0; channel < 3; ++channel) {
            double channel_value = base_color[channel];
            double final_value = channel_value * (ambient + total_diffuse)
                + 255.0 * (0.35 * total_specular);
            result[channel] = (unsigned char)std::min(255.0, final_value);
        }

        return { false, result };
    }
};

// ================================================================
//  Шейдер для глаз (глянцевый эффект)
// ================================================================
struct EyeShader : public IShader {
    const Model* model;
    vec3 key_light_dir_eye;
    vec3 rim_light_dir_eye;

    vec2 varying_uv[3];
    vec3 varying_position_eye[3];
    vec3 varying_normal_eye[3];

    EyeShader(const Model* model_ptr) : model(model_ptr) {}

    void initLightDirections(const vec3& key_dir_world, const vec3& rim_dir_world) {
        mat<3, 3> normal_matrix;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                normal_matrix[i][j] = ModelView[i][j];
            }
        }

        key_light_dir_eye = normalized(normal_matrix * key_dir_world);
        rim_light_dir_eye = normalized(normal_matrix * rim_dir_world);
    }

    vec4 vertex(int face_index, int vertex_index) override {
        vec3 vertex_position = model->vert(face_index, vertex_index);
        vec3 vertex_normal = model->normal(face_index, vertex_index);

        varying_uv[vertex_index] = model->uv(face_index, vertex_index);

        vec4 position_eye = ModelView * make_vec4(vertex_position[0],
            vertex_position[1],
            vertex_position[2],
            1.0);
        varying_position_eye[vertex_index] = position_eye.xyz();

        vec4 normal_eye = ModelView * make_vec4(vertex_normal[0],
            vertex_normal[1],
            vertex_normal[2],
            0.0);
        varying_normal_eye[vertex_index] = normal_eye.xyz();

        return Perspective * position_eye;
    }

    std::pair<bool, TGAColor> fragment(const vec3 barycentric) const override {
        vec3 position_eye = varying_position_eye[0] * barycentric[0]
            + varying_position_eye[1] * barycentric[1]
            + varying_position_eye[2] * barycentric[2];

        vec3 normal = normalized(varying_normal_eye[0] * barycentric[0]
            + varying_normal_eye[1] * barycentric[1]
            + varying_normal_eye[2] * barycentric[2]);

        vec2 uv = varying_uv[0] * barycentric[0]
            + varying_uv[1] * barycentric[1]
            + varying_uv[2] * barycentric[2];

        TGAColor base_color = model->diffuse(uv);
        vec3 view_direction = normalized(-position_eye);

        const double KEY_DIFFUSE_INTENSITY = 1.0;
        const double RIM_DIFFUSE_INTENSITY = 0.6;

        double key_diffuse = std::max(0.0, dot(normal, key_light_dir_eye))
            * KEY_DIFFUSE_INTENSITY;
        double rim_diffuse = std::max(0.0, dot(normal, rim_light_dir_eye))
            * RIM_DIFFUSE_INTENSITY;
        double total_diffuse = key_diffuse + rim_diffuse;

        // Усиленная зеркальная составляющая для глаз
        double specular_power = std::max(1.0, (double)model->specular(uv)) * 8.0;
        vec3 reflect_dir = normalized(normal * (2.0 * dot(normal, key_light_dir_eye))
            - key_light_dir_eye);
        double reflect_view_dot = std::max(0.0, dot(reflect_dir, view_direction));
        double specular = (reflect_view_dot > 0.0 ? std::pow(reflect_view_dot, specular_power) : 0.0);

        TGAColor result = base_color;
        for (int channel = 0; channel < 3; ++channel) {
            double channel_value = base_color[channel];
            double final_value = channel_value * (0.1 + total_diffuse)
                + 255.0 * (1.5 * specular);
            result[channel] = (unsigned char)std::min(255.0, final_value);
        }

        return { false, result };
    }
};

// ================================================================
//  Вспомогательные функции
// ================================================================

// Сохранение Z-буфера в виде изображения
static void save_zbuffer_image(const std::vector<double>& zbuffer,
    int width, int height,
    const char* filename) {
    TGAImage image(width, height, TGAImage::RGB);

    // Находим минимальное и максимальное значение глубины
    double min_depth = 1e9, max_depth = -1e9;
    for (double depth : zbuffer) {
        if (std::isfinite(depth)) {
            min_depth = std::min(min_depth, depth);
            max_depth = std::max(max_depth, depth);
        }
    }

    // Если буфер пуст
    if (!std::isfinite(min_depth)) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                image.set(x, y, TGAColor(255, 255, 255));
            }
        }
        image.write_tga_file(filename);
        return;
    }

    if (max_depth - min_depth < 1e-7) {
        max_depth = min_depth + 1e-7;
    }

    // Нормализация и сохранение
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double depth = zbuffer[x + y * width];
            unsigned char value = 255;

            if (std::isfinite(depth)) {
                double normalized = (depth - min_depth) / (max_depth - min_depth);
                value = (unsigned char)(255.0 * (1.0 - normalized));  // Ближе = темнее
            }

            image.set(x, y, TGAColor(value, value, value));
        }
    }

    image.write_tga_file(filename);
}

// Параметры SSAO (Screen Space Ambient Occlusion)
constexpr int AO_NUM_DIRECTIONS = 8;           // Количество направлений для сэмплинга
constexpr int AO_STEPS_PER_DIRECTION = 8;      // Шагов на каждое направление
constexpr double AO_SAMPLE_RADIUS = 16.0;      // Радиус сэмплинга в пикселях
constexpr double AO_OCCLUSION_THRESHOLD = 1e-3; // Порог для определения окклюзии
constexpr double AO_INTENSITY = 0.35;          // Интенсивность эффекта

// Вычисление SSAO в точке (x, y)
static double compute_ssao_at(const std::vector<double>& zbuffer,
    int width, int height,
    int pixel_x, int pixel_y) {
    double center_depth = zbuffer[pixel_x + pixel_y * width];
    if (!std::isfinite(center_depth)) return 1.0;

    int occluded_samples = 0, total_samples = 0;

    for (int direction = 0; direction < AO_NUM_DIRECTIONS; ++direction) {
        double angle = 2.0 * M_PI * direction / AO_NUM_DIRECTIONS;
        double dir_x = cos(angle), dir_y = sin(angle);

        for (int step = 1; step <= AO_STEPS_PER_DIRECTION; ++step) {
            double radius = (double)step / AO_STEPS_PER_DIRECTION * AO_SAMPLE_RADIUS;
            int sample_x = (int)round(pixel_x + dir_x * radius);
            int sample_y = (int)round(pixel_y + dir_y * radius);

            if (sample_x < 0 || sample_x >= width ||
                sample_y < 0 || sample_y >= height)
                continue;

            double sample_depth = zbuffer[sample_x + sample_y * width];
            if (!std::isfinite(sample_depth)) {
                total_samples++;
                continue;
            }

            if (sample_depth < center_depth - AO_OCCLUSION_THRESHOLD) {
                occluded_samples++;
            }
            total_samples++;
        }
    }

    if (total_samples == 0) return 1.0;

    double occlusion_factor = (double)occluded_samples / (double)total_samples;
    return 1.0 - occlusion_factor * AO_INTENSITY;
}

// Создание матрицы масштабирования
mat<4, 4> createScaleMatrix(double sx, double sy, double sz) {
    mat<4, 4> scale = mat<4, 4>::identity();
    scale[0][0] = sx;
    scale[1][1] = sy;
    scale[2][2] = sz;
    return scale;
}

// Создание матрицы переноса
mat<4, 4> createTranslationMatrix(double tx, double ty, double tz) {
    mat<4, 4> translation = mat<4, 4>::identity();
    translation[0][3] = tx;
    translation[1][3] = ty;
    translation[2][3] = tz;
    return translation;
}

// ================================================================
//  Основная функция
// ================================================================
int main(int argc, char** argv) {
    std::cout << "=== Renderer with ModelManager and Frustum Culling ===" << std::endl;

    // ================================================================
    // ИНИЦИАЛИЗАЦИЯ И ЗАГРУЗКА МОДЕЛЕЙ
    // ================================================================

    auto& modelManager = ModelManager::getInstance();

    const char* model_path = (argc > 1) ? argv[1] : DEFAULT_MODEL_PATH;
    const char* eyes_model_path = EYES_MODEL_PATH;
    const char* sponza_model_path = SPONZA_MODEL_PATH;

    std::cout << "Loading head model: " << model_path << std::endl;
    auto head_model = modelManager.loadModel(model_path);

    std::cout << "Loading eye model: " << eyes_model_path << std::endl;
    auto eye_model = modelManager.loadModel(eyes_model_path);

    std::cout << "Loading sponza model: " << sponza_model_path << std::endl;
    auto sponza_model = modelManager.loadModel(sponza_model_path);

    if (!head_model || !eye_model || !sponza_model) {
        std::cerr << "ERROR: Failed to load one or more models!" << std::endl;
        return 1;
    }

    modelManager.printStats();

    // ================================================================
// СОЗДАНИЕ МАТРИЦ МОДЕЛЕЙ
// ================================================================

// Матрица для Sponza (только смещение вниз)
    mat<4, 4> sponzaModelMatrix = createTranslationMatrix(0.0, -5.0, 0.0);

    // Матрица для головы (смещена вправо и вверх)
    mat<4, 4> headModelMatrix = createTranslationMatrix(15.0, 5.0, 0.0);
    mat<4, 4> eyeModelMatrix = headModelMatrix;  // Глаза там же где голова

    // ================================================================
    // АНАЛИЗ СЦЕНЫ
    // ================================================================

    std::cout << "\n=== Scene Analysis ===" << std::endl;

    // Локальные центры (до трансформаций)
    std::cout << "Local centers (before transformations):" << std::endl;
    std::cout << "  Sponza: (" << sponza_model->getCenter().x << ", "
        << sponza_model->getCenter().y << ", " << sponza_model->getCenter().z << ")" << std::endl;
    std::cout << "  Head: (" << head_model->getCenter().x << ", "
        << head_model->getCenter().y << ", " << head_model->getCenter().z << ")" << std::endl;

    // Мировые центры после трансформаций
    AABB sponzaWorldAABB = sponza_model->getWorldAABB(sponzaModelMatrix);
    AABB headWorldAABB = head_model->getWorldAABB(headModelMatrix);
    AABB eyeWorldAABB = eye_model->getWorldAABB(eyeModelMatrix);

    vec3 sponzaWorldCenter = sponzaWorldAABB.getCenter();
    vec3 headWorldCenter = headWorldAABB.getCenter();

    std::cout << "\nWorld centers (after transformations):" << std::endl;
    std::cout << "  Sponza: (" << sponzaWorldCenter.x << ", "
        << sponzaWorldCenter.y << ", " << sponzaWorldCenter.z << ")" << std::endl;
    std::cout << "  Head: (" << headWorldCenter.x << ", "
        << headWorldCenter.y << ", " << headWorldCenter.z << ")" << std::endl;

    // Расстояние между моделями
    vec3 distanceVec = headWorldCenter - sponzaWorldCenter;
    double distance = norm(distanceVec);

    std::cout << "\nDistance between Sponza and Head: " << distance << " units" << std::endl;

    // ================================================================
    // НАСТРОЙКА КАМЕРЫ
    // ================================================================

    Camera camera;

    // Позиция камеры ВНУТРИ Sponza (смещена от центра для лучшей композиции)
    vec3 cameraPosition;
    cameraPosition.x = 3.0;   // Смещение вправо от центрального прохода
    cameraPosition.y = 2.2;   // Высота примерно на уровне человеческого роста
    cameraPosition.z = 5.0;   // Расположена недалеко от входа, смотря вглубь

    // Смотрим вглубь атриума, слегка приподняв цель
    vec3 lookAtPoint;
    lookAtPoint.x = 0.0;      // Центр сцены по X
    lookAtPoint.y = 3.5;      // Цель выше камеры - взгляд слегка вверх на своды
    lookAtPoint.z = -15.0;    // Точка в глубине сцены

    camera.setEye(cameraPosition);
    camera.setTarget(lookAtPoint);
    camera.setUp(vec3{ 0, 1, 0 });
    camera.setFOV(70.0);      // Достаточно широкий угол для ощущения простора
    camera.setAspect((double)WIDTH / HEIGHT);
    // Дальнюю плоскость отсечения нужно увеличить для большой сцены
    camera.setClipping(0.5, 500.0);

    // ================================================================
    // ИНИЦИАЛИЗАЦИЯ РЕНДЕРИНГА
    // ================================================================

    TGAImage framebuffer(WIDTH, HEIGHT, TGAImage::RGB);
    init_zbuffer(WIDTH, HEIGHT);

    // Установка глобальных матриц
    ModelView = camera.getViewMatrix();
    Perspective = camera.getProjectionMatrix();
    init_viewport(0, 0, WIDTH, HEIGHT);

    // Направления источников света
    vec3 key_light_dir = normalized(make_vec3(1.0, 1.4, 1.0));
    vec3 fill_light_dir = normalized(make_vec3(-0.3, 0.5, 0.2));
    vec3 rim_light_dir = normalized(make_vec3(-1.0, 0.8, -1.5));

    // ================================================================
    // FRUSTUM CULLING
    // ================================================================

    mat<4, 4> viewProjection = Perspective * ModelView;
    Frustum frustum = Frustum::createFromMatrix(viewProjection);

    int models_culled = 0;
    int models_rendered = 0;
    int total_triangles = 0;
    int culled_triangles = 0;

    // ================================================================
    // РЕНДЕРИНГ SPONZA
    // ================================================================

    bool sponzaVisible = frustum.intersects(sponzaWorldAABB);
    if (sponzaVisible) {
        models_rendered++;
        std::cout << "\nRendering sponza model (VISIBLE)" << std::endl;

        mat<4, 4> originalModelView = ModelView;
        ModelView = ModelView * sponzaModelMatrix;

        PhongShader sponza_shader(sponza_model.get());
        sponza_shader.initLightDirections(key_light_dir, fill_light_dir, rim_light_dir);
        sponza_shader.normal_map_strength = 0.5;

        total_triangles += sponza_model->nfaces();
        for (int face = 0; face < sponza_model->nfaces(); ++face) {
            vec4 clip_space_triangle[3];
            for (int vertex = 0; vertex < 3; ++vertex) {
                clip_space_triangle[vertex] = sponza_shader.vertex(face, vertex);
            }
            rasterize(clip_space_triangle, sponza_shader, framebuffer);
        }

        ModelView = originalModelView;
    }
    else {
        models_culled++;
        culled_triangles += sponza_model->nfaces();
        std::cout << "\nSponza model CULLED by frustum" << std::endl;
    }

    // ================================================================
    // РЕНДЕРИНГ ГОЛОВЫ
    // ================================================================

    bool headVisible = frustum.intersects(headWorldAABB);
    if (headVisible) {
        models_rendered++;
        std::cout << "\nRendering head model (VISIBLE)" << std::endl;

        mat<4, 4> originalModelView = ModelView;
        ModelView = ModelView * headModelMatrix;

        PhongShader head_shader(head_model.get());
        head_shader.initLightDirections(key_light_dir, fill_light_dir, rim_light_dir);

        total_triangles += head_model->nfaces();
        for (int face = 0; face < head_model->nfaces(); ++face) {
            vec4 clip_space_triangle[3];
            for (int vertex = 0; vertex < 3; ++vertex) {
                clip_space_triangle[vertex] = head_shader.vertex(face, vertex);
            }
            rasterize(clip_space_triangle, head_shader, framebuffer);
        }

        std::vector<double> zbuffer_before_eyes = zbuffer;

        // ================================================================
        // РЕНДЕРИНГ ГЛАЗ
        // ================================================================

        bool eyesVisible = frustum.intersects(eyeWorldAABB);
        if (eyesVisible) {
            models_rendered++;
            std::cout << "\nRendering eye model (VISIBLE)" << std::endl;

            EyeShader eye_shader(eye_model.get());
            eye_shader.initLightDirections(key_light_dir, rim_light_dir);

            total_triangles += eye_model->nfaces();
            for (int face = 0; face < eye_model->nfaces(); ++face) {
                vec4 clip_space_triangle[3];
                for (int vertex = 0; vertex < 3; ++vertex) {
                    clip_space_triangle[vertex] = eye_shader.vertex(face, vertex);
                }
                rasterize(clip_space_triangle, eye_shader, framebuffer);
            }
        }
        else {
            models_culled++;
            culled_triangles += eye_model->nfaces();
            std::cout << "\nEye model CULLED by frustum" << std::endl;
        }

        ModelView = originalModelView;
        zbuffer = zbuffer_before_eyes;  // Для SSAO используем Z-буфер без глаз
    }
    else {
        models_culled++;
        culled_triangles += head_model->nfaces();
        std::cout << "\nHead model CULLED by frustum" << std::endl;
    }

    // ================================================================
    // СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
    // ================================================================

    if (models_rendered > 0) {
        if (framebuffer.write_tga_file("phong.tga")) {
            std::cout << "\nSaved: phong.tga" << std::endl;
        }
    }
    else {
        std::cout << "\nNo models rendered, skipping phong.tga" << std::endl;
    }

    save_zbuffer_image(zbuffer, WIDTH, HEIGHT, "zbuffer.tga");
    std::cout << "Saved: zbuffer.tga" << std::endl;

    // Вычисление SSAO
    std::cout << "Computing SSAO..." << std::endl;
    TGAImage ao_map(WIDTH, HEIGHT, TGAImage::RGB);
    for (int y = 0; y < HEIGHT; ++y) {
        for (int x = 0; x < WIDTH; ++x) {
            double ao_value = compute_ssao_at(zbuffer, WIDTH, HEIGHT, x, y);
            unsigned char intensity = (unsigned char)(255.0 * ao_value);
            ao_map.set(x, y, TGAColor(intensity, intensity, intensity));
        }
    }
    ao_map.write_tga_file("ao.tga");
    std::cout << "Saved: ao.tga" << std::endl;

    // Финальная композиция
    if (models_rendered > 0) {
        std::cout << "Compositing final image..." << std::endl;
        TGAImage final_result(WIDTH, HEIGHT, TGAImage::RGB);
        for (int y = 0; y < HEIGHT; ++y) {
            for (int x = 0; x < WIDTH; ++x) {
                TGAColor phong_color = framebuffer.get(x, y);
                TGAColor ao_color = ao_map.get(x, y);
                double ao_factor = ao_color[0] / 255.0;

                final_result.set(x, y, TGAColor(
                    (unsigned char)std::min(255.0, (double)phong_color[0] * ao_factor),
                    (unsigned char)std::min(255.0, (double)phong_color[1] * ao_factor),
                    (unsigned char)std::min(255.0, (double)phong_color[2] * ao_factor)
                ));
            }
        }
        final_result.write_tga_file("final.tga");
        std::cout << "Saved: final.tga" << std::endl;
    }

    // ================================================================
    // СТАТИСТИКА
    // ================================================================

    print_render_stats();

    std::cout << "\n=== Frustum Culling Statistics ===" << std::endl;
    std::cout << "  Total models: " << (models_rendered + models_culled) << std::endl;
    std::cout << "  Models rendered: " << models_rendered << std::endl;
    std::cout << "  Models culled: " << models_culled << std::endl;
    std::cout << "  Total triangles: " << total_triangles << std::endl;
    std::cout << "  Culled triangles: " << culled_triangles << std::endl;

    if ((total_triangles + culled_triangles) > 0) {
        std::cout << "  Triangle culling efficiency: "
            << (culled_triangles * 100.0 / (total_triangles + culled_triangles)) << "%" << std::endl;
    }

    std::cout << "\nCamera Information:" << std::endl;
    camera.printInfo();

    std::cout << "\nModel Positions:" << std::endl;
    std::cout << "  Sponza: scale 0.05, position (0, -5, 0)" << std::endl;
    std::cout << "  Head: position (15, 5, 0)" << std::endl;

    modelManager.printStats();

    std::cout << "\n=== Render Complete ===" << std::endl;
    std::cout << "Files saved:" << std::endl;
    std::cout << "  - zbuffer.tga (without eyes)" << std::endl;
    std::cout << "  - ao.tga (ambient occlusion)" << std::endl;

    if (models_rendered > 0) {
        std::cout << "  - phong.tga (with eyes)" << std::endl;
        std::cout << "  - final.tga (phong + ao)" << std::endl;
    }

    return 0;
}