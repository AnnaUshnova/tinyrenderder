// ================================================================
// camera.h - Камера и матрицы вида/проекции
// ================================================================
#pragma once
#include "geometry.h"

class Camera {
public:
    vec3 eye;       // Позиция камеры в мировых координатах
    vec3 center;    // Точка, на которую смотрит камера
    vec3 up;        // Вектор "вверх" в мировых координатах

    mat<4, 4> view; // Матрица вида (world -> camera space)
    mat<4, 4> proj; // Матрица проекции (camera -> clip space)

    Camera(const vec3& eye_pos, const vec3& target, const vec3& up_dir)
        : eye(eye_pos), center(target), up(up_dir) {
        updateView();
    }

    // Обновление матрицы вида
    void updateView() {
        vec3 z_axis = normalized(eye - center);  // Направление взгляда (от цели к камере)
        vec3 x_axis = normalized(cross(up, z_axis)); // Вектор "вправо"
        vec3 y_axis = cross(z_axis, x_axis);     // Вектор "вверх" в пространстве камеры

        view = mat<4, 4>::identity();

        // Заполняем матрицу поворота
        view[0][0] = x_axis[0]; view[0][1] = x_axis[1]; view[0][2] = x_axis[2];
        view[1][0] = y_axis[0]; view[1][1] = y_axis[1]; view[1][2] = y_axis[2];
        view[2][0] = z_axis[0]; view[2][1] = z_axis[1]; view[2][2] = z_axis[2];

        // Заполняем компоненты трансляции
        view[0][3] = -dot(x_axis, eye);
        view[1][3] = -dot(y_axis, eye);
        view[2][3] = -dot(z_axis, eye);
    }

    // Установка перспективной проекции
    void setPerspective(double fov_degrees, double aspect_ratio,
        double z_near, double z_far) {
        double fov_radians = fov_degrees * M_PI / 180.0;
        double tan_half_fov = std::tan(fov_radians / 2.0);

        proj = mat<4, 4>::identity();
        proj[0][0] = 1.0 / (aspect_ratio * tan_half_fov);
        proj[1][1] = 1.0 / tan_half_fov;
        proj[2][2] = (z_far + z_near) / (z_near - z_far);
        proj[2][3] = (2.0 * z_far * z_near) / (z_near - z_far);
        proj[3][2] = -1.0;  // Перенос z в w-компоненту для перспективного деления
        proj[3][3] = 0.0;
    }
};