// camera.h
#pragma once
#include "geometry.h"

class Camera {
public:
    vec3 eye;
    vec3 center;
    vec3 up;

    mat<4, 4> view;
    mat<4, 4> proj;

    Camera(const vec3& e, const vec3& c, const vec3& u)
        : eye(e), center(c), up(u)
    {
        updateView();
    }

    void updateView() {
        vec3 z = normalized(eye - center);
        vec3 x = normalized(cross(up, z));
        vec3 y = cross(z, x);

        view = mat<4, 4>::identity();
        view[0][0] = x[0]; view[0][1] = x[1]; view[0][2] = x[2];
        view[1][0] = y[0]; view[1][1] = y[1]; view[1][2] = y[2];
        view[2][0] = z[0]; view[2][1] = z[1]; view[2][2] = z[2];

        view[0][3] = -dot(x, eye);
        view[1][3] = -dot(y, eye);
        view[2][3] = -dot(z, eye);
    }

    void setPerspective(double fov_deg, double aspect, double znear, double zfar) {
        double fov = fov_deg * M_PI / 180.0;
        double t = std::tan(fov / 2.0);
        proj = mat<4, 4>::identity();

        proj[0][0] = 1.0 / (aspect * t);
        proj[1][1] = 1.0 / t;
        proj[2][2] = (zfar + znear) / (znear - zfar);
        proj[2][3] = (2.0 * zfar * znear) / (znear - zfar);
        proj[3][2] = -1.0;
        proj[3][3] = 0.0;
    }
};
