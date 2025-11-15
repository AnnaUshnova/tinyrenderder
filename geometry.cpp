#include "geometry.h"

// Correct specializations for type conversions
template <>
template <>
Vec3<int>::Vec3(const Vec3<float>& v) :
    x(static_cast<int>(v.x + 0.5f)),
    y(static_cast<int>(v.y + 0.5f)),
    z(static_cast<int>(v.z + 0.5f)) {
}

template <>
template <>
Vec3<float>::Vec3(const Vec3<int>& v) :
    x(static_cast<float>(v.x)),
    y(static_cast<float>(v.y)),
    z(static_cast<float>(v.z)) {
}