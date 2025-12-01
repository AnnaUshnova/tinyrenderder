#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <iostream>

#pragma pack(push,1)
struct TGAHeader {
    std::uint8_t  idlength = 0;
    std::uint8_t  colormaptype = 0;
    std::uint8_t  datatypecode = 2; // default = uncompressed true-color image
    std::uint16_t colormaporigin = 0;
    std::uint16_t colormaplength = 0;
    std::uint8_t  colormapdepth = 0;
    std::uint16_t x_origin = 0;
    std::uint16_t y_origin = 0;
    std::uint16_t width = 0;
    std::uint16_t height = 0;
    std::uint8_t  bitsperpixel = 24;
    std::uint8_t  imagedescriptor = 0;
};
#pragma pack(pop)

// -------------------------------------------------------------
// TGAColor Ч формат BGRA
struct TGAColor {
    std::uint8_t bgra[4];
    std::uint8_t bytespp;

    TGAColor() : bgra{ 0,0,0,255 }, bytespp(4) {}

    TGAColor(std::uint8_t R, std::uint8_t G, std::uint8_t B, std::uint8_t A = 255) {
        bgra[0] = B; bgra[1] = G; bgra[2] = R; bgra[3] = A;
        bytespp = 4;
    }

    // grayscale
    TGAColor(std::uint8_t v) {
        bgra[0] = v; bgra[1] = v; bgra[2] = v; bgra[3] = 255;
        bytespp = 1;
    }

    TGAColor(const std::uint8_t* p, std::uint8_t bpp) {
        bytespp = bpp;
        for (int i = 0; i < (int)bpp; ++i) bgra[i] = p[i];
        for (int i = bpp; i < 4; ++i) bgra[i] = 0;
    }

    std::uint8_t& operator[](int i) { return bgra[i]; }
    const std::uint8_t& operator[](int i) const { return bgra[i]; }

    TGAColor operator*(float intensity) const {
        TGAColor res = *this;
        if (intensity < 0.f) intensity = 0.f;
        if (intensity > 1.f) intensity = 1.f;
        for (int i = 0; i < 4; i++)
            res.bgra[i] = std::uint8_t(bgra[i] * intensity);
        return res;
    }
};

// -------------------------------------------------------------
// TGAImage (расширенна€ верси€ с RLE, scale, blur)
class TGAImage {
public:
    enum Format { GRAYSCALE = 1, RGB = 3, RGBA = 4 };

    TGAImage();
    TGAImage(const int width, const int height, const int bytespp, TGAColor clear = TGAColor());

    bool  read_tga_file(const std::string filename);
    bool write_tga_file(const std::string filename,
        const bool vflip = true,
        const bool rle = true) const;

    void flip_horizontally();
    void flip_vertically();

    TGAColor get(const int x, const int y) const;
    void set(const int x, const int y, const TGAColor& c);

    int width()  const;
    int height() const;

    // old TinyRenderer API kept
    bool scale(int new_w, int new_h);
    void gaussian_blur(const int radius);

    // direct access when needed
    std::uint8_t* buffer() { return data.empty() ? nullptr : data.data(); }

private:
    bool   load_rle_data(std::ifstream& in);
    bool unload_rle_data(std::ofstream& out) const;

private:
    int w = 0;
    int h = 0;
    std::uint8_t bpp = 0;
    std::vector<std::uint8_t> data;
};
