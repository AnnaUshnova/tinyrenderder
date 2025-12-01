#include "tgaimage.h"
#include <iostream>
#include <cstring>
#include <cmath>

TGAImage::TGAImage() : w(0), h(0), bpp(0) {}

TGAImage::TGAImage(const int width, const int height, const int bytespp, TGAColor clear)
    : w(width), h(height), bpp(bytespp) {
    data.resize(w * h * bpp, 0);
    // Fill with initial color
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            set(x, y, clear);
        }
    }
}

int TGAImage::width() const { return w; }
int TGAImage::height() const { return h; }

// ------------------- GET & SET -------------------

TGAColor TGAImage::get(const int x, const int y) const {
    if (!data.size() || x < 0 || y < 0 || x >= w || y >= h) {
        return TGAColor();
    }
    int idx = (x + y * w) * bpp;
    return TGAColor(&data[idx], bpp);
}

void TGAImage::set(const int x, const int y, const TGAColor& c) {
    if (!data.size() || x < 0 || y < 0 || x >= w || y >= h) return;

    int idx = (x + y * w) * bpp;
    for (int i = 0; i < bpp; i++) {
        data[idx + i] = c.bgra[i];
    }
}

// ------------------- FLIP -------------------

void TGAImage::flip_horizontally() {
    if (!data.size()) return;
    int half = w / 2;
    int bytes_per_line = w * bpp;

    for (int y = 0; y < h; y++) {
        int line = y * bytes_per_line;
        for (int x = 0; x < half; x++) {
            for (int i = 0; i < bpp; i++) {
                std::swap(data[line + x * bpp + i],
                    data[line + (w - 1 - x) * bpp + i]);
            }
        }
    }
}

void TGAImage::flip_vertically() {
    if (!data.size()) return;

    int bytes_per_line = w * bpp;
    std::vector<uint8_t> buf(bytes_per_line);

    for (int y = 0; y < h / 2; y++) {
        int line1 = y * bytes_per_line;
        int line2 = (h - 1 - y) * bytes_per_line;
        std::memcpy(buf.data(), &data[line1], bytes_per_line);
        std::memcpy(&data[line1], &data[line2], bytes_per_line);
        std::memcpy(&data[line2], buf.data(), bytes_per_line);
    }
}

// ------------------- READ TGA -------------------

bool TGAImage::read_tga_file(const std::string filename) {
    data.clear();

    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "can't open file " << filename << "\n";
        return false;
    }

    TGAHeader header;
    in.read((char*)&header, sizeof(header));
    if (!in.good()) {
        std::cerr << "can't read header\n";
        return false;
    }

    w = header.width;
    h = header.height;
    bpp = header.bitsperpixel >> 3;

    if (w <= 0 || h <= 0 || (bpp != 1 && bpp != 3 && bpp != 4)) {
        std::cerr << "invalid TGA format\n";
        return false;
    }

    data.resize(w * h * bpp);

    in.seekg(header.idlength, std::ios::cur);

    if (header.datatypecode == 2 || header.datatypecode == 3) {
        // uncompressed
        in.read((char*)data.data(), w * h * bpp);
    }
    else if (header.datatypecode == 10 || header.datatypecode == 11) {
        // RLE
        if (!load_rle_data(in)) return false;
    }
    else {
        std::cerr << "unknown TGA type\n";
        return false;
    }

    if (!(header.imagedescriptor & 0x20)) flip_vertically();
    if (header.imagedescriptor & 0x10)    flip_horizontally();

    return true;
}

bool TGAImage::load_rle_data(std::ifstream& in) {
    int pixelcount = w * h;
    int currentpixel = 0;
    int bytesread = 0;
    TGAColor c;

    while (currentpixel < pixelcount) {
        unsigned char chunkheader = 0;
        in.read((char*)&chunkheader, 1);

        if (chunkheader < 128) {
            int count = chunkheader + 1;
            for (int i = 0; i < count; i++) {
                in.read((char*)c.bgra, bpp);
                for (int t = 0; t < bpp; t++)
                    data[(currentpixel * bpp) + t] = c.bgra[t];
                currentpixel++;
                if (currentpixel > pixelcount) return false;
            }
        }
        else {
            int count = chunkheader - 127;
            in.read((char*)c.bgra, bpp);
            for (int i = 0; i < count; i++) {
                for (int t = 0; t < bpp; t++)
                    data[(currentpixel * bpp) + t] = c.bgra[t];
                currentpixel++;
                if (currentpixel > pixelcount) return false;
            }
        }
    }

    return true;
}

// ------------------- WRITE TGA -------------------

bool TGAImage::write_tga_file(const std::string filename, const bool vflip, const bool rle) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "can't open " << filename << "\n";
        return false;
    }

    TGAHeader header = {};
    header.bitsperpixel = bpp * 8;
    header.width = w;
    header.height = h;
    // set datatypecode according to whether we write RLE or RAW data
    // 2 = uncompressed true-color, 3 = uncompressed grayscale
    // 10 = RLE true-color,     11 = RLE grayscale
    header.datatypecode = (bpp == 1 ? (rle ? 11 : 3) : (rle ? 10 : 2));
    header.imagedescriptor = vflip ? 0x00 : 0x20;

    out.write((char*)&header, sizeof(header));

    if (rle) {
        if (!unload_rle_data(out)) {
            std::cerr << "RLE write fail\n";
            return false;
        }
    }
    else {
        out.write((char*)data.data(), w * h * bpp);
    }

    return true;
}

bool TGAImage::unload_rle_data(std::ofstream& out) const {
    const unsigned char max_chunk_length = 128;
    int npixels = w * h;
    int cur = 0;

    while (cur < npixels) {
        int chunkstart = cur * bpp;
        int run_length = 1;
        bool raw = true;

        while (cur + run_length < npixels && run_length < max_chunk_length) {
            bool equal = true;
            for (int i = 0; i < bpp; i++)
                if (data[(cur + run_length) * bpp + i] != data[chunkstart + i]) {
                    equal = false;
                    break;
                }
            if (!equal) break;
            run_length++;
            raw = false;
        }

        if (!raw) {
            // RLE chunk
            out.put(run_length - 1 + 128);
            out.write((char*)&data[chunkstart], bpp);
            cur += run_length;
        }
        else {
            // Raw chunk
            int rawstart = cur;
            run_length = 1;
            while (rawstart + run_length < npixels && run_length < max_chunk_length) {
                bool next_equal = true;
                for (int i = 0; i < bpp; i++)
                    if (data[(rawstart + run_length) * bpp + i] != data[(rawstart + run_length - 1) * bpp + i]) {
                        next_equal = false;
                        break;
                    }
                if (next_equal) break;
                run_length++;
            }
            out.put(run_length - 1);
            out.write((char*)&data[rawstart * bpp], run_length * bpp);
            cur += run_length;
        }
    }

    return true;
}

// ------------------- SCALE (из старой версии) -------------------

bool TGAImage::scale(int w2, int h2) {
    if (w2 <= 0 || h2 <= 0 || !data.size()) return false;

    std::vector<uint8_t> tdata(w2 * h2 * bpp);

    for (int y = 0; y < h2; y++) {
        for (int x = 0; x < w2; x++) {
            int srcx = x * w / w2;
            int srcy = y * h / h2;
            TGAColor c = get(srcx, srcy);

            int idx = (x + y * w2) * bpp;
            for (int i = 0; i < bpp; i++)
                tdata[idx + i] = c.bgra[i];
        }
    }

    w = w2;
    h = h2;
    data.swap(tdata);
    return true;
}

// ------------------- GAUSSIAN BLUR (из старой версии) -------------------

void TGAImage::gaussian_blur(const int radius) {
    if (radius <= 0 || !data.size()) return;

    // simple separable kernel
    std::vector<float> kernel(radius * 2 + 1);
    float sigma = radius / 2.0f;
    float sum = 0;

    for (int i = -radius; i <= radius; i++) {
        float v = std::exp(-(i * i) / (2 * sigma * sigma));
        kernel[i + radius] = v;
        sum += v;
    }
    for (auto& v : kernel) v /= sum;

    // TEMP buffer
    std::vector<uint8_t> tmp = data;

    // horizontal pass
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            float accum[4] = { 0,0,0,0 };
            for (int k = -radius; k <= radius; k++) {
                int xx = std::clamp(x + k, 0, w - 1);
                TGAColor c(&tmp[(xx + y * w) * bpp], bpp);
                float kv = kernel[k + radius];
                for (int i = 0; i < bpp; i++)
                    accum[i] += c.bgra[i] * kv;
            }
            int idx = (x + y * w) * bpp;
            for (int i = 0; i < bpp; i++)
                data[idx + i] = (uint8_t)accum[i];
        }
    }

    tmp = data;

    // vertical pass
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            float accum[4] = { 0,0,0,0 };
            for (int k = -radius; k <= radius; k++) {
                int yy = std::clamp(y + k, 0, h - 1);
                TGAColor c(&tmp[(x + yy * w) * bpp], bpp);
                float kv = kernel[k + radius];
                for (int i = 0; i < bpp; i++)
                    accum[i] += c.bgra[i] * kv;
            }
            int idx = (x + y * w) * bpp;
            for (int i = 0; i < bpp; i++)
                data[idx + i] = (uint8_t)accum[i];
        }
    }
}
