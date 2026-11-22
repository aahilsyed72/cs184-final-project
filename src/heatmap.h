#pragma once
// Step-count heatmap: maps per-pixel RK4 step counts through the 'inferno'
// colormap on a log scale, then appends a legend strip on the right with a
// color bar, tick marks and labels (drawn with a tiny built-in 5x7 bitmap font,
// so no font library is needed).

#include "colormap.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace heatmap {

struct RGB8 { uint8_t r, g, b; };

// Linear interpolation into the sampled inferno table, t in [0,1].
inline RGB8 inferno(double t) {
    t = std::clamp(t, 0.0, 1.0) * (INFERNO_N - 1);
    int i = std::min(int(t), INFERNO_N - 2);
    float f = float(t - i);
    auto ch = [&](int c) {
        float v = INFERNO[i][c] * (1 - f) + INFERNO[i + 1][c] * f;
        return (uint8_t)std::lround(v * 255.0f);
    };
    return {ch(0), ch(1), ch(2)};
}

// ---- 5x7 bitmap font: each glyph is 7 rows, low 5 bits of each byte = pixels.
inline const uint8_t* glyph(char c) {
    static const uint8_t blank[7] = {0};
    struct G { char c; uint8_t rows[7]; };
    static const G font[] = {
        {'0', {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}}, {'1', {0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}},
        {'2', {0x0E,0x11,0x01,0x02,0x04,0x08,0x1F}}, {'3', {0x1F,0x02,0x04,0x02,0x01,0x11,0x0E}},
        {'4', {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}}, {'5', {0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E}},
        {'6', {0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}}, {'7', {0x1F,0x01,0x02,0x04,0x08,0x08,0x08}},
        {'8', {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}}, {'9', {0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}},
        {'A', {0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}}, {'C', {0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}},
        {'D', {0x1C,0x12,0x11,0x11,0x11,0x12,0x1C}}, {'E', {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}},
        {'F', {0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}}, {'G', {0x0E,0x11,0x10,0x17,0x11,0x11,0x0F}},
        {'I', {0x0E,0x04,0x04,0x04,0x04,0x04,0x0E}}, {'K', {0x11,0x12,0x14,0x18,0x14,0x12,0x11}},
        {'L', {0x10,0x10,0x10,0x10,0x10,0x10,0x1F}}, {'M', {0x11,0x1B,0x15,0x15,0x11,0x11,0x11}},
        {'N', {0x11,0x19,0x15,0x13,0x11,0x11,0x11}}, {'O', {0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}},
        {'P', {0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}}, {'R', {0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}},
        {'S', {0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E}}, {'T', {0x1F,0x04,0x04,0x04,0x04,0x04,0x04}},
        {'V', {0x11,0x11,0x11,0x11,0x11,0x0A,0x04}}, {'X', {0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}},
        {'(', {0x02,0x04,0x08,0x08,0x08,0x04,0x02}}, {')', {0x08,0x04,0x02,0x02,0x02,0x04,0x08}},
        {'/', {0x01,0x02,0x02,0x04,0x08,0x08,0x10}}, {'.', {0x00,0x00,0x00,0x00,0x00,0x0C,0x0C}},
        {'=', {0x00,0x00,0x1F,0x00,0x1F,0x00,0x00}}, {'-', {0x00,0x00,0x00,0x1F,0x00,0x00,0x00}},
    };
    if (c >= 'a' && c <= 'z') c = char(c - 'a' + 'A');
    for (const G& g : font) if (g.c == c) return g.rows;
    return blank;
}

struct Canvas {
    int w, h;
    std::vector<RGB8> px;
    Canvas(int w_, int h_, RGB8 bg) : w(w_), h(h_), px(size_t(w_) * h_, bg) {}
    void set(int x, int y, RGB8 c) { if (x >= 0 && y >= 0 && x < w && y < h) px[size_t(y) * w + x] = c; }
    void rect(int x0, int y0, int x1, int y1, RGB8 c) {
        for (int y = y0; y < y1; ++y) for (int x = x0; x < x1; ++x) set(x, y, c);
    }
    // Draw text with its top-left at (x, y); each font pixel becomes s x s.
    void text(int x, int y, const std::string& s, int scale, RGB8 c) {
        for (char ch : s) {
            const uint8_t* g = glyph(ch);
            for (int row = 0; row < 7; ++row)
                for (int col = 0; col < 5; ++col)
                    if (g[row] & (0x10 >> col)) rect(x + col * scale, y + row * scale,
                                                      x + (col + 1) * scale, y + (row + 1) * scale, c);
            x += 6 * scale;
        }
    }
    int text_width(const std::string& s, int scale) const { return int(s.size()) * 6 * scale - scale; }
};

// steps: W*H per-pixel values (average RK4 steps per sample).
// lo/hi: color range; if lo <= 0, the range is taken from the data.
// title: printed at the bottom of the legend (e.g. "ADAPTIVE").
// Returns RGB8 pixels for an image of size (W + legend_w) x H.
inline std::vector<uint8_t> render(const std::vector<float>& steps, int W, int H, double lo,
                                   double hi, const std::string& title, int& out_w) {
    if (lo <= 0 || hi <= lo) {
        auto [mn, mx] = std::minmax_element(steps.begin(), steps.end());
        lo = std::max(1.0, double(*mn));
        hi = std::max(lo * 10.0, double(*mx));
    }
    const double llo = std::log10(lo), lhi = std::log10(hi);
    auto t_of = [&](double v) { return (std::log10(std::max(v, 1e-9)) - llo) / (lhi - llo); };

    const int scale = H >= 600 ? 2 : 1;  // font size
    const int legend_w = 110 * scale;
    out_w = W + legend_w;
    Canvas cv(out_w, H, {24, 24, 28});

    // Heatmap itself
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) cv.set(x, y, inferno(t_of(steps[size_t(y) * W + x])));

    // Legend: title, vertical color bar (high at top), ticks at 1-2-5 per decade.
    const RGB8 white{235, 235, 235}, grey{150, 150, 150};
    int tx = W + 10 * scale;
    cv.text(tx, 10 * scale, "STEPS/PX", scale, white);
    cv.text(tx, 21 * scale, "(LOG)", scale, grey);
    const int bar_x0 = W + 12 * scale, bar_x1 = bar_x0 + 16 * scale;
    const int bar_y0 = 40 * scale, bar_y1 = H - 30 * scale;
    for (int y = bar_y0; y < bar_y1; ++y) {
        RGB8 c = inferno(1.0 - double(y - bar_y0) / (bar_y1 - bar_y0 - 1));
        for (int x = bar_x0; x < bar_x1; ++x) cv.set(x, y, c);
    }
    auto y_of = [&](double v) { return int(std::lround(bar_y1 - 1 - t_of(v) * (bar_y1 - bar_y0 - 1))); };
    int last_label_y = 1 << 30;
    for (int k = int(std::floor(llo)); k <= int(std::ceil(lhi)); ++k)
        for (int m : {1, 2, 5}) {
            double v = m * std::pow(10.0, k);
            if (v < lo * 0.999 || v > hi * 1.001) continue;
            int y = y_of(v);
            cv.rect(bar_x1, y, bar_x1 + 5 * scale, y + scale, white);
            if (last_label_y - y >= 11 * scale) {  // skip labels that would overlap
                char buf[32];
                std::snprintf(buf, sizeof buf, "%.0f", v);
                cv.text(bar_x1 + 8 * scale, y - 3 * scale, buf, scale, white);
                last_label_y = y;
            }
        }
    cv.text(tx, H - 20 * scale, title, scale, grey);

    std::vector<uint8_t> out(size_t(out_w) * H * 3);
    for (size_t i = 0; i < cv.px.size(); ++i) {
        out[3 * i] = cv.px[i].r; out[3 * i + 1] = cv.px[i].g; out[3 * i + 2] = cv.px[i].b;
    }
    return out;
}

}  // namespace heatmap
