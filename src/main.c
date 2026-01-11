#include <notcurses/notcurses.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <locale.h>

int main([[maybe_unused]] int argc, char** argv)
{
    if (!setlocale(LC_ALL, "")) {
        fprintf(stderr, "%s: setlocale(LC_ALL, \"\") failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    struct notcurses_options opts = {};
    struct notcurses* nc = notcurses_core_init(&opts, stdout);
    if (!nc) {
        fprintf(stderr, "%s: notcurses_init() failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    uint32_t rows = 0;
    uint32_t cols = 0;
    struct ncplane* nstd = notcurses_stddim_yx(nc, &rows, &cols);
    if (!nstd) {
        fprintf(stderr, "%s: notcurses_stdplane() failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    cols *= 2;
    rows *= 2;

    size_t pixelBufferSize = rows * cols * sizeof(uint32_t);
    uint32_t* rgb = malloc(pixelBufferSize);
    if (!rgb) {
        fprintf(stderr, "%s: malloc() failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    for (uint32_t row = 0; row < rows; ++row) {
        for (uint32_t col = 0; col < cols; ++col) {
            float red = (float)col / (float)(cols - 1);
            float green = (float)row / (float)(rows - 1);
            float blue = 255.0f;

            uint32_t ired = (uint32_t)(red * 255.99f);
            uint32_t igreen = (uint32_t)(green * 255.99f);
            uint32_t iblue = (uint32_t)(blue * 255.99f);

            rgb[row * cols + col] = 0xFF000000 |
                                    (iblue << 16) |
                                    (igreen << 8) |
                                    ired;
        }
    }

    struct ncvisual_options vopts = {};
    vopts.n = nstd;
    vopts.leny = rows;
    vopts.lenx = cols;
    vopts.blitter = NCBLIT_2x2;
    vopts.scaling = NCSCALE_NONE;

    if (-1 == ncblit_rgba(rgb, (int32_t)(cols * sizeof(uint32_t)), &vopts)) {
        fprintf(stderr, "%s: ncblit_rgba() failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (-1 == notcurses_render(nc)) {
        fprintf(stderr, "%s: notcurses_render() failed\n", argv[0]);
        return EXIT_FAILURE;
    }

    for (;;) {
        struct ncinput input = {};
        notcurses_get(nc, NULL, &input);

        if (input.id == 'q' && input.evtype == NCTYPE_RELEASE) {
            break;
        }
    }
    free(rgb);
    return notcurses_stop(nc);
}
