#pragma once

#ifndef CLAMP
#define CLAMP(a, b, c) (MAX((a), MIN((b), (c))))
#endif

#ifndef ABS
#define ABS(a) ((((a) >= 0)) ? (a) : -(a))
#endif

#ifndef MAX
#define MAX(a, b) (((a) >= (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) >= (b)) ? (b) : (a))
#endif

#define IX2(i, j, nx) ((j) * (nx) + (i))

#define IX3(i, j, k, nx, ny) ((k) * (nx) * (ny) + (j) * (nx) + (i))

//#define THREE_D_TO_ONE_D(x, y, z, xMax, yMax)(((z * xMax * yMax) + (y * xMax) + x))
#define THREE_D_TO_ONE_D(b, x, y, width, height) ((b) * (width) * (height) + (y) * (width) + (x))

//radius theta phi in polar coord to x,y,z in cartesian coord
#define RTP2XYZ(r, t, p, x, y, z) \
    (z) = sin(t);                 \
    (x) = (r) * (z) * cos(p);     \
    (y) = (r) * (z) * sin(p);     \
    (z) = (r) * cos(t);

//x,y,z to radius theta phi
#define XYZ2RTP(x, y, z, r, t, p)                             \
    (r) = (x) * (x) + (y) * (y);                              \
    if ((r) < 0.00001)                                        \
    {                                                         \
        (r) = (float)sqrt((r) + (z) * (z));                   \
        (p) = 0.0;                                            \
        (t) = 0.0;                                            \
    }                                                         \
    else                                                      \
    {                                                         \
        (p) = (float)acos(CLAMP(-1.0, (x) / sqrt((r)), 1.0)); \
        if ((y) < 0.0)                                        \
            (p) = (float)(2.0 * M_PI) - (p);                  \
        (r) = (float)sqrt((r) + (z) * (z));                   \
        (t) = (float)acos(CLAMP(-1.0, (z) / (r), 1.0));       \
    }
