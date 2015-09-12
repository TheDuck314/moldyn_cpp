#ifndef VECTOR2_H
#define VECTOR2_H

#include <cmath>

struct Vector2
{
    double x;
    double y;

    void Zero() { x = 0; y = 0; }
};

inline Vector2 operator+(Vector2 a, Vector2 b)
{
    Vector2 ret{ a.x + b.x, a.y + b.y };
    return ret;
}

inline Vector2 operator-(Vector2 a, Vector2 b)
{
    Vector2 ret{ a.x - b.x, a.y - b.y };
    return ret;
}

inline void operator+=(Vector2 &v, Vector2 rhs)
{
    v = v + rhs;
}

inline void operator-=(Vector2 &v, Vector2 rhs)
{
    v = v - rhs;
}

inline Vector2 operator*(double c, Vector2 v)
{
    Vector2 ret{ c * v.x, c * v.y };
    return ret;
}

inline double norm2(Vector2 v)
{
    return v.x * v.x + v.y * v.y;
}

inline double norm(Vector2 v)
{
    return std::sqrt(norm2(v));
}

inline double dist(Vector2 a, Vector2 b)
{
    return norm(a - b);
}

#endif