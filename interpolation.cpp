#include <utility>
#include <stdexcept>
#include "interpolation.h"

constexpr double C_PI = 3.14159265358979323846;

static inline bool fuzzyIsNull(double d)
{
    return std::abs(d) < 0.000000000001;
}

static inline bool fuzzyIsNull(float f)
{
    return std::abs(f) <= 0.00001f;
}

Interpolation::Vector<Interpolation::data_type> SineInterp::interp(const Vector<data_type> &data, unsigned int interpNum) noexcept(false)
{
    // 理论上正弦插值每个周期至少需要三个点才能重现波形
    if(data.size() < 3)
    {
        throw std::invalid_argument("[SineInterp::interp] data.size() < 3;");
    }

    const int sz = static_cast<int>(data.size());
    int M = int(interpNum) + 1;

    Vector<data_type> ret((sz - 1) * interpNum + sz);
    for(int i = 0; i < sz; ++i)
    {
        ret[M * i] = data[i];
    }

    for(int k = 0; k < sz - 1; ++k)
    {
        for(int m = 1; m < M; ++m)
        {
            for(int n = 0; n < sz; ++n)
            {
                data_type x = (k - n + static_cast<data_type>(m) / M) * C_PI;
                ret[k * M + m] += data[n] * (std::sin(x) / x);
            }
        }
    }

    return std::move(ret);
}

Interpolation::Vector<Interpolation::data_type> CubicHermiteInterp::interp(const Vector<data_type> &data, unsigned int interNum) noexcept(false)
{
    // 只要x坐标是等间隔的,波形最终形状就与x无关,因此自己构造与y相等个数的x坐标
    const size_type sz = data.size();
    Vector<data_type> x(sz);
    for(size_type i = 0; i < sz; ++i)
    {
        x[i] = static_cast<data_type>(i * 10);
    }

    return std::move(interp(x, data, interNum).second);
}

std::pair<Interpolation::Vector<Interpolation::data_type>, Interpolation::Vector<Interpolation::data_type> >
CubicHermiteInterp::interp(const Vector<data_type> &x,
                           const Vector<data_type> &y,
                           unsigned int interpNum) noexcept(false)
{
    if(x.size() != y.size())
    {
        throw std::invalid_argument("[CubicHermiteInterp::interp] x.size() != y.size()");
    }

    const size_type sz = x.size();

    // 判断x值是否递增
    data_type deltaX = 0;
    for(size_type i = 0; i < sz - 1; ++i)
    {
        deltaX = x[i+1] - x[i];
        if(deltaX <= 0)
        {
            throw std::invalid_argument("[CubicHermiteInterp::interp] x vector not monotone increasing;");
        }
    }

    // 求取中间每个点处的斜率
    Vector<data_type> derivative(sz);
    for(size_type i = 1; i < sz - 1; ++i)
    {
        if((y[i] > y[i+1] && y[i] > y[i-1])
                || (y[i] < y[i+1] && y[i] < y[i-1]))
        {
            derivative[i] = 0;
            continue;
        }

        double deltaY1 = y[i] - y[i - 1];
        double deltaY2 = y[i+1] - y[i];
        double prevRatio = deltaY1 / (x[i] - x[i - 1]);
        double backRatio = deltaY2 / (x[i + 1] - x[i]);
        double deltaY = deltaY1 + deltaY2;
        // 非极值点,导数应该起到让通过该点的曲线能够单调的平滑过渡到极值点,这里存在优化空间
        if(fuzzyIsNull(deltaY1) || fuzzyIsNull(deltaY2))
        {
            derivative[i] = 0;
        }
        else
        {
            derivative[i] = prevRatio * (deltaY2 / deltaY) + backRatio * (deltaY1 / deltaY);
        }
    }

    // 求第一个点处的斜率
    double ratio0 = (y[0] - y[1]) / (x[0] - x[1]);
    if(fuzzyIsNull(derivative[1]))
    {
        derivative[0] = ratio0 * 2;  // 后一个点是极值点,曲线在第一段应该更陡
    }
    else
    {
        derivative[0] = ratio0 / 2;  // 后一个点不是极值点,曲线在第一段应该趋于平缓
    }

    // 求取最后一个点处的斜率
    const size_type iLast = sz - 1;
    double ratioN = (y[iLast] - y[iLast - 1]) / (x[iLast] - x[iLast - 1]);
    if(derivative[iLast-1] < -1 || derivative[iLast-1] > 1)
    {
        derivative[iLast] = ratioN / 2;    // 倒数第二个点处曲线较平滑,最后一段也应该趋于平缓
    }
    else
    {
        derivative[iLast] = ratioN * 2;    // 倒数第二个点处曲线较陡,最后一段也应该保持这个趋势
    }

    const size_type n = interpNum;
    const size_type M = interpNum + 1;
    const size_type retSz = (x.size() - 1) * interpNum + x.size();
    Vector<data_type> retY(retSz);
    Vector<data_type> retX(retSz);
    for(size_type k = 0; k < sz - 1; ++k)
    {
        retY[k * M] = y[k];
        retX[k * M] = x[k];

        for(size_type j = 1; j <= n; ++j)
        {
            auto currX = x[k] + (x[k+1] - x[k]) * j / (n + 1);
            auto div1 = (currX - x[k+1]) / (x[k] - x[k+1]);
            auto div2 = (currX - x[k]) / (x[k+1] - x[k]);

            auto currY = div1 * div1 * (1 + 2 * div2) * y[k]
                    + div2 * div2 * (1 + 2 * div1) * y[k+1]
                    + div1 * div1 * (currX - x[k]) * derivative[k]
                    + div2 * div2 * (currX - x[k+1]) * derivative[k+1];

            retY[k * M + j] = currY;
            retX[k * M + j] = currX;
        }
    }

    retY.back() = y.back();

    return std::move(std::make_pair(std::move(retX), std::move(retY)));
}


Interpolation::Vector<Interpolation::data_type> SplineInterp::interp(const Vector<data_type> &data,
                                                                     unsigned int interpNum) noexcept(false)
{
    const size_type sz = data.size();
    Vector<data_type> xData(sz);
    for(size_type i = 0; i < sz; ++i)
    {
        xData[i] = static_cast<data_type>(i * 10);
    }

    return std::move(interp(xData, data, interpNum).second);
}

std::pair<Interpolation::Vector<Interpolation::data_type>, Interpolation::Vector<Interpolation::data_type>>
SplineInterp::interp(const Vector<data_type> &x,
                     const Vector<data_type> &y,
                     unsigned int interpNum) noexcept(false)
{
    if(!setDatas(x, y))
    {
        return std::make_pair(Vector<data_type>(), Vector<data_type>());
    }

    const size_type sz = (x.size() - 1) * interpNum + x.size();
    const size_type M = interpNum + 1;

    Vector<data_type> xDatas(sz);
    Vector<data_type> yDatas(sz);
    for(size_type i = 0; i < x.size() - 1; ++i)
    {
        xDatas[M*i] = x[i];
        yDatas[M*i] = y[i];

        data_type delta = (x[i+1] - x[i]) / M;
        data_type xVal = 0;
        for(size_type j = 1; j <= interpNum; ++j)
        {
            xVal = x[i] + delta * j;
            xDatas[M*i+j] = xVal;
            yDatas[M*i+j] = value(xVal);
        }
    }

    xDatas.back() = x.back();
    yDatas.back() = y.back();

    return std::move(std::make_pair(std::move(xDatas), std::move(yDatas)));
}

bool SplineInterp::setDatas(const Vector<data_type> &x, const Vector<data_type> &y) noexcept(false)
{
    if(x.size() != y.size())
    {
        throw std::invalid_argument("[SplineInterp::setDatas] x.size() != y.size();");
    }

    if(x.size() <= 2)
    {
        reset();
        return false;
    }

    const size_type sz = x.size();
    m_x = x;
    m_y = y;

    m_a.resize(sz - 1);
    m_b.resize(sz - 1);
    m_c.resize(sz - 1);

    bool ok = buildSpline(x, y);

    if(!ok)
    {
        reset();
    }

    return ok;
}

Interpolation::data_type SplineInterp::value(data_type xVal) const
{
    if(m_a.size() == 0)
    {
        return 0.0;
    }

    const size_type i = lookup(xVal);

    const data_type delta = xVal - m_x[i];

    return ((((m_a[i] * delta) + m_b[i])
             * delta + m_c[i]) * delta + m_y[i]);
}

bool SplineInterp::isValid() const
{
    return m_a.size() != 0;
}

void SplineInterp::reset()
{
    m_a.resize(0);
    m_b.resize(0);
    m_c.resize(0);

    m_x.resize(0);
    m_y.resize(0);
}

bool SplineInterp::buildSpline(const Vector<data_type> &x, const Vector<data_type> &y)
{
    const size_type sz = x.size();

    const data_type *px = x.data();
    const data_type *py = y.data();

    data_type *a = m_a.data();
    data_type *b = m_b.data();
    data_type *c = m_c.data();

    // 判断x值是否递增
    Vector<data_type> h(sz - 1);
    for(size_type i = 0; i < sz - 1; ++i)
    {
        h[i] = px[i+1] - px[i];
        if(h[i] <= 0)
        {
            return false;
        }
    }

    Vector<data_type> d(sz - 1);
    data_type dy1 = (py[1] - py[0]) / h[0];
    for(size_type i = 1; i < sz - 1; ++i)
    {
        b[i] = c[i] = h[i];
        a[i] = 2.0 * (h[i-1] + h[i]);

        const data_type dy2 = (py[i+1] - py[i]) / h[i];
        d[i] = 6.0 * (dy1 - dy2);
        dy1 = dy2;
    }

    //
    // solve it
    //

    // L-U Factorization
    for(size_type i = 1; i < sz - 2; ++i)
    {
        c[i] /= a[i];
        a[i+1] -= b[i] * c[i];
    }

    // forward elimination
    Vector<data_type> s(sz);
    s[1] = d[1];
    for(size_type i = 2; i < sz - 1; ++i)
    {
        s[i] = d[i] - c[i-1] * s[i-1];
    }

    // backward elimination
    s[sz - 2] = - s[sz - 2] / a[sz - 2];
    for(size_type i = sz - 3; i > 0; --i)
    {
        s[i] = -(s[i] + b[i] * s[i+1]) / a[i];
    }
    s[sz - 1] = s[0] = 0.0;

    //
    // Finally, determine the spline coefficients
    //
    for(size_type i = 0; i < sz - 1; ++i)
    {
        a[i] = (s[i+1] - s[i]) / (6.0 * h[i]);
        b[i] = 0.5 * s[i];
        c[i] = (py[i+1] - py[i]) / h[i] - (s[i+1] + 2.0 * s[i]) * h[i] / 6.0;
    }

    return true;
}

Interpolation::size_type SplineInterp::lookup(data_type xVal) const
{
    size_type i1 = 0;
    const size_type sz = m_x.size();

    if(xVal <= m_x[0])
    {
        i1 = 0;
    }
    else if(xVal >= m_x[sz - 2])
    {
        i1 = sz - 2;
    }
    else
    {
        i1 = 0;
        size_type i2 = sz - 2;
        size_type i3 = 0;

        while(i2 - i1 > 1)
        {
            i3 = i1 + ((i2 - i1) >> 1);

            if(m_x[i3] > xVal)
            {
                i2 = i3;
            }
            else
            {
                i1 = i3;
            }
        }
    }

    return i1;
}
