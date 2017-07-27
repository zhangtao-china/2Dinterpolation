#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <cmath>


#ifdef _MSC_FULL_VER
    #if _MSC_VER < 1600
        #error "should use C++11 implementation"
    #else
        #pragma execution_character_set("utf-8")
    #endif
#elif __cplusplus < 201103L
    #error "should use C++11 implementation"
#endif

class Interpolation
{
public:
    using data_type = double;   
    template<typename T> using Vector = std::vector<T>;  
    using size_type = Vector<data_type>::size_type;

    Interpolation() = default;
    virtual ~Interpolation() = default;

    // 定义接口
    /// data: 需要插值的原始数据
    /// interNum: 每相邻两个点之间插入数据的个数
    virtual Vector<data_type> interp(const Vector<data_type> &data, unsigned int interNum) noexcept(false) = 0;
};

// 正弦插值
class SineInterp : public Interpolation
{
public:

    Vector<data_type> interp(const Vector<data_type> &data, unsigned int interpNum) noexcept(false) override;
};

// 三次埃尔米特插值
class CubicHermiteInterp : public Interpolation
{
public:

    /// x坐标等间隔时并且不关心x坐标,可调用此函数
    Vector<data_type> interp(const Vector<data_type> &data, unsigned int interNum) noexcept(false) override;

    /// x坐标非等间隔或同时需要x坐标, 调用此函数
    /// 返回值: pair--first为x坐标向量, second为y坐标向量
    std::pair<Vector<data_type>, Vector<data_type>>
    interp(const Vector<data_type> &x,
           const Vector<data_type> &y,
           unsigned int interpNum) noexcept(false);
};

// 三次样条曲线插值
class SplineInterp : public Interpolation
{
public:
    /// x坐标等间隔时并且不关心x坐标,可调用此函数
    Vector<data_type> interp(const Vector<data_type> &data, unsigned int interpNum)noexcept(false) override;

    /// x坐标非等间隔或同时需要x坐标, 调用此函数
    /// 返回值: pair--first为x坐标向量, second为y坐标向量
    std::pair<Vector<data_type>, Vector<data_type>>
    interp(const Vector<data_type> &x,
           const Vector<data_type> &y,
           unsigned int interpNum) noexcept(false);

protected:
    bool setDatas(const Vector<data_type> &x, const Vector<data_type> &y) noexcept(false);

    void reset();

    bool isValid() const;

    data_type value(data_type xVal) const;

    bool buildSpline(const Vector<data_type> &x, const Vector<data_type> &y);

    size_type lookup(data_type xVal) const;

private:
    // 系数
    Vector<data_type> m_a;
    Vector<data_type> m_b;
    Vector<data_type> m_c;

    // 控制点
    Vector<data_type> m_x;
    Vector<data_type> m_y;
};


#endif // INTERPOLATION_H
