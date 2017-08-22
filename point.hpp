#ifndef POINT_HPP
#define POINT_HPP

#include <algorithm>
#include <iterator>

// a class representing an n-dimensional point whose coordinates are of type value_type

template <std::size_t n, typename value_type>
struct point
{
    value_type values[n] ;

    point()
    {
    }

    point(const point &p)
    {
        std::copy(p.begin(), p.end(), begin()) ;
    }

    template <typename value_iterator_type>
    point(value_iterator_type first, value_iterator_type last)
    {
        std::copy(first, last, begin()) ;
    }

    point &operator = (const point &p)
    {
        std::copy(p.begin(), p.end(), begin()) ;
        return *this ;
    }

    value_type *begin()
    {
        return &values[0] ;
    }

    value_type *end()
    {
        return &values[n] ;
    }

    value_type const *begin() const
    {
        return &values[0] ;
    }

    value_type const *end() const
    {
        return &values[n] ;
    }
} ;

template <std::size_t n, typename value_type>
bool operator == (const point<n, value_type> &pa, const point<n, value_type> &pb)
{
    return std::equal(pa.begin(), pa.end(), pb.begin()) ;
}

template <std::size_t n, typename value_type>
bool operator != (const point<n, value_type> &pa, const point<n, value_type> &pb)
{
    return !std::equal(pa.begin(), pa.end(), pb.begin()) ;
}

template <std::size_t n, typename value_type>
std::ostream &operator << (std::ostream &stream, const point<n, value_type> &p)
{
    std::copy(p.begin(), p.end(), std::ostream_iterator<value_type>(stream, " ")) ;
    return stream ;
}

template <typename value_type>
point<3, value_type> make_point(const value_type &x, const value_type &y, const value_type &z)
{
    const value_type values[3] = {x, y, z} ;
    return point<3, value_type>(&values[0], &values[3]) ;
}

// todo(marke): make vector class by modifying a copy of point class, it would be nice if these classes included inplace operations such as +=
//              and the binary versions below eree defined in terms of these

// ADDITIONAL FUNCTIONS CREATED HERE

// todo(marke): dot should be on vectors
template <std::size_t n, typename value_type>
value_type dot(const point<n, value_type> &pa, const point<n, value_type> &pb) 
{
    // todo(marke): use standard algorithm beginning with i
    value_type sum=0;
    for(int i=0; i<n;i++)
    {
        sum+=(pa.values[i])*(pb.values[i]);
    }
    return sum;
}

// todo(marke): cross should be on vectors
template <typename value_type>
point<3, value_type> cross(const point<3, value_type> &pa, const point<3, value_type> &pb) 
{
    value_type x = (pa.values[1])*(pb.values[2])-(pa.values[2])*(pb.values[1]);
    value_type y = (pa.values[2])*(pb.values[0])-(pa.values[0])*(pb.values[2]);
    value_type z = (pa.values[0])*(pb.values[1])-(pa.values[1])*(pb.values[0]);

    // todo(marke): use make_point function instead of duplicating logic (actually use make_vector once vector class exists)
    const value_type values[3] = {x, y, z} ;
    return point<3, value_type>(&values[0], &values[3]) ;
}

// todo(marke): it would be nice if this was actually "modulus_sqr", avoiding the square root, and also to reuse the dot product here.
//              it is good to reuse code and minimise the duplication of logic
template <std::size_t n, typename value_type>
value_type modulus(const point<n, value_type> &pa) 
{
    value_type sum=0;
    for(int i=0; i<n;i++)
    {
        sum+=(pa.values[i])*(pa.values[i]);
    }
    // todo(marke): always use std versions of maths C runtime functions, such as sqrt; these will be overloaded for types other than doubles
    return sqrt(sum);
}

// todo(marke): normalise should be on vectors but would nice to avoid as it involved sqrt
template <std::size_t n, typename value_type>
const point<n, value_type> normalise(const point<n, value_type> &pa) 
{
    value_type mod=modulus(pa);

    // todo: use shared tolerance function?
    // check on modulus size
    value_type tol=1.0e-6;
    if(mod<=tol)
    {
        std::cout<<"!! W01: small modulus when normalising !!"<<std::endl;
    }

    return pa/mod;
}

// todo(marke): operator + should be on vectors not points
template <std::size_t n, typename value_type>
point<n, value_type> operator+(const point<n, value_type> &pa, const point<n, value_type> &pb)
{
    // todo(marke): use std algorithm starting with t
    point<n, value_type> out;
    for(int i=0;i<n;i++)
        out.values[i]=pa.values[i]+pb.values[i];
    return out ;
}

// todo(marke): result of operator - should be vector
template <std::size_t n, typename value_type>
point<n, value_type> operator-(const point<n, value_type> &pa, const point<n, value_type> &pb)
{
    // todo(marke): use std algorithm starting with t
    point<n, value_type> out;
    for(int i=0;i<n;i++)
        out.values[i]=pa.values[i]-pb.values[i];
    return out ;
}

// todo(marke): operator * should be on vector not point
template <std::size_t n, typename value_type>
point<n, value_type> operator*(const point<n, value_type> &pa, const value_type &s)
{
    // todo(marke): use std algorithm starting with t
    point<n, value_type> out;
    for(int i=0;i<n;i++)
        out.values[i]=pa.values[i]*s;
    return out ;
}

// todo(marke): operator / should be on vector not point
template <std::size_t n, typename value_type>
point<n, value_type> operator/(const point<n, value_type> &pa, const value_type &s)
{
    // todo(marke): use std algorithm starting with t
    point<n, value_type> out;
    for(int i=0;i<n;i++)
        out.values[i]=pa.values[i]/s;
    return out ;
}


#endif
