#ifndef INTERSECT_HPP
#define INTERSECT_HPP

#include "simplex.hpp"


// classes for 3d "line" and "plane" created here

// todo(marke): might be nice to template line on dimension as well as it is a dimension independent concept
// todo(marke): I would usually expect "struct line" to be underneath the template specification
template <typename value_type> struct line
{
    point<3,value_type> start,end;

    // todo (marke): it would be nice to avoid this default constructor and force clients to use the one which ensures we have a valid line
    line()
    {
    }

    // todo(marke): pass by reference (passing by value will copy at least 3*8=24 bytes per point onto call stack per point with double coordinates
    //              a reference may be eight bytes, so we'd save coping 24*2-8*2 bytes when constructing a line)
    // todo(marke): use initialiser lists (cleaner and more efficient as start and end are initialised with correct value in one step,
    //              rather than default constructed and then assigned to new value immediately)
    line(const point<3,value_type> istart, const point<3,value_type> iend)
    {
        start=istart;
        end=iend;
    }
    // using default assigment 
};

// todo(marke): apply same comments as line class to plane class
template <typename value_type> struct plane
{
    // todo(marke): plane normal should be a vector, not a point
    point<3,value_type> nor;    //plane equation: for X in plane, nor.(X)=dist
    value_type dist;            

    plane()
    {
    }

    plane(const value_type idist, const point<3,value_type> inor)
    {
        dist=idist;
        nor=inor;
    }
    // using default assignment
};

// subroutines used for the final intersect function
template <typename value_type>
plane<value_type> get_plane(const simplex<2,3,value_type> &triangle)
{
    // todo (marke): it is inefficint to declare a non-plain-old-data variable and then assign to it as the default constructor is called
    //               when it is declared and then we immediately overwrite that instance with a new one (see use of initialiser lists in
    //               line/plane classes. Use simultaneous declare and assign.
    // todo (marke): taking a copy of each point is wasteful, a const reference will be just fine
    point<3,value_type> a,b,c;
    a=triangle.vertices[0];
    b=triangle.vertices[1];
    c=triangle.vertices[2];

    // todo(marke): similar comments to the declaration/assignment of points
    // todo(marke): use constructor of plane class
    // todo(marke): make out const if you do not intend it to change
    plane<value_type> out;
    out.nor=cross(a-b,b-c);
    out.dist=dot(out.nor,a);

    // check for degenerate triangles (such as the initial input!)
    // todo (marke): it would be nice if the call to modulus avoided the use of a sqrt and we tested against tol * tol, or even tol_sqr here. 
    //               sqrt is much slower than a multiply, plus value_type could be some exotoc type which does not handle sqrt
    // todo (marke): if you reuse the 1.0e-6 tolerance value elsewhere then it would be good to define this in only one place, such as with a 
    //               function
    // todo (marke): you may want to look at throwing an exception instead of printing out a warning
    value_type tol=1.0e-6;
    if(modulus(out.nor)<tol)
        std::cout << "!! degenerate triangle input !! "<<std::endl;

    return out;
}

template <typename value_type>
line<value_type> intersect_plane_plane(const plane<value_type> &plane_1, const plane<value_type> &plane_2)
{
    
    // direction of line is perpendicular to both plane normals
    // todo(marke): apply comments above regarding declaring/assigning/using constructors
    point<3,value_type> dir;
    dir=cross(plane_1.nor,plane_2.nor);
    // alert for possible degeneracies

    // todo (marke): we're now looking at the square of the modulus and use a different tol
    // todo (marke): maybe throw an exception
    value_type tol=1.0e-8;
    if(dot(dir,dir)<tol)
    {
        std::cout << "!! degenerate plane plane intersect !! "<<std::endl;
    }
    // find a point on intersecting line:
    // intersect with additional plane, 
    plane<value_type> plane_3;
    // Either YZ,XZ, or XY plane is chosen. Decision based upon direction of line
    // Need to find largest coordinate (then encode in plane_3_flag respectively)
    int plane_3_flag=0; 
    if(fabs(dir.values[0])>fabs(dir.values[1])) plane_3_flag=0;
    else plane_3_flag=1;
    if(fabs(dir.values[2])>fabs(dir.values[plane_3_flag]))
        plane_3_flag=2;
    // building the third plane
    if        (plane_3_flag==0)    plane_3.nor=make_point(1.0,0.0,0.0);
    else if (plane_3_flag==1)    plane_3.nor=make_point(0.0,1.0,0.0);
    else                        plane_3.nor=make_point(0.0,0.0,1.0);
    plane_3.dist=0;
    // intersecting the three planes to get point
    // todo (marke): good to declare and assign and same time, nice if p was cost thoug as it does not change
    point<3,value_type> p=intersect_plane_plane_plane(plane_1,plane_2,plane_3);
    //dir=normalise(dir);
    // todo(amrke): use constructor of line
    line<value_type> out;
    out.start=p;
    out.end=p+dir;
    //std::cout << "modulus: " << modulus(dir) << std::endl;
    return out;
}

// todo (marke): it would be nice if the linear system solving was in a separate, reusable function
// todo (marke): should use vectors, not points here
template <typename value_type>
point<3,value_type> intersect_plane_plane_plane(const plane<value_type> &plane_1, const plane<value_type> &plane_2, const plane<value_type> &plane_3)
{
    // todo (marke): apply same comments above regarding declaration/assignment/const references/etc.
    point<3,value_type> N1,N2,N3;
    N1=plane_1.nor;
    N2=plane_2.nor;
    N3=plane_3.nor;
    value_type D1,D2,D3;
    D1=plane_1.dist;
    D2=plane_2.dist;
    D3=plane_3.dist;
    // must solve equation [A][X]=[c], where A=[N1 N2 N3], c=[D1 D2 D3]^T 
    // First find determinant: |A|=1/N1.(N2xN3)
    value_type det=dot(N1,cross(N2,N3));
    // Columns of A^-1 can be found from cross products A^-1=[N2xN3 N3xN1 N1xN2]^T / |A|
    point<3,value_type> C1,C2,C3;
    C1=cross(N2,N3);
    C2=cross(N3,N1);
    C3=cross(N1,N2);

    // todo(marke): would be nice to check determinant is not zero/very small here and perhaps throw an excepption if there's an issue
    point<3,value_type> out;
    out=C1*D1+C2*D2+C3*D3;
    out=out/det;
    return out;
}

template <typename value_type>
point<3,value_type> intersect_coplanar_lines(const line<value_type> &segment,const line<value_type> &line_2, value_type &convex_ratio)
{
    // todo (marke): apply same comments above regarding declaration/assignment/const references/etc.
    point<3,value_type> a0,a1,b0,b1,da,db,N;
    a0=segment.start;
    a1=segment.end;
    b0=line_2.start;
    b1=line_2.end;

    // need to solve a0+(a1-a0)*t=b0+(b1-b0)*s
    // cross with (b1-b0) to get a0x(b1-b0)+(a1-a0)x(b1-b0)*t=b0x(b1-b0)
    // rearranging gives (a1-a0)x(b1-b0)*t=(b0-a0)x(b1-b0)
    // label (a1-a0)x(b1-b0) as N (since its the plane normal)
    // t can thus be solved as t=(((b0-a0)x(b1-b0)).N)/(N.N)
    // t corresponds to the convex paramater for the segment

    da=a1-a0;
    db=b1-b0;
    N=cross(da,db);
    // todo(marke): might be nice to have a det function to go with your linear system solve function
    convex_ratio=dot(cross((b0-a0),db),N);
    // todo(marke): might be nice to have a norm_sqr function based on dot (would replace modulus and would not use sqrt)
    convex_ratio=convex_ratio/(dot(N,N));

    // pt is found by subs back into equation a0+(a1-a0)*t
    point<3,value_type> out;
    out=a0+da*convex_ratio;
    return out;
}

template <typename value_type>
bool intersect_line_triangle(const line<value_type> &iline, const simplex<2,3,value_type> &triangle, line<value_type> &segment_out)
{
    // todo (marke): apply same comments above regarding declaration/assignment/const references/etc.
    // intersect input line with edges of triangle
    line<value_type> edge1,edge2,edge3;
    edge1.start=triangle.vertices[0]; edge1.end=triangle.vertices[1];
    edge2.start=triangle.vertices[1]; edge2.end=triangle.vertices[2];
    edge3.start=triangle.vertices[2]; edge3.end=triangle.vertices[0];

    point<3,value_type> pts[3];
    value_type cr[3];
    
    pts[0]=intersect_coplanar_lines(edge1,iline,cr[0]);
    pts[1]=intersect_coplanar_lines(edge2,iline,cr[1]);
    pts[2]=intersect_coplanar_lines(edge3,iline,cr[2]);

    // determine which (if any) pair of edges are crossed by checking 0<cr<1
    bool convex[3];
    for(int i=0;i<3;i++)
    {
        // todo(marke): replace this with a single like "convex[i] = <boolean_expression>"
        if(cr[i]<0 || cr[i]>1)
            convex[i]=false;
        else
            convex[i]=true;
    }
    
    // check if all edges are not crossed by line
    // todo(marke): use std::find"
    if(!convex[0] && !convex[1] && !convex[2])
        return false;

    // find intersecting line segment before returning true
    // if any false, pick out the other two points corresponding to true
    if(!convex[0] || !convex[1] || !convex[2])    
    {
        // todo(marke): replace with ternary operator
        if(convex[0])
            segment_out.start=pts[0];
        else
            segment_out.start=pts[1];

        // todo(marke): replace with ternary operator
        if(convex[2])
            segment_out.end=pts[2];
        else
            segment_out.end=pts[1];
    }
    else // if all true, line must go through a vertex. The following 
    {
        segment_out.start=pts[0];
        segment_out.end=pts[2];
        // todo(marke): use shared tolerance value function..?
        value_type tol=1.0e-6;
        if(modulus(pts[0]-pts[2])<tol) //checks if pts[0]=pts[2]
        {
            segment_out.end=pts[1];
        }
    }

    return true;
}

template <typename value_type>
bool colinear_segments_overlap(const line<value_type> &seg1, const line<value_type> &seg2)
{
    //translate segments such that seg1.start is on the origin
    //then check seg2 end points are between origin and seg1.end;
    point<3,value_type> a=seg2.start-seg1.start,b=seg2.end-seg1.start,c=seg1.end-seg1.start;

    // todo (marke): it would be nice if these were const as they don't change
    value_type ch_dist=dot(c,c);  
    value_type dist1=dot(a,c); 
    value_type dist2=dot(b,c); 

    // todo (marke): could replace this with returning a single boolean expression
    if(dist1>=0 && dist1<=ch_dist)
        return true;
    else if(dist2>=0 && dist2<=ch_dist)
        return true;
    else if(dist1<=0 && dist2>=ch_dist)
        return true;
    else if(dist2<=0 && dist1>=ch_dist)
        return true;
    else
        return false;
    
    
    

}

template <typename value_type>
bool intersect(const simplex<2, 3, value_type> &triangle_a, const simplex<2, 3, value_type> &triangle_b)
{
    // STEP #1 intersect the two planes of the triangle to get a line

    // todo (marke): it would be nice if these were const as they don't change
    plane<value_type> plane_a=get_plane(triangle_a);
    plane<value_type> plane_b=get_plane(triangle_b);
    line<value_type> iline=intersect_plane_plane(plane_a,plane_b);

    // todo (marke): it would be nice if the below was a single line: exploit short-circuit evaluation of boolean expressions

    // STEP #2 intersect line with triangles to get segments
    
    line<value_type> seg1,seg2;
    bool check1=intersect_line_triangle(iline,triangle_a,seg1);
    // todo (marke): check1==false is just !check1
    if(check1==false) 
        return false;
    bool check2=intersect_line_triangle(iline,triangle_b,seg2);
    // todo (marke): check2==false is just !check2
    if(check2==false) 
        return false;
    
    // STEP #3 check segments overlap
    bool check3=colinear_segments_overlap(seg1,seg2);
    return check3;
}


#endif
