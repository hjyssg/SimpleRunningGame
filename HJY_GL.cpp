//
//
//  HJY_GL.cpp
//
//  Created by Junyang Huang on 11/6/13.
//  Copyright (c) 2013 . All rights reserved.
//




#include <vector>
#include <set>
#include <iostream>

#include <string>
#include <fstream>
#include <sstream>

#include <stdio.h>
#include <string.h>
#include "math.h"
#include <stdlib.h>

#if defined(_MSC_VER)
#include <Gl/glut.h>
#include <windows.h>
#elif defined(__APPLE__)
#include <GlUT/glut.h>
#include <unistd.h>
#endif


#define  CLOG(a)  std::cout<<a<<"\n"

#if DEBUG
#define  DLOG(a)  std::cout<<__PRETTY_FUNCTION__<<": "<<a<<"\n"
#else
#define DLOG(a)
#endif

#define PI 3.1415926
#define DEG_TO_RAD(a)  a/180*3.1415
#define RAD_TO_DEG(a)  a/3.1415*180

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


//return false if no root
inline bool solve_cube_equation(float a, float b, float c, float *t1, float *t2)
{
    float delta = b*b - 4*a*c;
    if(delta<0)
    {
        return false;
    }
    
    *t1 = (-b + sqrt(delta)) / (2*a);
    *t2 = (-b - sqrt(delta)) / (2*a);
    
    return true;
}


inline float dot_prod(float v1[], float v2[])
{
    float tot = 0;
    int i;
    for (i=0; i<4; i++)
        tot += v1[i]*v2[i];
    return tot;
}


inline void normalize(GLfloat *p)
{
    double d=0.0;
    int i;
    for(i=0; i<3; i++) d+=p[i]*p[i];
    d=sqrt(d);
    if(d > 0.0) for(i=0; i<3; i++) p[i]/=d;
}


inline  void cross_prod(GLfloat *a, GLfloat *b, GLfloat *c, GLfloat *d)
{
    d[0]=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]);
    d[1]=(b[2]-a[2])*(c[0]-a[0])-(b[0]-a[0])*(c[2]-a[2]);
    d[2]=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
    normalize(d);
}


inline void print_matrix(float my_matrix[])
{
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            printf ("%f ", my_matrix[i+j*4]);
        }
        printf ("\n");
    }
    printf ("\n");
}

inline int find_min_in_vector(std::vector<int>vect)
{
    int temp = vect.at(0);
    for (int ii = 0; ii<0; ii++)
    {
        if(vect.at(ii)<temp)
        {
            temp = vect.at(ii);
        }
    }
    return temp;
}

inline void draw_ground1()
{
    
    
    // Ground
    int minx = -20;
    int maxx = 20;
    int minz = -20;
    int maxz = 20;
    
    glLineWidth(3);
    glColor4f(0, 1, 0, 1);
    for (int x = minx; x <= maxx; x+=2)
    {
        glBegin(GL_LINES);
        glVertex3f(x,-1, minz);
        glVertex3f(x, -1, maxz);
        glEnd();
    }
    
    glColor4f(0, 0, 1, 1);
    for (int z = minz; z <= maxz; z+=2) {
        glBegin(GL_LINES);
        glVertex3f(minx,-1, z);
        glVertex3f(maxx,-1, z);
        glEnd();
    }
    
    
}

class Vertex
{
    
public:
    float x;
    float y;
    float z;
    float w;
    
public:
    Vertex()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
        this->w = 1;
    }
    
    Vertex(float x_p, float y_p, float z_p, float w_p)
    {
        this->x = x_p;
        this->y = y_p;
        this->z = z_p;
        this->w = w_p;
    }
    
    
    Vertex(float x_p, float y_p, float z_p)
    {
        this->x = x_p;
        this->y = y_p;
        this->z = z_p;
        this->w = 1;
    }
    
    Vertex(float * arr, int size)
    {
        if (size == 3)
        {
            x = arr[0];
            y = arr[1];
            z = arr[2];
            w = 1;
        }
        else if (size == 4)
        {
            x = arr[0];
            y = arr[1];
            z = arr[2];
            w = arr[4];
        }
        else
        {
            x = NAN;
            y = NAN;
            z = NAN;
            w = 1;
        }
    }
    
    void set_xyz(float x_p, float y_p, float z_p)
    {
        x = x_p;
        y = y_p;
        z = z_p;
    }
    
    Vertex get_copy( )
    {
        return Vertex(x,y,z,w);
    }
    
    
    Vertex  add(Vertex v2)
    {
        Vertex result;
        result.x = x + v2.x;
        result.y = y + v2.y;
        result.z = z + v2.z;
        result.w = 1;
        return result;
    }
    
    Vertex sub(Vertex v2)
    {
        Vertex result;
        result.x = x - v2.x;
        result.y = y - v2.y;
        result.z = z - v2.z;
        result.w = 1;
        return result;
    }
    
    Vertex crosss_product(Vertex v2)
    {
        Vertex result;
        result.x =   y * v2.z - z*v2.y;
        result.y = -(x * v2.z - z*v2.x);
        result.z =   x * v2.y - y*v2.x;
        result.w = 1;
        return result;
    }
    
    float dot_product(Vertex v2)
    {
        return  x * v2.x +  y * v2.y +  z * v2.z ;
    }
    
    
    Vertex  transfer( float x_p, float y_p, float z_p)
    {
        Vertex result;
        result.x = x + x_p;
        result.y = y + y_p;
        result.z = z + z_p;
        result.w = w;
        return result;
    }
    
    void  transfer_self( float x_p, float y_p, float z_p)
    {
        x = x + x_p;
        y = y + y_p;
        z = z + z_p;
    }
    
    Vertex scale( float x_s, float y_s, float z_s)
    {
        Vertex result;
        result.x = x * x_s;
        result.y = y * y_s;
        result.z = z * z_s;
        result.w = w;
        return result;
    }
    
    inline void scale_self( float x_s, float y_s, float z_s)
    {
        x = x * x_s;
        y = y * y_s;
        z = z * z_s;
    }
    
    Vertex rotate_x(float theta)
    {
        Vertex result;
        result.x = x;
        result.y = y*cos(theta) - sin(theta)*z;
        result.z = y*sin(theta) + cos(theta)*z;
        result.w = w;
        return result;
    }
    
    void rotate_x_self(float theta)
    {
        float tempy = y;
        float tempz = z;
        y = tempy*cos(theta) - sin(theta)*tempz;
        z = tempy*sin(theta) + cos(theta)*tempz;
    }
    
    Vertex rotate_y(float theta)
    {
        Vertex result;
        result.x = x*cos(theta) + sin(theta)*z;
        result.y = y;
        result.z = -x*sin(theta) + cos(theta)*z;
        result.w = w;
        return result;
    }
    
    void rotate_y_self(float theta)
    {
        float tempx = x;
        float tempz = z;
        
        x = tempx*cos(theta) + sin(theta)*tempz;
        z = -tempx*sin(theta) + cos(theta)*tempz;
    }
    
    Vertex rotate_z(float theta)
    {
        Vertex result;
        result.x = x*cos(theta) - sin(theta)*y;
        result.y = x*sin(theta) + cos(theta)*y;
        result.z = z;
        result.w = w;
        return result;
    }
    
    void rotate_z_self(float theta)
    {
        float tempx = x;
        float tempy = y;
        x = tempx*cos(theta) - sin(theta)*tempy;
        y = tempx*sin(theta) + cos(theta)*tempy;
    }
    
    //vtx is origined from the (0,0,0)
    Vertex rotate_around_axis(Vertex vtx, float theta)
    {
        //if angle is 0 no need to calculate
        if(theta==0)
        {
            return get_copy();
        }
        
        static Vertex cache_v = Vertex();
        static float cache_theta = 0.0;
        
        static Vertex r1 = Vertex();
        static Vertex r2 = Vertex();
        static Vertex r3 = Vertex();
        
        //only calculate matrix once for each vertex and degress combination
        if (!(vtx.isEqual(cache_v)&&cache_theta==theta))
        {
            float cc = cos(theta);
            float ss = sin(theta);
            
            float _cc = 1-cc;
            
            float xy = vtx.x * vtx.y;
            float zy = vtx.z * vtx.y;
            float zx = vtx.z * vtx.x;
            
            float xx = vtx.x * vtx.x;
            float yy = vtx.y * vtx.y;
            float zz = vtx.z * vtx.z;
            
            r1 = Vertex( cc+ xx* _cc,
                        xy*_cc - vtx.z*ss,
                        zx*_cc + vtx.y*ss);
            
            r2 = Vertex( xy*_cc + vtx.z*ss,
                        cc + yy *_cc,
                        zy*_cc- vtx.x*ss
                        );
            
            r3 = Vertex(zx *_cc - vtx.y*ss,
                        zy *_cc + vtx.x*ss,
                        cc +  zz*_cc
                        );
            
            cache_theta = theta;
            cache_v = vtx;
            
        }
        
        Vertex result;
        result.x = this->dot_product(r1);
        result.y = this->dot_product(r2);
        result.z = this->dot_product(r3);
        
        result.w = w;
        return result;
    }
    
    
    float length()
    {
        return  sqrtf(x*x + y*y + z*z);
    }
    
    Vertex normalize()
    {
        float len = this->length();
        float temp = 1/len;
        if (len==0)
        {
            temp = 1;
        }
        
        Vertex result;
        result.x = x * temp;
        result.y = y * temp;
        result.z = z * temp;
        result.w = 1;
        
        return result;
    }
    
    void normalize_self()
    {
        float len = this->length();
        if (len==0)
        {
            return;
        }
        
        float temp = 1/len;
        
        x = x * temp;
        y = y * temp;
        z = z * temp;
    }
    
    bool isZero()
    {
        return x==0&&y==0&&z==0;
    }
    
    bool isEqual(Vertex v)
    {
        return (this->x==v.x && this->y==v.y && this->z==v.z);
    }
    
    static Vertex get_x_axis()
    {
        return Vertex(1, 0, 0);
    }
    
    
    static Vertex get_y_axis()
    {
        return Vertex(0, 1, 0);
    }
    
    
    static Vertex get_z_axis()
    {
        return Vertex(0, 0, 1);
    }
    
    static Vertex get_null_vertex()
    {
        return Vertex(NAN, NAN, NAN);
    }
    
    bool isNull()
    {
        return  (isnan(x) || isnan(y) || isnan(z));
    }
    
    bool isParrell(Vertex v)
    {
        //zero vector is parrell to any vector
        if (v.isZero()||this->isZero())
        {
            return true;
        }
        
        float r1 = this->x / v.x;
        float r2 = this->y / v.y;
        float r3 = this->z / v.z;
        
        return (r1==r2 && r2 ==r3&&r3 == r1);
    }
    
    //returns in the range 0 to π radians.
    float angle_between_v(Vertex v)
    {
        float temp =  dot_product(v)/(length()*v.length());
        return  acosf(temp);
    }
    
    float distance_to_v(Vertex v)
    {
        return sub(v).length();
    }
    
    Vertex project_on_xy_plane()
    {
        return Vertex(this->x,this->y,0);
    }
    
    Vertex project_on_xz_plane()
    {
        return Vertex(this->x,0,this->z);
    }
    
    Vertex project_on_yz_plane()
    {
        return Vertex(0,this->y,this->z);
    }
    
    Vertex project_on_v(Vertex v)
    {
        float len = v.length();
        float temp = dot_product(v)/len*len;
        return v.scale(temp, temp, temp);
    }
    
    
    
    //use inline to solve duplicate simbol error
    inline static void draw_triangle_Vertex(int gl_mode_p, Vertex v1, Vertex v2, Vertex v3, int red, int green, int blue)
    {
        glBegin(gl_mode_p);
        {
            glColor3f(red, green, blue);
            glVertex4f(v1.x, v1.y, v1.z, v1.w);
            glVertex4f(v2.x, v2.y, v2.z, v2.w);
            glVertex4f(v3.x, v3.y, v3.z, v3.w);
        }
        glEnd();
    }
    
    inline static void draw_rectangle_Vertex(int gl_mode_p, Vertex v1, Vertex v2, Vertex v3, Vertex v4, int red, int green, int blue)
    {
        glBegin(gl_mode_p);
        {
            glColor3i(red, green, blue);
            
            glVertex4f(v1.x, v1.y, v1.z, v1.w);
            
            glVertex4f(v2.x, v2.y, v2.z, v2.w);
            
            glVertex4f(v3.x, v3.y, v3.z, v3.w);
            
            glVertex4f(v4.x, v4.y, v4.z, v4.w);
        }
        glEnd();
    }
    
    inline static void draw_two_vertex(int gl_mode_p, Vertex v1, Vertex v2, int red, int green, int blue)
    {
        glBegin(gl_mode_p);
        {
            glColor3f(red, green, blue);
            glVertex4f(v1.x, v1.y, v1.z, v1.w);
            glVertex4f(v2.x, v2.y, v2.z, v2.w);
        }
        glEnd();
    }
    
    std::string toString()
    {
        char tempStr[50];
        sprintf(tempStr, "[%f, %f, %f, %f]", x,y,z, w);
        std::string result(tempStr);
        return result;
    }
};

#define vertex_to_arr(v) {v.x, v.y, v.z, v.w}

//capture material info
class MaterialAttribute
{
    
public:
    //ambient intensity;
    Vertex amb_itn;
    
    //diffuse intensity
    Vertex dif_itn;
    
    //specular intensity
    Vertex spc_itn;
    
    
    //emission intensity
    Vertex ems_itn;
    
    //shiness intensity
    float shn_itn;
    
public:
    MaterialAttribute()
    {
        amb_itn = Vertex();
        dif_itn = Vertex();
        spc_itn = Vertex();
        ems_itn = Vertex();
        
        shn_itn = 0;
    }
    
    void set_up()
    {
        float amb[] = vertex_to_arr(amb_itn);
        float diff[] = vertex_to_arr(dif_itn);
        float spec[] = vertex_to_arr(spc_itn);
        float emi[] = vertex_to_arr(ems_itn);
        
        glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diff);
        glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
        glMaterialfv(GL_FRONT, GL_SHININESS, &shn_itn);
        glMaterialfv(GL_FRONT, GL_EMISSION, emi);
    }
    
    void print_debug_info()
    {
        CLOG("ambient:  "<<amb_itn.toString()<<"\n"<<
             "diffuse:  "<<dif_itn.toString()<<"\n"<<
             "specular: "<<spc_itn.toString()<<"\n"<<
             "shiness "<<shn_itn);
    }
};

#define GRAVIATIONAL_CONSTANT 0.01
class PointMass
{
public:
    
    double mass;
    Vertex pos;
    Vertex spd;
    Vertex acl;
    float radius;
    std::string name;
    
public:
    PointMass()
    {
        mass = 0;
        radius = 1;
        pos = Vertex();
        spd = Vertex();
        acl = Vertex();
    }
    
    void update_state(float time_step)
    {
        float ts = time_step; //syntax candy
        pos.x += spd.x*ts;
        pos.y += spd.y*ts;
        pos.z += spd.z*ts;
        
        spd.x += acl.x*ts;
        spd.y += acl.y*ts;
        spd.z += acl.z*ts;
    }
    
    
    //the gravitation that this poses on pm
    //this  <-attract--- pm
    Vertex gravitation_acceleration_on_pm(PointMass *pm)
    {
        Vertex direction = pos.sub(pm->pos);
        float distance = direction.length();
        
        float scalar_g = mass/distance/distance*GRAVIATIONAL_CONSTANT;
        
        // normalized and then times scalar
        float temp = 1/distance*scalar_g;
        direction.scale_self(temp, temp, temp);
        
        return direction;
    }
    
    void print_debug_info()
    {
        CLOG(name<<" position:  "<<pos.toString()<<"\n"<<
             "speed:  "<<spd.toString()<<"\n"<<
             "acceleration: "<<acl.toString()<<"\n"<<
             "mass: "<<mass);
    }
};

//a simulation to the real physics world where is not external force but the gravitation among point masses
class TheWorld
{
public:
    std::vector<PointMass *> pmArr;
    
public:
    TheWorld()
    {
        
    }
    
    void add_point_mass(PointMass * pm)
    {
        pmArr.push_back(pm);
    }
    
    //assert at least two point mass
    void update_state(float time_step)
    {
        for(int ii = 0; ii < pmArr.size();ii++)
        {   //clear acl
            PointMass *pm  = pmArr.at(ii);
            pm->acl.set_xyz(0, 0, 0);
        }
        
        //use n^2 algo now to firgure out graviation
        for(int ii = 0; ii < pmArr.size()-1;ii++)
        {
            PointMass *pm  = pmArr.at(ii);
            
            for(int jj = ii+1; jj < pmArr.size();jj++)
            {
                PointMass *pm_2  = pmArr.at(jj);
                
                Vertex g =  pm_2->gravitation_acceleration_on_pm(pm);
                
                pm->acl.transfer_self(g.x, g.y, g.z);
                
                float temp = -pm->mass/pm_2->mass;
                g.scale_self(temp,temp, temp);
                
                //pm_2->acl = pm_2->acl.add(g2);
                pm_2->acl.transfer_self(g.x, g.y, g.z);
                
                if(jj==pmArr.size()-1)
                {
                    pm_2->update_state(time_step);
                }
            }
            pm->update_state(time_step);
        }
    }
    
    std::set<std::pair<int, int> > detect_collision()
    {
        std::set<std::pair<int, int> >  result;
        
        for(int ii = 0; ii < pmArr.size()-1;ii++)
        {
            PointMass *pm  = pmArr.at(ii);
            
            for(int jj = ii+1; jj < pmArr.size();jj++)
            {
                PointMass *pm_2  = pmArr.at(jj);
                float distance = pm->pos.distance_to_v(pm_2->pos);
                
                if (distance<(pm_2->radius+pm->radius))
                {
                    std::pair<int, int> pm_pair(ii, jj);
                    result.insert(pm_pair);
                }
            }
        }
        return result;
    }
};



class ThreeDimensionObject
{
public:
    
    //the center position
    Vertex origin;
    
    //the bottom vertex and top vertex
    //used to define orientation of object
    Vertex bottom_vertex;
    Vertex top_vertex;
    
    float x_scaler;
    float y_scaler;
    float z_scaler;
    
    int num_stack;  //cut along y-axis
    int num_slice; //around the Y-axis
    
    int r;
    int g;
    int b;
    
    
    GLuint textName1;
    
    MaterialAttribute  material_attribute;
    
    //num_stack * num_stack * 4
    std::vector<std::vector<Vertex> > vertex_matrix;
    
    //their dimension are the same
    std::vector<std::vector<Vertex> > vertex_normal_matrix;
    
public:
    
    void init_default_parameters(float x_p, float y_p, float z_p)
    {
        this->origin = Vertex(x_p,y_p,z_p);
        
        this->x_scaler = 1;
        this->y_scaler = 1;
        this->z_scaler = 1;
        
        
        //use random color unless set by user
        r =  rand();
        g =  rand();
        b = rand();
    }
    
    //default constructor
    ThreeDimensionObject()
    {
        init_default_parameters(0,0,0);
        
        this->bottom_vertex = Vertex(0,-1,0);
        this->top_vertex = Vertex(0,1,0);
        
        this->num_stack = 0;
        this->num_slice = 0;
    }
    
    ThreeDimensionObject(float x_p, float y_p, float z_p)
    {
        init_default_parameters(x_p, y_p, z_p);
        
        this->bottom_vertex = Vertex(0,-1,0);
        this->top_vertex = Vertex(0,1,0);
        
        this->num_stack = 0;
        this->num_slice = 0;
        
    }
    
    //init from one dimension vertex array
    ThreeDimensionObject (float vertex_arr[][4], int arr_row_num)
    {
        ThreeDimensionObject();
        
        std::vector<Vertex> level_vertex;
        for (int ii = 0; ii < arr_row_num; ii++)
        {
            //float tempV[4] = {vertex_arr[ii][0],vertex_arr[ii][1],vertex_arr[ii][2],vertex_arr[ii][3]};
            Vertex tempV;
            tempV.x = vertex_arr[ii][0];
            tempV.y = vertex_arr[ii][1];
            tempV.w = vertex_arr[ii][2];
            tempV.w = vertex_arr[ii][3];
            
            level_vertex.push_back(tempV);
        }
        
        vertex_matrix.push_back(level_vertex);
    }
    
    
    //pretend vector's head connect its tail
    inline Vertex get_from_vector(std::vector<Vertex> arr, int index)
    {
        //DLOG("arr size "<<arr.size()<<"  index" <<index<<"\n";
        
        if(index>=(int)arr.size())
        {
            int ii = arr.size()%index;
            return arr.at(ii);
        }
        else if (index<0)
        {
            int ii = index;
            while (ii<0)
            {
                ii += arr.size();
            }
            return arr.at(ii);
        }
        else
        {
            return arr.at(index);
        }
    }
    
public:
    
    void gl_transfer_before_draw()
    {
        Vertex rv = get_rotation_axis();
        float ag = get_rotation_radian();
        
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        
        
        //remember!! in opengl, the transformation order is reversed
        glTranslatef(origin.x, origin.y, origin.z);
        glRotatef(RAD_TO_DEG(ag), rv.x, rv.y, rv.z);
        glScaled(x_scaler, y_scaler, z_scaler);
        
        
        //DLOG(getClassName()<<origin.toString());
    }
    
    
    void clean_after_draw()
    {
        glPopMatrix();
    }
    
#pragma mark Drawing Function
    //update the parameter and vertex accordingly
    //draw in the glut windows
    virtual void draw(int gl_mode)
    {
        gl_transfer_before_draw();
        
        for (int ii = 0; ii < vertex_matrix.size()-1; ii++)
        {
            std::vector<Vertex> level_vertex = vertex_matrix.at(ii);
            std::vector<Vertex> next_level_vertex = vertex_matrix.at((ii+1));
            
            for (int jj = 0; jj< level_vertex.size(); jj++)
            {
                Vertex v1 = get_from_vector(level_vertex, jj);
                Vertex v2 = get_from_vector(level_vertex, jj+1);
                Vertex v3 = get_from_vector(next_level_vertex, jj+1);//next_level_vertex.at(jj+1);
                Vertex v4 = get_from_vector(next_level_vertex, jj);
                
                Vertex::draw_rectangle_Vertex(gl_mode, v1, v2, v3, v4, r, g, b);
            }
        }
        
        clean_after_draw();
        
    }
    
    void draw_vertex_normal()
    {
        gl_transfer_before_draw();
        
        for (int ii = 0; ii < vertex_normal_matrix.size(); ii++)
        {
            std::vector<Vertex> level_normal_vertex = vertex_normal_matrix.at(ii);
            std::vector<Vertex> level_vertex = vertex_matrix.at(ii);
            for (int jj = 0; jj< level_normal_vertex.size(); jj++)
            {
                Vertex n1 = level_normal_vertex.at(jj);
                Vertex v1 = level_vertex.at(jj);
                
                Vertex v2 = v1.add(n1);
                
                glBegin(GL_LINES);
                float v[4] = vertex_to_arr(v1);
                float n[4] = vertex_to_arr(v2);
                glVertex4fv(v);
                glVertex4fv(n);
                
                glEnd();
            }
        }
        
        clean_after_draw();
    }
    
public:
    //after change geometry parameter of object
    //need to be called manually to update vertex
    virtual void updateVertex()
    {
        
    }
    
    //after change geometry parameter of object
    //always follow after updateVertex()
    //need to be called manually to update vertex normal
    virtual void updateVertexNormal()
    {
        
    }
    
#pragma mark modeling function
    
    virtual void transfer(float x_p, float y_p, float z_p)
    {
        //update one by one
        this->origin.transfer_self(x_p, y_p, z_p);
    }
    
    
    virtual void scale(float xs, float ys, float zs)
    {
        this->x_scaler *= xs;
        this->y_scaler *= ys;
        this->z_scaler *= zs;
    }
    
    void rotate_x(float theta)
    {
        this->origin.rotate_x_self(theta);
        this->bottom_vertex.rotate_x_self(theta);
        this->top_vertex.rotate_x_self(theta);
    }
    
    void rotate_y(float theta)
    {
        //update one by one
        this->origin.rotate_y_self(theta);
        this->bottom_vertex.rotate_y_self(theta);
        this->top_vertex.rotate_y_self(theta);
    }
    
    void rotate_z(float theta)
    {
        //update one by one
        this->origin.rotate_z_self(theta);
        this->bottom_vertex.rotate_z_self(theta);
        this->top_vertex.rotate_z_self(theta);
    }
    
    
    
    //vtx is origined from the (0,0,0)
    void rotate_around_axis(Vertex vtx, float theta)
    {
        this->origin.rotate_around_axis(vtx, theta);
        this->bottom_vertex.rotate_around_axis(vtx, theta);
        this->top_vertex.rotate_around_axis(vtx, theta);
    }
    
    
#pragma mark local cordinate functions
    
    virtual Vertex get_center_axis()
    {
        Vertex  center_axis = top_vertex.sub(this->bottom_vertex).normalize();
        return center_axis;
    }
    
    //tet the vector to
    //only rotate once to the current orientation
    virtual Vertex get_rotation_axis()
    {
        return Vertex::get_y_axis().crosss_product(get_center_axis()).normalize();
    }
    
    virtual float get_rotation_radian()
    {
        float ag = Vertex::get_y_axis().angle_between_v(get_center_axis());
        return fabsf(ag);
    }
    
    
    virtual void draw_center_axis()
    {
        gl_transfer_before_draw();
        
        Vertex  center_axis = get_center_axis();
        
        glBegin(GL_LINES);
        float v[4] = vertex_to_arr(top_vertex.add(center_axis));
        float n[4] = vertex_to_arr(bottom_vertex.sub(center_axis));
        glVertex4fv(v);
        glVertex4fv(n);
        glEnd();
        
        clean_after_draw();
    }
    
    virtual Vertex get_intersect_vertex(Vertex ray_start_point, Vertex ray_direction)
    {
        return  Vertex::get_null_vertex();
    }
    
    
    virtual void rotate_along_self_center_axis(float theta)
    {
        Vertex center_axis = get_center_axis();
        
        for (int ii = 0; ii < vertex_matrix.size(); ii++)
        {
            std::vector<Vertex> level_vertex = vertex_matrix.at(ii);
            
            for (int jj = 0; jj< level_vertex.size(); jj++)
            {
                Vertex v1 = level_vertex.at(jj);
                
                v1.transfer_self(-(this->origin.x), -(this->origin.y), -(this->origin.z));
                v1 =  v1.rotate_around_axis(center_axis, theta);
                v1.transfer_self(this->origin.x, this->origin.y, this->origin.z);
                
                level_vertex.at(jj) = v1;
            }
            vertex_matrix.at(ii) = level_vertex;
        }
    }
    
    //checkimage dimension will [64][64][4]
    void generate_check_image_and_bind_texttue(GLuint * text_name)
    {
        
        static GLubyte checkImage[64][64][4];
        
        int i,j, c;
        for (i=0; i<64; i++)
        {
            for (j=0; j<64; j++)
            {
                c = (((i&0x8)==0)^((j&0x8)==0))*255;
                checkImage[i][j][0] = (GLubyte)c;
                checkImage[i][j][1] = (GLubyte)c;
                checkImage[i][j][2] = (GLubyte)c;
                checkImage[i][j][3] = (GLubyte)255;
            }
        }
        
        glGenTextures(1, text_name);
        glBindTexture(GL_TEXTURE_2D, *text_name);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                        GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                        GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 64,
                     64, 0, GL_RGBA, GL_UNSIGNED_BYTE,
                     checkImage);
    }
    
    virtual std::string getClassName()
    {
        return std::string("3D Object");
    }
    
    virtual void printDebugInfo()
    {
        DLOG("vertex:");
        
        for (int ii = 0; ii < vertex_matrix.size(); ii++)
        {
            std::vector<Vertex> level = vertex_matrix.at(ii);
            
            for (int jj = 0; jj < level.size(); jj++)
            {
                Vertex v = level.at(jj);
                
                DLOG(v.toString());
            }
        }
        
        //DLOG("\v\v\v\vvertex normal");
        for (int ii = 0; ii < vertex_normal_matrix.size(); ii++)
        {
            std::vector<Vertex> level = vertex_normal_matrix.at(ii);
            
            for (int jj = 0; jj < level.size(); jj++)
            {
                Vertex v = level.at(jj);
                
                DLOG(v.toString());
            }
        }
    }
};


class CADObject: public ThreeDimensionObject {
    
public:
    //each vector contain usually 3,4, 5 index to vertex and texture
    //e.g 1//2 2//3 3//3. before // is the vertex index.behind // is the texture index
    //inspired by .obj file
    std::vector<std::vector<int> > face_list;
    
    int face_vertex_number;
    
    
public:
    inline std::vector<int> get_int_from_face_string(std::string s)
    {
        std::vector<int> result;
        std::string tempStr;
        
        for (int ii = 0; ii< s.size(); ii++)
        {
            char tempC = s.at(ii);
            
            
            if (isdigit(tempC))
            {
                tempStr.push_back(tempC);
            }
            else
            {
                if (tempStr.size()==0){
                    continue;
                }
                int tempI = atoi(tempStr.c_str());
                
                tempStr.erase();
                
                if (tempI!=0)
                {
                    result.push_back(tempI);
                }
            }
        }
        
        int tempI = atoi(tempStr.c_str());
        
        if (tempI!=0&&tempStr.size()!=0)
        {
            result.push_back(tempI);
        }
        return result;
    }
    
    //loading from a object
    //very imformal implememtation. no io robustness is guaranteed
    CADObject (const char * object_file_path)
    {
        this->origin = Vertex(0,0,0);
        this->bottom_vertex = Vertex(0,-1,0);
        this->top_vertex = Vertex(0,1,0);
        
        //3 is the default
        this->face_vertex_number = 3;
        
        this->num_stack = 0;
        this->num_slice = 0;
        
        
        this->x_scaler = 1;
        this->y_scaler = 1;
        this->z_scaler = 1;
        
        std::vector<Vertex> level_vertex;
        std::vector<Vertex> level_vertex_normal;
        
        std::ifstream obj_file(object_file_path);
        
        std::string line;
        while (std::getline(obj_file, line))
        {
            //ignore empty line
            if(line.size()==0){continue;}
            
            std::istringstream iss(line);
            
            static std::string comment("#");
            static std::string v("v");
            static std::string vt("vt");
            static std::string vn("vn");
            static std::string f("f");
            static std::string s("s");
            
            //spilt by space
            std::string prefix;
            getline( iss, prefix, ' ' );
            
            
            if (prefix.compare(comment)==0)
            {
                //ignore comment
                continue;
            }
            else if (prefix.compare(v)==0)
            {
                std::string xS;
                std::string yS;
                std::string zS;
                getline( iss, xS, ' ' );
                getline( iss, yS, ' ' );
                getline( iss, zS, ' ' );
                
                float x = atof(xS.c_str());
                float y = atof(yS.c_str());
                float z = atof(zS.c_str());
                
                Vertex tempV = Vertex(x,y,z,1);
                
                level_vertex.push_back(tempV);
                //DLOG(tempV.toString());
                
            }
            else if (prefix.compare(vt)==0)
            {
                std::string xS;
                std::string yS;
                std::string zS;
                getline( iss, xS, ' ' );
                getline( iss, yS, ' ' );
                getline( iss, zS, ' ' );
                
                float x = atof(xS.c_str());
                float y = atof(yS.c_str());
                float z = atof(zS.c_str());
                
                Vertex tempV = Vertex(x,y,z,1);
                
                level_vertex_normal.push_back(tempV);
            }
            else if (prefix.compare(vn)==0)
            {
                
            }
            else if (prefix.compare(f)==0)
            {
                std::vector<int> face =  get_int_from_face_string(line);
                
                if (face.size()!=0)
                {
                    for (int ii = 0; ii<face.size(); ii ++ )
                    {
                        face.at(ii) = face.at(ii)--;
                    }
                    face_list.push_back(face);
                }
            }
            else if (prefix.compare(s)==0)
            {
                
            }
            else
            {
                
            }
        }
        vertex_matrix.push_back(level_vertex);
        vertex_normal_matrix.push_back(level_vertex_normal);
    }
    
    
    virtual void draw(int gl_mode)
    {
        gl_transfer_before_draw();
        
        std::vector<Vertex> level_vertex = vertex_matrix.at(0); //for object, there is only one-dimension vertex array
        for (int ii = 0; ii < face_list.size(); ii++)
        {
            std::vector<int> face = face_list.at(ii);
            
            glBegin(gl_mode);
            
            int f_size = (int)face.size();
            int itv = (int)face.size()/face_vertex_number;
            
            //DLOG(f_size<<" "<<itv);
            
            for (int jj = 0; jj<f_size; jj += itv )
            {
                Vertex v;
                
                int tempI = face.at(jj) -1;
                if (tempI<level_vertex.size()&&tempI>=0)
                {
                    v = level_vertex.at(tempI);
                }
                else
                {
                    //DLOG(tempI);
                }
                
                if (v.isZero())
                {
                    // DLOG(f_size<<" "<<itv);
                    //DLOG(tempI);
                    continue;
                }
                
                float v_arr[4] = vertex_to_arr(v);
                
                glVertex4fv(v_arr);
            }
            glEnd();
        }
        
        clean_after_draw();
    }
    
    virtual std::string getClassName()
    {
        return std::string("off file Object");
    }
};

class Cube: public ThreeDimensionObject
{
public:
    float height;
    float width;
    float depth;
    
    
public:
    Cube (float x_p, float y_p, float z_p, float height, float width, float depth)
    {
        this->num_slice = 50;
        this->num_stack = 50;
        
        init_default_parameters(x_p, y_p, z_p);
        
        this->bottom_vertex = Vertex(0,-height/2, 0).transfer(x_p, y_p, z_p);
        this->top_vertex = Vertex(0,height/2, 0).transfer(x_p, y_p, z_p);
        
        this->height = height;
        this->width = width;
        this->depth = depth;
        
        this->updateVertex();
    }
    
    
    void updateVertex()
    {
        vertex_matrix.clear();
        std::vector<Vertex> level_vertex;
        
        float half_h = height/2;
        float half_w = width/2;
        float half_d = depth/2;
        
        
        for (int ii = -1; ii<2; ii = ii+2)
        {
            for (int jj = -1; jj<2; jj = jj+2)
            {
                for (int kk = -1; kk<2; kk = kk+2)
                {
                    Vertex  tempVertex = Vertex(ii * half_h, jj * half_w, kk * half_d, 1);
                    
                    level_vertex.push_back(tempVertex);
                }
            }
        }
        
        vertex_matrix.push_back(level_vertex);
    }
    
    void draw(int gl_mode)
    {
        
        gl_transfer_before_draw();
        
        std::vector<Vertex> level_vertex = vertex_matrix.at(0);
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(0), level_vertex.at(1), level_vertex.at(2), 0, 1, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(1), level_vertex.at(3), level_vertex.at(2), 0, 1, 0);
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(7), level_vertex.at(6), level_vertex.at(3), 1, 0, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(3), level_vertex.at(6), level_vertex.at(2), 1, 0, 0);
        
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(4), level_vertex.at(5), level_vertex.at(0), 0, 0, 1);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(1), level_vertex.at(0), level_vertex.at(5), 0, 0, 1);
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(1), level_vertex.at(7), level_vertex.at(3), 0, 0, 1);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(1), level_vertex.at(5), level_vertex.at(7), 0, 0, 1);
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(4), level_vertex.at(5), level_vertex.at(6), 1, 0, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(5), level_vertex.at(7), level_vertex.at(6), 1, 0, 0);
        
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(6), level_vertex.at(0), level_vertex.at(2), 1, 0, 1);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(6), level_vertex.at(4), level_vertex.at(0), 1, 0, 1);
        
        clean_after_draw();
    }
    
    
    virtual void scale(float xs, float ys, float zs)
    {
        ThreeDimensionObject::scale(xs, ys, zs);
        height *= ys;
        width *= xs;
        depth *= zs;
    }
    
private:Vertex get_intersect_with_origin_cube(Vertex ray_start_point, Vertex ray_direction)
    {
        float dx = ray_direction.x;
        float dy = ray_direction.y;
        float dz = ray_direction.z;
        float px = ray_start_point.x;
        float py = ray_start_point.y;
        float pz = ray_start_point.z;
        
        float t[6];
        
        t[0] = (-px - width/2)/dx;  //left face
        t[1] = (-px + width/2)/dx;  //right face
        t[2] = (-py - height/2)/dy;  //bottom face
        t[3] = (-py + height/2)/dy;  //top face
        t[4] = (-pz - depth/2)/dz;   //back face
        t[5] = (-pz + depth/2)/dz;   //front face
        
        
        std::vector<float> arr;
        
        for(int ii= 0; ii<6; ii++)
        {
            //printf("%i %f\n",ii,t[ii]);
            if(t[ii]>0)
            {
                arr.push_back(t[ii]);
            }
        }
        
        
        if(arr.size()==0)
        {
            return Vertex::get_null_vertex();
        }
        
        //sort vectot from min positive to mac positive
        std::sort(arr.begin(),arr.end());
        
        for (int ii = 0; ii<arr.size(); ii++)
        {
            float s = arr.at(ii);
            
            Vertex result = ray_start_point.add(ray_direction.scale(s, s, s));
            //DLOG("result "<<result.toString());
            
            //check if the point is inside/on the cube
            if(   result.y <= height/2 && result.y >= -height/2
               && result.x <= width/2  && result.x >= -width/2
               && result.z <= depth/2  && result.z >= -depth/2)
            {
                return result;
            }
        }
        return Vertex::get_null_vertex();
    }
    
public:    virtual Vertex get_intersect_vertex(Vertex ray_start_point, Vertex ray_direction)
    {
        Vertex rv = get_rotation_axis();
        float ag = get_rotation_radian();
        
        ray_start_point.transfer_self(-origin.x, -origin.y, -origin.z);
        ray_start_point = ray_start_point.rotate_around_axis(rv, -ag);
        
        ray_direction.normalize_self();
        ray_direction = ray_direction.rotate_around_axis(rv, -ag);
        
        Vertex v = get_intersect_with_origin_cube(ray_start_point, ray_direction);
        if(v.isNull())
        {
            return v;
        }
        v.transfer_self(origin.x, origin.y, origin.z);
        v = v.rotate_around_axis(rv, ag);
        return v;
    }
    
    virtual std::string getClassName()
    {
        return std::string("Cube Object");
    }
};

class House: public ThreeDimensionObject {
    
    float height;
    float width;
    float depth;
    float roof;
public:
    House (float x_p, float y_p, float z_p, float height, float width, float depth, float roof)
    {
        this->num_slice = 50;
        this->num_stack = 50;
        
        init_default_parameters(x_p, y_p, z_p);
        
        this->bottom_vertex = Vertex(0,-height/2, 0);
        this->top_vertex = Vertex(0,height/2, 0);
        
        this->height = height;
        this->width = width;
        this->depth = depth;
        
        this->updateVertex();
    }
    
    void updateVertex()
    {
        vertex_matrix.clear();
        std::vector<Vertex> level_vertex;
        
        float half_h = height/2;
        float half_w = width/2;
        float half_d = depth/2;
        
        
        for (int ii = -1; ii<2; ii = ii+2)
        {
            for (int jj = -1; jj<2; jj = jj+2)
            {
                for (int kk = -1; kk<2; kk = kk+2)
                {
                    Vertex  tempVertex = Vertex(ii * half_w, jj * half_h, kk * half_d, 1);
                    
                    //vertex is generate based on the case when the object's center is in origin
                    //transfer to correct center
                    level_vertex.push_back(tempVertex);
                }
            }
        }
        
        Vertex roof_top = Vertex(0,height+roof,0,1);
        level_vertex.push_back(roof_top);
        
        vertex_matrix.push_back(level_vertex);
    }
    
    void draw(int gl_mode)
    {
        gl_transfer_before_draw();
        
        std::vector<Vertex> level_vertex = vertex_matrix.at(0);
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(0), level_vertex.at(1), level_vertex.at(2), 0, 1, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(1), level_vertex.at(3), level_vertex.at(2), 0, 1, 0);
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(7), level_vertex.at(6), level_vertex.at(3), 1, 0, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(3), level_vertex.at(6), level_vertex.at(2), 1, 0, 0);
        
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(4), level_vertex.at(5), level_vertex.at(0), 0, 0, 1);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(1), level_vertex.at(0), level_vertex.at(5), 0, 0, 1);
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(1), level_vertex.at(7), level_vertex.at(3), 0, 0, 1);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(1), level_vertex.at(5), level_vertex.at(7), 0, 0, 1);
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(4), level_vertex.at(5), level_vertex.at(6), 1, 0, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(5), level_vertex.at(7), level_vertex.at(6), 1, 0, 0);
        
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(6), level_vertex.at(0), level_vertex.at(2), 1, 0, 1);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(6), level_vertex.at(4), level_vertex.at(0), 1, 0, 1);
        
        
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(8), level_vertex.at(3), level_vertex.at(7), 1, 0, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(8), level_vertex.at(7), level_vertex.at(6), 1, 0, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(8), level_vertex.at(6), level_vertex.at(3), 1, 0, 0);
        Vertex::draw_triangle_Vertex(gl_mode, level_vertex.at(8), level_vertex.at(2), level_vertex.at(3), 1, 0, 0);
        
        clean_after_draw();
    }
    
    virtual void scale(float xs, float ys, float zs)
    {
        ThreeDimensionObject::scale(xs, ys, zs);
        height *= ys;
        width *= xs;
        depth *= zs;
        roof *=ys;
    }
    
    virtual std::string getClassName()
    {
        return std::string("House Object");
    }
    
    //i reuse code from cube wihtout changing. it is wrong, but only serve as an approximation
private:Vertex get_intersect_with_origin_house(Vertex ray_start_point, Vertex ray_direction)
    {
        float dx = ray_direction.x;
        float dy = ray_direction.y;
        float dz = ray_direction.z;
        float px = ray_start_point.x;
        float py = ray_start_point.y;
        float pz = ray_start_point.z;
        
        float t[6];
        
        t[0] = (-px - width/2)/dx;  //left face
        t[1] = (-px + width/2)/dx;  //right face
        t[2] = (-py - height/2)/dy;  //bottom face
        t[3] = (-py + height/2)/dy;  //top face
        t[4] = (-pz - depth/2)/dz;   //back face
        t[5] = (-pz + depth/2)/dz;   //front face
        
        
        std::vector<float> arr;
        
        for(int ii= 0; ii<6; ii++)
        {
            //printf("%i %f\n",ii,t[ii]);
            if(t[ii]>0)
            {
                arr.push_back(t[ii]);
            }
        }
        
        
        if(arr.size()==0)
        {
            return Vertex::get_null_vertex();
        }
        
        //sort vectot from min positive to mac positive
        std::sort(arr.begin(),arr.end());
        
        for (int ii = 0; ii<arr.size(); ii++)
        {
            float s = arr.at(ii);
            
            Vertex result = ray_start_point.add(ray_direction.scale(s, s, s));
            
            //check if the point is inside/on the cube
            if(   result.y <= height/2 && result.y >= -height/2
               && result.x <= width/2  && result.x >= -width/2
               && result.z <= depth/2  && result.z >= -depth/2)
            {
                return result;
            }
        }
        return Vertex::get_null_vertex();
    }
    
public:    virtual Vertex get_intersect_vertex(Vertex ray_start_point, Vertex ray_direction)
    {
        Vertex rv = get_rotation_axis();
        float ag = get_rotation_radian();
        
        ray_start_point.transfer_self(-origin.x, -origin.y, -origin.z);
        ray_start_point = ray_start_point.rotate_around_axis(rv, -ag);
        
        ray_direction.normalize_self();
        ray_direction = ray_direction.rotate_around_axis(rv, -ag);
        
        Vertex v = get_intersect_with_origin_house(ray_start_point, ray_direction);
        if(v.isNull())
        {
            return v;
        }
        v.transfer_self(origin.x, origin.y, origin.z);
        v = v.rotate_around_axis(rv, ag);
        return v;
    }
};



class Sphere : public ThreeDimensionObject
{
    public :
    float radius;
    
public:Sphere (float x_p, float y_p, float z_p, float radius_p)
    {
        init_default_parameters(x_p, y_p, z_p);
        
        this->bottom_vertex = Vertex(0,-radius_p/2, 0);
        this->top_vertex = Vertex(0,radius_p/2, 0);
        
        this->num_slice = 30;
        this->num_stack = 30;
        this->radius = radius_p;
        this->updateVertex();
    }
    
    void updateVertex()
    {
        vertex_matrix.clear();
        
        float vertical_interval = radius*2.0/num_stack;  // in length
        float slice_interval = M_PI*2.0/num_slice;   // in radian
        
        for (int ii = 0; ii <= num_stack; ii++)
        {
            float temp_y = -radius + ii * vertical_interval;
            
            std::vector<Vertex> level_vertex;
            for (int jj = 0; jj< num_slice; jj++)
            {
                float delta = slice_interval * jj;
                
                float level_radius = sqrtf(radius * radius - temp_y * temp_y);
                
                float temp_x = level_radius * cos(delta);
                float temp_z = level_radius * sin(delta);
                
                Vertex  tempVertex = Vertex(temp_x, temp_y, temp_z, 1);
                //vertex is generate based on the case when the object's center is in origin
                //transfer to correct center
                level_vertex.push_back(tempVertex);
            }
            
            vertex_matrix.push_back(level_vertex);
        }
    }
    
    void updateVertexNormal()
    {
        vertex_normal_matrix.clear();
        
        std::vector<Vertex> bottom_vertex_normals;
        std::vector<Vertex> bottom_vertex = vertex_matrix.at(0);
        for (int ii = 0; ii < bottom_vertex.size();ii++)
        {
            bottom_vertex_normals.push_back(Vertex(0, -radius/radius, 0, 1));
        }
        
        vertex_normal_matrix.push_back(bottom_vertex_normals);
        
        for (int ii = 1; ii < vertex_matrix.size()-1; ii++)
        {
            std::vector<Vertex> previous_level_vertex = vertex_matrix.at(ii-1);
            std::vector<Vertex> level_vertex = vertex_matrix.at(ii);
            std::vector<Vertex> next_level_vertex = vertex_matrix.at((ii+1));
            
            std::vector<Vertex> level_normal ;
            
            for (int jj = 0; jj< level_vertex.size(); jj++)
            {
                
                //get up, down, left, right vertex get surfuce normal and then vertex normal
                Vertex c_v = level_vertex.at(jj);
                Vertex up_v = next_level_vertex.at(jj);
                Vertex down_v = previous_level_vertex.at(jj);
                Vertex left_v = get_from_vector(level_vertex, jj-1);
                Vertex right_v = get_from_vector(level_vertex, jj+1);
                
                Vertex up_vector = up_v.sub(c_v);
                Vertex down_vector = down_v.sub(c_v);
                Vertex right_vector = right_v.sub(c_v);
                Vertex left_vector = left_v.sub(c_v);
                
                
                //get surface normal
                Vertex face_normal1 =  up_vector.crosss_product(left_vector).normalize();
                Vertex face_normal2 = left_vector.crosss_product( down_vector).normalize();
                Vertex face_normal3 = down_vector.crosss_product(right_vector).normalize();
                Vertex face_normal4 = right_vector.crosss_product(up_vector).normalize();
                
                
                //sum them up
                Vertex temp_vertex_normal =  face_normal1.add(face_normal2).add(face_normal3).add(face_normal4);                //...then divide
                Vertex temp_vertex_normal2 = temp_vertex_normal.scale(-0.25, -0.25, -0.25);
                
                Vertex temp_vertex_normal3 = temp_vertex_normal2.normalize();
                
                level_normal.push_back(temp_vertex_normal3);
            }
            
            vertex_normal_matrix.push_back(level_normal);
        }
        
        std::vector<Vertex> ceil_vertex_normals;
        std::vector<Vertex> ceil_vertex = vertex_matrix.at(vertex_matrix.size()-1);
        for (int ii = 0; ii < ceil_vertex.size();ii++)
        {
            //get_vertex(0, radius, 0, 1)
            ceil_vertex_normals.push_back(Vertex(0, radius/radius, 0, 1) );
        }
        
        vertex_normal_matrix.push_back(ceil_vertex_normals);
    }
    
private:
    Vertex get_intersect_with_origin_sphere(Vertex ray_start_point, Vertex ray_direction)
    {
        
        //DLOG(radius);
        ray_direction.normalize_self();
        
        float dx = ray_direction.x;
        float dy = ray_direction.y;
        float dz = ray_direction.z;
        float px = ray_start_point.x;
        float py = ray_start_point.y;
        float pz = ray_start_point.z;
        
        float a = dx*dx + dy*dy + dz*dz;
        float b = 2*(px*dx + py*dy + pz*dz);
        float c = px*px + py*py + pz*pz - radius*radius;
        
        float t1;
        float t2;
        
        if(!solve_cube_equation(a, b, c, &t1, &t2))
        {
            return Vertex::get_null_vertex();
        }
        
        float s = MIN(t1, t2);
        
        Vertex result = ray_start_point.add(ray_direction.scale(s, s, s));
        return result;
    }
    
public:
    virtual Vertex get_intersect_vertex(Vertex ray_start_point, Vertex ray_direction)
    {
        
        ray_start_point.transfer_self(-origin.x, -origin.y, -origin.z);
        
        Vertex v = get_intersect_with_origin_sphere(ray_start_point, ray_direction);
        if(v.isNull())
        {
            return v;
        }
        
        v.transfer_self(origin.x, origin.y, origin.z);
        
        return v;
    }
    
    virtual void scale(float xs, float ys, float zs)
    {
        ThreeDimensionObject::scale(xs, ys, zs);
        radius *= xs;
        
    }
    
    //update the parameter and vertex accordingly
    //draw in the glut windows
    virtual void draw(int gl_mode)
    {
        ThreeDimensionObject::draw(gl_mode);
    }
    
    virtual std::string getClassName()
    {
        return std::string("Sphere Object");
    }
    
};


class Cylinder : public ThreeDimensionObject
{
public:
    float radius;
    float height;
    
public:Cylinder(float x_p, float y_p, float z_p, float radius_p, float height_p)
    {
        this->num_slice = 30;
        this->num_stack = 30;
        
        init_default_parameters(x_p, y_p, z_p);
        
        this->bottom_vertex = Vertex(0,-height_p/2, 0);
        this->top_vertex = Vertex(0,height_p/2, 0);
        this->radius = radius_p;
        this->height = height_p;
        
        this->updateVertex();
    }
    
    
    void updateVertex()
    {
        vertex_matrix.clear();
        
        float vertical_interval = height/num_stack;  // in length
        float slice_interval = M_PI*2.0/num_slice;   // in radian
        
        for (int ii = 0; ii <= num_stack; ii++)
        {
            float temp_y = -height/2.0 + ii * vertical_interval;
            
            std::vector<Vertex> level_vertex;
            for (int jj = 0; jj< num_slice; jj++)
            {
                float delta = slice_interval * jj;
                
                float temp_x = radius * cos(delta);
                float temp_z = radius * sin(delta);
                
                Vertex  tempVertex = Vertex(temp_x, temp_y, temp_z, 1);
                level_vertex.push_back(tempVertex);
            }
            
            vertex_matrix.push_back(level_vertex);
        }
    }
    
    void updateVertexNormal()
    {
        vertex_normal_matrix.clear();
        
        std::vector<Vertex> bottom_vertex_normals;
        std::vector<Vertex> bottom_vertex = vertex_matrix.at(0);
        for (int ii = 0; ii < bottom_vertex.size();ii++)
        {
            
            bottom_vertex_normals.push_back(Vertex(0, -radius, 0, 1).normalize());
        }
        
        vertex_normal_matrix.push_back(bottom_vertex_normals);
        
        
        for (int ii = 1; ii < vertex_matrix.size()-1; ii++)
        {
            std::vector<Vertex> previous_level_vertex = vertex_matrix.at(ii-1);
            std::vector<Vertex> level_vertex = vertex_matrix.at(ii);
            std::vector<Vertex> next_level_vertex = vertex_matrix.at((ii+1));
            
            std::vector<Vertex> level_normal ;
            
            for (int jj = 0; jj< level_vertex.size(); jj++)
            {
                //get up, down, left, right vertex get surfuce normal and then vertex normal
                Vertex c_v = level_vertex.at(jj);
                Vertex up_v = next_level_vertex.at(jj);
                Vertex down_v = previous_level_vertex.at(jj);
                Vertex left_v = get_from_vector(level_vertex, jj-1);
                Vertex right_v = get_from_vector(level_vertex, jj+1);
                
                Vertex up_vector = up_v.sub(c_v);
                Vertex down_vector = down_v.sub(c_v);
                Vertex right_vector = right_v.sub(c_v);
                Vertex left_vector = left_v.sub(c_v);
                
                //get surface normal
                Vertex face_normal1 =  up_vector.crosss_product(left_vector).normalize();
                Vertex face_normal2 = left_vector.crosss_product( down_vector).normalize();
                Vertex face_normal3 = down_vector.crosss_product(right_vector).normalize();
                Vertex face_normal4 = right_vector.crosss_product(up_vector).normalize();
                
                //sum them up
                Vertex temp_vertex_normal =  face_normal1.add(face_normal2).add(face_normal3).add(face_normal4);                //...then divide
                Vertex temp_vertex_normal2 = temp_vertex_normal.scale(-0.25, -0.25, -0.25);
                
                Vertex temp_vertex_normal3 = temp_vertex_normal2.normalize();
                
                level_normal.push_back(temp_vertex_normal3);
            }
            
            vertex_normal_matrix.push_back(level_normal);
        }
        
        std::vector<Vertex> ceil_vertex_normals;
        std::vector<Vertex> ceil_vertex = vertex_matrix.at(vertex_matrix.size()-1);
        for (int ii = 0; ii < ceil_vertex.size();ii++)
        {
            ceil_vertex_normals.push_back(Vertex(0, radius, 0, 1).normalize() );
        }
        
        vertex_normal_matrix.push_back(ceil_vertex_normals);
    }
    
    
    void draw(int gl_mode)
    {
        ThreeDimensionObject::draw(gl_mode);
        
        //default did not draw the bottom and up
        gl_transfer_before_draw();
        
        std::vector<Vertex> bottom_vertex = vertex_matrix.at(0);
        std::vector<Vertex> ceil_vertex = vertex_matrix.at(vertex_matrix.size()-1);
        
        for (int ii = 0; ii < bottom_vertex.size()-1;ii++)
        {
            Vertex v1 = bottom_vertex.at(ii);
            Vertex v2 = bottom_vertex.at(ii+1);
            
            Vertex::draw_triangle_Vertex(gl_mode, this->bottom_vertex, v1, v2, 1, 1, 0);
        }
        
        for (int ii = 0; ii < ceil_vertex.size()-1;ii++)
        {
            Vertex v1 = ceil_vertex.at(ii);
            Vertex v2 = ceil_vertex.at(ii+1);
            
            Vertex::draw_triangle_Vertex(gl_mode, this->top_vertex, v1, v2, 1, 1, 0);
        }
        
        clean_after_draw();
    }
    
    virtual void scale(float xs, float ys, float zs)
    {
        ThreeDimensionObject::scale(xs, ys, zs);
        height *= ys;
        radius *=xs;
    }
    
private:
    Vertex get_intersect_with_origin_cynlinder(Vertex ray_start_point, Vertex ray_direction)
    {
        ray_direction.normalize_self();
        
        float dx = ray_direction.x;
        float dz = ray_direction.z;
        float px = ray_start_point.x;
        float pz = ray_start_point.z;
        
        float a = dx*dx  + dz*dz;
        float b = 2*(px*dx + pz*dz);
        float c = px*px + pz*pz - radius*radius;
        
        float t1;
        float t2;
        
        if(!solve_cube_equation(a, b, c, &t1, &t2))
        {
            return Vertex::get_null_vertex();
        }
        
        float s = MIN(t1, t2);
        
        Vertex result = ray_start_point.add(ray_direction.scale(s, s, s));
        
        if (result.y<-height/2||result.y>height/2)
        {
            return Vertex::get_null_vertex();
        }
        
        //DLOG("result"<< result.toString()<<"\n");
        return result;
    }
    
public:
    virtual Vertex get_intersect_vertex(Vertex ray_start_point, Vertex ray_direction)
    {
        Vertex rv = get_rotation_axis();
        float ag = get_rotation_radian();
        
        ray_start_point.transfer_self(-origin.x, -origin.y, -origin.z);
        ray_start_point = ray_start_point.rotate_around_axis(rv, -ag);
        
        ray_direction = ray_direction.rotate_around_axis(rv, -ag);
        
        Vertex v = get_intersect_with_origin_cynlinder(ray_start_point, ray_direction);
        
        if(v.isNull())
        {
            return v;
        }
        
        v.transfer_self(origin.x, origin.y, origin.z);
        v = v.rotate_around_axis(rv, ag);
        return v;
    }
    
    virtual std::string getClassName()
    {
        return std::string("Cylinder Object");
    }
};

class Cone: public ThreeDimensionObject {
    
    float radius;   //the bottom
    float height;
    
public:Cone(float x_p, float y_p, float z_p, float radius_p, float height_p)
    {
        this->num_slice = 30;
        this->num_stack = 30;
        
        init_default_parameters(x_p, y_p, z_p);
        
        this->bottom_vertex = Vertex(0,-height_p/2, 0);
        this->top_vertex = Vertex(0,height_p/2, 0);
        this->radius = radius_p;
        this->height = height_p;
        
        this->updateVertex();
    }
    
    
    void updateVertex()
    {
        vertex_matrix.clear();
        
        float vertical_interval = height/num_stack;  // in length
        float slice_interval = M_PI*2.0/num_slice;   // in radian
        
        for (int ii = 0; ii <= num_stack; ii++)
        {
            float temp_y = -height/2.0 + ii * vertical_interval;
            float temp_h = height - ii * vertical_interval;
            
            std::vector<Vertex> level_vertex;
            for (int jj = 0; jj< num_slice; jj++)
            {
                float delta = slice_interval * jj;
                
                float temp_radius = radius * temp_h / height;
                
                float temp_x = temp_radius * cos(delta);
                float temp_z = temp_radius * sin(delta);
                
                Vertex  tempVertex = Vertex(temp_x, temp_y, temp_z, 1);
                level_vertex.push_back(tempVertex);
            }
            vertex_matrix.push_back(level_vertex);
        }
    }
    
    void updateVertexNormal()
    {
        vertex_normal_matrix.clear();
        
        std::vector<Vertex> bottom_vertex_normals;
        std::vector<Vertex> bottom_vertex = vertex_matrix.at(0);
        for (int ii = 0; ii < bottom_vertex.size();ii++)
        {
            //get_vertex(0, -radius, 0, 1)
            
            bottom_vertex_normals.push_back(Vertex(0, -radius, 0, 1).normalize());
        }
        
        vertex_normal_matrix.push_back(bottom_vertex_normals);
        
        
        for (int ii = 1; ii < vertex_matrix.size()-1; ii++)
        {
            std::vector<Vertex> previous_level_vertex = vertex_matrix.at(ii-1);
            std::vector<Vertex> level_vertex = vertex_matrix.at(ii);
            std::vector<Vertex> next_level_vertex = vertex_matrix.at((ii+1));
            
            std::vector<Vertex> level_normal ;
            
            for (int jj = 0; jj< level_vertex.size(); jj++)
            {
                //get up, down, left, right vertex get surfuce normal and then vertex normal
                Vertex c_v = level_vertex.at(jj);
                Vertex up_v = next_level_vertex.at(jj);
                Vertex down_v = previous_level_vertex.at(jj);
                Vertex left_v = get_from_vector(level_vertex, jj-1);
                Vertex right_v = get_from_vector(level_vertex, jj+1);
                
                Vertex up_vector = up_v.sub(c_v);
                Vertex down_vector = down_v.sub(c_v);
                Vertex right_vector = right_v.sub(c_v);
                Vertex left_vector = left_v.sub(c_v);
                
                //get surface normal
                Vertex face_normal1 =  up_vector.crosss_product(left_vector).normalize();
                Vertex face_normal2 = left_vector.crosss_product( down_vector).normalize();
                Vertex face_normal3 = down_vector.crosss_product(right_vector).normalize();
                Vertex face_normal4 = right_vector.crosss_product(up_vector).normalize();
                
                //sum them up
                Vertex temp_vertex_normal =  face_normal1.add(face_normal2).add(face_normal3).add(face_normal4);                //...then divide
                Vertex temp_vertex_normal2 = temp_vertex_normal.scale(-0.25, -0.25, -0.25);
                
                Vertex temp_vertex_normal3 = temp_vertex_normal2.normalize();
                
                level_normal.push_back(temp_vertex_normal3);
            }
            
            vertex_normal_matrix.push_back(level_normal);
        }
        
        std::vector<Vertex> ceil_vertex_normals;
        std::vector<Vertex> ceil_vertex = vertex_matrix.at(vertex_matrix.size()-1);
        for (int ii = 0; ii < ceil_vertex.size();ii++)
        {
            //get_vertex(0, radius, 0, 1)
            ceil_vertex_normals.push_back(Vertex(0, radius, 0, 1).normalize());
        }
        
        vertex_normal_matrix.push_back(ceil_vertex_normals);
    }
    
    
    void draw(int gl_mode)
    {
        ThreeDimensionObject::draw(gl_mode);
        
        gl_transfer_before_draw();
        
        //default did not draw the bottom and up
        std::vector<Vertex> bottom_vertex = vertex_matrix.at(0);
        
        
        for (int ii = 0; ii < bottom_vertex.size()-1;ii++)
        {
            Vertex v1 = bottom_vertex.at(ii);
            Vertex v2 = bottom_vertex.at(ii+1);
            
            Vertex::draw_triangle_Vertex(gl_mode, this->bottom_vertex, v1, v2, 1, 1, 0);
        }
        glPopMatrix();
    }
    
    virtual void scale(float xs, float ys, float zs)
    {
        ThreeDimensionObject::scale(xs, ys, zs);
        height *= ys;
        radius *=xs;
    }
    
private:
    Vertex get_intersect_with_origin_cone(Vertex ray_start_point, Vertex ray_direction)
    {
        ray_direction.normalize_self();
        
        float dx = ray_direction.x;
        float dy = ray_direction.y;
        float dz = ray_direction.z;
        float px = ray_start_point.x;
        float py = ray_start_point.y;
        float pz = ray_start_point.z;
        
        float r2 = radius * radius;
        float h2 = height * height;
        float a = dx*dx  + dz*dz - r2*dy*dy/h2;
        float b = 2*(px*dx + pz*dz)-r2/h2*2*py*dy+2*r2/height*dy;
        float c = px*px + pz*pz  - r2 - r2/h2*py*py+r2/h2*2*height*py;
        
        float t1;
        float t2;
        
        if(!solve_cube_equation(a, b, c, &t1, &t2))
        {
            return Vertex::get_null_vertex();
        }
        
        float s = MIN(t1, t2);
        
        Vertex result = ray_start_point.add(ray_direction.scale(s, s, s));
        
        if (result.y<-height/2||result.y>height/2)
        {
            return Vertex::get_null_vertex();
        }
        return result;
    }
    
public:
    virtual Vertex get_intersect_vertex(Vertex ray_start_point, Vertex ray_direction)
    {
        Vertex rv = get_rotation_axis();
        float ag = get_rotation_radian();
        
        ray_start_point.transfer_self(-origin.x, -origin.y, -origin.z);
        ray_start_point = ray_start_point.rotate_around_axis(rv, -ag);
        
        ray_direction = ray_direction.rotate_around_axis(rv, -ag);
        
        Vertex v = get_intersect_with_origin_cone(ray_start_point, ray_direction);
        if(v.isNull())
        {
            return v;
        }
        
        v.transfer_self(origin.x, origin.y, origin.z);
        v = v.rotate_around_axis(rv, ag);
        return v;
    }
    
    virtual std::string getClassName()
    {
        return std::string("Cone Object");
    }
};


class Torus: public ThreeDimensionObject
{
    
    //continue to use the varible inheritate from parent but imply differenrt meaning
    //    int num_stack;  //Number of circles across the tube
    //    int num_slice; //Number of circles along the tube
    //
    //    //num_stack * num_stack * 4
    //    std::vector<std::vector<Vertex>> vertex_matrix;
    
public:
    float big_radius;
    float tube_radius;
    
public:
    Torus(float x_p, float y_p, float z_p,  float big_radius, float tube_radius)
    {
        this->num_slice = 50;
        this->num_stack = 50;
        
        init_default_parameters(x_p, y_p, z_p);
        
        this->bottom_vertex = Vertex(0,-tube_radius/2, 0);
        this->top_vertex = Vertex(0,tube_radius/2, 0);
        
        this->big_radius = big_radius;
        this->tube_radius = tube_radius;
        
        this->updateVertex();
    }
    
    
    void updateVertex()
    {
        vertex_matrix.clear();
        
        float tube_interval = M_PI*2.0/num_stack;
        float big_interval = M_PI*2.0/num_slice;
        
        for (int ii = 0; ii< num_slice; ii++)
        {
            float delta = big_interval * ii;
            
            std::vector<Vertex> tube_vertex;
            
            for (int jj = 0; jj < num_stack; jj++)
            {
                float theta = tube_interval * jj;
                
                float tempf = (big_radius + tube_radius * cos(theta));
                
                float x = tempf *sin(delta);
                float y = tube_radius * sin(theta);
                float z = tempf *cos(delta);
                
                Vertex tempV  = Vertex(x,y,z, 1);
                
                tube_vertex.push_back(tempV);
            }
            
            vertex_matrix.push_back(tube_vertex);
        }
    }
    
    void updateVertexNormal()
    {
        vertex_normal_matrix.clear();
        
        for (int ii = 0; ii < vertex_matrix.size(); ii++)
        {
            std::vector<Vertex> previous_level_vertex = vertex_matrix.at((ii==0?vertex_matrix.size()-1:ii-1));
            
            std::vector<Vertex> level_vertex = vertex_matrix.at(ii);
            
            std::vector<Vertex> next_level_vertex = vertex_matrix.at((ii == vertex_matrix.size()-1?0:ii+1));
            
            std::vector<Vertex> level_normal ;
            
            for (int jj = 0; jj< level_vertex.size(); jj++)
            {
                //get up, down, left, right vertex get surfuce normal and then vertex normal
                Vertex c_v = level_vertex.at(jj);
                Vertex up_v = next_level_vertex.at(jj);
                Vertex down_v = previous_level_vertex.at(jj);
                Vertex left_v = get_from_vector(level_vertex, jj-1);
                Vertex right_v = get_from_vector(level_vertex, jj+1);
                
                Vertex up_vector = up_v.sub(c_v);
                Vertex down_vector = down_v.sub(c_v);
                Vertex right_vector = right_v.sub(c_v);
                Vertex left_vector = left_v.sub(c_v);
                
                //get surface normal
                Vertex face_normal1 =  up_vector.crosss_product(left_vector).normalize();
                Vertex face_normal2 = left_vector.crosss_product( down_vector).normalize();
                Vertex face_normal3 = down_vector.crosss_product(right_vector).normalize();
                Vertex face_normal4 = right_vector.crosss_product(up_vector).normalize();
                
                //sum them up
                Vertex temp_vertex_normal =  face_normal1.add(face_normal2).add(face_normal3).add(face_normal4);
                //...then divide
                Vertex temp_vertex_normal2 = temp_vertex_normal.scale(-0.25, -0.25, -0.25);
                
                Vertex temp_vertex_normal3 = temp_vertex_normal2.normalize();
                
                level_normal.push_back(temp_vertex_normal3);
            }
            
            vertex_normal_matrix.push_back(level_normal);
        }
    }
    
    void draw(int gl_mode)
    {
        ThreeDimensionObject::draw(gl_mode);
        
        gl_transfer_before_draw();
        
        std::vector<Vertex> head_vertex = vertex_matrix.at(0);
        std::vector<Vertex> tail_vertex = vertex_matrix.at(vertex_matrix.size()-1);
        
        for (int ii = 0; ii < head_vertex.size(); ii++)
        {
            Vertex v1 = get_from_vector(tail_vertex, ii);
            Vertex v2 = get_from_vector(head_vertex, ii);
            Vertex v3 = get_from_vector(head_vertex, ii+1);//next_level_vertex.at(jj+1);
            Vertex v4 = get_from_vector(tail_vertex, ii+1);
            
            Vertex::draw_rectangle_Vertex(gl_mode, v1, v2, v3, v4, 1, 0, 0);
        }
        
        clean_after_draw();
    }
    
    virtual void scale(float xs, float ys, float zs)
    {
        ThreeDimensionObject::scale(xs, ys, zs);
        big_radius *= ys;
        tube_radius *= xs;
    }
    
    virtual std::string getClassName()
    {
        return std::string("Torus Object");
    }
};



class HJYCamera {
    
public:
    Vertex position ;
    Vertex lookAt;
    Vertex up;
    float left;
    float right ;
    float near;
    float far;
    float top;
    float bottom;
    float fovy;
    float apsect;
    
private:
	HJYCamera()
   	{
        position = Vertex(0, 0, 50);
        lookAt = Vertex(0, 0, 0);
        up = Vertex(0, 1, 0);
        left = -10;
        right = 10;
        near = 0.1;
        far = 100;
        top = 10;
        bottom = -10;
        fovy = 50;
        apsect = 1;
    }
    
public:
    //singleton
    static HJYCamera* getCamera()
    {
        static HJYCamera* instance = new HJYCamera();
        return instance;
    }
    
    Vertex get_lookVector()
    {
        return lookAt.sub(position).normalize();
    }
    
    Vertex get_w_axis()
    {
        return get_lookVector().scale(-1, -1, -1);
    }
    
    Vertex get_v_axis()
    {
        Vertex w = get_w_axis();
        Vertex v;
        Vertex temp_up = up;
        if (up.isParrell(w))
        {
            temp_up = Vertex(0, 1, 0);
        }
        
        float dr = temp_up.dot_product(w);
        v = temp_up.sub( w.scale(dr,dr,dr));
        v.normalize_self();
        
        return v;
    }
    
    Vertex get_u_axis()
    {
        Vertex w = get_w_axis();
        Vertex v = get_v_axis();
        return v.crosss_product(w);
    }
    
    //the same function glutLookAt
    void look()
    {
        gluLookAt(position.x, position.y, position.z,lookAt.x, lookAt.y, lookAt.z, up.x, up.y, up.z);
    }
    
    //the samt function as gluortho
    void ortho()
    {
        glOrtho(left, right, bottom, top, near, far);
    }
    
    void perspective()
    {
        gluPerspective(fovy, apsect, near, far);
    }
    
    void move_near_and_far_plane(float near_dist, float far_dist)
    {
        //DLOG(near_dist<<" "<<far_dist);
        //near plane
        float tempNear = near + near_dist;
        if (tempNear< far)
        {
            near = tempNear;
        }
        
        //far plane
        float tempFar = far + far_dist;
        if (tempFar>near)
        {
            far = tempFar;
        }
    }
    
    void zoom(float percentage)
    {
        fovy *=percentage;
        
        left *=percentage;
        right *=percentage;
        top *= percentage;
        bottom *=percentage;
    }
    
    void move (float x_p, float y_p, float z_p)
    {
        position.transfer_self(x_p, y_p, z_p);
        lookAt.transfer_self(x_p, y_p, z_p);
        
        look();
    }
    
    void move_local_cord(float u_p, float v_p, float w_p)
    {
        Vertex w = get_w_axis();
        Vertex v = get_v_axis();
        Vertex u = get_u_axis();
        
        position.transfer_self(u_p*u.x, u_p*u.y, u_p*u.z);
        position.transfer_self(v_p*v.x, v_p*v.y, v_p*v.z);
        position.transfer_self(w_p*w.x, w_p*w.y, w_p*w.z);
        
        lookAt.transfer_self(u_p*u.x, u_p*u.y, u_p*u.z);
        lookAt.transfer_self(v_p*v.x, v_p*v.y, v_p*v.z);
        lookAt.transfer_self(w_p*w.x, w_p*w.y, w_p*w.z);
        
        look();
    }
    
    void rotate_aound_w(float radian)
    {
        Vertex w = get_w_axis();
        
        up.transfer_self(-position.x, -position.y, -position.z);
        up = up.rotate_around_axis(w, radian);
        up.transfer_self(position.x, position.y, position.z);
        
        look();
    }
    
    void rotate_aound_v(float radian)
    {
        Vertex v = get_v_axis();
        
        lookAt.transfer_self(-position.x, -position.y, -position.z);
        lookAt = lookAt.rotate_around_axis(v, radian);
        lookAt.transfer_self(position.x, position.y, position.z);
        
        look();
    }
    
    void rotate_aound_u(float radian)
    {
        Vertex u = get_u_axis();
        
        //may still have bug
        
        lookAt.transfer_self(-position.x, -position.y, -position.z);
        lookAt = lookAt.rotate_around_axis(u, radian);
        //lookAt.rotate_x_self(radian);
        lookAt.transfer_self(position.x, position.y, position.z);
        
        up.transfer_self(-position.x, -position.y, -position.z);
        up = up.rotate_around_axis(u, radian);
        //up.rotate_x_self(radian);
        up.transfer_self(position.x, position.y, position.z);
        
        if (up.y<0)
        {
            up.y = -up.y;
        }
        
        look();
    }
    
    void reset()
    {
        position = Vertex(4, 5, 30);
        lookAt = Vertex(0, 0, 0);
        up = Vertex(0, 1, 0);
    }
    
    Vertex get_ray_direction(int mouse_x, int mouse_y)
    {
        double modelViewMatrix[16];
        double projMatrix[16];
        int viewport[4];
        
        Vertex result = Vertex();
        
        // first we need to get the modelview matrix, the projection matrix, and the viewport
        glGetDoublev(GL_MODELVIEW_MATRIX, modelViewMatrix);
        glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
        glGetIntegerv(GL_VIEWPORT, viewport);
        
        //the windows y cord is top-down
        //we need to down-top
        mouse_y = viewport[3]-mouse_y;
        
        // gluUnProject with a Z value of 1 will find the point on the far clipping plane
        // corresponding the the mouse click. This is not the same as the vector
        double clickPoint[3];
        gluUnProject(mouse_x, mouse_y, 1.0, modelViewMatrix, projMatrix, viewport, &clickPoint[0], &clickPoint[1], &clickPoint[2]);
        
        result.x = clickPoint[0];
        result.y = clickPoint[1];
        result.z = clickPoint[2];
        
        result =  result.sub(position).normalize();
        return result;
    }
    
    void change_aspect_ratio(float percentage)
    {
        //for perspective
        fovy *= percentage;
        
        //for orth
        left *= percentage;
        right*= percentage;
    }
    
    void change_height_angle(float percentage)
    {
        //for orth
        apsect*= percentage;
        
        //for perspective
        bottom *= percentage;
        top *= percentage;
    }
    
    void print()
    {
        std::cout<<"\nCam Position: "<<position.toString()<<
        "\nLookAt:       "<<lookAt.toString()<<
        "\nlook vector:  "<<get_lookVector().toString()<<
        "\nUp vector:    "<<up.toString()<<
        "\nU axis:       "<<get_u_axis().toString()<<
        "\nV axis:       "<<get_v_axis().toString()<<
        "\nW axis:       "<<get_w_axis().toString()<<
        "\nfar, near :  ["<<far  << ", "<<near<<"]"<<
        "\nleft,right:  ["<<left << ", "<<right<<"]"<<
        "\ntop,bottom:  ["<<top  << ", "<<bottom<<"]\n";
    }
};

//capture info for an light source
class LightSource
{
public:
    //ambient intensity;
    Vertex amb_itn;
    
    //diffuse intensity
    Vertex dif_itn;
    
    //specular intensity
    Vertex spc_itn;
    
    //position
    Vertex position;
    
public:
    LightSource ()
    {
        amb_itn = Vertex();
        dif_itn = Vertex();
        spc_itn = Vertex();
        position = Vertex();
    }
    
    //set up light in the scene
    //light index will be the number of light in openGl
    void set_up(int light_index)
    {
        float amb[] = vertex_to_arr(amb_itn);
        float diff[] = vertex_to_arr(amb_itn);
        float spec[] = vertex_to_arr(amb_itn);
        float pos[] = vertex_to_arr(amb_itn);
        
        glLightfv(GL_LIGHT0+light_index, GL_AMBIENT, amb);
        glLightfv(GL_LIGHT0+light_index, GL_DIFFUSE, diff);
        glLightfv(GL_LIGHT0+light_index, GL_SPECULAR, spec);
        glLightfv(GL_LIGHT0+light_index, GL_POSITION, pos);
    }
    
    void print_debug_info()
    {
        CLOG("ambient:  "<<amb_itn.toString()<<"\n"<<
             "diffuse:  "<<dif_itn.toString()<<"\n"<<
             "specular: "<<spc_itn.toString()<<"\n"<<
             "position: "<<position.toString());
    }
};



#define x_i 0
#define y_i 1
#define z_i 2



enum Direction_t
{
    
LEFT, RIGHT
};

const float leg_scl[] = {2.0, 0.4, 0.4};
const float body_scl[] = {0.5, leg_scl[x_i]*1.3, 1.5};
const float hand_scl[] = {1.6, leg_scl[y_i], leg_scl[z_i]};
const float foot_scl[] = {1.0, 0.3, leg_scl[z_i]};

class Human
{
private:
    
    int left_high_leg_angle;
    int right_high_leg_angle;
    
    int left_lower_leg_angle;
    int right_lowe_leg_angle;
    
    int left_higher_arm_angle;
    int right_higher_arm_angle;
    
    int left_lower_arm_angle;
    int right_lower_arm_angle;
    
    int left_foot_angle;
    int right_foor_angle;
    
    Direction_t motion_dir = LEFT;
    
    float human_overall_scale;

    float head_radius;

    
    
public:
    
    Human()
    {
        left_high_leg_angle = -60, right_high_leg_angle = -120;
        left_lower_leg_angle = 0, right_lowe_leg_angle = 0;
        left_higher_arm_angle = -120, right_higher_arm_angle = -60;
        left_lower_arm_angle = -60, right_lower_arm_angle = 0;
        left_foot_angle = 0, right_foor_angle = 0;
        
         human_overall_scale =  0.5;
    
         head_radius =  human_overall_scale * body_scl[x_i];
    }
    
    
    void draw_arm(int higherHandAngle, int lowerHandAngle,int leftOrRight)
    {
        // Left hand
        glColor3f(0.5, 0.5, 0.5);
        glPushMatrix();
        {
            glTranslatef(0, human_overall_scale * body_scl[y_i], 0.0);
            glRotatef((GLfloat) higherHandAngle, 0.0, 0.0, 1.0);
            if (leftOrRight == LEFT)
            {
                glTranslatef(human_overall_scale * hand_scl[x_i] / 2, 0.0, human_overall_scale *
                             (body_scl[z_i] + hand_scl[z_i]) / 2);
            }else
            {
                glTranslatef(human_overall_scale * hand_scl[x_i] / 2, 0.0, -human_overall_scale * (body_scl[z_i] + hand_scl[z_i]) / 2);
            }
            
            glPushMatrix();
            {
                glScalef(hand_scl[x_i], hand_scl[y_i], hand_scl[z_i]);
                glutSolidCube(human_overall_scale);
            }
            glPopMatrix();
            glTranslatef(human_overall_scale * hand_scl[x_i] / 2, 0.0, 0.0);
            
            glRotatef((GLfloat) lowerHandAngle, 0.0, 0.0, 1.0);
            glTranslatef(human_overall_scale * hand_scl[x_i] / 2, 0.0, 0.0);
            glPushMatrix();
            {
                glScalef(hand_scl[x_i], hand_scl[y_i], hand_scl[z_i]);
                glutSolidCube(human_overall_scale);
                
            }
            glPopMatrix();
        }
        glPopMatrix();
        
    }
    
    void draw_Leg(int higherAngle, int loweAngle, int footAngle, int leftOrRight)
    {
        glPushMatrix();
        {
            glRotatef((GLfloat) higherAngle, 0.0, 0.0, 1.0);
            
            if (leftOrRight == LEFT)
            {
                glTranslatef(human_overall_scale * leg_scl[x_i]  / 2, 0.0, human_overall_scale * leg_scl[z_i]);
            }else
            {
                glTranslatef(human_overall_scale * leg_scl[x_i] / 2, 0.0, -human_overall_scale * leg_scl[z_i]);
            }
            
            
            //higher leg
            glPushMatrix();
            {
                glScalef(leg_scl[x_i], leg_scl[y_i], leg_scl[z_i]);
                glutSolidCube(human_overall_scale);
            }
            glPopMatrix();
            
            //lower leg
            glTranslatef(human_overall_scale * leg_scl[z_i], 0.0, 0.0);
            glRotatef((GLfloat) loweAngle, 0.0, 0.0, 1.0);
            glTranslatef(human_overall_scale *leg_scl[x_i] / 2, 0.0, 0.0);
            glPushMatrix();
            {
                glScalef(leg_scl[x_i], leg_scl[y_i], leg_scl[z_i]);
                glutSolidCube(human_overall_scale);
            }
            glPopMatrix();
            
            // Foot
            glTranslatef(human_overall_scale * leg_scl[x_i] / 2, 0.0, 0.0);
            glRotatef(90.0, 0.0, 0.0, 1.0);
            glTranslatef(-human_overall_scale * foot_scl[x_i] / 2 + human_overall_scale * leg_scl[y_i] / 2, 0.0, 0.0);
            glRotatef((GLfloat) footAngle, 0.0, 0.0, 1.0);
            glPushMatrix();
            {
                glScalef(foot_scl[x_i], foot_scl[y_i], foot_scl[z_i]);
                glutSolidCube(human_overall_scale);
            }
            glPopMatrix();
        }
        glPopMatrix();
    }
    
    
    void draw()
    {
        static bool single = false;
        if (!single)
        {
            static GLUquadricObj *cube = gluNewQuadric();
            gluQuadricNormals(cube, GLU_SMOOTH);
            gluQuadricTexture(cube, GL_TRUE);
            single = true;
        }
        
        glRotatef(-15, 0, 0, 1);
        glRotatef(180, 0, 1, 0);
        
        glPushMatrix();
        {
            draw_Leg(left_high_leg_angle, left_lower_leg_angle, left_foot_angle, LEFT);
            draw_Leg(right_high_leg_angle, right_lowe_leg_angle, right_foor_angle, RIGHT);
            
            draw_arm(left_higher_arm_angle, left_lower_arm_angle, LEFT);
            draw_arm(right_higher_arm_angle, right_lower_arm_angle, RIGHT);
            
            // Body
            glColor3f(1.0, 0.7, 1.0);
            glPushMatrix();
            {
                glTranslatef(0, human_overall_scale * body_scl[y_i] / 2, 0.0);
                glScalef(body_scl[x_i],body_scl[y_i], body_scl[z_i]);
                glutSolidCube(human_overall_scale);
            }
            glPopMatrix();
            
            // Head
            glPushMatrix();
            {
                glTranslatef(0, human_overall_scale * body_scl[y_i] + head_radius, 0.0);
                glutSolidSphere(head_radius, 10, 8);
            }
            glPopMatrix();
        }
        glPopMatrix();
    }
    
    void update_human_motion()
    {
        const int angle1 = -75;
        const int angle2= -105;
        const int lower_leg_change_spd = 6;
        const int right_angle = 90;
        const int left_leg_max_angle = -60;
        const int left_leg_min_angle = -120;
        
        
        if (motion_dir == LEFT)
        {
            //right leg swing ahead
            //left leg swing backaward
            left_high_leg_angle --;
            right_high_leg_angle ++;
            
            if (left_high_leg_angle >= angle1)
            {
                left_lower_leg_angle += lower_leg_change_spd;
                left_foot_angle = 0;
            }
            else if (left_high_leg_angle < angle2){
                left_lower_leg_angle -= lower_leg_change_spd;
            }
            
            if (left_high_leg_angle <= left_leg_min_angle) {
                motion_dir = RIGHT;
            }
            
            //foot angle
            right_foor_angle = -right_high_leg_angle - right_angle;
            
            //update arm angle
            left_higher_arm_angle ++;
            left_lower_arm_angle ++;
            right_higher_arm_angle --;
            right_lower_arm_angle --;
         
        }
        else
        {
            //left leg swing ahead
            //right leg swing backaward
            left_high_leg_angle ++;
            right_high_leg_angle --;
            
            if (right_high_leg_angle >= angle1)
            {
                right_lowe_leg_angle += lower_leg_change_spd;
                right_foor_angle = 0;
            }
            else if (right_high_leg_angle < angle2)
            {
                right_lowe_leg_angle -= lower_leg_change_spd;
            }
            
            if (left_high_leg_angle >= left_leg_max_angle)
            {
                motion_dir = LEFT;
            }
            
            //foot angle
            left_foot_angle = -left_high_leg_angle - right_angle;
            
             //update arm angle
            left_higher_arm_angle --;
            left_lower_arm_angle --;
            right_higher_arm_angle ++;
            right_lower_arm_angle ++;
        }
    }
};

