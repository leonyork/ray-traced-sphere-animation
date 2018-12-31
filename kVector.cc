#include "kVector.h"
#include <math.h>


namespace kBase {
    
    /***********************************
     * v3D: 3D Vector
     ***********************************/
    
    v3D::v3D(real x, real y, real z) {
        set(x, y, z);
    }
    
    v3D::v3D() {
        set();
    }
    
    void v3D::set(real x, real y, real z) {
        this->x=x;
        this->y=y;
        this->z=z;
    }
    
    void v3D::set(const v3D &v) {
        this->x=v.x;
        this->y=v.y;
        this->z=v.z;
    }
    
    void v3D::set() {
        this->set(0, 0, 0);
    }
    
    real v3D::getSqrLength() const{
        return this->x*this->x+this->y*this->y+this->z*this->z;
    }
    
    real v3D::getLength() const{
        return sqrt(this->getSqrLength());
    }
    
    bool v3D::isZero() const
    {
        return (getSqrLength()<EPSILON);
    }
    
    
    // Operator Functions
    
    
    v3D v3D::operator+ (const v3D &v) const
    {
        return add(*this,v);
    }
    
    void v3D::operator+= (const v3D &v)
    {
        this->x+=v.x;
        this->y+=v.y;
        this->z+=v.z;
    }
    
    v3D v3D::operator- (const v3D &v) const
    {
        return v3D(this->x-v.x,this->y-v.y, this->z-v.z);
    }
    
    void v3D::operator-= (const v3D &v)
    {
        this->x-=v.x;        
        this->y-=v.y;
        this->z-=v.z;
    }
        
    
    void v3D::operator*= (real a)
    {
        this->x*=a;
        this->y*=a;
        this->z*=a;
    }
    
    bool v3D::operator== (v3D &v) const
    {
        v3D diff = v-(*this);
        return (diff.getSqrLength()<EPSILON);
    }    
    
    
    // Global Functions
     
    
    v3D add(const v3D &u, const v3D &v) {
        return v3D(u.x+v.x, u.y+v.y, u.z+v.z);
    }
    
    v3D scalarMultiply(real a, const v3D &v) 
    {
        return v3D(a*v.x, a*v.y, a*v.z);
    }
    
    real dotProduct(const v3D &u, const v3D &v) {
        return u.x*v.x+u.y*v.y+u.z*v.z;
    }
    
    v3D crossProduct(const v3D &u, const v3D &v)
    {
        return v3D(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
    }
    
    v3D componentProduct(const v3D &u, const v3D &v)
    {
        return v3D(u.x*v.x, u.y*v.y, u.z*v.z);
    }
    
    v3D operator* (real a, const v3D &v)
    {
        return scalarMultiply(a,v);
    }
    
    v3D operator* (const v3D &v, real a)
    {
        return scalarMultiply(a,v);
    }
    
    std::ostream &operator<< (std::ostream &os, const v3D &v)
    {
        os.precision(4);
        os << std::scientific << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
    
    std::istream &operator>> (std::istream &is, v3D &v)
    {
        is >> v.x >> v.y >> v.z;
        return is;
    }
    
    /***********************************
     * v2D: 2D Vector
     ***********************************/  
    
    v2D::v2D(real x, real y) {
        set(x, y);
    }
    
    v2D::v2D() {
        set();
    }
    
    void v2D::set(real x, real y) {
        this->x=x;
        this->y=y;
    }
    
    void v2D::set(const v2D &v) {
        this->x=v.x;
        this->y=v.y;
    }
    
    void v2D::set() {
        this->set(0, 0);
    }
    
    real v2D::getSqrLength() const{
        return this->x*this->x+this->y*this->y;
    }
    
    real v2D::getLength() const{
        return sqrt(this->getSqrLength());
    }
    
    bool v2D::isZero() const
    {
        return (getSqrLength()<EPSILON);
    }    
    
    
    // Operator Functions
    
    
    v2D v2D::operator+ (const v2D &v) const
    {
        return add(*this,v);
    }
    
    void v2D::operator+= (const v2D &v)
    {
        this->x+=v.x;
        this->y+=v.y;
    }
    
    v2D v2D::operator- (const v2D &v) const
    {
        return v2D(this->x-v.x, this->y-v.y);
    }
    
    void v2D::operator-= (const v2D &v)
    {
        this->x-=v.x;
        this->y-=v.y;
    }    
    
    void v2D::operator*= (real a)
    {
        this->x*=a;
        this->y*=a;
    }
    
    bool v2D::operator== (v2D &v) const
    {
        v2D diff = v-(*this);
        return (diff.getSqrLength()<EPSILON);
    }
    
    
    // Global Functions
     
    
    v2D add(const v2D &u, const v2D &v) {
        return v2D(u.x+v.x, u.y+v.y);
    }
    
    v2D scalarMultiply(real a, const v2D &v) 
    {
        return v2D(a*v.x, a*v.y);
    }
    
    real dotProduct(const v2D &u, const v2D &v) {
        return u.x*v.x+u.y*v.y;
    }
    
    real determinant(const v2D &u, const v2D &v)
    {
        return u.x*v.y-u.y*v.x;
    }
    
    v2D componentProduct(const v2D &u, const v2D &v)
    {
        return v2D(u.x*v.x, u.y*v.y);
    }
    
    v2D operator* (real a, const v2D &v)
    {
        return scalarMultiply(a,v);
    }
    
    v2D operator* (const v2D &v, real a)
    {
        return scalarMultiply(a,v);
    } 
    
    std::ostream &operator<< (std::ostream &os, const v2D &v)
    {
        os.precision(4);
        os << std::scientific << "(" << v.x << ", " << v.y << ")";
        return os;
    }
    
    std::istream &operator>> (std::istream &is, v2D &v)
    {
        is >> v.x >> v.y;
        return is;
    }    
    
    /**********************************
     * r3D: 3D ray
     **********************************/
    
    r3D::r3D()
    {
        this->set();
    } 
    
    r3D::r3D(v3D x, v3D v)
    {
        this->set(x,v);
    }
    
    void r3D::set()
    {
        this->x = v3D();
        this->v = v3D();
    }
    
    void r3D::set(const v3D &x, const v3D &v)
    {
        this->x=x;
        this->v=v;
    }
    
    void r3D::set(const r3D &r)
    {
        this->x = r.x;
        this->v = r.v;
    }
    
    r3D r3D::operator+ (real t) const
    {
        return r3D((*this).x+(*this).v*t, (*this).v);
    }
    
    void r3D::operator+= (real t)
    {
        this->x += (this->v)*t;
    }
    
    std::ostream &operator<< (std::ostream &os, const r3D &r)
    {
        os << "[" << r.x << " -> " << r.v << "]";
        return os;
    }
    std::istream &operator>> (std::istream &is, r3D &r)
    {
        is >> r.x >> r.v;
        return is;
    }
    
    
    /******************************************************
     * kPointStack : Stack of 3D Points
     ******************************************************/
    
        kPointStack::kPointStack(unsigned char maxLength)
        {
            this->maxLength = maxLength;
            this->length = 0;
            this->curPos = 0;
            this->points = new v3D[maxLength];
        }
        
        kPointStack::kPointStack(unsigned char maxLength, const v3D &firstPoint)
        {
            this->maxLength = maxLength;
            this->length = 1;
            this->curPos = 1;
            this->points = new v3D[maxLength]; 
            this->points[0] = firstPoint;
        }
        kPointStack::kPointStack(unsigned char maxLength, const v3D *points, unsigned char numOfPoints)
        {
            this->maxLength = maxLength;
            this->points = new v3D[maxLength]; 
            if (numOfPoints > maxLength)
            {
                this->length = maxLength;  
                this->curPos = 0;
                unsigned char difference = numOfPoints-maxLength;
                for (unsigned char i=0; i<maxLength; i++)
                {
                    this->points[i] = points[i+difference];                     
                }
            }
            else
            {
                this->length = numOfPoints;
                this->curPos = numOfPoints;
                for (unsigned char i=0; i<numOfPoints; i++)
                {
                    this->points[i] = points[i];                     
                }        
            }
        }
        
        void kPointStack::newPoint(const v3D &point)
        {
            if (curPos == maxLength)
                curPos=0;
            
            this->points[curPos] = point;
            curPos++; 
            
            if (curPos>length)
                length=curPos;
        }
        
        unsigned char kPointStack::getLength() const
        {
            return length;
        }
        
        unsigned char kPointStack::getMaxLength() const
        {
            return maxLength;
        }
        
        v3D kPointStack::getPoint(unsigned char n) const
        {
            unsigned char toGet = this->curPos + n;
            if (toGet > maxLength)
                toGet-=maxLength;
            return this->points[toGet];
        }
        
    /******************************************************
     * kPointStack : Stack of 2D Points
     ******************************************************/
    
        kPointStack2D::kPointStack2D(unsigned char maxLength)
        {
            this->maxLength = maxLength;
            this->length = 0;
            this->curPos = 0;
            this->points = new v2D[maxLength];
        }
        
        kPointStack2D::kPointStack2D(unsigned char maxLength, const v2D &firstPoint)
        {
            this->maxLength = maxLength;
            this->length = 1;
            this->curPos = 1;
            this->points = new v2D[maxLength]; 
            this->points[0] = firstPoint;
        }
        kPointStack2D::kPointStack2D(unsigned char maxLength, const v2D *points, unsigned char numOfPoints)
        {
            this->maxLength = maxLength;
            this->points = new v2D[maxLength]; 
            if (numOfPoints > maxLength)
            {
                this->length = maxLength;  
                this->curPos = 0;
                unsigned char difference = numOfPoints-maxLength;
                for (unsigned char i=0; i<maxLength; i++)
                {
                    this->points[i] = points[i+difference];                     
                }
            }
            else
            {
                this->length = numOfPoints;
                this->curPos = numOfPoints;
                for (unsigned char i=0; i<numOfPoints; i++)
                {
                    this->points[i] = points[i];                     
                }        
            }
        }
        
        void kPointStack2D::newPoint(const v2D &point)
        {
            if (curPos == maxLength)
                curPos=0;
            
            this->points[curPos] = point;
            curPos++; 
            
            if (curPos>length)
                length=curPos;
        }
        
        unsigned char kPointStack2D::getLength() const
        {
            return length;
        }
        
        unsigned char kPointStack2D::getMaxLength() const
        {
            return maxLength;
        }
        
        v2D kPointStack2D::getPoint(unsigned char n) const
        {
            unsigned char toGet = this->curPos + n;
            if (toGet > maxLength)
                toGet-=maxLength;
            return this->points[toGet];
        } 
        
        /**************************************************
         * kParticle: A particle with constant acceleration
         ***************************************************/

        kParticle::kParticle():x(), v(), a(){};    //!< Constructs Particle with initial displacement, velocity and acceleration zero
        kParticle::kParticle(const v3D &x):x(x), v(), a(){};   //!< Constructs Particle with initial displacement x, velocity and acceleration zero
        kParticle::kParticle(const v3D &x, const v3D &v):x(x), v(v), a(){}; //!< Constructs Particle with initial displacement x, velocity v and acceleration zero
        kParticle::kParticle(const v3D &x, const v3D &v, const v3D &a):x(x), v(v), a(a){};//!< Constructs Particle with initial displacement x, velocity v and acceleration a
        
        kParticle kParticle::operator+(real t) const
        {
            return kParticle(this->x+this->v*t+0.5*t*t*this->a, this->v+t*this->a, this->a);
        }
        
        void kParticle::operator+=(real t)   
        {
            this->x += this->v*t+0.5*t*t*this->a;
            this->v += t*this->a;
        }
}

