#include<iostream>
#include<cmath>
#include<cstdio>
#include<algorithm>
#include<vector>
const long double eps=1e-10;
const long double PI=acos(-1);
 
using namespace std;
 
 
struct Point{
    long double x;
    long double y;
    Point(long double x=0,long double y=0):x(x),y(y){}
    void operator<<(Point &A) {cout<<A.x<<' '<<A.y<<endl;}
};
 
int dcmp(long double x)  {return (x>eps)-(x<-eps); }
 
typedef  Point  Vector;
 
Vector  operator +(Vector A,Vector B) { return Vector(A.x+B.x,A.y+B.y);}
 
Vector  operator -(Vector A,Vector B) { return Vector(A.x-B.x,A.y-B.y); }
 
Vector  operator *(Vector A,long double p) { return Vector(A.x*p,A.y*p);  }
 
Vector  operator /(Vector A,long double p) {return Vector(A.x/p,A.y/p);}
 
 
 
ostream &operator<<(ostream & out,Point & P) { out<<P.x<<' '<<P.y<<endl; return out;}
//
bool  operator< (const Point &A,const Point &B) { return dcmp(A.x-B.x)<0||(dcmp(A.x-B.x)==0&&dcmp(A.y-B.y)<0); }
 
bool  operator== ( const Point &A,const Point &B) { return dcmp(A.x-B.x)==0&&dcmp(A.y-B.y)==0;}
 
 
long double  Dot(Vector A,Vector B) {return A.x*B.x+A.y*B.y;}
 
long double  Cross(Vector A,Vector B)  {return A.x*B.y-B.x*A.y; }
 
long double  Length(Vector A)  { return sqrt(Dot(A, A));}
 
 
long double  Angle(Vector A,Vector B) {return acos(Dot(A,B)/Length(A)/Length(B));}
 
long double  Area2(Point A,Point B,Point C ) {return Cross(B-A, C-A);}
 
Vector Rotate(Vector A,long double rad) { return Vector(A.x*cos(rad)-A.y*sin(rad),A.x*sin(rad)+A.y*cos(rad));}
Vector Normal(Vector A) {long double L=Length(A);return Vector(-A.y/L,A.x/L);}
 
Point GetLineIntersection(Point P,Vector v,Point Q,Vector w)
{
    Vector u=P-Q;
    long double t=Cross(w, u)/Cross(v,w);
    return P+v*t;
    
}
 
long double DistanceToLine(Point P,Point A,Point B)
{
    Vector v1=P-A; Vector v2=B-A;
    return fabs(Cross(v1,v2))/Length(v2);
    
}
 
long double DistanceToSegment(Point P,Point A,Point B)
{
    if(A==B)  return Length(P-A);
    
    Vector v1=B-A;
    Vector v2=P-A;
    Vector v3=P-B;
    
    if(dcmp(Dot(v1,v2))==-1)    return  Length(v2);
    else if(Dot(v1,v3)>0)    return Length(v3);
    
    else return DistanceToLine(P, A, B);
    
}
 
Point GetLineProjection(Point P,Point A,Point B)
{
    Vector v=B-A;
    Vector v1=P-A;
    long double t=Dot(v,v1)/Dot(v,v);
    
    return  A+v*t;
}
 
bool  SegmentProperIntersection(Point a1,Point a2,Point b1,Point b2)
{
    long double c1=Cross(b1-a1, a2-a1);
    long double c2=Cross(b2-a1, a2-a1);
    long double c3=Cross(a1-b1, b2-b1);
    long double c4=Cross(a2-b1, b2-b1);
    
    return dcmp(c1)*dcmp(c2)<0&&dcmp(c3)*dcmp(c4)<0 ;
    
}
 
bool  OnSegment(Point P,Point A,Point B)
{
    return dcmp(Cross(P-A, P-B))==0&&dcmp(Dot(P-A,P-B))<0;
}
 
long double PolygonArea(Point *p,int n)
{
    long double area=0;
    
    for(int i=1;i<n-1;i++)
    {
        area+=Cross(p[i]-p[0], p[i+1]-p[0]);
        
    }
    return area/2;
    
}
 
Point  read_point()
{
    Point P;
    scanf("%Lf%Lf",&P.x,&P.y);
    return  P;
}
 
// ---------------与圆有关的--------
 
struct Circle
{
    Point c;
    long double r;
    
    Circle(Point c=Point(0,0),long double r=0):c(c),r(r) {}
    
    Point point(long double a)
    {
        return Point(c.x+r*cos(a),c.y+r*sin(a));
    }
    
    
};
 
struct  Line
{
    Point p;
    Vector v;
    Line(Point p=Point(0,0),Vector v=Vector(0,1)):p(p),v(v) {}
    
    Point point(long double t)
    {
        return Point(p+v*t);
    }
    
};
 
int getLineCircleIntersection(Line L,Circle C,long double &t1,long double &t2,vector<Point> &sol)
{
    long double a=L.v.x;
    long double b=L.p.x-C.c.x;
    long double c=L.v.y;
    long double d=L.p.y-C.c.y;
    
    long double e=a*a+c*c;
    long double f=2*(a*b+c*d);
    long double g=b*b+d*d-C.r*C.r;
    
    long double delta=f*f-4*e*g;
    
    if(dcmp(delta)<0) return 0;
    
    if(dcmp(delta)==0)
    {
        t1=t2=-f/(2*e);
        sol.push_back(L.point(t1));
        return 1;
    }
    
    else
    {
        t1=(-f-sqrt(delta))/(2*e);
        t2=(-f+sqrt(delta))/(2*e);
        
        sol.push_back(L.point(t1));
        sol.push_back(L.point(t2));
        
        return 2;
    }
    
}
 
// 向量极角公式
 
long double angle(Vector v)  {return atan2(v.y,v.x);}
 
int getCircleCircleIntersection(Circle C1,Circle C2,vector<Point> &sol)
{
    long double d=Length(C1.c-C2.c);
    
     if(dcmp(d)==0)
     {
         if(dcmp(C1.r-C2.r)==0)  return -1;  // 重合
         else return 0;    //  内含  0 个公共点
     }
    
    if(dcmp(C1.r+C2.r-d)<0)  return 0;  // 外离
    if(dcmp(fabs(C1.r-C2.r)-d)>0)  return 0;  // 内含
    
    long double a=angle(C2.c-C1.c);
    long double da=acos((C1.r*C1.r+d*d-C2.r*C2.r)/(2*C1.r*d));
    
    Point p1=C1.point(a-da);
    Point p2=C1.point(a+da);
    
    sol.push_back(p1);
    
    if(p1==p2)  return 1; // 相切
    else
    {
        sol.push_back(p2);
        return 2;
    }
}
 
 
//  求点到圆的切线
 
int getTangents(Point p,Circle C,Vector *v)
{
    Vector u=C.c-p;
    
    long double dist=Length(u);
    
    if(dcmp(dist-C.r)<0)  return 0;
    
    else if(dcmp(dist-C.r)==0)
    {
        v[0]=Rotate(u,PI/2);
        return 1;
    }
    
    else
    {
        
        long double ang=asin(C.r/dist);
        v[0]=Rotate(u,-ang);
        v[1]=Rotate(u,+ang);
        return 2;
    }
    
}
 
//  求两圆公切线
 
// 本题需要  加一个flag
bool flag=0;
 
int getTangents(Circle A,Circle B,Point *a,Point *b)
{
    int cnt=0;
    
    if(A.r<B.r)
    {
        flag=1;
        swap(A,B); swap(a, b);  //  有时需标记
    }
    
    long double d=Length(A.c-B.c);
   
    long double rdiff=A.r-B.r;
    long double rsum=A.r+B.r;
    
    if(dcmp(d-rdiff)<0)  return 0;   // 内含
    
    long double base=angle(B.c-A.c);
    
    if(dcmp(d)==0&&dcmp(rdiff)==0)   return -1 ;  // 重合 无穷多条切线
    
    if(dcmp(d-rdiff)==0)             // 内切   外公切线
    {
        a[cnt]=A.point(base);
        b[cnt]=B.point(base);
        cnt++;
        return 1;
    }
    
     // 有外公切线的情形
    
    long double ang=acos(rdiff/d);
    a[cnt]=A.point(base+ang);
    b[cnt]=B.point(base+ang);
    cnt++;
    a[cnt]=A.point(base-ang);
    b[cnt]=B.point(base-ang);
    cnt++;
    
    if(dcmp(d-rsum)==0)     // 外切 有内公切线
    {
        a[cnt]=A.point(base);
        b[cnt]=B.point(base+PI);
        cnt++;
    }
    
    else  if(dcmp(d-rsum)>0)   // 外离   又有两条外公切线
    {
        long double  ang_in=acos(rsum/d);
        a[cnt]=A.point(base+ang_in);
        b[cnt]=B.point(base+ang_in+PI);
        cnt++;
        a[cnt]=A.point(base-ang_in);
        b[cnt]=B.point(base-ang_in+PI);
        cnt++;
    }
    
    return cnt;
}
 
int main()
{
    Circle A,B;
    Point Zero_p;
    while(1)
    {
        flag=0;
        A.c=read_point();
        cin>>A.r;
        B.c=read_point();
        cin>>B.r;
        if(A.c==Zero_p&&B.c==Zero_p&&A.r==0&&B.r==0)  break;
        Point a[5],b[5];
        int ans=getTangents(A,B,a,b);
        for(int i=0;i<ans;i++)
            for(int j=i+1;j<ans;j++)
            {
                if(a[j]<a[i])
                {
                    swap(a[i],a[j]);
                    swap(b[i],b[j]);
                }
            }
        
        printf("%d\n",ans);
        for(int i=0;i<ans;i++)
            printf("%.5Lf %.5Lf %.5Lf %.5Lf %.5Lf\n",a[i].x,a[i].y,b[i].x,b[i].y,Length((a[i]-b[i])));
    }
}
