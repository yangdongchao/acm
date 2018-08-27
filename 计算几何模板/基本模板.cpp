#include<bits/stdc++.h>
using namespace std;
const double eps=1e-8;
const double inf=1e20;
const double PI=acos(-1.0);
const int maxn=1e5+3;
int sgn(double x)
{
    if(fabs(x)<eps) return 0;
    if(x<0) return -1;
    else return 1;
}
inline double square(double x) {return x*x;}
struct Point
{
    double x,y;
    Point()=default;
    Point(double x_, double y_):x(x_),y(y_){}
    void input(){scanf("%lf%lf",&x,&y);}
    void output(){printf("%.2lf %.2lf\n",x,y);}
    Point operator - (const Point &b) const
    {
        return Point(x-b.x,y-b.y);
    }
    bool operator==(Point b) const
    {
        return sgn(x-b.x)==0&&sgn(y-b.y)==0;
    }
    bool operator <(Point b) const
    {
        return sgn(x-b.x)==0?sgn(y-b.y)<0:x<b.x;
    }
    Point operator + (const Point &b) const
    {
        return Point(x+b.x,y+b.y);
    }
    double operator ^ (const Point &b) const
    {
        return x*b.y-y*b.x;
    }
    double operator * (const Point &b) const
    {
        return x*b.x+y*b.y;
    }
    Point operator/(const double &k)
    {
        return Point(x/k,y/k);
    }
    double len(){return hypot(x,y);}
    double len2(){return x*x+y*y;}
    double distance(Point b){return hypot(x-b.x,y-b.y);}
    double rad(Point a,Point b)
    {//求p与a,b的夹角
        Point p=*this;
        return fabs(atan2(fabs((a-p)^(b-p)),(a-p)*(b-p)));
    }
    Point rotate(Point p,double angle)
    {
        Point v=(*this)-p;
        double c=cos(angle),s=sin(angle);
        return Point(p.x+v.x*c-v.y*s, p.y+v.x*s + v.y*c);
    }
    Point rotleft()
    {
        return Point(-y,x);//逆90
    }
    Point rotright() {return Point(y,-x);}//顺90
    //化为长度为 r 的向量
    Point trunc(double r)//化为长度为r的向量
    {
        double l = len();
        if(!sgn(l))return *this;
        r /= l;
        return Point(x*r,y*r);
    }
};
struct Line
{
    Point s,e;
    Line()=default;
    Line(Point s_,Point e_):s(s_),e(e_){}
    bool operator == (Line v){return (s==v.s)&&(e==v.e);}
    Line(Point p,double angle)
    {
        s=p;
        if(sgn(angle-PI/2)==0)
        {
            e=s+Point(0,1);
        }
        else
        {
            e=s+Point(1,tan(angle));
        }
    }
    double angle()
    {
        double k=atan2(e.y-s.y,e.x-s.x);
        if(sgn(k)<0) k+=PI;
        if(sgn(k-PI)==0) k-=PI;
        return k;
    }
    int relation(Point p)
    {
        int c=sgn((p-s)^(e-s));
        if(c<0) return 1;//左
        else if(c>0) return 2;//右
        else return 3;//线上
    }
    bool pointSeg(Point p)
    {
        return sgn((p-s)^(e-s))==0&&sgn((p-s)*(p-e))<=0;
    }
    int segcrossseg(Line v)//2 规范
    {//1非规  0 不交
        int d1 = sgn((e-s)^(v.s-s));
        int d2 = sgn((e-s)^(v.e-s));
        int d3 = sgn((v.e-v.s)^(s-v.s)); int d4 = sgn((v.e-v.s)^(e-v.s));
        if( (d1^d2)==-2 && (d3^d4)==-2 ) return 2;
        return (d1==0 && sgn((v.s-s)*(v.s-e))<=0) ||
        (d2==0 && sgn((v.e-s)*(v.e-e))<=0) ||
        (d3==0 && sgn((s-v.s)*(s-v.e))<=0) ||
        (d4==0 && sgn((e-v.s)*(e-v.e))<=0);
    }
    int linecrossseg(Line v)
    {
        int d1 = sgn((e-s)^(v.s-s));
        int d2 = sgn((e-s)^(v.e-s));
        if((d1^d2)==-2) return 2;
        return (d1==0||d2==0);//0 不交
    }
    bool parallel(Line v)
    {
        return sgn((e-s)^(v.e-v.s))==0;
    }
    int linecrossline(Line v)
    {
        if((*this).parallel(v))
        {
            return v.relation(s)==3;
        }
        return 2;//0 平行,1重和,2相交
    }
    Point crossPoint(Line v)
    {
        double a1=(v.e-v.s)^(s-v.s);
        double a2=(v.e-v.s)^(e-v.s);
        return Point((s.x*a2-e.x*a1)/(a2-a1),(s.y*a2-e.y*a1)/(a2-a1));
    }
    double length()
    {
        return s.distance(e);
    }
    double disPointLine(Point p)
    {
        return fabs(((p-s)^(e-s))/length());
    }
    double disPointSeg(Point p)
    {
        if(sgn((p-s)*(e-s))<0||sgn((p-e)*(s-e))<0)
        {
            return min(p.distance(s),p.distance(e));
        }
        return disPointLine(p);
    }
    //返回点 p 在直线上的投影
};
struct circle
{
    Point p;//圆心
    double r;
    circle()=default;
    circle(Point a,double b):p(a),r(b){}
    circle(Point a,Point b,Point c)//三角行外接圆
    {
        Line u=Line((a+b)/2,( (a+b)/2+(b-a).rotleft()));
        Line v=Line( (b+c)/2 ,(b+c)/2+ ( (c-b).rotleft()));
        p=u.crossPoint(v);
        r=p.distance(a);
    }
    circle(Point a,Point b,Point c,bool t)//内接圆
    {
        Line u,v;
        double m = atan2(b.y-a.y,b.x-a.x), n = atan2(c.y-a.y,c.x-a.
        x);
        u.s = a;
        u.e = u.s + Point(cos((n+m)/2),sin((n+m)/2));
        v.s = b;
        m = atan2(a.y-b.y,a.x-b.x) , n = atan2(c.y-b.y,c.x-b.x);
        v.e = v.s + Point(cos((n+m)/2),sin((n+m)/2));
        p = u.crossPoint(v);
        r = Line(a,b).disPointSeg(p);
    }
    void input()
    {
        p.input();
        scanf("%lf",&r);
    }
    void output(){printf("%.2lf %.2lf %.2lf\n",p.x,p.y,r);}
    bool operator == (circle v){return (p==v.p)&&sgn(r-v.r)==0;}
    double area(){return PI*r;}//面积
    double circumference(){return 2*PI*r;}//周长
    bool operator < (circle v)const{ return ((p<v.p)||((p==v.p)&&sgn(r-v.r)<0));}
    int relation(Point b)//点与圆的位置关系 
    {
        double dst = b.distance(p);//0在圆外
        if(sgn(dst-r) < 0)return 2;//1 在圆上
        else if(sgn(dst-r)==0)return 1; //2在圆内
        return 0;
    }
    int relationseg(Line v)//线段和圆
    {
        double dst = v.disPointSeg(p);
        if(sgn(dst-r) < 0)return 2;
        else if(sgn(dst-r) == 0)return 1;
        return 0;
    }
    int relationline(Line v)//直线和圆
    {
        double dst = v.disPointLine(p);
        if(sgn(dst-r) < 0)return 2;
        else if(sgn(dst-r) == 0)return 1;
        return 0;
    }
    int relationcircle(circle v)//圆与圆的位置关系
    {
        double d = p.distance(v.p);     //5 相离
        if(sgn(d-r-v.r) > 0)return 5;   //4 外切
        if(sgn(d-r-v.r) == 0)return 4;  //3 相交
        double l = fabs(r-v.r); //2 内切
        if(sgn(d-r-v.r)<0 && sgn(d-l)>0)return 3;   //1 内含
        if(sgn(d-l)==0)return 2;
        if(sgn(d-l)<0)return 1;
    }


//求两个圆的交点，返回 0 表示没有交点，返回 1 是一个交点，2 是两个交点
    int pointcrosscircle(circle v,Point &p1,Point &p2)
    {
        int rel = relationcircle(v);
        if(rel == 1 || rel == 5)return 0;
        double d = p.distance(v.p);
        double l = (d*d+r*r-v.r*v.r)/(2*d);
        double h = sqrt(r*r-l*l);
        Point tmp = p + (v.p-p).trunc(l);
        p1 = tmp + ((v.p-p).rotleft().trunc(h));
        p2 = tmp + ((v.p-p).rotright().trunc(h));
        if(rel == 2 || rel == 4)
            return 1;
        return 2;
    }

    //得到过 a,b 两点，半径为 r1 的两个圆
    int gercircle(Point a,Point b,double r1,circle &c1,circle &c2)
    {
        circle x(a,r1),y(b,r1);
        int t = x.pointcrosscircle(y,c1.p,c2.p);
        if(!t)return 0;
        c1.r = c2.r = r;
        return t;
    }   
};

