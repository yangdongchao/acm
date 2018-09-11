#include<cstdio>
#include<cmath>
#include<cstring>
#include<iostream>
#include<algorithm>
#include<cstdlib>
#include<queue>
#include<map>
#include<stack>
#include<set>
#define e exp(1.0); //2.718281828
#define mod 1000000007
#define INF 0x7fffffff
#define inf 0x3f3f3f3f
typedef long long LL;
using namespace std;

#define zero(x) (((x)>0?(x):(-x))<eps)
const double eps=1e-8;
const double pi=acos(-1.0);

int sgn(double x) {
    if(fabs(x)<eps) return 0;
    if(x>0) return 1;
    return -1;
}
inline double sqr(double x) {
    return x*x;
}
struct point {
    double x,y;
    point() {};
    point(double a,double b):x(a),y(b) {};
    void input() {
        scanf("%lf %lf",&x,&y);
    }
    friend point operator + (const point &a,const point &b) {
        return point(a.x+b.x,a.y+b.y);
    }
    friend point operator - (const point &a,const point &b) {
        return point(a.x-b.x,a.y-b.y);
    }
    friend bool operator == (const point &a,const point &b) {
        return sgn(a.x-b.x)==0&&sgn(a.y-b.y)==0;
    }
    friend point operator * (const point &a,const double &b) {
        return point(a.x*b,a.y*b);
    }
    friend point operator * (const double &a,const point &b) {
        return point(a*b.x,a*b.y);
    }
    friend point operator / (const point &a,const double &b) {
        return point(a.x/b,a.y/b);
    }
    friend bool operator < (const point &a, const point &b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    }
    double norm() {
        return sqrt(sqr(x)+sqr(y));
    }
};
//计算两个向量的叉积
double cross(const point &a,const point &b) {
    return a.x*b.y-a.y*b.x;
}
double cross3(point A,point B,point C) { //叉乘
    return (B.x-A.x)*(C.y-A.y)-(B.y-A.y)*(C.x-A.x);
}
//计算两个点的点积
double dot(const point &a,const point &b) {
    return a.x*b.x+a.y*b.y;
}
double dot3(point A,point B,point C) { //点乘
    return (C.x-A.x)*(B.x-A.x)+(C.y-A.y)*(B.y-A.y);
}
//求平面点集的凸包 点集p,个数cnt,凸包点集 res，函数返回值凸包点集个数
int ConvexHull(point* P, int cnt, point* res) {  
    sort(P, P + cnt);
    cnt = unique(P, P + cnt) - P;
    int m = 0;
    for (int i = 0; i < cnt; i++) {
        while (m > 1 && cross(res[m - 1] - res[m - 2], P[i] - res[m - 2]) <= 0)
            m--;
        res[m++] = P[i];
    }
    int k = m;
    for (int i = cnt - 2; i >= 0; i--) {
        while (m > k && cross(res[m - 1] - res[m - 2], P[i] - res[m - 2]) <= 0)
            m--;
        res[m++] = P[i];
    }
    if (cnt > 1) m--;
    return m;
}

double len(point A,point B) { //返回向量AB的模平方
    return (A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y);
}
double minRetangleCover(point *res,int n) { //最小矩形面积覆盖(旋转卡壳)
    if(n < 3) return 0.0;
    res[n] = res [0];
    double ans = -1;
    int r = 1, p = 1,q;
    for(int i = 0; i < n; ++ i) {
        //卡出离边 res[i]-res[i+1]最远的点
        while(sgn(cross3(res[i],res[i+1],res[r+1])-cross3(res[i],res[i+1],res[r]))>=0)
            r = (r+1)% n;
        //卡出res[i]-res[i+1]方向上正向n最远的点
        while(sgn(dot3(res[i],res[i+1],res[p+1])-dot3(res[i],res[i+1],res[p]))>=0)
            p = (p+1)% n;
        if(i == 0) q = p;
        //卡出res[i]-res[i+1]方向上负向最远的点
        while(sgn(dot3(res[i],res[i+1],res[q+1])-dot3(res[i],res[i+1],res[q]))<=0)
            q = (q+1)% n;
        double d = len(res[i],res[i+1]);
        double temp = cross3(res[i],res[i+1],res[r])*(dot3(res[i],res[i+1],res[p])-dot3(res[i],res[i+1],res[q]))/d;
        if(ans < 0 || ans > temp) ans = temp;
    }
    return ans;
}
const int N=1005;
point a[N];
point cha[N];
int main() {
    int n;
    while(scanf("%d",&n)==1 && n) {
        for(int i = 0; i < n; ++ i) {
            scanf("%lf %lf",&a[i].x,&a[i].y);
        }
        n = ConvexHull(a,n,cha);//构造凸包
        printf("%.4f\n",minRetangleCover(cha,n));
    }
    return 0;
}