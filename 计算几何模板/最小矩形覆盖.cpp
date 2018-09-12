#pragma GCC optimize(2)
#pragma G++ optimize(2)
#include<cstring>
#include<cmath>
#include<iostream>
#include<algorithm>
#include<cstdio>

#define eps 0.00000001
#define N 50007
using namespace std;
inline int read()
{
    int x=0,f=1;char ch=getchar();
    while(!isdigit(ch)){if(ch=='-')f=-1;ch=getchar();}
    while(isdigit(ch)){x=(x<<1)+(x<<3)+ch-'0';ch=getchar();}
    return x*f;
}

int n,tot;
double ans=1e60;
struct P
{
    double x,y;
    P(){}
    P(double _x,double _y):x(_x),y(_y){}
    friend bool operator<(P a,P b){return fabs(a.y-b.y)<eps?a.x<b.x:a.y<b.y;}
    friend bool operator==(P a,P b){return fabs(a.x-b.x)<eps&&fabs(a.y-b.y)<eps;}
    friend bool operator!=(P a,P b){return !(a==b);}
    friend P operator+(P a,P b){return P(a.x+b.x,a.y+b.y);}
    friend P operator-(P a,P b){return P(a.x-b.x,a.y-b.y);}
    friend double operator*(P a,P b){return a.x*b.y-a.y*b.x;}
    friend P operator*(P a,double b){return P(a.x*b,a.y*b);}
    friend double operator/(P a,P b){return a.x*b.x+a.y*b.y;}
    friend double dis(P a){return sqrt(a.x*a.x+a.y*a.y);}
}p[N],q[N],t[5];

bool cmp(P a,P b)
{
    double t=(a-p[1])*(b-p[1]);
    if(fabs(t)<eps)return dis(p[1]-a)-dis(p[1]-b)<0;
    return t>0;
}
void Graham()
{
    for (int i=2;i<=n;i++)
        if(p[i]<p[1])swap(p[i],p[1]);
    sort(p+2,p+n+1,cmp);
    q[++tot]=p[1];
    for (int i=2;i<=n;i++)
    {
        while(tot>1&&(q[tot]-q[tot-1])*(p[i]-q[tot])<eps)tot--;
        q[++tot]=p[i];
    }
    q[0]=q[tot];//凸包是一个回路。
}
void RC()
{
    int l=1,r=1,p=1;
    double L,R,D,H;
    for (int i=0;i<tot;i++)
    {
        D=dis(q[i]-q[i+1]);
        while((q[i+1]-q[i])*(q[p+1]-q[i])-(q[i+1]-q[i])*(q[p]-q[i])>-eps)p=(p+1)%tot;
        while((q[i+1]-q[i])/(q[r+1]-q[i])-(q[i+1]-q[i])/(q[r]-q[i])>-eps)r=(r+1)%tot;
        if(i==0)l=r;
        while((q[i+1]-q[i])/(q[l+1]-q[i])-(q[i+1]-q[i])/(q[l]-q[i])<eps)l=(l+1)%tot;
        L=(q[i+1]-q[i])/(q[l]-q[i])/D,R=(q[i+1]-q[i])/(q[r]-q[i])/D;
        H=(q[i+1]-q[i])*(q[p]-q[i])/D;
        if(H<0)H=-H;
        double tmp=(R-L)*H;
        if(tmp<ans)
        {
            ans=tmp;
            t[0]=q[i]+(q[i+1]-q[i])*(R/D);
            t[1]=t[0]+(q[r]-t[0])*(H/dis(t[0]-q[r]));
            t[2]=t[1]-(t[0]-q[i])*((R-L)/dis(q[i]-t[0]));
            t[3]=t[2]-(t[1]-t[0]);
        }
    }
}
int main()
{
    n=read();
    for (int i=1;i<=n;i++)
        scanf("%lf%lf",&p[i].x,&p[i].y);
    Graham();
    RC();
    printf("%.5lf\n",ans);
    int fir=0;
    for (int i=1;i<=3;i++)
        if(t[i]<t[fir])fir=i;
    for (int i=0;i<=3;i++)
        printf("%.5lf %.5lf\n",t[(i+fir)%4].x,t[(i+fir)%4].y);
}