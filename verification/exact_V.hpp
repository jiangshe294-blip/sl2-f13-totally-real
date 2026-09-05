// Public source-only copy of the independently audited V arithmetic.
// Archive main removed; no archive files, native binaries or absolute paths.
#pragma once
#include <gmpxx.h>
#include <vector>
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <unistd.h>
using Z=mpz_class;using Q=mpq_class;using Vec=std::vector<Z>;using Mat=std::vector<Vec>;
struct V{Mat a=Mat(14,Vec(13));Z d=1;};using K=std::vector<V>;
auto start=std::chrono::steady_clock::now();Z G,G38;Vec gi(15);int reductions=0;long long scalar_products=0;
void cap(){if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>35)throw std::runtime_error("35 second native reserve");}
Z rz(){std::string s;if(!(std::cin>>s))throw std::runtime_error("integer input");return Z(s);}Q rq(){std::string s;if(!(std::cin>>s))throw std::runtime_error("rational input");Q q(s);q.canonicalize();return q;}
void norm(V&v){Z c=v.d;for(auto&row:v.a)for(auto&z:row)mpz_gcd(c.get_mpz_t(),c.get_mpz_t(),z.get_mpz_t());if(c!=1){v.d/=c;for(auto&row:v.a)for(auto&z:row)z/=c;}}
V readv(){V v;v.d=rz();if(v.d<=0)throw std::runtime_error("denominator");for(auto&row:v.a)for(auto&z:row)z=rz();norm(v);return v;}
V add(const V&a,const V&b,int sign=1){V c;mpz_lcm(c.d.get_mpz_t(),a.d.get_mpz_t(),b.d.get_mpz_t());Z x=c.d/a.d,y=sign*(c.d/b.d);for(int i=0;i<14;i++)for(int j=0;j<13;j++)c.a[i][j]=x*a.a[i][j]+y*b.a[i][j];norm(c);return c;}
void conv(Mat&p,const V&a,const V&b){for(int i=0;i<14;i++){for(int j=0;j<13;j++)if(a.a[i][j]!=0)for(int k=0;k<14;k++)for(int l=0;l<13;l++)if(b.a[k][l]!=0){mpz_addmul(p[i+k][j+l].get_mpz_t(),a.a[i][j].get_mpz_t(),b.a[k][l].get_mpz_t());scalar_products++;}cap();}}
V reduce(Mat p,const Z&den){
 for(int k=26;k>=14;k--){Vec top=p[k];for(auto&row:p)for(auto&z:row)z*=G;for(int i=0;i<15;i++)for(int j=0;j<25;j++)p[k-14+i][j]-=top[j]*gi[i];}
 for(int k=24;k>=14;k--){Vec top(14);for(int i=0;i<14;i++)top[i]=p[i][k];for(int i=0;i<14;i++)for(auto&z:p[i])z*=G;for(int j=0;j<15;j++)for(int i=0;i<14;i++)p[i][k-14+j]-=top[i]*gi[j];}
 V c;c.d=den*G38;
 for(int j=0;j<13;j++){Vec a(27);for(int i=0;i<14;i++)a[i]=G*p[i][j];for(int i=0;i<14;i++)for(int k=0;k<=13-j;k++)a[i+k]-=p[i][13]*gi[j+k+1];
  for(int k=26;k>=14;k--){Z top=a[k];for(auto&z:a)z*=G;for(int i=0;i<15;i++)a[k-14+i]-=top*gi[i];}
  for(int i=14;i<27;i++)if(a[i]!=0)throw std::runtime_error("high reduced coordinate");for(int i=0;i<14;i++)c.a[i][j]=a[i];}
 norm(c);reductions++;cap();return c;
}
V mul(const V&a,const V&b){Mat p(27,Vec(25));conv(p,a,b);return reduce(p,a.d*b.d);}
