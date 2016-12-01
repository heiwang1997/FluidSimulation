#pragma once

#include <cassert>
#include <cmath>
#include <iostream>

template<class T>
struct Tensor22
{
	T t[2][2];

	Tensor22(void){}

	Tensor22(T value_for_all)
	{t[0][0]=t[0][1]=t[1][0]=t[1][1]=value_for_all;}

	template<class S>
	Tensor22(const Tensor22<S> source)
	{t[0][0]=(T)source.t[0][0];t[0][1]=(T)source.t[0][1];t[1][0]=(T)source.t[1][0];t[1][1]=(T)source.t[1][1];}

	template<class S>
	Tensor22(const S source[4])
	{t[0][0]=(T)source[0];t[0][1]=(T)source[1];t[1][0]=(T)source[2];t[1][1]=(T)souce[3];}

	template<class S>
	Tensor22(const S source[2][2])
	{t[0][0]=(T)source[0][0];t[0][1]=(T)source[0][1];t[1][0]=(T)source[1][0];t[1][1]=(T)souce[1][1];}

	Tensor22(T t00, T t01, T t10, T t11)
	{t[0][0]=t00;t[0][1]=t01;t[1][0]=t10;t[1][1]=t11;}

	T &operator[](int index)
	{
		assert(0<=index && (unsigned int)index<2);
		return t[index];
	}

	const T &operator[](int index) const
	{
		assert(0<=index && (unsigned int)index<2);
		return t[index];
	}

	// +=,-=,*=,/=

	inline Tensor22<T> transpose()
	{return Tensor22<T>(t[0][0],t[1][0],t[0][1],t[1][1]);}
};

template<class T>
struct Vec2
{
   T v[2];

   Vec2(void)
   {}

   Vec2(T value_for_all)
   { v[0]=v[1]=value_for_all; }

   template<class S>
   Vec2(const Vec2<S> source)
   { v[0]=(T)source.v[0]; v[1]=(T)source.v[1];}

   template<class S>
   Vec2(const S source[2])
   { v[0]=(T)source[0]; v[1]=(T)source[1];}

   Vec2(T v0, T v1)
   { v[0]=v0; v[1]=v1;}

   T &operator[](int index)
   {
      assert(0<=index && (unsigned int)index<2);
      return v[index];
   }

   const T &operator[](int index) const
   {
      assert(0<=index && (unsigned int)index<2);
      return v[index];
   }

   Vec2<T> operator+=(const Vec2<T> &w)
   { v[0]+=w.v[0]; v[1]+=w.v[1]; return *this; }

   Vec2<T> operator-=(const Vec2<T> &w)
   { v[0]-=w.v[0]; v[1]-=w.v[1]; return *this; }

   template<class S>
   Vec2<T> operator*=(S scalar)
   { v[0]*=scalar; v[1]*=scalar; return *this; }

   template<class S>
   Vec2<T> operator/=(S scalar)
   { v[0]/=scalar; v[1]/=scalar; return *this; }
};

template<class T> inline T mag2(const Vec2<T> &a)
{ return a.v[0]*a.v[0] + a.v[1]*a.v[1]; }

template<class T> inline T mag(const Vec2<T> &a)
{ return std::sqrt(mag2(a)); }

template<class T> inline T dist2(const Vec2<T> &a, const Vec2<T> &b)
{ return sqr(a.v[0]-b.v[0])+sqr(a.v[1]-b.v[1]); }

template<class T> inline T dist(const Vec2<T> &a, const Vec2<T> &b)
{ return std::sqrt(dist2(a,b)); }

template<class T> inline bool operator==(const Vec2<T> &a, const Vec2<T> &b)
{ return a.v[0]==b.v[0] && a.v[1]==b.v[1]; }

template<class T> inline bool operator!=(const Vec2<T> &a, const Vec2<T> &b)
{ return a.v[0]!=b.v[0] || a.v[1]!=b.v[1]; }

template<class T> inline Vec2<T> operator-(const Vec2<T> &a)
{ return Vec2<T>(-a.v[0], -a.v[1]); }

template<class T> inline Vec2<T> operator+(const Vec2<T> &a, const Vec2<T> &b)
{ return Vec2<T>(a.v[0]+b.v[0], a.v[1]+b.v[1]); }

template<class T> inline Vec2<T> operator-(const Vec2<T> &a, const Vec2<T> &b)
{ return Vec2<T>(a.v[0]-b.v[0], a.v[1]-b.v[1]); }

template<class S, class T>
inline Vec2<T> operator*(const Vec2<T> &a, S scalar)
{ return Vec2<T>(scalar*a.v[0], scalar*a.v[1]); }

template<class S, class T>
inline Vec2<T> operator*(S scalar, const Vec2<T> &a)
{ return Vec2<T>(scalar*a.v[0], scalar*a.v[1]); }

template<class S, class T>
inline Vec2<T> operator/(const Vec2<T> &a, S scalar)
{ return Vec2<T>(a.v[0]/scalar, a.v[1]/scalar); }

template<class T> inline T dot(const Vec2<T> &a, const Vec2<T> &b)
{ return a.v[0]*b.v[0] + a.v[1]*b.v[1]; }

//template<class T> inline Vec2<T> cross(const Vec2<T> &a, const Vec2<T> &b)
//{ return Vec2<T>(a.v[1]*b.v[2]-a.v[2]*b.v[1],
//                 a.v[2]*b.v[0]-a.v[0]*b.v[2],
//                 a.v[0]*b.v[1]-a.v[1]*b.v[0]); }

template<class T> inline void normalize(Vec2<T> &a)
{
	float l = mag(a);
	if (l<1e-6) {
		a.v[0] = a.v[1] = 0;
	} else a/=l; 
}

template<class T> inline Vec2<T> normalized(const Vec2<T> &a)
{ return a/mag(a); }

template<class T>
inline std::ostream &operator<<(std::ostream &out, const Vec2<T> &a)
{ return out<<a.v[0]<<' '<<a.v[1]; }

template<class T>
inline std::istream &operator>>(std::istream &in, Vec2<T> &a)
{ return in>>a.v[0]>>a.v[1]; }

template<class T>
inline unsigned int hash(const Vec2<T> &a)
{ return HashTable::hash(a.v[0]) ^ a.v[1]; }

template<class T>
inline void assign(const Vec2<T> &a, T &a0, T &a1, T &a2)
{ a0=a.v[0]; a1=a.v[1]; }

template<class T>
inline void clamp(Vec2<T> &a, T lower, T upper) {
	for (int i=0; i<2; i++) {
		if (a.v[i]<lower) a.v[i] = lower;
		else if (a.v[i]>upper) a.v[i] = upper;
	}
}

// common types of vectors ===================================================
typedef Vec2<float> Vec2f;
typedef Vec2<double> Vec2d;
typedef Vec2<int> Vec2i;

// type-specific operations ==================================================

template<class T> inline Vec2i round(const Vec2<T> &a)
{ return Vec2i(int(a.v[0]+0.5f), int(a.v[1]+0.5f)); }
template<class T> inline Vec2i floor(const Vec2<T> &a)
{ return Vec2i(int(a.v[0]), int(a.v[1]));}

//interact with Tensor22
template<class S, class T>
inline Vec2<T> operator*(Tensor22<S> tensor, const Vec2<T> &a)
{ return Vec2<T>(tensor.t[0][0]*a.v[0]+tensor.t[0][1]*a.v[1] , tensor.t[1][0]*a.v[0]+tensor.t[1][1]*a.v[1]); }


template<class T>
struct Vec3
{
   T v[3];

   Vec3(void)
   {}

   Vec3(T value_for_all)
   { v[0]=v[1]=v[2]=value_for_all; }

   template<class S>
   Vec3(const Vec3<S> source)
   { v[0]=(T)source.v[0]; v[1]=(T)source.v[1]; v[2]=(T)source.v[2]; }

   template<class S>
   Vec3(const S source[3])
   { v[0]=(T)source[0]; v[1]=(T)source[1]; v[2]=(T)source[2]; }

   Vec3(T v0, T v1, T v2)
   { v[0]=v0; v[1]=v1; v[2]=v2; }

   T &operator[](int index)
   {
      assert(0<=index && (unsigned int)index<3);
      return v[index];
   }

   const T &operator[](int index) const
   {
      assert(0<=index && (unsigned int)index<3);
      return v[index];
   }

   Vec3<T> operator+=(const Vec3<T> &w)
   { v[0]+=w.v[0]; v[1]+=w.v[1]; v[2]+=w.v[2]; return *this; }

   Vec3<T> operator-=(const Vec3<T> &w)
   { v[0]-=w.v[0]; v[1]-=w.v[1]; v[2]-=w.v[2]; return *this; }

   template<class S>
   Vec3<T> operator*=(S scalar)
   { v[0]*=scalar; v[1]*=scalar; v[2]*=scalar; return *this; }

   template<class S>
   Vec3<T> operator/=(S scalar)
   { v[0]/=scalar; v[1]/=scalar; v[2]/=scalar; return *this; }
};

template<class T> inline T mag2(const Vec3<T> &a)
{ return a.v[0]*a.v[0] + a.v[1]*a.v[1] + a.v[2]*a.v[2]; }

template<class T> inline T mag(const Vec3<T> &a)
{ return std::sqrt(mag2(a)); }

template<class T> inline T dist2(const Vec3<T> &a, const Vec3<T> &b)
{ return sqr(a.v[0]-b.v[0])+sqr(a.v[1]-b.v[1])+sqr(a.v[2]-b.v[2]); }

template<class T> inline T dist(const Vec3<T> &a, const Vec3<T> &b)
{ return std::sqrt(dist2(a,b)); }

template<class T> inline bool operator==(const Vec3<T> &a, const Vec3<T> &b)
{ return a.v[0]==b.v[0] && a.v[1]==b.v[1] && a.v[2]==b.v[2]; }

template<class T> inline bool operator!=(const Vec3<T> &a, const Vec3<T> &b)
{ return a.v[0]!=b.v[0] || a.v[1]!=b.v[1] || a.v[2]!=b.v[2]; }

template<class T> inline Vec3<T> operator-(const Vec3<T> &a)
{ return Vec3<T>(-a.v[0], -a.v[1], -a.v[2]); }

template<class T> inline Vec3<T> operator+(const Vec3<T> &a, const Vec3<T> &b)
{ return Vec3<T>(a.v[0]+b.v[0], a.v[1]+b.v[1], a.v[2]+b.v[2]); }

template<class T> inline Vec3<T> operator-(const Vec3<T> &a, const Vec3<T> &b)
{ return Vec3<T>(a.v[0]-b.v[0], a.v[1]-b.v[1], a.v[2]-b.v[2]); }

template<class S, class T>
inline Vec3<T> operator*(const Vec3<T> &a, S scalar)
{ return Vec3<T>(scalar*a.v[0], scalar*a.v[1], scalar*a.v[2]); }

template<class S, class T>
inline Vec3<T> operator*(S scalar, const Vec3<T> &a)
{ return Vec3<T>(scalar*a.v[0], scalar*a.v[1], scalar*a.v[2]); }

template<class S, class T>
inline Vec3<T> operator/(const Vec3<T> &a, S scalar)
{ return Vec3<T>(a.v[0]/scalar, a.v[1]/scalar, a.v[2]/scalar); }

template<class T> inline T dot(const Vec3<T> &a, const Vec3<T> &b)
{ return a.v[0]*b.v[0] + a.v[1]*b.v[1] + a.v[2]*b.v[2]; }

template<class T> inline Vec3<T> cross(const Vec3<T> &a, const Vec3<T> &b)
{ return Vec3<T>(a.v[1]*b.v[2]-a.v[2]*b.v[1],
                 a.v[2]*b.v[0]-a.v[0]*b.v[2],
                 a.v[0]*b.v[1]-a.v[1]*b.v[0]); }

template<class T> inline void normalize(Vec3<T> &a)
{
	float l = mag(a);
	if (l<1e-6) {
		a.v[0] = a.v[1] = a.v[2] = 0;
	} else a/=l; 
}

template<class T> inline Vec3<T> normalized(const Vec3<T> &a)
{ return a/mag(a); }

template<class T>
inline std::ostream &operator<<(std::ostream &out, const Vec3<T> &a)
{ return out<<a.v[0]<<' '<<a.v[1]<<' '<<a.v[2]; }

template<class T>
inline std::istream &operator>>(std::istream &in, Vec3<T> &a)
{ return in>>a.v[0]>>a.v[1]>>a.v[2]; }

template<class T>
inline unsigned int hash(const Vec3<T> &a)
{ return HashTable::hash(HashTable::hash(a.v[0]) ^ a.v[1]) ^ a.v[2]; }

template<class T>
inline void assign(const Vec3<T> &a, T &a0, T &a1, T &a2)
{ a0=a.v[0]; a1=a.v[1]; a2=a.v[2]; }

template<class T>
inline void clamp(Vec3<T> &a, T lower, T upper) {
	for (int i=0; i<3; i++) {
		if (a.v[i]<lower) a.v[i] = lower;
		else if (a.v[i]>upper) a.v[i] = upper;
	}
}

// common types of vectors ===================================================
typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;

// type-specific operations ==================================================

template<class T> inline Vec3i round(const Vec3<T> &a)
{ return Vec3i(int(a.v[0]+0.5f), int(a.v[1]+0.5f), int(a.v[2]+0.5f)); }
template<class T> inline Vec3i floor(const Vec3<T> &a)
{ return Vec3i(int(a.v[0]), int(a.v[1]), int(a.v[2]));}




