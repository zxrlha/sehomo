#ifndef _WAMAM_FINALSTATE_HPP
#define _WAMAM_FINALSTATE_HPP 11

#include <cassert>
#include <cmath>
#include <iostream>
#include <exception>
#include "momentum.hpp"

namespace wamam
{
	namespace _w_detail
	{
		template<typename T>
		void kkboost(const momentum<T>& k1, const momentum<T>& k2, const momentum<T>& p1, momentum<T>& p2)
		{
			T ksq = k1 * k1;
			momentum<T> kksum = k1 + k2;
			p2 = p1 - 2. * kksum * (kksum * p1) * (1./ (kksum * kksum)) + 2. * k2 * (k1 * p1) * (1./ ksq);
		}

		template<typename T>
		T rcm_pmax(const momentum<T>& Q, T m1, T m2)
		{
			T E2 = Q * Q;
			T t1 = m1 * m1 - m2 * m2 - E2;
			T t2 = t1 * t1 / 4 / E2 - m2 * m2;
			assert(t2 >= 0);
			T p = sqrt(t2);
			return p;
		}

		template<typename T>
		void rcm_2(const momentum<T>& Q, T x[], T m[], momentum<T> p[], T& jkb)
		{
			T theta = x[0] * M_PI;
			T phi = x[1] * 2 * M_PI;
			T pmax = rcm_pmax(Q, m[0], m[1]);
			T E1CM = std::sqrt(pmax * pmax + m[0] * m[0]);
			T E2CM = std::sqrt(pmax * pmax + m[1] * m[1]);
			T ECM = sqrt(Q * Q);
			jkb = pmax * std::sin(theta) / 8 / ECM;
			momentum<T> p1cm;
			p1cm.px() = pmax * std::sin(theta) * std::cos(phi);
			p1cm.py() = pmax * std::sin(theta) * std::sin(phi);
			p1cm.pz() = pmax * std::cos(theta);
			p1cm.e() = E1CM;
			momentum<T> QCM;
			QCM.px() = 0;
			QCM.py() = 0;
			QCM.pz() = 0;
			QCM.e() = ECM;
			kkboost(QCM, Q, p1cm, p[0]);
			p[1] = Q - p[0];
		}

		template<typename T>
		void rcm_n(int n, const momentum<T>& Q, T x[], T m[], momentum<T> p[], T& jkb)
		{
			if (n > 2)
			{
				//determine particle 1
				T mremain = T();
				for (int i = 1; i < n; ++i)
				{
					mremain += m[i];
				}
				T theta = M_PI * x[0];
				T phi = 2 * M_PI * x[1];
				T pmax = rcm_pmax(Q, m[0], mremain);
				T p1 = pmax * x[2];
				momentum<T> p1cm;
				p1cm.px() = p1 * std::sin(theta) * std::cos(phi);
				p1cm.py() = p1 * std::sin(theta) * std::sin(phi);
				p1cm.pz() = p1 * std::cos(theta);
				p1cm.e() = std::sqrt(p1 * p1 + m[0] * m[0]);
				momentum<T> QCM;
				QCM.px() = 0;
				QCM.py() = 0;
				QCM.pz() = 0;
				QCM.e() = std::sqrt(Q * Q);
				kkboost(QCM, Q, p1cm, p[0]);
				momentum<T> Qremain = Q - p[0];
				rcm_n(n - 1, Qremain, x + 3, m + 1, p + 1, jkb);
				jkb *= pmax * pmax * pmax * x[2] * x[2] * sin(theta) / 8 / M_PI / p1cm.e();
			}
			else if (n == 2)
			{
				return rcm_2(Q, x, m, p, jkb);
			}
			else
			{
				std::cerr << "ERROR:" << "wamam:rcm" << "particle number should >= 2" << std::endl;
				throw std::invalid_argument("n");
			}
		}
	}

	//this first method: recursive method
	//suitable for massive state
	template<typename T>
	void rcm(int n, const momentum<T>& Q, T x[], T mass[], momentum<T> mom[], T& jkb)
	{
		_w_detail::rcm_n(n, Q, x, mass, mom, jkb);
	}
}

#endif
