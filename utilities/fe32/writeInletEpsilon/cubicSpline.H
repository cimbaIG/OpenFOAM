#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#define BOOST_UBLAS_TYPE_CHECK 0

#include <cmath>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

class CubicSpline
{
   private:
      struct SplineCoeffs { double a,b,c,d; };
      bool valid_;
      std::vector<std::pair<double, double> > points_;
      std::vector<SplineCoeffs> data_;
      double lowBC_;
      double highBC_;
      std::vector<double> ddy_;

   public:
      CubicSpline() : valid_(false), lowBC_(0), highBC_(0) { }
      ~CubicSpline() {}

      unsigned int size() const
      {
         return points_.size();
      }

      void addPoint(double x, double y)
      {
         points_.push_back(std::make_pair(x, y));
      }

      double evaluate(double x_param, unsigned int derivative = 0)
      {
         prepare();
         if (!valid_)
            return 0.;

         // extrapolation below
         if (x_param < x(0))
         {
            const double t = x_param - x(0);
            const double d = (y(1) - y(0)) / h(0) - 2 * h(0) * (data_[0].b + 2*data_[1].b) / 6.0;
            if (derivative == 0)
               return t * t * lowBC_ + d * t + y(0);
            if (derivative == 1)
               return 2*t * lowBC_ + d;
         }

         // extrapolation above
         if (x_param >= x(size()-1))
         {
            const double t = x_param - x(size()-1);
            const double d = 2*h(size() - 2) * (ddy_[size() - 2] + 2*ddy_[size() - 1]) / 6 + (y(size() - 1) - y(size() - 2)) / h(size() - 2);
            if (derivative == 0)
               return t * t * highBC_ + d * t + y(size() - 1);
            if (derivative == 1)
               return 2*t * highBC_ + d;
         }

         // normal evaluation
         for (unsigned int i = 0; i < size()-1; ++i)
            if (x_param >= x(i) && x_param < x(i+1))
            {
               const double t = x_param - x(i);
               const SplineCoeffs &coef = data_[i];

               if (derivative == 0)
                  return ((coef.a * t + coef.b) * t + coef.c) * t + coef.d;
               if (derivative == 1)
                  return (3*coef.a * t + 2*coef.b) * t + coef.c;
            }

         return 0.;
      }

   private:

      double x(unsigned int i) const
      {
         return points_[i].first;
      }

      double y(unsigned int i) const
      {
         return points_[i].second;
      }

      double h(unsigned int i) const
      {
         return x(i+1) - x(i);
      }

      void prepare()
      {
         using namespace boost::numeric;
         if (valid_ || size() < 2)
            return;

         // sort the points based on x
         std::sort(points_.begin(), points_.end()); // identical points not taken into account

         const unsigned int e = size() - 1;

         // prepare matrix coefficients
         ublas::matrix<double> A = ublas::identity_matrix<double>(size(), size());

         for (unsigned int i = 1; i < e; ++i)
         {
            A(i-1,i) = h(i-1);
            A( i, i) = (h(i-1) + h(i)) * 2;
            A(i+1,i) = h(i);
         }

         // prepare RHS
         ublas::vector<double> C(size());
         C(0) = lowBC_;
         C(e) = highBC_;
         for (unsigned int i = 1; i < e; ++i)
            C(i) = 6 * ((y(i+1) - y(i)) / h(i) - (y(i) - y(i-1)) / h(i-1));

         // create a permutation matrix for the LU-factorization
         ublas::permutation_matrix<std::size_t> pm(A.size1());
         int res = ublas::lu_factorize(A, pm);
         if (res)
            return;

         // create identity matrix of "inverse"
         ublas::matrix<double> AInv(size(), size());
         AInv.assign(ublas::identity_matrix<double>(A.size1()));

         // backsubstitute to get the inverse
         ublas::lu_substitute(A, pm, AInv);

         // get the 2nd derivatives
         ublas::vector<double> ddy = ublas::prod(C, AInv);
         ddy_.resize(ddy.size());
         std::copy(ddy.begin(), ddy.end(), ddy_.begin());

         data_.resize(size() - 1);
         for (unsigned int i = 0; i < e; ++i)
         {
            SplineCoeffs &coef = data_[i];
            coef.a = (ddy(i+1) - ddy(i)) / (6.0 * h(i));
            coef.b = ddy(i) * 0.5;
            coef.c = (y(i+1) - y(i)) / h(i) - ddy(i+1) * h(i) / 6.0 - ddy(i) * h(i) / 3.0;
            coef.d = y(i);
         }

         valid_ = true;
      }
};

#endif // CUBICSPLINE_H
