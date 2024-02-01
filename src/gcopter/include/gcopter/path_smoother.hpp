#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include "cubic_spline.hpp"
#include "lbfgs.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>

namespace path_smoother
{

    class PathSmoother
    {
    private:
        cubic_spline::CubicSpline cubSpline;

        int pieceN;
        int obstaN;
        Eigen::Matrix3Xd diskObstacles;
        double penaltyWeight;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        Eigen::Matrix2Xd points;
        Eigen::Matrix2Xd gradByPoints;

        lbfgs::lbfgs_parameter_t lbfgs_params;

    private:
        static inline double costFunction(void *ptr,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &g)
        {
           //TODO
           /*
                需要从PathSmoother中调用函数以及变量 PathSmoother &obj = *(PathSmoother *)ptr; 
                例如代价函数的调用 
           */
            PathSmoother &obj = *(PathSmoother *)ptr;  
            for(int i = 0; i < obj.pieceN-1;i++)
            {
                obj.points(0,i)=x(2*i+0);
                obj.points(1,i)=x(2*i+1);
            }

            Eigen::Matrix2Xd gE;
            double energy = obj.costEnergy(gE);

            Eigen::Matrix2Xd gP;
            double potential = obj.constPotential(gP);

            obj.gradByPoints = gP + gE;

            for(int i = 0; i < obj.pieceN-1; i++)
            {
                g(2*i +0) = obj.gradByPoints(0,i);
                g(2*i +1) = obj.gradByPoints(1,i);
            }
            double cost = energy + potential;
            return cost;
        }

        inline double costEnergy(Eigen::Matrix2Xd &gEnergy)
        {
            gEnergy.resize(2,pieceN-1);
            gEnergy.setZero();
            double energy = 0;
          
            cubSpline.setInnerPoints(points); //要计算getGrad里面的coffs
        
            cubSpline.getStretchEnergy(energy);
   
            cubSpline.getGrad(gEnergy);
       
            return energy;
        }

        inline double constPotential(Eigen::Matrix2Xd &gPotential)
        {
            gPotential.resize(2,pieceN-1);
            gPotential.setZero();
            double potential = 0, overlap = 0;
            for(int i = 0; i < pieceN-1; i++)
            {
                for(int j = 0; j < obstaN; j++)
                {
                    overlap = diskObstacles(2,j) - std::sqrt(std::pow(diskObstacles(0,j)-points(0,i),2)
                                + std::pow(diskObstacles(1,j)-points(1,i),2));
                    if(overlap>0)
                    {
                        gPotential(0,i)+=(diskObstacles(0,j)-points(0,i))/std::sqrt(std::pow(diskObstacles(0,j)-points(0,i),2)
                                + std::pow(diskObstacles(1,j)-points(1,i),2));
                        
                        gPotential(1,i)+=(diskObstacles(1,j)-points(1,i))/std::sqrt(std::pow(diskObstacles(0,j)-points(0,i),2)
                                + std::pow(diskObstacles(1,j)-points(1,i),2));
                    }
                    potential += std::max(overlap,0.0);
                }
            }
           
            return potential;
        }

    public:
        inline bool setup(const Eigen::Vector2d &initialP,
                          const Eigen::Vector2d &terminalP,
                          const int &pieceNum,
                          const Eigen::Matrix3Xd &diskObs,
                          const double penaWeight)
        {
            pieceN = pieceNum;
            diskObstacles = diskObs;
            obstaN = diskObs.cols();
            penaltyWeight = penaWeight;
            headP = initialP;
            tailP = terminalP;
     
            cubSpline.setConditions(headP, tailP, pieceN);

            points.resize(2, pieceN - 1);
            gradByPoints.resize(2, pieceN - 1);

            return true;
        }

        inline double optimize(CubicCurve &curve,
                               const Eigen::Matrix2Xd &iniInPs,
                               const double &relCostTol)
        {
          //TODO
            Eigen::VectorXd x(2*pieceN-2);
            for(int i = 0; i < pieceN-1; i++)
            {
                x(2*i+0) = iniInPs(0,i);
                x(2*i+1) = iniInPs(1,i);
            }
            double minCost;
            lbfgs_params.mem_size=128;
            lbfgs_params.delta = relCostTol;

            int ret = lbfgs::lbfgs_optimize(x,minCost,&PathSmoother::costFunction,nullptr,this,lbfgs_params);

            for(int i = 0; i < pieceN-1;i++)  //将优化后的值传给Points
            {
                points(0, i) = x(2*i+0);
                points(1, i) = x(2*i+1);
            }
            if(ret >=0 )
            {
                cubSpline.setInnerPoints(points);//将优化后的点生成三次样条曲线
                cubSpline.getCurve(curve);
            }
            else
            {
                minCost = INFINITY;
                std::cout << "Optimization Failed: "
                          << lbfgs::lbfgs_strerror(ret)
                          << std::endl;
            }
            return minCost;
        }
    };

}

#endif

