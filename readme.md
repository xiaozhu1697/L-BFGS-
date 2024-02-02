# 一、L-BFGS的具体代码实现使用（进阶）

## 1.涉及的主要文件

程序用到的L-BFGS主要在include/gcopter/cubic_spline.hpp    path_smother.hpp

## 2.具体代码分析

### 2.1 lbfgs_evaluate_t 的代码

```cpp
/*
	lbfgs_evaluate_t对应的是path_smoother.hpp里面的costFunction
	返回值 : 代价值
	 参数  : 需要计算梯度值
*/
static inline double costFunction(void *ptr,
                                  const Eigen::VectorXd &x, //通过lbfgs_optimize传进来
                                  Eigen::VectorXd &g)       //将梯度保存在g中
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
            double energy = obj.costEnergy(gE);  //代价函数1 将梯度信息保存在gE中 返回代价值

            Eigen::Matrix2Xd gP;
            double potential = obj.constPotential(gP); //代价函数2  将梯度信息保存在gP中

            obj.gradByPoints = gP + gE;             //将两个梯度加起来

            for(int i = 0; i < obj.pieceN-1; i++)
            {
                g(2*i +0) = obj.gradByPoints(0,i);   //每个点的梯度信息保存在g中，一个x对应一个g
                g(2*i +1) = obj.gradByPoints(1,i);
            }
            double cost = energy + potential;
            return cost;
        }

        inline double costEnergy(Eigen::Matrix2Xd &gEnergy) //返回代价值，梯度信息保存在参数中
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
```

### 2.2 lbfgs_optimize 代码使用

```cpp
/*
	可以设置参数值
*/
lbfgs_params.mem_size=128;
lbfgs_params.delta = relCostTol;

int ret = lbfgs::lbfgs_optimize(x,minCost,&PathSmoother::costFunction,nullptr,this,lbfgs_params);

```

# 二 、这个程序的代码框架

从src/curve_gen.cpp中看出，先从move_base_simple/goal 中获取起始点和目标点，根据步长，获取中间点的位置，因为这个相当于分段函数，我们需要知道每个分段函数的初始点位置和梯度，才能进行L-BFGS优化然后计算L-BFGS所需要的代价函数，将点经过L-BFGS优化后，更新x，然后进行三次样条插值，让曲线平滑
