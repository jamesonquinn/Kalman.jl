module load julia/1.0.0-fasrc01
module load gcc/7.1.0-fasrc01 arpack/96-fasrc01
cd Kalman.jl
git pull
cd src
#julia installPackages.jl
cd ../test
julia repeatableTest.jl
