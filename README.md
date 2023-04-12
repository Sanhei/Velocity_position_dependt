# Velocity_position_dependt
The program is based on C++, for 5GB traject runs for 1 mins (Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz).
The trajectory should be one dimension.
```
mkdir build & cd build
cmake ..
make -j4
cd src
./D_trans --filename $Your_trajectory$
```
You need to change your time interval.

This program is design to get rid of transition path, so there is only the velocity position distribution.
Program will create three files: position and velocity ticks and heatmap, you can use python contour3D to plot the counts.
Likely,
![alt text](https://github.com/Sanhei/Velocity_position_dependt/blob/main/Without_Trans.png)
