This folder contains the codes and data used for papers:
Shi, Yuanyuan, Bolun Xu, Yushi Tan, Daniel Kirschen, and Baosen Zhang. 
["Optimal Battery Control Under Cycle Aging Mechanisms in Pay for Performance Settings"](https://arxiv.org/pdf/1709.05715.pdf), IEEE Transactions on Automatic Control (2018).

### System requirement:
Install Matlab; 
Matlab rainflow cycle counting package (already included in the folder): https://www.mathworks.com/matlabcentral/fileexchange/3026-rainflow-counting-algorithm

### Data:
["PJM_Reg_Signal_2013_06-201405.mat"](https://drive.google.com/file/d/1zd2RGHvgBwo47X5TnjqpBYpL3fyCNLLu/view?usp=sharing) PJM fast frequency regulation signal (RegD) from 2013/06/01 to 2014/05/31 with 2 seconds resolution.

Code:

Run the online_controller.m file to check out the proposed battery online control strategy under cycle aging mechanisms. It returns the online battery power output trajectory, SoC trajectory and cycle identification result as output.

 

The main function online_controller.m will call functions init_opts.m (for parameter initialization and setting), sig2ext.m (to extract local extreme points from SoC profile), rainflow.m (to count cycle depth via rainflow algorithm), cal_cost.m (to calculate the battery degradation cost, regulation mismatch penalty & energy cost), and cal_soc.m (to calculate battery SoC using battery power traces) during its execution.


