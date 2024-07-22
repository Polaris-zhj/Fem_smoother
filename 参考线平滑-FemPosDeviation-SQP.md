`FemPosDeviation`参考线平滑方法是离散点平滑方法，`Fem`是`Finite element estimate`的意思。

## 1. 优化目标

![在这里插入图片描述](https://img-blog.csdnimg.cn/5a83fe13b9c44ab6aa1a6cea962fbbb4.png#pic_center)

## 2 代价函数

**平滑考虑因素**

- 平滑性
- 几相似何性
- 均匀性

### 2.1 平滑性

参考线平滑的首要目标当然是平滑性，以三个点为例，使用向量的模$\overrightarrow{|P_2P_2^{\prime}}|$来表示，显然$|\overrightarrow{P_2P_2^{\prime}}|$越小，三个点$P_1,P_2,P_3$越接近一条直线，越平滑。其中$p_1= (x_1, y_2)，p_2 = (x_2, y_2)，p_3 = (x_3, y_3)$。

![image-20240313140912635](https://image-1314148267.cos.ap-nanjing.myqcloud.com/typora-img/image-20240313140912635.png)
$$
J_{smooth}=|\vec{P_{2}P_{2}^{\prime}}|^{2}=|\vec{P_{2}P_{1}}\:+\vec{P_{2}P_{3}}|^{2}=(x_{1}\:-x_{2}\:,y_{1}\:-y_{2}\:)^{2}+(x_{3}\:-x_{2}\:,y_{3}\:-y_{2}\:)^{2}
$$

$$
(x_{1}+x_{3}-2x_{2},y_{1}+y_{3}-2y_{2})=(x_{1},y_{1},x_{2},y_{2},x_{3},y_{3})\left.\left(\begin{matrix}{1}&{0}\\{0}&{1}\\{-2}&{0}\\{0}&{-2}\\{1}&{0}\\{0}&{1}\\\end{matrix}\right.\right) \\
记： \left.x = \left(\begin{array}{c}x_1\\y_1\\x_2\\y_2\\x_3\\y_3\end{array}\right.\right) , A_1=\begin{pmatrix}1&0&-2&0&1&0\\0&1&0&-2&0&1\end{pmatrix} \\
则 ： (x_{1}+x_{1}-2x_{2})^{2}+(y_{1}+y_{3}-2y_{2})^{2}=x^TA_{1}^TA_{1}x
$$

当有n个点时 
$$
\sum_{i=1}^{n-2}(x_{i}+x_{i+2}-2x_{i+1})^{2}+(y_{i}+y_{i+2}-2y_{i+1})^{2} \\
=(x_{1}+x_{3}-2x_{1},y_{1}+y_{3}-2y_{2},x_{2}+x_{4}-2x_{3},y_{2}+y_{4}-2y_{3},\cdots)\\(x_{1}+x_{3}-2x_{2},y_{1}+y_{3}-2y_{2},x_{2}+x_{4}-2x_{3},y_{2}+y_{4}-2y_{3},\cdots)^{T} \\
= (x_1,y_1,...) \left.\left(\begin{array}{rrrrrrrrr}1&0&&&&&&&\\0&1&&&&&&\\-2&0&1&0&&&&\\0&-2&0&1&&&&\\1&0&-2&0&1&&&\\0&1&0&-2&0&&&\\&&1&0&-2&0&&\\&&0&1&0&-2&&\\&&&&1&0&\ddots&\\&&&&0&1&&\ddots\\\\\end{array}\right.\right) \\ 
= X^TA_1^TA_1X
$$
其中：X 为 （1 × 2n ）矩阵， $A_1$为（2n ×( 2n - 4)）的矩阵。（注：$A_1$有2n行不难想象，列数由于3个点是有2列，n个点是只有最后两个点组不成3点，所以能组成（n-2）个三点组合，所以有2*（n-2）列）

平滑性代价 ： $W_{\mathrm{cost-smooth}}\cdot X^TA_1^TA_1X$。

### 2.2 几何相似性

平滑后的参考线，希望能够保留原始道路的几何信息，不会把弯道的处的参考线平滑成一条直线。使用平滑后点与原始点的距离来表示。

![image-20240313141208649](https://image-1314148267.cos.ap-nanjing.myqcloud.com/typora-img/image-20240313141208649.png)
$$
J_{deviation}=|\vec{P_{r,1}P_1}|^2+|\vec{P_{r,2}P_2}|^2+|\vec{P_{r,3}P_3}|^2= (x_1-x_{1,r})^2+(y_1-y_{1,r})^2+(x_2-x_{2,r})^2+(y_2-y_{2,r})^2+(x_3-x_{3,r})^2+(y_3-y_{3,r})^2
$$

$$
\begin{aligned}\mathrm{cost}&=\sum_{i=1}^{n}(x_{i}-x_{ir})^{2}+(y_{i}-y_{ir})^{2} \\&
=\sum_{i=1}^{n}(x_{i}^{2}+y_{i}^{2})+\sum_{i=1}^{n}(-2x_{ir}x_{i}+-2y_{ir}y_{i})+    \sum_{i=1}^{n}(x_{ir}^{2}+y_{ir}^{2})(该项为固定值可以省略)  \\&
=  (x_1,y_1,\cdots) \left.\left(\begin{array}{ccc}1&&\\&1&\\&&\ddots\end{array}\right.\right)  (x_1,y_1,\cdots)^T + (-2)(x_{1r},y_{1r},\cdots) \left.\left(\begin{array}{l}x_1\\y_1 \\ .\\.\\x_n\\y_n\end{array}\right.\right)
\end{aligned}
$$

几何相似性代价：$w_{\mathrm{cost-ref}}\cdot(X^T A_3^TA_3X+h^TX)$。  $A_3$ 为 2n×2n

其中 $h = \left.\left(\begin{matrix}-2x_{1r}\\-2y_{1r}\\\vdots\\-2x_{nr}\\-2y_{nr}\\\end{matrix}\right.\right)$。

### 2.3 均匀性

平滑后的参考线的每两个相邻点之间的长度尽量均匀一直。以三个点为例：

![image-20240313165505744](https://image-1314148267.cos.ap-nanjing.myqcloud.com/typora-img/image-20240313165505744.png)
$$
J_{length}=|\vec{P_1P_2}|^2+|\vec{P_2P_3}|^2=(x_2-x_1)^2+(y_2-y_1)^2+(x_3-x_2)^2+(y_3-y_2)^2
$$
当有n个点时
$$
\sum_{i=1}^{n-1}(x_{i}-x_{i+1})^{2}+(y_{i}-y_{i+1})^{2} \\
= (x_1-x_2,y_1-y_2,\quad x_2-x_3,\quad y_2-y_3,\cdots)\cdot(x_1-x_2,y_1-y_2,\quad x_2-x_3,\quad y_1-y_3\cdots)^T \\
= (x_1,y_1,\cdots)\left.\left(\begin{array}{rrrrrr}1&0&&&&\\0&1&&&\\-1&0&1&0&\\0&-1&0&1&\\&&-1&0&\\&&o&-1&\\&&&&\ddots&\\&&&&&\ddots\end{array}\right.\right) \\
=  X^{T}A_{2}^{T}A_{2}X
$$
均匀性代价：$W_{cost-length}\cdot x^{1}A_{2}^{T}A_{2}x$。 $A_2$ 为2n×2n

**总代价函数为**：$cost_{all} = W_{\mathrm{cost-smooth}}\cdot X^TA_1^TA_1X + W_{cost-length}\cdot x^{1}A_{2}^{T}A_{2}x + w_{\mathrm{cost-ref}}\cdot(X^T A_3^TA_3X+h^TX)$ 

## 3 约束

3.1 离散点可以在一定范围内变化

即 $ x_i - lb < x_i < x_i+ ub， y_i - lb < y_i < y_i+ ub$

3.2 前一个离散点x坐标小于后一个离散点x坐标

即 $x_i <= x_{i+1}$，前提是原路径的x不会递减， 如果出现180度拐弯的路径不要加这个约束。

**代码连接**

## 参考资料：

[自动驾驶决策规划算法第二章第二节(上) 参考线模块](https://www.cnblogs.com/zhjblogs/p/16174366.html) 

[参考线平滑-FemPosDeviation-Ipopt](https://blog.csdn.net/mpt0816/article/details/127650709?ops_request_misc=%257B%2522request%255Fid%2522%253A%2522171025161316800184123186%2522%252C%2522scm%2522%253A%252220140713.130102334.pc%255Fall.%2522%257D&request_id=171025161316800184123186&biz_id=0&utm_medium=distribute.pc_search_result.none-task-blog-2~all~first_rank_ecpm_v1~rank_v31_ecpm-6-127650709-null-null.142^v99^pc_search_result_base8&utm_term=%E7%A6%BB%E6%95%A3%E7%82%B9%E5%B9%B3%E6%BB%91fem&spm=1018.2226.3001.4187)